#!/usr/bin/perl

use Modern::Perl '2013';
use autodie;
use Statistics::Descriptive;
use List::Util qw/min max reduce/;
use Text::CSV;
use File::Copy qw(cp);
use Cwd 'abs_path';
use YAML;
use Getopt::Long;


my ($namesfile, $outfile, $help);
GetOptions('outfile=s' => \$outfile,
           'help|h|?' => \$help,
           'names=s' => \$namesfile,
           
           );

my %names;
if ($namesfile) {
    eval{
        %names = %{YAML::LoadFile($namesfile)};
    };
    if ($@) {
        warn "Something is wrong with namesfile '$namesfile': $@";
    }
}
$outfile ||= 'summary.txt';
$outfile .= '.txt' unless $outfile =~ /\./;
           
my @par_list = qw/r_wp r_wp_dash r_p r_bragg gof penalties_weighting_K1/;
my @par_regexes;

foreach my $par (@par_list) {
    push @par_regexes, qr/($par) \s+ ([-+\d.]+)/x;
}

my %summaries;
my $startdir ||= '.';
$startdir = abs_path($startdir);
# die $startdir;
chdir $startdir;

opendir( my $startdirh, $startdir);
foreach my $subdir (grep {-d "$startdir/$_"} readdir $startdirh) {
    opendir( my $subdirh, $subdir);
    # chdir $subdir;
    my @files = readdir $subdirh;
    my $inp = (grep {/inp$/i} @files)[0];
    my $txt = (grep {/toparunning_result.*\.txt$/i} @files)[0];
    next unless $inp and $txt;
    next unless -s "$subdir/$txt" > 200;
    # say "$inp -- $txt" if $inp and $txt;
    unless (-d "$subdir/zero" and -e "$subdir/zero_bonds_summary.txt") {
        eval {mkdir "$subdir/zero"};
        cp "$subdir/$inp", "$subdir/zero/$inp";
        chdir "$subdir/zero/";
        system("toparun.pl -f -s 8,2 $inp");
        chdir '..';
        system("spaghetti zero");
        chdir $startdir;
    }
    my @properties;
    opendir( my $zero_runh, "$subdir/zero");
    foreach my $out (grep {/out$/} readdir $zero_runh) {
        open my $outh, "<", "$subdir/zero/$out";
        push @properties, get_properties($outh, \@par_regexes);
    }
    $summaries{$subdir} = reduce {abs($a->{rms_delta_d} - 0.01) < abs($b->{rms_delta_d} - 0.01) ? $a : $b} @properties;  #choose point with RMS delta d closest to 0.01 for comparison
    @properties = sort {$a->{rms_delta_d} <=> $b->{rms_delta_d}} @properties;
    $summaries{$subdir}{min_rms_delta_d} = $properties[0]->{rms_delta_d};
    $summaries{$subdir}{max_rms_delta_d} = $properties[-1]->{rms_delta_d};

    say $subdir;
    my $csv = Text::CSV->new({eol => undef, sep_char => "\t", allow_whitespace => 1});
    open my $tableh, "<", "$subdir/$txt" or die "Cannot open $txt for reading: $!\n";
    
    $csv->column_names(map {lc $_} @{$csv->getline($tableh)});  #use lc'd Mogul headers and do not care about order
    my $st = Statistics::Descriptive::Full->new();
    while (my $row = $csv->getline_hr($tableh)) {
        $st->add_data(($row->{upper} - $row->{lower})/2);
    }
    $summaries{$subdir}{HUW_mean} = $st->mean;
    $summaries{$subdir}{HUW_sd} = $st->standard_deviation;
    $summaries{$subdir}{HUW_nice} = out_prec($st->mean, $st->standard_deviation);
    $summaries{$subdir}{Total_bonds} = $st->count;
    warn "Toparunning not finished for $subdir" if $summaries{$subdir}{Total_bonds} != $summaries{$subdir}{restr_count}
}
open my $summh, ">", $outfile;

my @headers = (@par_list, qw/HUW_mean HUW_sd HUW_nice TI Total_bonds rms_delta_d min_rms_delta_d max_rms_delta_d/);
say $summh join "\t", 'Name', @headers;
say $summh join "\t", $names{$_} || $_, @{$summaries{$_}}{@headers} foreach sort keys %summaries;


sub get_properties {
    my ($inph, $regexes) = @_;
    my %factors;
    my %harmonics;
    my $search_harm;
    my ($msqd_d, $d_count, $max) = (0,0,0);
    while (<$inph>) {
        s/'.*//;
        foreach my $regex (@$regexes) {
            if (/$regex/) {
                $factors{$1} = $2;
            }
        }
        $search_harm = 1 if /PO_Spherical_Harmonics/;
        while ($search_harm and m{y(?<name>(?<order>\d)[0-9pm]+)
              \s+
            !? (?<var>\w+)_c\g{name}
               \s+
            (?<value>[-+\d.]+)}xg) {
            $harmonics{$+{var}}{$+{order}}{$+{name}} = $+{value};
        }        
        $search_harm = 0 if /\)/;


        if (/Distance_Restrain(?:_Morse|_Breakable)?\([^,]+,\s*([0-9.]+)\s*,\s*([0-9.]+)/i) {
            my $delta = ($1-$2);
            $max = abs($delta) if abs($delta) > $max;
            $msqd_d += $delta**2;
            $d_count++;
        }
    }
    $factors{restr_count} = $d_count;
    $factors{rms_delta_d} = sqrt($msqd_d/$d_count);
    my @ti;
    foreach my $var (keys %harmonics) {
        no warnings 'once';
		my $ti = 0;
		foreach my $order ( keys %{$harmonics{$var}}) {
			$ti += (reduce {$a + $b**2} 0, values %{$harmonics{$var}{$order}})/($order*2+1);
		}
		push @ti, $ti;
	}
    warn "More than one PO_Spherical_Harmonics found!\n" if @ti > 1;
    $factors{TI} = $ti[0] || 0;
    return \%factors;
}

sub out_prec {
	my ($num, $error) = @_;
	unless ($error) {
		return $num;
	}
	$error < 1 or die "Cannot display numbers with error > 1" ;
	
	my $out_prec = sprintf "%.0f", (-log10($error) + 1);
	my $error_str = sprintf "%.${out_prec}f", $error;
	my ($out_error) = ($error_str =~ /([1-9]\d*)$/);
	my $out_num = sprintf "%.${out_prec}f", $num;
	return "$out_num($out_error)";
}

sub log10 { log($_[0])/log(10)}