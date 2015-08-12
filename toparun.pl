#!/usr/bin/perl

use warnings;
use strict;
use autodie;
use Cwd qw();
use File::Spec;
use Getopt::Long;
# use Win32;
use Statistics::Descriptive;
# use Modern::Perl '2013';
use YAML::Any;
my $tc = $ENV{TOPAS_COMMANDLINE} || 'C:\Topas4-2\tc.exe';
my $topasdir = $ENV{TOPAS_DIR} || 'C:\Topas4-2';


my @points;
my @stepsizes;
my $force;
my $clever;
my $limit = 15;
my $nextgen;
GetOptions( "points=s", \@points,
            "stepsizes=s", \@stepsizes,
			"clever", \$clever,
            "nextgen|X" => \$nextgen,
			"limit=i", => \$limit,
			"force", \$force);

if (!$force and $^O =~ /Win32/ ) {		
    require Win32;
    Win32::GetSystemMetrics(0x1000) #SM_REMOTESESSION
       and die "You are in the remote session! Use TightVNC or --force to run TOPAS!\n";
} 
			
@points = map {split /,/} @points;
@stepsizes = map {split /,/} @stepsizes;
@points = qw/100 15 0.25/ unless @points;
@stepsizes = qw/1 0.25/ unless @stepsizes;
			
unless (scalar @points == 1+@stepsizes) {
	my $points = scalar @points;
	my $stepsizes = scalar @stepsizes;
	die "Combination $points points(@points) and $stepsizes stepsizes (@stepsizes) is impossible!"
}

@ARGV or die <<USAGE;
Usage:
   toparun.pl [options] start_inp_filename
   
   Runs topas from command line, changing the setting of K1 in process. 
   If CIF output is specified in the input file, writes CIF on each step.
   Command to run TOPAS and TOPAS working directory are defined by 
   environment variables TOPAS_COMMANDLINE and TOPAS_DIR respectively.
 
Options:
   --points, -p
     Comma-separated array of points defining the ranges with different
	 stepsizes. 
	 Default: 100,15,0.25
   --stepsizes, -s
     Self-explanatory. Must be one POSITIVE stepsize for each range.
	 Default: 1,0.25
   --clever, -c 
     Stop after finding %limit% outlying K1 points.
   --limit, -l
     Limit for --clever. Default 15.
   --force, -f
	 Allows running from RDP session
USAGE

my $filename = shift @ARGV;
$filename = Cwd::abs_path($filename);
(my $file_base = $filename) =~ s/[\d.]*\.(inp|out)$//i;


my @steps;
foreach my $i (0..$#points) {
	last if !defined $points[$i+1];
	push @steps, range($points[$i], $points[$i+1], $stepsizes[$i]);
}
push @steps, $points[-1];
#print join ", ", @steps, "\n";
my $current_file = $filename;
my $cwd = Cwd::cwd;
chdir $topasdir;
my $outlier_count = 0;
my %outlying_bonds;
my %fencevals = split /\s+/, <<TABLE;
5 5.6261925136995
6 3.5968976710598
7 2.9968395038194
8 2.9141936236358
9 3.3723204143719
10 2.931963136169
11 2.727715854628
12 2.704176821725
13 2.865768968530
14 2.669211023776
15 2.570927052216
16 2.558814438419
17 2.649954414974
18 2.530756187825
19 2.477919093629
20 2.474182867954
21 2.533845464945
22 2.459833793193
23 2.421228663830
24 2.415942917107
25 2.457848258947
26 2.408792375900
27 2.380904433451
28 2.382168350389
29 2.410482010937
30 2.369893425707
31 2.352150223843
32 2.352469518201
33 2.377720356491
34 2.346434067225
35 2.332695588445
36 2.333793147463
37 2.352135811869
38 2.329312688899
39 2.315859738387
40 2.322118347288
41 2.337192462580
42 2.318222192086
43 2.308369256524
44 2.309286475482
45 2.320533558833
46 2.308809999170
47 2.299839558111
48 2.300816351901
49 2.314275570300
50 2.298270129066
TABLE
# die Dump(\%fencevals);
K1:
foreach my $k1 (@steps) {
	open my $currenth, "<", $current_file;
	if ($clever && $k1 != $steps[0]) {
		my $restraints = read_restraints($currenth);
		$_->{error} = $_->{real} - $_->{ideal} foreach @$restraints;
		# die Dump($restraints)."\n";
		my $stats = Statistics::Descriptive::Full->new();
		$stats->add_data(map {$_->{error}} @$restraints);
		my $IQR = $stats->quantile(3) - $stats->quantile(1);
        my @outliers;
        my $multiplier = $nextgen ? $fencevals{$stats->count()} : 1.5;
        
        my $lower_fence = $stats->quantile(1) - $multiplier*$IQR;
        my $upper_fence = $stats->quantile(3) + $multiplier*$IQR;
        unless ($nextgen) {    
            # die "IQR $IQR, lower fence $lower_fence, upper fence $upper_fence\n";
            @outliers = grep {$_->{error} < $lower_fence-$_->{sigma}*2 or $_->{error} > $upper_fence+$_->{sigma}*2} @$restraints;
        } else {
            @outliers = grep {$_->{error} < $lower_fence or $_->{error} > $upper_fence} @$restraints;
        }
		# $outlier_count += @outliers;
		$outlier_count++ if @outliers;		#Incompatible change: now count instances of bad K1
		$outlying_bonds{$_->{name}}++ foreach @outliers;
		if ($outlier_count > $limit) {
			print "Stopped due to exceeded outlier limit.\nOutlying bonds are ", join(", ", keys %outlying_bonds), "\n";
			last K1;
		}
	}
	
	seek $currenth, 0, 0;
	change_k1($currenth, "$file_base$k1.inp", $k1);
	system("$tc $file_base$k1");
	$current_file = "$file_base$k1.out";
}
chdir $cwd;

sub change_k1 {
	my ($fromh, $to, $k1) =@_;
	my $result = '';
	while (<$fromh>) {
		s/(penalties_weighting_K1\s+)[\d.]+/$1$k1/;
		s/[\d.]*(\.cif)/sprintf('%.4g', $k1).$1/eg;
        s{(Out_CIF_IUCR\([^)]+?)
            [\d.]*
            \s*\)}
        {
            my ($before, $k1_f) = ($1, sprintf('%.4g', $k1));
            $k1_f =~ tr/./_/;
            $before.$k1_f.")";
        }xeg;
		$result .= $_ unless /C_matrix_normalized/../}/;
	}
	open my $toh, '>', $to;
	print $toh $result;
	return 1;
}

sub range  {		#generates steps for a range with a given stepsize. endpoint not included.
	my ($start, $finish, $stepsize) = @_;
	my $step = $start;
	my @steps;
	if ($start < $finish) {
		while (sprintf('%.3f', $step) < sprintf('%.3f', $finish)) {
			push @steps, sprintf '%.6g', $step;
			$step += $stepsize;
		}
	}
	elsif ($start > $finish) {
		while (sprintf('%.3f', $step) > sprintf('%.3f', $finish)) {
			push @steps, , sprintf '%.6g', $step;
			$step -= $stepsize;
		}
	}
	return @steps;
}

sub read_restraints {        #generates a @AoH of {name, ideal, real} from filehandle; returns referecne to it}
    my $handle = shift;
	my @bondlist;
    while (<$handle>) {
        if (m{Distance_Restrain\w*\s*\(
                 \s* ([^,]+?)        #name
                 \s*,\s*
                 ([\d.]+) (?:`[^,]*)?   #ideal
                 \s*,\s*
                 ([\d.]+) (?:`_([-.\d]+))?   #real and sigma
            }xi) {
            my ($name, $ideal, $real, $sigma) = ($1, $2, $3, $4);
            $name =~ s/\s+/ /;
            $name =~ s/\s+$//;
            push @bondlist, {name => $name, ideal => $ideal, real => $real, sigma => $sigma || 0};
        }
    }
    return \@bondlist;
}
