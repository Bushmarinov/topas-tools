#!/usr/bin/perl

use Modern::Perl '2013';
use autodie;
use Text::CSV;
use Getopt::Long;
use YAML;

my ($outfile, $help, $namesfile);

GetOptions('outfile=s' => \$outfile,
           'help|h|?' => \$help,
           'names=s' => \$namesfile,
           
           );
die <<HELP if $help or !@ARGV;
Usage: melt_toparun.pl [-o outfile] [-n names.yml] <dirnames>

Prits out toparunning result from given directories in molten style.
Accepts wildcards in dirnames even when shell does not (i.e. on Win32).

names.yml allows renaming the folders to prettify output.
HELP

my %names;
if ($namesfile) {
    eval{
        %names = %{YAML::LoadFile($namesfile)};
    };
    if ($@) {
        warn "Something is wrong with namesfile '$namesfile': $@";
    }
}
$outfile ||= 'molten_toparun.txt';
$outfile .= '.txt' unless $outfile =~ /\./;
@ARGV = grep {-d $_} map { glob $_ } @ARGV;

my @molten = [qw/Name Bond Ideal Lower Upper/];

DIR:
foreach my $dir (@ARGV) {
    my $tableh;
    eval {
        -e "$dir/toparunning_result.txt" or die "toparunning_result.txt does not exist!\n";
        -s "$dir/toparunning_result.txt" or die "toparunning_result.txt is empty!\n";
        open $tableh, "$dir/toparunning_result.txt";
        my $csv = Text::CSV->new({eol => undef, sep_char => "\t", allow_whitespace => 1});
        $csv->column_names(map {lc $_} @{$csv->getline($tableh)});  #use lc'd headers and do not care about order
        ROW:
        while (my $row = $csv->getline_hr($tableh)) {
            push @molten, [$names{$dir} || $dir, @$row{qw/bond ideal lower upper/}];
        }  
    };
    if ($@) {
        say "Problems with opening toparunning_result.txt in `$dir':\n$@Moving on...\n";
        next DIR;
    }    
}

open my $outh, '>', $outfile;

say $outh join "\t", @$_ foreach @molten;

say "Output written to $outfile";