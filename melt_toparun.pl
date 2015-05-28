#!/usr/bin/perl

use Modern::Perl '2013';
use autodie;
use Text::CSV;

die <<HELP unless @ARGV;
Usage: melt_toparun.pl <dirnames>

Prits out toparunning result from given directories in molten style.
Accepts wildcards in dirnames even when shell does not (i.e. on Win32).
HELP

@ARGV = grep {-d $_} map { glob $_ } @ARGV;

my @molten = [qw/Name Bond Ideal Lower Upper/];

DIR:
foreach my $dir (@ARGV) {
    my $tableh;
    eval {
        open $tableh, "$dir/toparunning_result.txt";
        my $csv = Text::CSV->new({eol => undef, sep_char => "\t", allow_whitespace => 1});
        $csv->column_names(map {lc $_} @{$csv->getline($tableh)});  #use lc'd headers and do not care about order
        ROW:
        while (my $row = $csv->getline_hr($tableh)) {
            push @molten, [$dir, @$row{qw/bond ideal lower upper/}];
        }  
    };
    if ($@) {
        say "Problems with opening toparunning_result.txt\nin `$dir':\n$@Moving on...";
        next DIR;
    }    
}

open my $outh, '>', 'molten_toparun.txt';

say $outh join "\t", @$_ foreach @molten;

say "Output written to molten_toparun.txt";