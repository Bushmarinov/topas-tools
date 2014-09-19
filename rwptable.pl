#!/usr/bin/perl

use warnings;
use strict;
use autodie;

my @par_list = qw/penalties_weighting_K1  r_wp r_wp_dash r_p r_p_dash r_bragg gof/;
my @par_regexes;

foreach my $par (@par_list) {
    push @par_regexes, qr/($par) \s+ ([-+\d.]+)/x;
}
my %par_values;

foreach my $file (@ARGV) {
    my $filename = $file;
    open my $fileh, "<", $file;
    while (<$fileh>) {
        foreach my $regex (@par_regexes) {
            if (/$regex/) {
                $par_values{$filename}{$1} = $2;
            }
        }
    }
}

my @output;    #table; AoA
@output = map {[sprintf '%30s', $_]} '', @par_list;     #"pretty" formatting of the parameters names and of the corner above

foreach my $file (@ARGV) {
    my $filename = $file;
    push @{$output[0]}, $filename;
    foreach my $i (0..$#par_list) {
        if (defined $par_values{$filename}{$par_list[$i]}) {
            push @{$output[$i+1]}, sprintf '%.2f', $par_values{$filename}{$par_list[$i]};
        } else {
            push @{$output[$i+1]}, 'none'
        }
        
    }
}

open my $table_file, '>par_table.txt';
foreach my $row (@output) {
    print $table_file join "\t", @$row;
    print $table_file "\n";
}
print "Table output to par_table.txt\n";