#!/usr/bin/perl -n

use warnings;
use strict;

BEGIN {
	our $msqd_d = 0;
	our $d_count = 0;
	our $max = 0;
}

#Distance_Restrain_Morse(C7		C6,		1.4007,	1.41863`,	4, 1000)
#Distance_Restrain_Morse(C5		C6,		1.3908,	1.37090`,	4, 1000)

#Angle_Restrain(C16	N15	C14,	119.97,	115.85862`,	1.0, 10)

our $msqd_d;
our $d_count;
our $max;

if (/Distance_Restrain(?:_Morse|_Breakable)?\([^,]+,\s*([0-9.]+)\s*,\s*([0-9.]+)/i) {
	my $delta = ($1-$2);
	$max = abs($delta) if abs($delta) > $max;
	$msqd_d += $delta**2;
	$d_count++;
}

END {
	printf "RMS |delta d|: %.4f\n", sqrt($msqd_d/$d_count);
    print "Total bonds: $d_count\n";
	printf "Max |delta d|: %.5g\n", $max;
}