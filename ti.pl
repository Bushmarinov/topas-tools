#!/usr/bin/perl

use warnings;
use strict;
use List::Util qw/reduce/;
use YAML;
use v5.14;

my %harmonics;

die "Please specify a TOPAS *.out file!\n" unless @ARGV;

while (<ARGV>) {
	while (m{y(?<name>(?<order>\d)[0-9pm]+)
		  \s+
		!? (?<var>\w+)_c\g{name}
		   \s+
		(?<value>[-+\d.]+)}xg) {
		$harmonics{$+{var}}{$+{order}}{$+{name}} = $+{value};
	}
}
print YAML::Dump(\%harmonics);
{
no warnings 'once';
	foreach my $var (keys %harmonics) {
		my $ti = 0;
		foreach my $order ( keys %{$harmonics{$var}}) {
			$ti += (reduce {$a + $b**2} 0, values %{$harmonics{$var}{$order}})/($order*2+1);
		}
		say "For $var TI = $ti";
	}
}
__END__
		PO_Spherical_Harmonics(, 4 load sh_Cij_prm {
			y00   !m51c97912_3_c00  1.00000
			y20   m51c97912_3_c20   0.07364`
			y22p  m51c97912_3_c22p  0.16495`
			y22m  m51c97912_3_c22m -0.09070`
			y40   m51c97912_3_c40   0.04427`
			y42p  m51c97912_3_c42p  0.00914`
			y42m  m51c97912_3_c42m -0.06180`
			y44p  m51c97912_3_c44p  0.02603`
			y44m  m51c97912_3_c44m -0.00153`
			} )