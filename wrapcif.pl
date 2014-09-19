#!/usr/bin/perl

use warnings;
use strict;
use Text::Wrap qw(wrap fill $columns);

$columns = 80;
my $cif_file = $ARGV[0] || 'a.cif';

(my $cif_wrapped = $cif_file) =~ s/\./_wrapped./;
open CIF, "<", $cif_file or die "Cannot open file $cif_file: $!";
open WRAPPED, ">", $cif_wrapped  or die "Cannot open file $cif_wrapped: $!";

my @buffer;
my $is_text = 0;

while (<CIF>) {
	if (/^;/) {
		if ($is_text) {
			print WRAPPED fill("", "", @buffer), "\n";
			@buffer = ();
			$is_text = 0;;
		} else {
			$is_text = 1;
		}
		print WRAPPED $_;
		next;
	}
	if ($is_text) {
		s/^[ \t]+//;
		#chomp if /\S/;
		push @buffer, $_;
		next;
	}
	print WRAPPED $_;
}
	
