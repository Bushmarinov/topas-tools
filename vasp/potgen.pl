#!/usr/bin/perl

my $command = join ' ', 'cat', (map {"/opt/vasp/POT_PBE/$_/POTCAR"} @ARGV), '>', 'POTCAR';
print "$command\n";
print `$command`;