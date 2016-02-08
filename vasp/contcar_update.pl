#!/usr/bin/perl

use Modern::Perl '2015';
use lib 'C:\SAXI\SXTL';

use File::Copy qw/move/;
use YAML qw();
use Chemistry::Crystal 1.02;



my @coords;
my $resfile = shift @ARGV or die "Usage: contcar_update.pl <res-file>\n";

open my $resh, '<', $resfile 			or die "cannot open file `$resfile': $!\n";
my $tmpfile = $resfile.'.tmp';
open my $resouth, '>', $tmpfile or die "cannot open file `$tmpfile': $!\n";
open my $contcarh, 'CONTCAR' 			or die "cannot open CONTCAR file: $!\n";



my @matrix;
while (<$contcarh>) {
	push @matrix, [split " ", $_] if $. >= 3 and $. <= 5;
	if (/Direct/i..!/\S/ and !/Direct/i) {
		chomp;
		s/(-?\d*\.\d+)/sprintf '%.6f', $1/eg;
		push @coords, $_;
	}
}
close $contcarh;

my $crystal = Chemistry::Crystal->new();
$crystal->cell_matrix(\@matrix);
my @parameters = $crystal->parameters;
my $celline = join " ", ((map {sprintf '%.4f', $_} @parameters[0..2]), (map {sprintf '%.3f', $_} @parameters[3..5]));
#print YAML::Dump(\@matrix);
=for old_times_sake
my ($aa, $ba, $ca) = map {sprintf '%.4f', $_} map {Modulus($_)} @matrix;
my ($al, $be, $ga) = map {sprintf '%.3f', $_} map {180*$_/pi} map {VecAngle($matrix[$_->[0]], $matrix[$_->[1]])} ([1,2], [0,2], [0,1]);
my $celline = "$aa $ba $ca $al $be $ga";
=cut 

my $atom_count=0;
while (<$resh>) {
	s{^(cell \s+ [\d.]+  \s+)			#start
	   ((?:[\d.]+\s+){5}[\d.]+)			#parameters
	}
	{$2 ne $celline ? "$1$celline\nREM OLD $1$2\n": "$1$2\n"}ixe;
	s/^((?:[A-Za-z]{1,2})(?![A-Za-z])[\w']*    	#symbol
		\s+ \d+	 								#type
	   \s+) ((?:[-0-9.]+\s+){2}[-0-9.]+)		# $2: coords
	/$1.$coords[$atom_count++]/xei unless /END/..eof;
	print $resouth $_;
}
close $resouth;
close $resh;

unlink $resfile.'.old';

#let's do some reasonable 'old' generation
# 'old' file is a copy, isn't it?
# obviously, we want to have one true 'old' fil
opendir my $dirh, '.';
(my $copy_pattern = $resfile) =~ s{^([^.]+)(\..*)}
                                  {'('.quotemeta($1).'_old)(\d*)('.quotemeta($2).')'}e;
$copy_pattern = qr/$copy_pattern/i;
my @file_list = readdir $dirh;
#say $basename;
#say join "\n", @file_list;
my $newest_copy = (sort {$b cmp $a} grep {$_ =~ /$copy_pattern/} @file_list)[0];
if ($newest_copy) {
    # die "HERE: $newest_copy";
    $newest_copy =~ s/$copy_pattern/$1.($2 ? $2+1 : 1).$3/e;
} else {
    # die "THERE";
    ($newest_copy = $resfile) =~ s/^([^.]+)(\..*)/${1}_old$2/;
}

move($resfile, $newest_copy);
move($tmpfile, $resfile);
say "File $resfile updated. Copy of the old file is now in $newest_copy";
