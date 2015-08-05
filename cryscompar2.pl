#!/usr/bin/perl

use lib "C:/SAXI/SXTL";
use YAML;
use Modern::Perl '2013';
use Chemistry::Crystal 1.07;
use List::Util 'reduce';
use Chemistry::File::CIF;
use autodie;

die <<HELP unless @ARGV ==2;
Usage: cryscompar2.pl <cif1> <cif2>

Calculates RMSd by van der Streek and Neumann.
Works best if both CIF files are converted to P1 with pack cell in OLEX2.
May provide spurious results otherwise.
HELP

my ($file1, $file2) = @ARGV;
my ($crys1, $crys2) = map {Chemistry::Crystal->read($_)} ($file1, $file2);

printf "RMS diff is %.4f A\n", crys_compar($crys1, $crys2);

sub crys_compar {
    my ($c1, $c2) = map {$_->clone} @_;
    foreach my $cr ($c1, $c2) {
        $cr->delete_atom($_) foreach grep {$_->Z == 1} $cr->atoms;
    }
    die "Unequal number of non-hydrogen atoms!\n" unless scalar (@{[$c1->atoms]}) == scalar (@{[$c2->atoms]});
    my $sumsq = 0;
    my $count = 0;
    foreach my $atom1 ($c1->atoms) {
         my @closest = sort { $a->{dist} <=> $b->{dist} }
                        map { {%$_, dist => $c2->frac2cart($_->{diff})->length} }
                        map {
                              { 
                                atom => $_,
                                diff =>  min_transl_diff($_->attr("crystal/fract_coords"), $atom1->attr("crystal/fract_coords")),
                              }
                            } $c2->atoms;
        # die "nooo\n";
        my ($atom2, $diff) = @{$closest[0]}{qw/atom diff/};
        printf "%s %s %.4f\n", $atom1->name, $atom2->name, $closest[0]->{dist};
        my $avdist = ($c1->frac2cart($diff)->length + $c2->frac2cart($diff)->length)/2;
        
        $sumsq += $avdist**2;
        $count++;
        $c1->delete_atom($atom1);
        $c2->delete_atom($atom2);
    }
    return sqrt($sumsq)/$count;
}

sub min_transl_diff {
	my ($v1, $v2) = @_;
	my $result = $v1 - $v2;
    my @abc = $result->array;
	foreach (my @abc) {
		$_ -= 1 while $_ > 0.5;
		$_ += 1 while $_ <= -0.5 ;
	}
    $result = Math::VectorReal->new(@abc);
    # die Dump($result) if rand() < 0.1;
	return $result;
}

