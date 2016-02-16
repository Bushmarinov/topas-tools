#!/usr/bin/perl

use warnings;
use strict;

use Chemistry::Elements qw(get_Z get_symbol);
use Math::Trig;
use YAML qw();

@ARGV > 1 or die <<HELP;
CRYSCOMPAR.PL

   Performs a comparison of two structures 
   according to rules defined by the "DFT verification"
   Acta B article.
   
   Accepts two .res files.
HELP

my @files = @ARGV;
my @atom_lists;  #AoAoH
my @matrices;
foreach my $file (@files) {
	my @cell_matrix;
	my @atom_list;
	open INPUT, '<', $file or die "Cannot open file `$file': $!\n";

	while (<INPUT>) {
		my $el_num;
		my $el_coord;
		#convert unit cell to lattice vectors
		# CELL 0.71073 11.25746 8.15014 6.22333 81.2526 122.7034 107.8732
		#ZERR 1        0.00039 0.00022 0.00022  0.002    0.0016   0.0022

		if (/^CELL \s+ [0-9.]+  ((?:\s+[0-9.]+){6})/xi) {
			my ($aa, $ba, $ca, $al, $be, $ga) = split " ", $1;
			foreach my $angle ($al, $be, $ga) {
				$angle = $angle/180*pi;
			}
			my $v = sqrt(1-cos($al)**2-cos($be)**2-cos($ga)**2+2*cos($al)*cos($be)*cos($ga));
			@cell_matrix = ( [$aa, $ba*cos($ga), $ca*cos($be) ],
							 [  0, $ba*sin($ga), $ca*(cos($al)-cos($be)*cos($ga))/sin($ga)],
							 [  0,            0, $ca*$v/sin($ga)] );
			push @matrices, \@cell_matrix;
		}
		next if /AFIX/i;
        next if /^\s*H\d+/;
		if (/^([A-Za-z]{1,2})(?![A-Za-z])\w* \s+ \d+ \s+ ((?:[-0-9.]+\s+){3})/xi) {
			(my $el_symbol, $el_coord) = ($1, $2);
			$el_num = get_Z($el_symbol);	
			#next unless $el_num > 1;
			
			my $abc = [split " ", $el_coord];
			#print "$el_symbol @$abc\n";
			#print YAML::Dump($abc, [$abc], transpose([$abc]), \@cell_matrix);
			
			
			push @atom_list, {symbol => $el_symbol, Z => $el_num, abc => $abc};
			#print "$el_symbol @$xyz\n";
		} 
		
	}
	push @atom_lists, \@atom_list;
}

#rework. Cartesian diff is currently wrong if matrices are different
my $sum;
my $count = 0;
foreach my $id (0..$#{$atom_lists[0]}){
	die "lists are incompatible\n" unless $atom_lists[0][$id] 
											and $atom_lists[1][$id] 
											and $atom_lists[0][$id]->{Z} == $atom_lists[1][$id]->{Z};
	next unless $atom_lists[0][$id]->{Z} > 1;
	$count++;
	my $diff = 0;
	my $abc_diff = MinTranslDiff($atom_lists[0][$id]{abc}, $atom_lists[1][$id]{abc});
	$diff += Modulus(FracToCart($abc_diff,$_))/2 foreach @matrices[0,1];
	$sum += $diff**2;
}
my $rms_diff = sqrt($sum/$count);
printf "RMS diff is %.3f\n", $rms_diff;

sub MinTranslDiff {
	my ($v1, $v2) = @_;
	my $result = VecDiff($v1,$v2);
	foreach (@$result) {
		$_ -= 1 while $_ > 0.5;
		$_ += 1 while $_ <= -0.5 ;
	}
	return $result;
}

sub FracToCart {
	my ($v, $m) = @_;
	return transpose(MatrixMult($m, transpose([$v])))->[0];
}

sub MatrixMult {
	my ($m1, $m2) = @_;
	my $result;
	my $dim1 = scalar @{$m1->[0]};
	my $dim2 = scalar @$m2;
	scalar @{$m1->[0]} == scalar @$m2 or die "cannot multiply matrices of different dimensionality $dim1 and $dim2";
	foreach my $i (0..$#$m1) {
		foreach my $j (0..$#{$m2->[0]}) {
			foreach my $k (0..$#$m2) {
				$result->[$i][$j] ||= 0;
				$result->[$i][$j] += $m1->[$i][$k]*$m2->[$k][$j];
			}
		}
	}
	return $result;
}

sub ScalarMult {
	my ($v1, $v2) = @_;
	scalar @$v1 == scalar @$v2 or die "cannot multiply vectors of different dimensionality";
	my $sum =  0;
	foreach my $i (0..$#$v1) {
		$sum += $v1->[$i]*$v2->[$i];
	}
	return $sum;
}

sub Modulus {
	my $vec = shift;
	return sqrt(ScalarMult($vec, $vec));
}

sub VecByScalar {
	my ($vec, $scalar) = @_;
	my @vector = @$vec;	#dereference
	foreach my $i (0..$#vector) {
		$vector[$i] *= $scalar;
	}
	return \@vector;
}

sub VecAdd {
	my ($v1, $v2) = @_;
	my @vector;
	scalar @$v1 == scalar @$v2 or die "cannot add vectors of different dimensionality";
	foreach my $i (0..$#$v1) {
		$vector[$i] = $v1->[$i] + $v2->[$i];
	}
	return \@vector;
}

sub VecDiff {
	my ($v1, $v2) = @_;
	return VecAdd($v1, VecByScalar($v2, -1));
}

sub BisectralScaled {
	my ($vec1, $vec2) = @_;
	my $bisectral = VecAdd(VecByScalar($vec1, Modulus($vec2)),VecByScalar($vec2, Modulus($vec1)));
	$bisectral = VecByScalar($bisectral, 1/Modulus($bisectral));
	return $bisectral;
}

sub transpose {
	my $matrix = shift;
	my $new_matrix;
	foreach my $i (0..$#$matrix) {
		foreach my $j (0..$#{$matrix->[$i]}) {
			$new_matrix->[$j][$i] = $matrix->[$i][$j]
		}
	}
	return $new_matrix;
}