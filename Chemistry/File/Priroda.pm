#!/usr/bin/perl

#writes a PRIRODA input file
# 1.04 rotation vector


package Chemistry::File::Priroda;
use base qw(Chemistry::File);
use warnings;
use strict;

our $VERSION = 1.04;

use Carp;
use Math::Trig;

use Chemistry::Mol;
use Chemistry::Atom;
Chemistry::Mol->register_format("Priroda", __PACKAGE__);

our $W = ' ' x 8;       #default whitespace

sub name_is {
    my ($class, $fname, %options) = @_;
    return $fname =~ /\.(?:in)$/i ? 1 : 0;
}

sub file_is {
    my ($class, $fname, %options) = @_;
    return $fname =~ /\.(?:in)$/i ? 1 : 0;
}

sub write_header { 
	my $self  = shift;
	my ($fh, %opts) = ($self->fh, %{$self->{opts}});
	my ($basis, $four);
	if ($opts{relativistic}) {
		$four = 1;
		$basis = 'basis4';
	} else {
		$four = 0;
		$basis = 'basis';
	}
	print $fh <<HEADER;
\$system
memory=3000
disk=-3500
path=/home/_ib/tmp
\$end
\$control
 task=optimize
 theory=DFT
 basis=/home/_ib/lib/${basis}.in
 four=$four
\$end
\$dft
functional=PBE
\$end
HEADER
}

sub write_mol {
	my ($self, $fh, $mol, %opts)  = @_;
	my $basis = $opts{basis} // "L1";
	my $charge = $opts{charge} // 0;
	my $mult = $opts{multiplicity} // 1;
	print $fh <<MOL_HEADER;
\$molecule
 charge=$charge
 mult=$mult
 cartesian
 set=$basis
MOL_HEADER
	foreach my $atom ($mol->atoms) {
        printf $fh "%-3s   %10.6f %10.6f %10.6f\n", 
                  $atom->Z, $atom->coords->array;
    }
    print $fh "\$end\n";
}

1;