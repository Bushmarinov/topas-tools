#!/usr/bin/perl

package Chemistry::File::CIF;
use base qw(Chemistry::File);
use warnings;
use strict;

our $VERSION = 1.03;

#1.02 - moved all math to Chemistry::Crystal

use Chemistry::Mol;
use Chemistry::Atom;
use Chemistry::Crystal;
use STAR::Parser;
use Math::MatrixReal qw();
use Math::VectorReal qw();
use Math::Trig;
use File::Temp;
use Carp;
use File::Copy qw(cp);

Chemistry::Mol->register_format("cif", __PACKAGE__);




# override the read_mol method
# universal (numeric) crystal parameters (coords, unitcell parameters, matrix etc) are written to the crystal/ namespace
# those in cif-specific format (spacegroup, symops) - to the cif/ one

sub read_mol {
    my ($self, $fh, %opts) = @_;
    my $cryst = Chemistry::Crystal->new;
    my $mol_class  = $opts{mol_class}  || 'Chemistry::Mol';
    my $mol = $mol_class->new;
    unless ($self->{datablocks}) {
        my $tmp = File::Temp->new(SUFFIX => '.cif');
        cp($fh, $tmp);
        $self->{datablocks} = [STAR::Parser->parse($tmp->filename)];
        $self->{block_counter} = 0;    
    }
    
    my $cif = $self->{datablocks}->[$self->{block_counter}++];
    return unless $cif;
    
    my $spacegroup = ($cif->get_item_data(-item=>"_space_group_name_H-M_alt"))[0]                #OLEX2 style
                      ||  ($cif->get_item_data(-item=>"_symmetry_space_group_name_H-M"))[0];     #SHELX style
    $cryst->attr('cif/spacegroup', $spacegroup);
    $mol->attr('cif/spacegroup', $spacegroup);
    
    my $name = $cif->title;
    $cryst->name($name);    
    $mol->name($name);    
    my @cell_param_names = qw/_cell_length_a    _cell_length_b    _cell_length_c    _cell_angle_alpha _cell_angle_beta  _cell_angle_gamma/;
    my @parameters = map {remove_esd($_)} map {$cif->get_item_data($_)} @cell_param_names;
    croak "The CIF file contains no unit cell parameters!" if grep {!$_} @parameters;
    $cryst->parameters(@parameters);
    $mol->attr('crystal/parameters', $cryst->attr('crystal/parameters'));
    my @symops = $cif->get_item_data('_symmetry_equiv_pos_as_xyz');        #SHELX style
    @symops or @symops = $cif->get_item_data('_space_group_symop_operation_xyz'); # OLEX2 style
    $cryst->attr('cif/symops', \@symops);
    $mol->attr('cif/symops', \@symops);
    my @cif_atoms;  # AoH; {$label, $type_symbol, $frac_x, $frac_y, $frac_z}

    #read atom data
    
    foreach my $item (qw/label type_symbol fract_x fract_y fract_z occupancy B_iso_or_equiv U_iso_or_equiv/) {
        my @data = $cif->get_item_data("_atom_site_$item");
        $cif_atoms[$_]->{$item} = remove_esd($data[$_]) foreach 0..$#data;
    }
    
    # final cleanup of atoms before adding to molecule
    
    foreach my $atom (@cif_atoms) {
        if (!$atom->{U_iso_or_equiv} and $atom->{B_iso_or_equiv}) {
            $atom->{U_iso_or_equiv} = $atom->{B_iso_or_equiv}/(8*pi**2);
        }
        $atom->{U_iso_or_equiv} ||= 0.05; # no Uiso found, give up. But 0.05 is better then 0 for visualizing

        $atom->{fract_vector} = Math::VectorReal->new(@$atom{qw/fract_x fract_y fract_z/});
        #print $atom->{fract_vector};
        $atom->{cart_vector} = $cryst->frac2cart($atom->{fract_vector});
    }
    
    foreach (@cif_atoms) {
        $mol->new_atom(
            name   => $_->{label},
            coords => $_->{cart_vector},
            symbol => $_->{type_symbol},
            attr   => {
                "crystal/fract_coords" => $_->{fract_vector},
                "crystal/occupancy"    => defined($_->{occupancy}) ? $_->{occupancy} : 1,
                "crystal/Uiso"         => $_->{U_iso_or_equiv}
            }
        );
    }
    
    return $mol;
}

sub name_is {
    my ($class, $fname, %options) = @_;
    return $fname =~ /\.(?:cif|acc)$/i ? 1 : 0;
}

sub file_is {
    my ($class, $fname, %options) = @_;
    return $fname =~ /\.(?:cif|acc)$/i ? 1 : 0;
}

# methods I need to override:

# name_is      => name check
# parse_string =>  returns @ of mols - don't need if implement read_mol correctly
# read_mol     => return next molecule 
# write_mol - maybe... write minimal cif readable from mercury.




sub remove_esd {
    my $number = shift;
    $number =~ s/^\s*([-+.0-9]+)\([0-9]+\)\s*$/$1/;
    return $number;
}




1;
