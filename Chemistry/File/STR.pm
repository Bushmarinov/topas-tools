#!/usr/bin/perl

#writes a TOPAS str file
# 1.04 rotation vector
#1.05 breakable
#1.06 free RB
#1.08 loaded atoms

package Chemistry::File::STR;
use base qw(Chemistry::File);
use warnings;
use strict;

our $VERSION = 1.08;
our $loaded_atoms = 1;
use Carp;
use Math::Trig;

use Chemistry::Mol;
use Chemistry::Atom;
Chemistry::Mol->register_format("str", __PACKAGE__);

our $W = ' ' x 8;       #default whitespace

sub name_is {
    my ($class, $fname, %options) = @_;
    return $fname =~ /\.(?:str)$/i ? 1 : 0;
}

sub file_is {
    my ($class, $fname, %options) = @_;
    return $fname =~ /\.(?:str)$/i ? 1 : 0;
}

sub write_mol {
    my ($self, $fh, $mol, %opts) = @_;
    $mol->isa("Chemistry::Crystal") or croak "Cannot print non-crystal as STR!\n";
    print $fh topas_header($mol); 
    foreach my $atom ($mol->atoms) {
        print $fh $loaded_atoms ? topas_loaded_atom($atom) : topas_atom($atom);
    }    
    print $fh "${W}}\n" if $loaded_atoms;
    #process restraints
    print $fh "\n${W}'Bond restraints\n";
    my %bond_by_type;
    foreach my $bond_r ($mol->bond_restraints) {
        push @{$bond_by_type{$bond_r->{attrib}{type}}}, $bond_r;
    }
    my @internal_types = grep {/internal/i} keys %bond_by_type;
    my @hydrogen_types = grep {/hydrogen/i} keys %bond_by_type;
    my @other_types = grep {!(/internal/i || /hydrogen/i)} keys %bond_by_type;
    foreach my $type (@internal_types, @other_types, @hydrogen_types) {
        print $fh "${W}'$type\n";
        foreach my $bond_r (@{$bond_by_type{$type}}) {
            my @atoms = @{$bond_r->{atoms}};
            my @names = map {$_->name} @atoms;
            #Distance_Restrain_Morse(O3        C16,    1.2685,    1.33334`,    4, 1000)
            if ($opts{parabolic}) {
                printf $fh "${W}Distance_Restrain(%4s  %4s, %6.3f, %8.5f`, 0.00001, 16000)\n", @names, $bond_r->{value}, ($atoms[0]->distance($atoms[1]))[0];
            } elsif ($opts{breakable}) {
				printf $fh "${W}Distance_Restrain_Breakable(%4s  %4s, %6.3f, %8.5f`, 4, 1000)\n", @names, $bond_r->{value}, ($atoms[0]->distance($atoms[1]))[0];
			}
			else {
                printf $fh "${W}Distance_Restrain_Morse(%4s  %4s, %6.3f, %8.5f`, 4, 1000)\n", @names, $bond_r->{value}, ($atoms[0]->distance($atoms[1]))[0];
            }
        }
        print $fh "\n";
    }
    
    print $fh "\n${W}'Angle restraints\n";
    my %angle_by_type;
    foreach my $angle_r ($mol->angle_restraints) {
        push @{$angle_by_type{$angle_r->{attrib}{type}}}, $angle_r;
    }
    @internal_types = grep {/internal/i} keys %angle_by_type;
    @hydrogen_types = grep {/hydrogen/i} keys %angle_by_type;
    @other_types = grep {!(/internal/i || /hydrogen/i)} keys %angle_by_type;
    foreach my $type (@internal_types, @other_types, @hydrogen_types) {
        print $fh "${W}'$type\n";
        foreach my $angle_r (@{$angle_by_type{$type}}) {
            my @atoms = @{$angle_r->{atoms}};
            my @names = map {$_->name} @atoms;
            #Angle_Restrain(C16    N15    C14,    119.97,    115.85862`,    1.0, 10)
            printf $fh "${W}Angle_Restrain(%4s  %4s  %4s, %8.3f, %10.5f`, 1.0, %.1f)\n", 
			                                      @names, 
														  $angle_r->{value}, 
																 ($atoms[0]->angle_deg(@atoms[1,2]))[0], 
																			   $angle_r->{attrib}{ahardness} || 10;
        }
        print $fh "\n";
    }
    print $fh "\n${W}'Rigid bodies\n";
=for comment
           rigid 
               load point_for_site ua ub uc {
                  C2 -0.20073   0.06725   0.65563 
                  H2 -0.198361  0.102509  0.565456
                }
          Rotate_about_axies(@  4.33057`, @  1.93286`, @ -2.51469`)
          Translate(@ -0.03452`, @  0.03983`, @  0.01914`)
=cut
    my $rb_rot_limits = '';
    my $ref_H = '';
	if ($opts{free_rigid_bodies}) {
		$rb_rot_limits = '';
        $ref_H = '@';
	}
	
    foreach my $rigid ($mol->rigid_bodies) {                        #write rigids in Cartesian: we have these anyway
																	#substract shift_origin
        my @atoms = @{$rigid->{atoms}};								
        print $fh "${W}rigid\n";
        print $fh "${W}    load point_for_site ux uy uz {\n";
        foreach (@atoms) {
            printf $fh "${W}${W}%4s %7.4f %7.4f %7.4f\n", 
			         $_->name, $mol->frac2cart($_->attr('crystal/fract_coords') - $rigid->{shift_origin})->array;
        }
        print $fh "${W}    }\n";
        if ($rigid->{rotation_vector}) {
            printf $fh "${W} rotate @ 0 qx %8.5f qy %8.5f qz %8.5f\n", $rigid->{rotation_vector}->array;
        }
        printf $fh "${W}Rotate_about_axies($ref_H  %8.5f $rb_rot_limits, $ref_H  %8.5f $rb_rot_limits, $ref_H %8.5f $rb_rot_limits)\n", $rigid->{rotate}->array;
        printf $fh "${W}Translate(@  %8.5f, @  %8.5f, @ %8.5f)\n", ($rigid->{translate} + $rigid->{shift_origin})->array;
        print $fh "\n";
    }
}

sub topas_header {
    my $mol = shift;
    my $name = $mol->name;
    (my $spacegroup = $mol->attr('cif/spacegroup') || '?') =~ s/\s//g;
    my ($aa, $ba, $ca, $al, $be, $ga) = $mol->parameters;
    $_ = sprintf "%.4f", $_ foreach $aa, $ba, $ca;
    $_ = sprintf "%.3f", $_ foreach $al, $be, $ga;
    my $vol = sprintf '%.2f', $mol->volume;
	my $string = <<HEADER;
    str 
        CS_L(@, 300)
        penalties_weighting_K1 1000
        'TCHZ_Peak_Type(@,-0.03374`,@, 0.03592`,@,-0.00515`,, 0,@, 0.14017`,, 0)
        r_bragg  10
        scale @  0.0002
        phase_name "$name"
        space_group $spacegroup
        MVW( 1738.158, $vol, 100.000`)
        
        Out_CIF_STR("$name-topas.cif")
        Phase_LAC_1_on_cm( 6.70350`)
        Phase_Density_g_on_cm3( 1.24456`)
        a  $aa
        b  $ba
        c  $ca
        al $al
        be $be
        ga $ga
        
HEADER
	my %seen; 
	my @beqtypes = grep {!$seen{$_}++} map {lc $_->symbol} sort {$b->Z <=> $a->Z} $mol->atoms;
	foreach (@beqtypes) {
		$string .= "${W}prm beq$_ 1.000\n";
	}
	$string .= "\n";
	foreach (@beqtypes) {
		$string .= "${W}prm =beq$_/(8*Pi^2); : 0.05\n";
	}
	$string .= "\n";
    $string .= "${W}load site num_posns x y z occ beq {\n" if $loaded_atoms;
	return $string;
	
}

#sub topas_atom {
#    my ($atom) = @_;
#    croak "No fractional coordinates in atom $atom" unless $atom->attr('crystal/fract_coords');
#    my $fmt= "        site %-4s  x \@ % .5f y \@ % .5f z \@ % .5f occ %s %.3f beq %.3f\n";
#    $fmt =~ s/\@/ /g if $atom->attr('topas/norefine');
#    return sprintf $fmt, $atom->name, $atom->attr('crystal/fract_coords')->array, $atom->symbol, 
#                          $atom->attr('crystal/occupancy'), $atom->attr('crystal/Uiso')*(8*pi**2);
#}

sub topas_atom {
    my ($atom) = @_;
    croak "No fractional coordinates in atom $atom" unless $atom->attr('crystal/fract_coords');
    my $fmt= "        site %-4s num_posns 0 x \@ % .5f y \@ % .5f z \@ % .5f occ %s %.3f beq %s\n";
    $fmt =~ s/\@/ /g if $atom->attr('topas/norefine');
    return sprintf $fmt, $atom->name, $atom->attr('crystal/fract_coords')->array, $atom->symbol, 
                          $atom->attr('crystal/occupancy'), "=beq".lc($atom->symbol).";";
}

sub topas_loaded_atom {
    my ($atom) = @_;
    croak "No fractional coordinates in atom $atom" unless $atom->attr('crystal/fract_coords');
    my $fmt= "${W}    %-4s 0 \@ % .5f \@ % .5f \@ % .5f %s %.3f %s\n";
    $fmt =~ s/\@/ /g if $atom->attr('topas/norefine');
    return sprintf $fmt, $atom->name, $atom->attr('crystal/fract_coords')->array, $atom->symbol, 
                          $atom->attr('crystal/occupancy'), "=beq".lc($atom->symbol).";";
}

1;
