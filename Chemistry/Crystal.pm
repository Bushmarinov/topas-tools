#!/usr/bin/perl

package Chemistry::Crystal;

use warnings;
use strict;

our $VERSION = 1.07;

# 1.02 cell_matrix became a setter
# 1.03 added frac2cart
# 1.04 shift_origin for rigid bodies; use it in fix hydrogens
# 1.05 better restraints for CH2
# 1.06 rotation vector
# 1.07 no hydrogen restraints by default

use Math::Trig;
use Math::MatrixReal qw();
use Math::VectorReal qw();
use Chemistry::Mol;
use Scalar::Util qw(blessed);
use Carp;

require Exporter;
our @ISA = qw(Chemistry::Mol Exporter);

our @EXPORT_OK = qw(sorted_atom_pairs sorted_atom_triples mnl gen_cell_matrix);

sub mnl ($$);

foreach my $func (qw(new read parse)) {                 #auto-generate inheritance from Chemistry::Mol
    no strict 'refs';
    
    *{__PACKAGE__."::$func"}  = sub {
        my $class = shift;
        my $mol = Chemistry::Mol->$func(@_);
        bless($mol, $class);
        return $mol;
    }
}

foreach my $func (qw/bond_restraints angle_restraints rigid_bodies/){       #auto-generate functions for table listings
    no strict 'refs';
    *{__PACKAGE__."::$func"}  = sub {
        my $mol = shift;
        $mol->attr('topas/'.$func, []) unless $mol->attr('topas/'.$func);
        return @{$mol->attr('topas/'.$func)};
    }
}

sub set_restraint {                                      #sets a restrain; ovverides one previously set, if any
    my ($mol, $atoms, $value, $attrib) = @_;            # $atoms is arrayref to atom_ids or something which stringifyes to them
                                                        #attrib is a hashref with any content 
    my $table_name;
    if (@$atoms == 2) {
        $table_name = 'topas/bond_restraints';
    } elsif (@$atoms ==3) {
        $table_name = 'topas/angle_restraints';
    } else {
        croak "Please specify exaclty two or three atoms to restrain!";
    }
    my @atoms = map {$mol->by_id($_)} @$atoms;                                   #allows to use atom_ids
    @atoms[0,-1] = sort mnl @atoms[0,-1]; #swap atoms in restraint if in wrong order, allows for unique keys;
    $attrib ||= {};
    my $restraint = {atoms => \@atoms, value => $value, attrib => $attrib};
    my $registry_name = 'topas/restraint_registry';
    $mol->attr($table_name) or $mol->attr($table_name, []);
    $mol->attr($registry_name) or $mol->attr($registry_name, {});
    my $table = $mol->attr($table_name);
    my $registry = $mol->attr($registry_name);
    my $key = join "|", map {$_->id} @atoms;
    if ($registry->{$key}) {
        #my $names = join " ", map {$_->name} @atoms;
        #printf "attempt to update restraint %s (key %s) by value %.3f (previous value %.3f)\n", $names, $key, $value, $registry->{$key}{value};
        %{$registry->{$key}} = %$restraint;  
    } else {
        push @$table, $restraint;
        $registry->{$key} = $restraint;
    }
    return $restraint;
}

sub get_restraint {
    my $mol = shift;
    my @atoms;
    eval {
        if (@_ == 3 or @_ == 2) {
            @atoms = @_;
        } elsif (@_==1) {
            @atoms = @{$_[0]};
        } else {
            die "Please specify one to three arguments for get_restraint";
        }
    };
    if ($@) {
        croak "Wrong arguments for get_restraint: $@";
    }
    @atoms == 3 or @atoms == 2 or croak "Please specify 2 or 3 atoms";
    my $key = join "|", map {$_->id} @atoms;
    return $mol->attr("topas/restraint_registry")->{$key};
}



sub sorted_atom_pairs { #returns bonded atom pairs, sorted cif-style
    my ($mol, %params) = @_;
    my @ordered_atoms = sort mnl $mol->atoms;
    my %atom_order = map {; $ordered_atoms[$_] => $_} 0..$#ordered_atoms;   #store ordering in a hash to speedup and clarify code
    my @atom_pairs = grep {@$_ == 2 } map {[$_->atoms]} $mol->bonds;
    
    if ($params{exclude_hydrogens}) {
        foreach my $i (0..1) {
            @atom_pairs = grep {$_->[$i]->symbol ne 'H'} @atom_pairs;
        }
    }
    @$_ = sort { $atom_order{$a} <=> $atom_order{$b}  } @$_ foreach @atom_pairs;
    @atom_pairs = sort { $atom_order{$a->[0]} <=> $atom_order{$b->[0]}
                                            ||
                         $atom_order{$a->[1]} <=> $atom_order{$b->[1]}} @atom_pairs;
    return \@atom_pairs;
}

sub sorted_atom_triples {               # for angle A B C the triple is [A B C], as per Chemistry::Atom::angle syntax
    my ($mol, %params) = @_;
    my @ordered_atoms = sort mnl $mol->atoms;
    my %atom_order = map {; $ordered_atoms[$_] => $_} 0..$#ordered_atoms;   #store ordering in a hash to speedup and clarify code
    my @atom_pairs = @{sorted_atom_pairs($mol, %params)};
    my @atom_triples;
    foreach my $i (0..$#atom_pairs) {
        foreach my $j ($i+1..$#atom_pairs) {
            my @triple = _form_angle($atom_pairs[$i], $atom_pairs[$j]);
            push @atom_triples, \@triple if @triple;
        }
    }
    @$_[0,2] = sort {$atom_order{$a} <=> $atom_order{$b}} @$_[0,2] foreach @atom_triples;
    
    @atom_triples =  sort { $atom_order{$a->[1]} <=> $atom_order{$b->[1]}       #mimic CIF sorting
                                            ||
                            $atom_order{$a->[0]} <=> $atom_order{$b->[0]}
                                            ||
                            $atom_order{$a->[2]} <=> $atom_order{$b->[2]}
                            } @atom_triples;
    return \@atom_triples;
}

sub add_rigid_body {                                    #adds a rigid body and sets a norefine flag to corresponding atoms
                                                        # accepts [atoms], {params}
    my ($mol, $atoms, $params) = @_;
    my @atoms = map {$mol->by_id($_)} @$atoms; 
    my ($rotate, $translate, $shift_origin, $rotation_vector) = @{$params}{qw/rotate translate shift_origin rotation_vector/};
    foreach ($rotate, $translate, $shift_origin, $rotation_vector) {
        $_ ||= [0,0,0];
        unless (blessed($_) and $_->isa('Math::VectorReal')) {
            $_ = Math::VectorReal->new(@$_);
        }
    }
    
    my $rigid_body = {atoms => \@atoms, rotate => $rotate, translate => $translate, shift_origin => $shift_origin};
    if ($rotation_vector->length > 0) {
        $rigid_body->{rotation_vector} = $rotation_vector;
    }
    my $table_name = 'topas/rigid_bodies';
    $mol->attr($table_name, []) unless $mol->attr($table_name);
    push @{$mol->attr($table_name)}, $rigid_body;
    $_->attr('topas/norefine', 1) foreach @atoms;
    return $rigid_body;
}

sub fix_hydrogens {
    my $mol = shift;
	my %options = @_;
    my @atoms = sort mnl grep {$_->total_hydrogens > 0} $mol->atoms;
    foreach my $atom (@atoms) {
        my @neighbors = sort mnl $atom->neighbors;
        my @h = grep {$_->symbol eq 'H'} @neighbors;
        my @non_h = grep {$_->symbol ne 'H'} @neighbors;
        my @dummy;
        # my $shift_origin;
		if ($options{set_restraints}) {
			if (@h == 3 and @non_h == 1) {
				$mol->set_restraint([$non_h[0], $atom, $h[0]], 109.5, {type => 'hydrogen', ahardness => 1});
				$mol->set_restraint([$non_h[0], $atom, $h[1]], 109.5, {type => 'hydrogen', ahardness => 1});
				# (my $dummy_label = $atom->name) =~ s/^[A-Za-z]+/X/;
				# my $dummy_coords = ($h[0]->coords+$h[1]->coords+$h[2]->coords)/3;
				# my $dummy = $mol->new_atom(
					# name   => $dummy_label,
					# coords => $dummy_coords,
					# symbol => 'C',
					# attr   => {
						# "crystal/fract_coords" => $mol->cart2frac($dummy_coords),
						# "crystal/occupancy"    => 0,
						# "crystal/Uiso"         => 1
					# }
				# );
				# push @dummy, $dummy;
			}                                
			elsif (@h == 2 and @non_h == 2) {
				$mol->set_restraint([$non_h[0], $atom, $h[0]], 109.5, {type => 'hydrogen', ahardness => 1});
				$mol->set_restraint([$non_h[1], $atom, $h[0]], 109.5, {type => 'hydrogen', ahardness => 1});
				$mol->set_restraint([$non_h[0], $atom, $h[1]], 109.5, {type => 'hydrogen', ahardness => 1});
				$mol->set_restraint([$non_h[1], $atom, $h[1]], 109.5, {type => 'hydrogen', ahardness => 1})
			}                                                                            
			elsif (@h == 1 and @non_h == 3) {                                            
				$mol->set_restraint([$non_h[0], $atom, $h[0]], 109.5, {type => 'hydrogen', ahardness => 1});
				$mol->set_restraint([$non_h[1], $atom, $h[0]], 109.5, {type => 'hydrogen', ahardness => 1});
				$mol->set_restraint([$non_h[2], $atom, $h[0]], 109.5, {type => 'hydrogen', ahardness => 1});
			}
			elsif (@h == 1 and @non_h == 2) {
				$mol->set_restraint([$non_h[0], $atom, $h[0]], 120, {type => 'hydrogen', ahardness => 1});
				$mol->set_restraint([$non_h[1], $atom, $h[0]], 120, {type => 'hydrogen', ahardness => 1});
			} else {
				my $message = sprintf "Don't know how to generate restraints for atom %s with %d hydrogens and %d non-h\n", $atom->name, scalar @h, scalar @non_h;
				warn $message;
			}
		}
        my %param;
        croak "No fractional coordinates in atom $atom" unless $atom->attr('crystal/fract_coords');
        $param{shift_origin} = $atom->attr('crystal/fract_coords');
        if (@non_h == 1) {
            $param{rotation_vector} =  ($atom->coords - $non_h[0]->coords)->norm; 
        }        
        $mol->add_rigid_body([$atom, @h, @dummy], \%param);
    }
    return scalar @atoms;
}


sub mnl ($$) {           # sort a list of atoms by Mass, Number, Label - takes into account the fact that name and label can be unknown
    $_[1]->Z <=> $_[0]->Z
            ||
    _label_number($_[0]->name || '') <=> _label_number($_[1]->name || '')
            ||
    ($_[0]->name || '') cmp ($_[1]->name || '')
}

sub _label_number {          #returns a number for the atomic label; e.g. 12 for H12B
    my $label = shift;
    return ($label =~ /(\d+)/)[0] || 0;
}

sub _form_angle {                                      #check if the atoms in two atom pairs form an angle
    my @pairs = @_;
    @pairs == 2 or croak "Please specify exactly 2 atom pairs!";
    my @atoms;
    if (ref $pairs[0] eq 'ARRAY' and ref $pairs[1] eq 'ARRAY') {
        @atoms = map {@$_} @pairs;
    } 
    elsif ($pairs[0]->isa("Chemistry::Bond") and $pairs[1]->isa("Chemistry::Bond")) {
        @atoms = map {$_->atoms} @pairs;
    } else {
        croak "Received @pairs, which are neither two pairs nor two bonds";
    }
    
    my %atoms = map {; "$_" => $_} @atoms;            #stringification is overloaded to id for Chemistry::Atom
                                                      # see http://blogs.perl.org/users/tom_wyant/2012/01/the-case-of-the-overloaded-curlys.html
                                                      #for explanation of map {; } @array
    my %seen;
    my @unique = grep {!$seen{$_}++} @atoms;
    my @common = grep {$seen{$_} == 2} @unique;
    @common == 1 or return;                         # Bonds @bonds do not have exactly one common atom
    my @different = grep {$seen{$_} == 1} @unique;
    @different == 2 or return;                      #Bonds @bonds do not have exactly two different atoms
    return wantarray ? ($different[0], $common[0], $different[1]) : 1; 
}

# creates a matrix which will convert a frac vector to cart by $v*$cell_matrix
# by coincidence, it is just the matrix needed by VASP at the start of a CONTCAR file
sub gen_cell_matrix {                            
    my ($aa, $ba, $ca, $al, $be, $ga) = @_;
    foreach my $angle ($al, $be, $ga) {
        $angle = $angle/180*pi;
    }
    # from http://en.wikipedia.org/wiki/Fractional_coordinates
    my $v = sqrt(1-cos($al)**2-cos($be)**2-cos($ga)**2+2*cos($al)*cos($be)*cos($ga));
    my @cell_matrix = ( [$aa, $ba*cos($ga), $ca*cos($be) ],
                        [  0, $ba*sin($ga), $ca*(cos($al)-cos($be)*cos($ga))/sin($ga)],
                        [  0,            0, $ca*$v/sin($ga)] );
    my $matrix_obj = Math::MatrixReal->new_from_rows(\@cell_matrix);
    $matrix_obj = ~$matrix_obj;            #transpose, as we will use the vector in row form 
    return $matrix_obj;
}

sub parameters_from_matrix {
    my $matrix = shift;
    my @vecs = map {$matrix->matrix_row2vector($_)} 0..2;
    my ($aa, $ba, $ca) = map {$_->length} @vecs;
    my ($al, $be, $ga) = map {180*$_/pi} map {_vec_angle($vecs[$_->[0]], $vecs[$_->[1]])} ([1,2], [0,2], [0,1]);
    return {a => $aa, b => $ba, c => $ca, alpha => $al, beta => $be, gamma => $ga};
    
}

sub _vec_angle {
	my ($v1, $v2) = @_;
	return acos($v1->norm . $v2->norm);
}

sub parameters {
    my $mol = shift;
    if (@_ == 0) {
        croak "No unit cell parameters in molecule $mol" unless $mol->attr('crystal/parameters') || $mol->attr('crystal/matrix');
        unless ($mol->attr('crystal/parameters')) {
            $mol->attr('crystal/parameters', parameters_from_matrix($mol->attr('crystal/matrix')));
        }
    } elsif (@_ == 1) {
        my $param_hashref = shift;
        croak "Bad hashref supplied to parameters!" if grep {!$param_hashref->{$_}} qw/a b c alpha beta gamma/;
        $mol->attr('crystal/parameters', $param_hashref);
    } elsif (@_ == 6) {
        my @parameters = @_;
        my @names = qw/a b c alpha beta gamma/;
        my %param_hash = map {; $names[$_] => $parameters[$_]} 0..$#names;
        $mol->attr('crystal/parameters', \%param_hash);
    } else {
        croak "Wrong number of parameters for parameters method!";
    }
    
    return @{$mol->attr('crystal/parameters')}{qw/a b c alpha beta gamma/}
}

sub cell_matrix {
    my $mol = shift;
	if (@_ == 1) {
		$mol->attr('crystal/matrix', Math::MatrixReal->new_from_rows($_[0]));
	}
    unless ($mol->attr('crystal/matrix')) {
        $mol->attr('crystal/matrix', gen_cell_matrix($mol->parameters));
    }
    return $mol->attr('crystal/matrix');
}

sub volume {
    my $mol = shift;
    my $matrix = $mol->cell_matrix;
    my @vectors = map {$matrix->matrix_row2vector($_)} 0..2;
    return $vectors[0] . ($vectors[1] x $vectors[2]);
}

sub frac2cart {		#returns a Math::VectorReal object
	my $mol = shift;
	my $vec;
	if (@_ == 1) {
		$vec = shift;
		unless (blessed($vec) and $vec->isa('Math::VectorReal')) {
            $vec = Math::VectorReal->new(@$vec);
        }
	}
	elsif (@_ == 3) {
		$vec = Math::VectorReal->new(@$vec);
	}
	else {
		croak "Wrong number of arguments!\n";
	}
	return $vec * $mol->cell_matrix;
}

sub cart2frac {
	my $mol = shift;
	my $vec;
	if (@_ == 1) {
		$vec = shift;
		unless (blessed($vec) and $vec->isa('Math::VectorReal')) {
            $vec = Math::VectorReal->new(@$vec);
        }
	}
	elsif (@_ == 3) {
		$vec = Math::VectorReal->new(@$vec);
	}
	else {
		croak "Wrong number of arguments!\n";
	}
	return $vec * $mol->cell_matrix->inverse;
}

1;

