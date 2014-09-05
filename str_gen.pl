#!/usr/bin/perl

#v 1.11 Distangle restraints
#v 1.10 MOGUL input
#v 1.09 free RB
#v 1.08 no hydrogen restraints by default
#v 1.07 help, Breakable
#v 1.07
#template tolerance

#v 1.06
#multiple matches!

use lib "C:/SAXI/SXTL";
# use autodie;

use Getopt::Long;
use STAR::Parser;
use YAML::Any;
use File::Basename;
use Chemistry::Mol;
use Chemistry::Pattern;
use Chemistry::Atom;
use Chemistry::File::XYZ;
use Text::CSV;

use Chemistry::Crystal 1.01 qw(sorted_atom_pairs sorted_atom_triples); #my module
use Chemistry::File::CIF 1.03;   #my module
use Chemistry::File::STR 1.06;   #my module
use Chemistry::File::MDLMol;
use Carp;
use Chemistry::Bond::Find qw(find_bonds);

use Modern::Perl '2013';

# my $default = "schekina.cif";

my $scale_cutoff = 1.1;
my $template_scale_cutoff = 1.1;
my ($use_morse, $write_parabolic, $need_help, $restrain_hydrogens, $distangle, $listatoms);
my $gotopt = GetOptions("tolerance|t=f", \$scale_cutoff,
						"templatetolerance|e=f", \$template_scale_cutoff,
						"parabolic|p", \$write_parabolic,
						"morse|m" => \$use_morse,
						"restrain-hydrogens|r" => \$restrain_hydrogens,
                        "distangle" => \$distangle,
                        "listatoms" =>\$listatoms,
						"help|h|?" => \$need_help
);
my $usage = <<USAGE;
Usage: str_gen.pl [options] cif_to_restrain.cif [template1] [template2]..

Template can be CIF, XYZ or MOGUL-generated TSV.
Without a template restrains bonds and angles to current values.
MOGUL templates override molecular ones.
TSVs can be created by Mogul->Load CIF->All Fragments->Export.

Options:
 --morse, -m     
    force use of Morse restraints (deprecated)
 --parabolic, -p 
    print additional file with parabolic restraints
 --restrain-hydrogens, -r
    add restraints for hydrogens in addition to rigid bodies
	(also removes limits on RB rotation if used; use with care on OH)
 --tolerance, -t <cutoff>
    Cutoff for bond-finding in the model (default 1.1)
 --templatetolerance, -e <number>
    Cutoff for bond-finding in the templates (default 1.1)
 --distangle, -d
    Set angle restraints as 1,3 distances
 --help, -h, -?  
    Displays this message
 --listatoms
    Lists which atoms from template map to where (makes sence for CIF)
USAGE

die $usage if $need_help;
my $input = (shift @ARGV) or die $usage;
my @templates = @ARGV;


#print "Scale cutoff: $scale_cutoff\n";
die $usage unless $gotopt;
(my $output = $input) =~ s/(\.cif)?$/_mol.str/;
(my $output_par = $input) =~ s/(\.cif)?$/_mol_parab.str/;
open my $logh, ">str_gen.log" or die "Cannot open file str_gen.log for writing: $!\n";
my @molecules = Chemistry::Crystal->read($input);
#print "Read successful\n";
print "The first datablock from file $input will be used.\n" if @molecules > 1;
my $mol = $molecules[0];

my $name = $mol->name;
unless ($mol->name) {
    ($name = basename($input)) =~ s/\..*//g;
    $mol->name($name);
}
#print "Name set successful\n";

# diagnostics
# my $out_xyz = $name."_test_mol.xyz";
# $mol->write($out_xyz);


find_bonds($mol, tolerance => $scale_cutoff);
#print $mol->bonds(1)->order;

my $pairs_to_restrain = $mol->sorted_atom_pairs(exclude_hydrogens => 1);


foreach (@$pairs_to_restrain) {
    my @atoms = @$_;
    $mol->set_restraint($_, ($atoms[0]->distance($atoms[1]))[0], {type => "internal"});
    #printf $logh "        %s %s %.3f\n", $_->[0]->name, $_->[1]->name, $_->[0]->distance($_->[1]);
}

my $angles_to_restrain = $mol->sorted_atom_triples(exclude_hydrogens => 1);
foreach (@$angles_to_restrain) {
    my @atoms = @$_;
    unless ($distangle) {
        $mol->set_restraint($_, ($atoms[0]->angle_deg(@atoms[1,2]))[0], {type => "internal"});
    } else {
        $mol->set_restraint([$atoms[0], $atoms[2]], ($atoms[0]->distance($atoms[2]))[0], {type => "internal distangle"});
    }
    #printf $logh "        %s %s %s %.3f\n", $_->[0]->name, $_->[1]->name, $_->[2]->name, $_->[0]->angle_deg($_->[1], $_->[2]);
}

$mol->fix_hydrogens(set_restraints => $restrain_hydrogens);

my $mogul_ext = qr/(?i:toparunning\w*\.txt|\.tsv)$/i;       #also processes toparunner output now
my (@mogul_templates, @struc_templates);
foreach (@templates) {
    if (/$mogul_ext/) {
        push @mogul_templates, $_;
    } else {
        push @struc_templates, $_;
    }
}
# push /$mogul_ext/ ? \@mogul_templates : \@struc_templates, $_ foreach @templates; #works, but strange
# say "Debug:\n".Dump([\@mogul_templates, \@struc_templates]);

foreach my $template (@struc_templates) {
    my $patt = Chemistry::Pattern->read($template);
    foreach my $h (grep {$_->symbol eq 'H'} $patt->atoms) {
        $patt->delete_atom($h);
    }
	foreach my $bad ( grep {!($_->Z)}  $patt->atoms) {
		$patt->delete_atom($bad);
	}
	foreach my $bond ($patt->bonds) {
		$patt->delete_bond($bond);
	}
    find_bonds($patt, tolerance => $template_scale_cutoff);
    my $matched =0;
    while ($patt->match($mol)) {
        $matched++;
        #print "Matched!\n";
        if ($listatoms) {
            foreach ($patt->atoms) {
                print $_->name, "=>", $_->map_to->name, "\n";
            }
        }
        my $template_pairs = sorted_atom_pairs($patt, exclude_hydrogens => 1);
        my $template_triples = sorted_atom_triples($patt, exclude_hydrogens => 1);
        foreach (@$template_pairs) {
            my @t_atoms = @$_;
            my @mol_ids = map {$_->map_to} @t_atoms;
            $mol->set_restraint(\@mol_ids, ($t_atoms[0]->distance($t_atoms[1]))[0], {type => "external from $template"});
            #$mol->set_restraint(\@mol_ids, 0, {type => "external"});
        }
        foreach (@$template_triples) {
            my @t_atoms = @$_;
            my @mol_ids = map {$_->map_to} @t_atoms;
            unless ($distangle) {
                $mol->set_restraint(\@mol_ids, ($t_atoms[0]->angle_deg(@t_atoms[1,2]))[0], {type => "external from $template"});
            } else {
                $mol->set_restraint([$mol_ids[0], $mol_ids[2]], ($t_atoms[0]->distance($t_atoms[2]))[0], {type => "external distangle from $template"});
            }
        }        
    } 
    unless ($matched) {
        print "Structure from template $template did not match!\n";
    }
}

foreach my $tablefile (@mogul_templates) {
    my $csv = Text::CSV->new({eol => undef, sep_char => "\t", allow_whitespace => 1});
    open my $tableh, "<", $tablefile or die "Cannot open $tablefile for reading: $!\n";
    $csv->column_names(map {lc $_} @{$csv->getline($tableh)});  #use lc'd Mogul headers and do not care about order
    ROW:
    while (my $row = $csv->getline_hr($tableh)) {
        # say $row->{number};
        next ROW if $row->{number} =~ /no\s+hits/i or $row->{number} =~ /excluded/i;
        if ($row->{middle}) {                       #patch for processing toparunner output
            $row->{fragment} = $row->{bond};
            $row->{median} = $row->{middle};
        }
        die "Fragment column not found in $tablefile! Are you sure that it is in correct format?\n" unless $row->{fragment};
        my @atoms = split " ", $row->{fragment};
        die "Only bond or and angle restraints supported, but found '$row->{fragment}' as Fragment in $tablefile\n" unless @atoms == 2 or @atoms == 3;
        $row->{median} or die "Median column not found in $tablefile. Are you sure that it is in correct format?\n\n";
        my @atom_ids;
        foreach my $atomname (@atoms) {
            my @found = $mol->atoms_by_name($atomname);
            unless ($found[0]) {
                print "Atom $atomname from table $tablefile not found in molecule $input\n";
                next ROW;
            }
            unless (@found == 1) {
                print "Two or more atoms with name $atomname found in molecule $input! Cannot generate restraints from $tablefile.\n";
                next ROW;
            }
            push @atom_ids, @found;
        }
        $mol->set_restraint(\@atom_ids, $row->{median}, {type => "external from Mogul table $tablefile"});
    }
}

if ($use_morse) {
	$mol->write($output, free_rigid_bodies => $restrain_hydrogens);
} else {
	$mol->write($output, breakable => 1, free_rigid_bodies => $restrain_hydrogens);
}
$mol->write($output_par, parabolic => 1, free_rigid_bodies => $restrain_hydrogens) if $write_parabolic;

printf "Structure with %s restraints output to $output\n", $use_morse ? 'Morse' : 'breakable';
print "Structure with parabolic restraints output to $output_par\n" if $write_parabolic;



