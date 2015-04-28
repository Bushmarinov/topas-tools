#!/usr/bin/perl

use Modern::Perl '2013';
use lib "C:/SAXI/SXTL";

use autodie;
use Chemistry::Mol;
use Chemistry::File::XYZ;
use Chemistry::File::SUM;
use Chemistry::File::SMILES;
use Chemistry::Bond::Find qw(find_bonds assign_bond_orders);
use Chemistry::Pattern;
use Getopt::Long;
use List::Util qw(sum);
use File::Basename;
use YAML::Any;

my $scale_cutoff = 1.1;
my ($need_help, $append, $findh, $findc, $lencut, $ecut, $kremer, $mode, $xyz);
$mode //= '';
my $gotopt = GetOptions("tolerance|t=f", \$scale_cutoff,
                        "append" => \$append,
                        "hydrogen|h" => \$findh,
                        "carbon5" => \$findc,
                        "energy=f" => \$ecut,
                        "length=f" => \$lencut,
                        "kremer" => \$kremer,
						"help|?" => \$need_help,
                        "mode=s" => \$mode,
                        "xyz" => \$xyz,
);

die <<HELP if $need_help or !@ARGV;
This program reads one or more  AIMAll .sum files and outputs charges of found fragments
in Excel-ready tab-separated format into file frags.txt, and properties of weak bonds 
between these fragments into file weakbonds.txt. Understands filesystem globs. 
Criteria for weak bonds can be applied simultaneously.

Example usage:
   chargetransfer.pl  file1.sum file2.sum
   chargetransfer.pl *.sum > frags.txt

Options:
  --tolerance, -t <cutoff>
    Cutoff for bond-finding in the model (default 1.1)
  --hydrogen, -h
    Find hydrogen bonds and split by them (default)
  --carbon5, -c
    find 'fifth' C, for CH2 contacts, and split by them (default)
  --energy, -e <E>
    Define bonds with |Econt| < E as weak
  --length, l <L>
    Define bonds with length > L as weak
  --append, -a
    Append the results to output files instead of overwriting them
  --help, -?
    Show this messsage
  --mode <mode>
    analyzes the files according to some hard-coded idea.
    e.q.: COOH
  --xyz
    Output an *.xyz file based on the *.sum one
  
HELP

unless ($findh || $lencut || $ecut || $findc || $kremer) {
    $findh = 1;
    $findc = 1;
}

my @files = map { glob $_ } @ARGV;
my $opensymb = '>';
$opensymb = '>>' if $append;
open my $weakh, $opensymb, 'weakbonds.txt';
say $weakh join "\t", qw/File Atom1 Atom2 Length Rho DelSqRho V G K Ellipticity Econt/ unless $append;
open my $fragsh, $opensymb, 'frags.txt'; 

my %COOH_report;

foreach my $file (@files) {
    # say $file;
    my $mol = Chemistry::Mol->read($file, read_cps => 1);
    
    (my $filename = basename($file)) =~ s/\..*//g;
    unless ($mol->name) {
        $mol->name($filename);
    }
    # say Dump($mol);
    # die;
    #print "Name set successful\n";

    # diagnostics
    if ($xyz) {
        my $out_xyz = $filename."_test_mol.xyz";
        $mol->write($out_xyz);
    }
    my @weak_bonds;
    if ($findh) {
        foreach my $h (grep {$_->symbol =~ /^[H|D|T]$/} $mol->atoms) {
            my @bonds = $h->bonds;
            next if @bonds <= 1;
            @bonds = sort {$a->length <=> $b->length} @bonds;
            push @weak_bonds, @bonds[1..$#bonds];
        }
    }
    if ($findc) {
        foreach my $c (grep {$_->symbol =~ /^C$/} $mol->atoms) {
            my @bonds = $c->bonds;
            next if @bonds <= 4;
            @bonds = sort {$a->length <=> $b->length} @bonds;
            push @weak_bonds, @bonds[4..$#bonds];
        }
    }
    if ($kremer) {
        push @weak_bonds, grep {-($_->attr('bader/G') +$_->attr('bader/V')) < 0 and  -4*$_->attr('bader/L') > 0} $mol->bonds;
    }
    if ($ecut) {
        push @weak_bonds, grep { abs($_->attr('bader/Econt')) < $ecut} $mol->bonds;
    }
    if ($lencut) {
        push @weak_bonds, grep { $_->length > $lencut} $mol->bonds;
    }
    { #make weak_bonds unique and sort'em by E
        my %seen;
        @weak_bonds = sort {$a->attr('bader/Econt') <=> $b->attr('bader/Econt') } grep {!$seen{$_->id}++} @weak_bonds;
    }
    if ($mode =~ /cooh/) {
        my $alt;
        $alt = 1 if $filename =~ /alt/i;
        my ($calc_name, $angle) = split /_/, $filename, 3;
        my $subkey = $angle. ($alt ? 'alt' : '');
        my %entry;
        my $hbond = $weak_bonds[0];
        $entry{hbond}{Rho} = $hbond ? $hbond->attr('bader/Rho') : 0;
        $entry{hbond}{Econt} = $hbond ? $hbond->attr('bader/Econt') : 0;
        $entry{'Econt-total'} = $hbond ? sum map {$_->attr('bader/Econt')} @weak_bonds : 0;
        $COOH_report{$calc_name}{$subkey} //= {};
        $COOH_report{$calc_name}{$subkey} = {%{$COOH_report{$calc_name}{$subkey}}, %entry};
    }
    foreach my $weak (@weak_bonds) {
        say $weakh join "\t", 
            $filename,
            (map {$_->id} $weak->atoms),
            $weak->length,
            $weak->attr('bader/Rho'),
            -4*$weak->attr('bader/L'),
            $weak->attr('bader/V'),
            $weak->attr('bader/G'),
            $weak->attr('bader/K'),
            $weak->attr('bader/Bond Ellipticity'),
            $weak->attr('bader/Econt');
        $mol->delete_bond($weak);
    }
    # find_bonds($mol, tolerance => $scale_cutoff);
    my @frags = $mol->separate;
    # say Dump(\@frags);
    my @output = ($filename);
    foreach my $frag (@frags) {
        $frag->attr("bader/q", sum map {$_->attr("bader/q")} $frag->atoms);
        assign_bond_orders($frag);
        push @output, $frag->formula, $frag->print(format => 'smiles'), $frag->attr('bader/q');
    }
    say $fragsh join "\t", @output;
    
    if ($mode =~ /cooh/) {
        my $alt = 1 if $filename =~ /alt/i;
        my ($calc_name, $angle) = split /_/, $filename, 3;
        my $subkey = $angle. ($alt ? 'alt' : '');
        my $patt_str = '[H]OC(O)C([H])([H])[H]';
        my $patt = Chemistry::Pattern->parse($patt_str, format => 'smiles');
        die "Not found acid!" unless $patt->match($mol);
        my @atoms = $patt->atom_map;
        my @bonds = $patt->bond_map;
        # say join "\t", map {sprintf '%.3f', $_->length} @bonds;
        my %entry;     # hash to output in report
        @{$entry{q}}{qw/H O1 C O2/} = map {$_->attr('bader/q')} @atoms[0..3];
        $entry{q}{CH3} = sum map {$_->attr('bader/q')} @atoms[4..7];
        $entry{energy} = $mol->attr('wfn/energy');
        $entry{'O-H'}{len} = $bonds[0]->length;
        $entry{'CO'}{len} = $bonds[2]->length;
        $COOH_report{$calc_name}{$subkey} //= {};
        $COOH_report{$calc_name}{$subkey} = {%{$COOH_report{$calc_name}{$subkey}}, %entry};
        # say join "\t", $_->id, $_->attr('bader/q');
        # die "in peace";
    }
}

if ($mode =~ /cooh/) {
    open my $reporth, '>report.txt';
    
    my @headers = qw/energy q_H q_O1 q_C q_O2 q_CH3 O-H_len CO_len hbond_Rho hbond_Econt Econt-total/;
    say $reporth join "\t", 'Name', (map {"0_$_"} @headers), (map {"180_$_"} @headers);
    foreach my $calcname (sort keys %COOH_report) {
        my @line;
        push @line, $calcname;
        foreach my $angle (0, 180) {
            foreach my $header (@headers) {
                my @hierarchy = split /_/, $header;
                my $record = $COOH_report{$calcname}{$angle};
                $record = $record->{$_} foreach @hierarchy;  # process subkeys of any depth
                push @line, $record;
            }
        }
        say $reporth join "\t", @line;
    }
}


say "Output written to frags.txt and weakbonds.txt";