#!/usr/bin/perl

package Chemistry::File::SUM;
use base qw(Chemistry::File);
use warnings;
use strict;

our $VERSION = 1.00;


use Chemistry::Mol;
use Chemistry::Atom;
use Carp;
use File::Copy qw(cp);

Chemistry::Mol->register_format("sum", __PACKAGE__);

use constant BOHR_RADIUS => 0.529177249;
our $scinum = qr/[+-]?\d\.\d{10}E[-+]\d{2}/i;

# override the read_mol method
# universal (numeric) crystal parameters (coords, unitcell parameters, matrix etc) are written to the crystal/ namespace
# those in cif-specific format (spacegroup, symops) - to the cif/ one

sub read_mol {
    my ($self, $fh, %opts) = @_;
    return if $fh->eof;
    my $mol_class  = $opts{mol_class}  || "Chemistry::Mol";
    my $mol = $mol_class->new;
    my %atoms;
    my @cps;
    while (<$fh>) {
        if (/^Wf[xn] Title:\s+(.*)$/i) {
            (my $name = $1) =~ s/\s*$//;
            $mol->name($name);
        }
        if (/molecular/i and /E\(Mol\)/i) {
            $mol->attr('wfn/energy', (split)[-1]);
        }
        if (/^Nuclear Charges and Cartesian Coordinates/../^Some Atomic Properties/) {
            if (/^(?<name>(?<elem>[A-Za-z]+)\d*) \s+ [\d.]+ \s+ (?<X>$scinum) \s+ (?<Y>$scinum) \s+ (?<Z>$scinum)\s*$/x) {
                $atoms{$+{name}} = {%+};
            }
        }
        if (/^Some Atomic Properties/../^Total/) {
            # die "horribly";
            if (/^(?<name>[A-Za-z]+\d*) \s+ (?<q>$scinum) \s+ (?<L>$scinum) \s+ (?<K>$scinum) \s+ (?<K_Scaled>$scinum) \s+ (?<Mu_Intra>$scinum)\s*$/x) {
                # die "horribly";
                $atoms{$+{name}} = {%{$atoms{$+{name}}}, %+};
            }
        }
        if ($opts{read_cps} and /Critical Point Analysis/i) {
        #TODO:  implement OPTIONAL reading of bonding information from CPs
            @cps = _read_CPs($fh);  #we are close to the end of file anyway
            last;
        }
    }
    foreach my $name (keys %atoms) {
        my %custom_attr;
        my %exclude = map {$_ => 1} qw/X Y Z name elem/;
        foreach my $key (grep {!$exclude{$_}} keys %{$atoms{$name}}) {
            $custom_attr{"bader/$key"} = $atoms{$name}{$key};
        }
        $mol->new_atom(
            id => $name,
            symbol => $atoms{$name}{elem},
            coords => [map {$_*BOHR_RADIUS} @{$atoms{$name}}{qw/X Y Z/}],
            attr => \%custom_attr
        );
    }
    foreach my $cp (grep {defined $_} @cps) {
        next unless $cp->{Name} =~ /BCP/i;
        $cp->{Econt} = $cp->{V}*0.5*627.5095;
        foreach my $len (grep {/GBL/i} keys %$cp) {
            $cp->{$len} *= BOHR_RADIUS;
        }
        # $_->{Rho} *= 6.7483;
        # $_->{L}   *= -4*24.1;
        my @bcp_atoms = split " ", $cp->{Atoms};
        my %custom_attr;
        my %exclude = map {$_ => 1} qw/Atoms/;
        foreach my $key (grep {!$exclude{$_}} keys %$cp) {
            $custom_attr{"bader/$key"} = $cp->{$key};
        }
        $mol->new_bond(
            atoms => [map {$mol->by_id($_)} @bcp_atoms],
            attr => \%custom_attr
            );
    }

    return $mol;
}

sub name_is {
    my ($class, $fname) = @_;
    $fname =~ /\.sum$/i;
}

sub file_is {
    my ($class, $fname) = @_;
    $fname =~ /\.sum$/i;
}

sub _read_CPs {
    my $fh = shift;
    my $id;
    my $type_regex = qr/^\s*Type\s+=\s+(\(3,[+-]\d\))\s+(BCP|RCP|CCP|NACP|NNACP)\s+(.*)$/i;
    my @datatypes = (qw/Rho V G K L GBL_[A-Z]+/, 'Bond Ellipticity');
    my $datatypes = join '|', @datatypes;
    my $data_regex = qr/^\s*($datatypes)\s+=\s+($scinum|NA)\b/i;
    my @cps;
    
    while (<$fh>) {
        if (/CP# (\d+)/) {
            $id = $1;
        }
        if (/$type_regex/) {
            @{$cps[$id]}{qw/Type Name Atoms/} = ($1, $2, $3);
        } 
        if (/$data_regex/) {
            $cps[$id]{$1} = $2;
        }
    }
    return @cps;
}


# methods I need to override:

# name_is      => name check
# parse_string =>  returns @ of mols - don't need if implement read_mol correctly







1;
