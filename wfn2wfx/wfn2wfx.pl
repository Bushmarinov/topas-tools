#!/usr/bin/perl

use Modern::Perl '2015';
use Text::Xslate;
use YAML::Any;
use List::Util qw /sum/;
use Getopt::Long;
# use File::Handle;
use autodie;



my ($title, $orbitals_count, $primitives_count, $nuclei_count, $energy, $virial, 
    @MOs, @centers, @types, @exponents, @atoms, @EDF_centers, @EDF_coefs, @EDF_exponents);
my ($outf, $help);
my ($pp_line, $r_line, %PPs, %radii);
my $multip;


GetOptions(
    "outfile=s" => \$outf,
    "pseudopotentials=s" => \$pp_line,
    "multiplicity=i" => \$multip,
    "radii=s" => \$r_line,
    "help" => \$help,
);

my $wfnf = (shift @ARGV);

die <<HELP if $help || !(defined($wfnf) && ($wfnf =~/\.wfn$/));
Syntax: wfn2wfx.pl [options] <*.wfn>

Options:
  --pseudopotentials, -p "AtomType1 ElNumber1 AtomType2 ElNumber2 ..."
      Specify pairs of atom types and electron numbers to insert a 
      pseudopotential section into the resulting wfx. Atom types are 
      case-insensitive.
  --radii, -r "AtomType1 RadiusMultiplier1 ..."
      Specify pairs of atom types and pseudopotential radius multiplier 
      for the resulting wfx. Usually not needed, good value for Cl in
      the SR2 basis set is 2.
  --multiplicity, -m <number>
      Provide a multiplicity (recommended unless it is 1 for 
      closed-shell or 2 for open-shell)
  --outfile, -o <name>
      Custom output file name.
  --help, -h, -?
      Display this message.
HELP

open my $wfnh, "<", $wfnf;

unless (defined $outf) {
    ($outf = $wfnf) =~ s/.[^.]*$/_converted.wfx/;
}

if ($pp_line) {
    for ($pp_line) {
        s/[,= ]+/ /g;
        my @items = split;
        die "Odd number of elements in pseudopotential line!\n" if @items % 2;
        %PPs = map {lc $_} @items;
    }
}

if ($r_line) {
    for ($r_line) {
        s/[,= ]+/ /g;
        my @items = split;
        die "Odd number of elements in pseudopotential line!\n" if @items % 2;
        %radii = map {lc $_} @items;
    }
}


# die Dump(\%PPs);

my $readmo;

while (<$wfnh>){
    chomp;
    if ($. == 1) {
        $title = $_;
        next;
    }
    #GAUSSIAN              5 MOL ORBITALS    126 PRIMITIVES        3 NUCLEI
    if (/GAUSSIAN \s* (\d+) \s* MOL\s*ORBITALS \s* (\d+) \s* PRIMITIVES \s*(\d+) \s* NUCLEI/x) {
        ($orbitals_count, $primitives_count, $nuclei_count) = ($1, $2, $3);
        next;
    }
    #CENTRE ASSIGNMENTS    1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
    if (/CENTRE \s* ASSIGNMENTS ((?:\s+\d+)+)/x) {
        push @centers, split ' ', $1;
        next;
    }
    #TYPE ASSIGNMENTS      1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  2  2  2  3
    if (/TYPE \s* ASSIGNMENTS ((?:\s+\d+)+)/x) {
        push @types, split ' ', $1;
        next;
    }
    #EXPONENTS  0.1533000D+05 0.2299000D+04 0.5224000D+03 0.1473000D+03 0.4755000D+02
    if (/EXPONENTS ((?:\s+[\d.DdEe+-]+)+)/x) {
        my $expline = $1;
        $expline =~ s/D/e/gi;
        # say $expline;
        push @exponents, split ' ', $expline;
        next;
    }
    #MO    1     MO 0.0        OCC NO =    2.0000000  ORB. ENERGY =    0.000000
    if (/MO \s* (?<number>\d+) \s* (?:MO \s* (?<mo>[\d.+-]+))? \s* OCC\s*NO\s*=\s* (?<occupancy>[\d.+-]+) \s* ORB.\s*ENERGY\s*=\s* (?<energy>[\d.+-]+)/x) {
        push @MOs, {%+};
        $readmo = 1;
        next;
    }
    # -0.24775821D-12 -0.46178057D-12 -0.78301385D-12 -0.11850082D-11 -0.14782043D-11
    if (/^  ((?:\s*[\d.DdEe+-]+)+) \s* $/x and $readmo) {
        my $mo_line = $1;
        $mo_line =~ s/D/e/gi;
        $mo_line =~ s/(\d)-/$1e-/g;
        push @{$MOs[-1]{coefs}}, split ' ', $mo_line;
        next;
    }
    #  O    1    (CENTRE  1)   0.00000000  0.00000000  0.22527161  CHARGE =  8.0
    if (/\(CENTRE/) {
        if (/^\s* (?<symbol>[A-Za-z]+) \s* (?<id>\d+) \s* \(CENTRE \s* (?<id>\d+)\) (?<coords>(?:\s+[\d.+-]+)+) \s* CHARGE \s* = \s* (?<charge>[\d.]+)/x) {
            my %atom = %+;
            $atom{coords} = [split ' ', $atom{coords}];
            push @atoms, \%atom;
            next;
        } else {
            die "Cannot parse atom string:\n$_\n";
        }
    }
    # TOTAL ENERGY =    -76.429980506591 THE VIRIAL(-V/T)=   2.00504062
    if (m{\w+ \s* ENERGY \s* = \s* ([\d.+-]+) \s* (?:THE)? \s* \QVIRIAL(-V/T)\E \s* = \s*  ([\d.+-]+)}x) {
        ($energy, $virial) = ($1, $2);
        next;
    }
    if (/^\s*END \s* DATA \s*$/x) {
        next;
    }
    next unless /\S/;
    die "Unparsed string $. in file $wfnf:\n$_\n";
    
}
die "Nuclei count mismatch\n" if $nuclei_count != @atoms;
die "MO count mismatch\n" if $orbitals_count != @MOs;
die "Primitive count mismatch\n" if grep {$primitives_count != scalar @$_} \@centers, \@types, \@exponents, map {$_->{coefs}} @MOs;

foreach my $atom (@atoms) {
    if ($PPs{lc $atom->{symbol}}) {
        my $r = $radii{lc $atom->{symbol}} // 1;
        push @EDF_centers, $atom->{id};
        push @EDF_coefs, $PPs{lc $atom->{symbol}}*8/($r*sqrt($r));
        push @EDF_exponents, 12.56637061436/$r;
    }
}
my $electrons = sum(map {$_->{occupancy}} @MOs);
$multip //= $electrons % 2 ? 2 : 1;
my $alpha = int($electrons/2) + ($multip-1),
my $beta = int($electrons/2)- ($multip-1) + $electrons % 2,
my $charge = -sum(map {$_->{occupancy}} @MOs) + sum(map {$_->{charge}} @atoms) - sum(map {$PPs{lc $_->{symbol}} || 0} @atoms);
say $charge;
# OUTPUT
my $tmpl_path = './';
my $tx = Text::Xslate->new(
    function => {
        wrap => \&wrapnum,
        formexp => \&formexp,
        formint => \&formint,
        count => sub {
            my $ary = shift;
            return scalar @$ary;
        }
    },
    type => 'text',
    path => $tmpl_path,
    verbose => 2,
    # cache => 0,
);


open my $outh, ">:raw", $outf;

print $outh $tx->render('wfx.tx', {
    atoms => \@atoms,
    MOs => \@MOs,
    exponents => \@exponents,
    types => \@types,
    centers => \@centers,
    charge => $charge,
    electrons => $electrons,
    alpha => $alpha,
    beta => $beta,
    energy => $energy,
    virial => $virial,
    title => $title,
    EDF_centers => \@EDF_centers,
    EDF_coefs => \@EDF_coefs,
    EDF_exponents => \@EDF_exponents,
});

say "Converted wfx written to $outf";

sub formexp {
    if (@_ > 1) {
        return join " ", map {formexp($_)} @_;
    }
    elsif (@_ == 1) {
        if (ref $_[0] eq 'ARRAY') {
            return join " ", map {formexp($_)} @{$_[0]};
        } else {
            return sprintf('% .14E', $_[0])
        }
    } else {
        return;
    }    
}

sub formint {
    if (@_ > 1) {
        return join " ", map {formint($_)} @_;
    }
    elsif (@_ == 1) {
        if (ref $_[0] eq 'ARRAY') {
            return join " ", map {formint($_)} @{$_[0]};
        } else {
            return sprintf('% 6i', $_[0])
        }
    } else {
        return;
    }    
}

sub wrapnum {
    my $numline = shift;
    $numline =~ s/(.{70,100}\d)(?=\s)/$1\n/g;
    return " ".$numline;
}