#!/usr/bin/perl

use warnings;
use strict;

#print out_prec(1.54927, 0.0256);
my $help = <<HELP;
Usage: fix_atom_prec <TOPAS out> <cif file>
HELP

my $inp = shift @ARGV or die $help;
my $cif = shift @ARGV or die $help;
open my $inph, "<", $inp or die "Cannot open file $inp: $!\n";

#read "good" atoms
my @atoms;
my %atoms;
my $cur_rigid_pivot;
while (<$inph>) {
    next if /^\s*'/;
	if (m{site \s+ (\w+) \s+  
			x  [@\s]+ (-?[.\d]+)(?:`_([.\d]+))? \s+
			y  [@\s]+ (-?[.\d]+)(?:`_([.\d]+))? \s+
			z  [@\s]+ (-?[.\d]+)(?:`_([.\d]+))? \s+
	    }ix) 
	{
		my $atom;
		$atom->{name} = $1;
		$atom->{coords} = [$2, $4, $6];
		$atom->{errors} = [$3, $5, $7];
		push @atoms, $atom;
		$atoms{$atom->{name}} = $atom;
	}
	if (/rigid/i../Translate/i) {
		if (/{/) {
			$cur_rigid_pivot = '';
		}
		if (/{/../}/) {
			if (/{/) {
				$cur_rigid_pivot = '';
				next;
			}
			if (/^\s+([^\sH]\w+)/i) {
				$cur_rigid_pivot ||= $1;
			}
		}
		if (m{Translate \s* \(   
			[@\s]* (?:-?[.\d]+)`_([.\d]+) \s*,\s*
			[@\s]* (?:-?[.\d]+)`_([.\d]+) \s*,\s*
			[@\s]* (?:-?[.\d]+)`_([.\d]+) 
	    }ix)
		{
			#die "Caught!";
			my @errors = ($1, $2, $3);
			$atoms{$cur_rigid_pivot}->{errors} = \@errors;
			#print "set errors for $cur_rigid_pivot!\n";
		}
	}
	
}

 foreach my $atom (@atoms) {
	 # print out_atom_prec($_), "\n";
	 $atom->{outcoords} = [map {out_prec($atom->{coords}[$_], $atom->{errors}[$_])} (0..2)];
 }

# my $cif = shift @ARGV or die "Please specify a CIF file\n";
(my $fixed = $cif) =~ s/(\.[^\.]+)$/_fixed$1/;

open my $cifh, "<", $cif;
open my $fixedh, ">", $fixed;

my @name_and_coords = qw/_atom_site_label  _atom_site_fract_x _atom_site_fract_y _atom_site_fract_z/;



my %header;
my $read_header;
my $fix_coords;
my $header_counter;
while (<$cifh>) {
    unless (/\S/) {
        print $fixedh $_;
        next;
    }
    if (/^\s*#/) {
        print $fixedh $_;
        next;
    }
    if (/^\s*loop_/) {
        $read_header = 1;
        $header_counter = 0;
        %header = ();
        print $fixedh $_;
        next;
    }
    if ($read_header and !/^\s*_/) {
        $read_header = 0;
        #print keys %header, "\n";
        my @good = grep {defined($header{$_})} @name_and_coords;
        if (@good == @name_and_coords) {            
            $fix_coords = 1;
        }
    }
    if ($fix_coords and (/^\s*_/ or /^\s*loop_/)) {
        $fix_coords = 0;
    }
    if ($read_header) {
        if (/^\s*(_\w+)/) {
            $header{$1} = $header_counter++;
        }
    }
    if ($fix_coords) {
        chomp;
        my $copy = $_;
        my @quoted;
        $copy =~ s/("[^"]*")/push @quoted, $1; '@@@'.$#quoted.'@@@'/ge;
        $copy =~ s/('[^']*')/push @quoted, $1; '@@@'.$#quoted.'@@@'/ge;
        
        my @line = split " ", $copy;
        foreach my $elem (@line) {
            s/@@@(\d+)@@@/$quoted[$1]/ge;
            s/@@@(\d+)@@@/$quoted[$1]/ge;
        }
		my $name = $line[$header{_atom_site_label}];
        $line[$header{_atom_site_fract_x}] = $atoms{$name}{outcoords}[0];
		$line[$header{_atom_site_fract_y}] = $atoms{$name}{outcoords}[1];
		$line[$header{_atom_site_fract_z}] = $atoms{$name}{outcoords}[2];
        print $fixedh join " ", @line;
        print $fixedh "\n";
        next;
    }
    print $fixedh $_;
}

print "Corrected cif file output to $fixed\n";

sub out_atom_prec {
	my $atom = shift;
	my $coord_line = join '  ', map {sprintf '%-14s', $_} map {out_prec($atom->{coords}[$_], $atom->{errors}[$_])} (0..2);
	return sprintf '%-6s  %s', $atom->{name}, $coord_line;
}

sub out_prec {
	my ($num, $error) = @_;
	unless ($error) {
		return $num;
	}
	$error < 1 or die "Cannot display numbers with error > 1" ;
	
	my $out_prec = sprintf "%.0f", (-log10($error) + 1);
	my $error_str = sprintf "%.${out_prec}f", $error;
	my ($out_error) = ($error_str =~ /([1-9]\d*)$/);
	my $out_num = sprintf "%.${out_prec}f", $num;
	return "$out_num($out_error)";
}

sub log10 { log($_[0])/log(10)}

__END__
 site C1    x    0.50323 y   -0.14556 z    0.66880 occ C 0.250 beq =beqc;
        site C2    x    0.66076 y   -0.30397 z    0.70107 occ C 0.250 beq =beqc;
        site C3    x    0.87087 y   -0.24772 z    0.69634 occ C 0.250 beq =beqc;
        site C4    x    0.92101 y   -0.09412 z    0.65556 occ C 0.250 beq =beqc;
        site C5    x @  0.76556`_0.00066 y @  0.15765`_0.00223 z @  0.62075`_0.00020 occ C 0.250 beq =beqc;
        site C9    x @  0.55529`_0.00067 y @  0.14253`_0.00225 z @  0.62713`_0.00022 occ C 0.250 beq =beqc;
        site C6    x @  0.65999`_0.00064 y @  0.33823`_0.00203 z @  0.53546`_0.00022 occ C 0.250 beq =beqc;
        site C7    x @  0.45102`_0.00069 y @  0.24821`_0.00220 z @  0.53743`_0.00017 occ C 0.250 beq =beqc;
        site C8    x @  0.39327`_0.00053 y @  0.23597`_0.00286 z @  0.58357`_0.00015 occ C 0.250 beq =beqc;
        site C10   x    0.30048 y    0.34558 z    0.49701 occ C 0.250 beq =beqc;
        site N1    x    0.80918 y    0.23752 z    0.57484 occ N 0.250 beq =beqc;
        site O1    x @  0.21521`_0.00199 y @  0.25979`_0.00695 z @  0.58688`_0.00028 occ O 0.250 beq =beqc;
        site H1    x    0.37590 y   -0.24758 z    0.66649 occ H 0.250 beq =beqh;
        site H2    x    0.63420 y   -0.34144 z    0.73274 occ H 0.250 beq =beqh;
        site H3    x    0.96290 y   -0.39114 z    0.71631 occ H 0.250 beq =beqh;
        site H4    x    1.05542 y   -0.10813 z    0.65109 occ H 0.250 beq =beqh;
        site H1a   x    0.93331 y    0.22477 z    0.57057 occ H 0.250 beq =beqh;
        site H10   x    0.17632 y    0.23100 z    0.48993 occ H 0.250 beq =beqh;
		
		        'Rigid bodies
        rigid
            load point_for_site ux uy uz {
                 N13 14.0932  3.5691  0.0882
                 H13 13.7053  3.3621 -0.6509
            }
        Rotate_about_axies(0, 0, 0)
        Translate( @  0.00000,  @  0.00000,  @ 0.00000)