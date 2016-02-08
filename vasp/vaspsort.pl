#/usr/bin/perl

use Modern::Perl '2015';
use Chemistry::Elements qw(get_Z get_symbol);
use Math::Trig;
use Getopt::Long;
use autodie;

my @lines; # AoA [line, el_num]
my %types; #holds atomic numbers
my @cell_matrix;
my $poscarf = "POSCAR";
my $outfile;
my ($force, $help);
GetOptions( "force" => \$force,
            "poscar=s" => \$poscarf,
            "resfile=s" => \$outfile,
            "help|?" => \$help);

(my $input = shift @ARGV) && !$help or die <<HELP;
VASPSORT.PL [options] <p1.res>

   Sorts a .res file according to VASP rules, removes AFIX.
   Writes a sorted .res and a complete POSCAR file 
   assumping P1 space group.
   To get such a nice file, HIMP the hydrogens and run the 
   following commands in OLEX2:
      ChangeSg P1
	  pack cell
	  file vasp.res
      
   Options:
     --force, -f     
        force overwrite of output files
     --poscar, -p    
        name of output POSCAR file (default POSCAR)
     --resfile, -r   
        name of output RES file (append _sorted by default)
     --help, -h, -?  
        show this message 
HELP

(my $title = $input) =~ s/\..*//;
#$title = uc $title;
$outfile //= $title.'_sorted.res';
open my $inputh, '<', $input or die "Cannot open file `$input': $!\n";
$force || !(-e $poscarf) or die "File $poscarf exists! Delete it or use -f to overwrite.\n";
$force || !(-e $outfile) or die "File $outfile exists! Delete it or use -f to overwrite.\n";
open my $outfileh, ">", $outfile;
open my $poscarh, ">", $poscarf;

while (<$inputh>) {
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
		@cell_matrix = @{transpose([ [$aa, $ba*cos($ga), $ca*cos($be) ],
						 [  0, $ba*sin($ga), $ca*(cos($al)-cos($be)*cos($ga))/sin($ga)],
						 [  0,            0, $ca*$v/sin($ga)] ])};
	}
	next if /AFIX/i;
	if (/^([A-Za-z]{1,2})(?![A-Za-z])[\w']* \s+ \d+ \s+ ((?:[-0-9.]+\s+){2}[-0-9.]+)/xi) {
		(my $el_symbol, $el_coord) = ($1, $2);
		$el_num = get_Z($el_symbol);	
		if ($el_num) {
			$types{$el_num}++;
		} else {
			$el_num = 999;
		}
        my $newcoord = join " ", map {sprintf '%.5f', $_} 
                                map {$_ > 9 ? $_-10 : $_} 
                                split " ", $el_coord;       # fix an OLEX2 error of assigning "constant" coods
        s/\Q$el_coord\E/$newcoord/;
        $el_coord = $newcoord;
	} elsif (/HKLF/i || /END/i) {
		$el_num = 0;
	} else {
		$el_num = 999
	}
	push @lines, [$_, $el_num, $el_coord];
}
my $cell_matrix = join '', map {('    ', (join " ",  map { sprintf '%.6f', $_} @$_), "\n")} @cell_matrix;

@lines = sort {$b->[1] <=> $a->[1]} @lines;

print $outfileh shrinkspaces(map {$_->[0]} @lines);

print $outfileh "\n-----------------------------------\n\n";
print $outfileh map {get_symbol($_)." ".$types{$_}."\n"} sort {$b <=> $a} keys %types;
print $outfileh "\n\n-----------------------------------\n\n";
print $outfileh $cell_matrix;
print $outfileh "\n\n------------POSCAR---------------\n\n";
chomp $cell_matrix;
foreach my $outh ($outfileh, $poscarh) {
    print $outh <<POSCAR;
$title CRYSTAL
   1.0000
$cell_matrix
POSCAR
    print $outh '   ', (join ' ', map {get_symbol($_)} sort {$b <=> $a} keys %types), "\n";
    print $outh '   ', (join ' ', map {$types{$_}} sort {$b <=> $a} keys %types), "\n";
    print $outh "Direct\n";
    print $outh map {" $_->[2]\n"} grep {$_->[2]} @lines;
}

say "RES file with additional info written to $outfile";
say "POSCAR file written to POSCAR";

sub transpose {
	my $matrix = shift;
	my $new_matrix;
	foreach my $i (0..$#$matrix) {
		foreach my $j (0..$#{$matrix->[$i]}) {
			$new_matrix->[$i][$j] = $matrix->[$j][$i]
		}
	}
	return $new_matrix;
}

sub shrinkspaces {				#collapse all empty strings in an array to one
	my $seen_space = 0;
	my @result;
	foreach my $line (@_) {
		my $is_space = !($line =~ /\S/);
		next if $is_space and $seen_space;
		$seen_space = $is_space;
		push @result, $line;
	}
	return @result;
}
