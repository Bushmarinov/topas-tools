#!/usr/bin/perl

use Modern::Perl '2013';
use autodie;
use threads;
use Thread::Queue;
# use threads::shared;
use Getopt::Long;
use File::Path;
use Carp;
use Sys::Info;
use Fcntl;
use YAML;
use File::Basename;
use Chemistry::Elements qw/get_Z/;

sub ml ($$);
my $ncpu;
$ncpu = 1;
{
    my $safe = async {
        #Sys::Info unfortunately is VERY thread-unsafe: crashes Perl upon joining >1 thread
        my $info = Sys::Info->new;          # so I load it inside a separate thread and destroy with that thread
        my $cpu  = $info->device( 'CPU' );
        my $count = 0;
        $count += $_->{number_of_cores} foreach $cpu->identify;
        return $count;
    };
    $ncpu = $safe->join();
    $ncpu or die "No cpu count returned!\n";
}

if ($^O =~ /win32/i) {
    my $proc;
    require Win32::Process;
    Win32::Process::Open($proc, $$, 0);
    $proc->SetPriorityClass(Win32::Process::IDLE_PRIORITY_CLASS());
    # $format = 'C:/AIMALL/aimint.exe %s %s >nul';
}

# my $dir_pattern = '^\w+_\w+_[-+\d._]+$';
# my $file_pattern = '\.inp$';
my $help;
my $keep;
my $precision = 0.005;
my $delta = 0.2;
my $multiplier = 2;
my $nextgen;
GetOptions( "help|?" => \$help, 
            "keep" => \$keep,
            "precision=f" => \$precision,
            "delta=f" => \$delta,
            "number-of-cpus|n=i" => \$ncpu,
            "multiplier=f" => \$multiplier,
            "nextgen|X" => \$nextgen,
);

die <<HELP if $help or !@ARGV;
Usage: toparunner3.pl [options] <inpfile>

Performs a search for optimal values of all bond restraints
in inpfile using restraint stabilization.

Options:
   --help, -h, ?
     Displays this message
   --keep, -k
     do not delete toparun files adter running (space-consuming!)
   --precision
     expected error (A) on finding the bond stability range
     (default 0.005)
   --delta
     starting delta (A) for searching bond stability range
     assuming that it outlies up there (default 0.2).
   --number-of-cpus
     number of CPUs ("cores") to use. By default all available,
     but with IDLE priority, so no worry.
   --multiplier, -m
     factor to multiply toparun default steps "1,0.25".
     Default 2.0
   --nextgen, -X
     use "next generation" IQR multiplier
HELP

say "I have $ncpu cores";

my $startfile = shift @ARGV;
$startfile =~ /\.inp$/i or die "The starting model should be in an INP file!\n";

{
    my $nsteps = int log($delta/$precision)/log(2);
    printf "Will require %i steps per bond, with final error of %g\n", $nsteps, $delta*2**(-$nsteps-1);
}
my $bond_limits = Thread::Queue->new();
my $startname = (fileparse($startfile))[0];

my $restraints;
my $model_text;
{
    open my $inph, "<", $startfile;
    $model_text = join "", <$inph>;
    seek $inph, 0, 0;
    $restraints = read_restraints($inph);
}
my $q = Thread::Queue->new();
my %ideal = map {$_->{name} => $_->{ideal}} @$restraints;
# die Dump(\%ideal);
# my %bond_ordering;

sysopen my $resulth, "toparunning_result.txt",  O_RDWR | O_CREAT;
my %already;
while (<$resulth>) {
    chomp;
    my $bond = (split /\t+/, $_)[0];
    next unless $bond =~ /\w+\s+\w+/;
    say "Bond '$bond' already processed, will not repeat";
    $already{$bond}++;
}

$q->enqueue(map {$_->{bond}} 
            sort ml 
            map {
                m/^([A-Za-z]+)(\d+)\w* \s+ ([A-Za-z]+)(\d+)\w*\b.*$/x or die "Bond '$_' has unusual format!\n";
                {Z1 => get_Z($1), L1 => $2, Z2 => get_Z($3), L2 => $4, bond => $_};
            } 
            grep {!$already{$_}}
            keys %ideal);
$q->end;



my @workers;
foreach (1..$ncpu) {
    push @workers, async {
        my $thread_id = $_;
        while (my $bond = $q->dequeue) {
            my $pending = $q->pending || 0;
            printf "Thread %i searches for limits for bond %s, %i bond%s left\n", $thread_id, $bond, $pending, $pending == 1 ? '' : 's';
            (my $dir = "./$bond") =~ s/\s/_/g;
            -e $dir or mkdir $dir;
            my $func = sub {
               run_toparun($model_text, $bond, $_[0], $dir, $startname, $keep, $multiplier, $nextgen);
            };
            my $upper_limit = bisection([$ideal{$bond}, 1], [$ideal{$bond}+$delta, -1], $precision, $func);
            my $lower_limit = bisection([$ideal{$bond}, 1], [$ideal{$bond}-$delta, -1], $precision, $func);
            $bond_limits->enqueue([$bond, $lower_limit, $upper_limit]);
            printf "Thread %i done with bond %s\n", $thread_id, $bond;
        }
        # $bond_limits->end();
        return;
    }
    sleep 1;
}

$resulth->autoflush(1);

while (grep {$_->is_running} @workers or $bond_limits->pending) {
    my $line = $bond_limits->dequeue_timed(5);         #to consider small racing condition when the last thread finishes job
    next unless $line;
    say $resulth "Bond\tLower\tUpper\tMiddle\tIdeal\tDiff" unless keys %already;
    my ($bond, $lower, $upper) = @$line;
    my $middle = ($lower+$upper)/2;
    say $resulth join "\t", $bond,  $lower, $upper, $middle, $ideal{$bond}, $middle-$ideal{$bond};
    $already{$bond}++;
}
$bond_limits->end;
# say $resulth "Finished";

foreach (@workers) {
    my ($error) = $_->join();
    # say "Joining worker $_";
    if ($error) {
        print $error;
    }
}
# print $resulth Dump(\%bond_limits);

say "Toparunning finished";


sub read_restraints {        #generates a @AoH of {name, ideal, real} from filehandle; returns referecne to it}
    my $handle = shift;
    my @bondlist;
    while (<$handle>) {
        if (m{Distance_Restrain\w*\s*\(
                 \s* ([^,]+?)        #name
                 \s*,\s*
                 ([\d.]+) (?:`[^,]*)?   #ideal
                 \s*,\s*
                 ([\d.]+) (?:`_([-.\d]+))?   #real and sigma
            }xi) {
            my ($name, $ideal, $real, $sigma) = ($1, $2, $3, $4);
            $name =~ s/\s+/ /;
            $name =~ s/\s+$//;
            push @bondlist, {name => $name, ideal => $ideal, real => $real, sigma => $sigma || 0};
        }
    }
    return \@bondlist;
}

sub bisection {
    my (undef, undef, $precision, $func) = @_;
    my ($l, $r) = sort {$a->[0] <=> $b->[0]} map {[$_->[0], $_->[1] // $func->($_->[0])]} ($_[0], $_[1]);
    croak "No zero evident" if $l->[1]*$r->[1] > 0;
    my $step = ($r->[0] - $l->[0])/2;
    while ($step >= $precision) {
        my $next = $l->[0] + $step;
        my $result = $func->($next);
        if ($l->[1]*$result > 0) {
            $l = [$next, $result];
        } else {
            $r = [$next, $result];
        }
        $step /= 2;
    }
    return $l->[0] + $step;
}


sub run_toparun {
    my ($model, $bond, $value, $dir, $filename, $keep, $multiplier, $nextgen) = @_;
    # die Dump([$model, $bond, $value, $dir, $filename]);
    (my $bondr = $bond) =~ s/\s+/\\s*/g;
    my @steps = (1, 0.25);
    my $X = $nextgen ? '-X' : '';
    @steps = map {$_*$multiplier} @steps;
    my $limit = sprintf "%.0f", 16/$multiplier;
    my $steps = join ",", @steps;
    $bondr = qr/$bondr/i;
    $model =~ s{(Distance_Restrain\w*\s*\(
                 \s* ($bondr)        #name
                 \s*,\s*)
                 [\d.]+
               }
               {$1$value}x or die "Bond $bondr not found in model!\n";
    (my $subdir_name = $value) =~ s/\./_/g;
    {
        -e "$dir/$subdir_name" or mkdir "$dir/$subdir_name";
        open my $modelh, ">", "$dir/$subdir_name/$filename";
        print $modelh $model;
    }
    my $result = 1;
    {    
        open my $runnerh, "toparun.pl -f -c $X -l $limit -s $steps $dir/$subdir_name/$filename |";
        $runnerh->autoflush(1);
        while (<$runnerh>) {
            if (/Stopped/i and /limit/i) {
                $result = -1;
            }
            #print;
        }
    }
    system "spaghetti $dir/$subdir_name >nul";
    rmtree("$dir/$subdir_name") unless $keep;
    return $result;
}

sub ml ($$) {
    return $_[1]->{Z1} <=> $_[0]->{Z1}
                      ||
           $_[0]->{L1} <=> $_[1]->{L1}
                      ||
           $_[1]->{Z2} <=> $_[0]->{Z2}
                      ||
           $_[0]->{L2} <=> $_[1]->{L2}
}