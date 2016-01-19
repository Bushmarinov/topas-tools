#!/usr/bin/perl
# finds all uninstalled modeules for all scripts (*.pl) in a directory
# with -i, installs them

use Modern::Perl '2015';
use Module::Extract::Use;
use autodie;

my (%seen, @modules);
my $extor = Module::Extract::Use->new;

foreach (glob('*.pl')) {
    push @modules, grep {!$seen{$_}++} grep {defined $_} $extor->get_modules($_);
}
@modules = 
    grep { !eval "require $_"} 
    grep {/^[A-Z]/}       #pragmas start lowercase
    @modules;
say join "\n", @modules;
if ($ARGV[0] and $ARGV[0] eq '-i') {
    say system "cpan -i ".join " ", @modules;
}