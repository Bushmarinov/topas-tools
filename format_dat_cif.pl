#!/usr/bin/perl

use warnings;
use strict;
use autodie;
use File::Basename;

my $dat = shift @ARGV or die "Please specify a a dataset file in Out_Graphs format!\n";
(my $name = basename($dat)) =~ s/\..*//g;
(my $outname = $dat) =~ s/(\.[^.]*)?$/_dataset.cif/;
open my $outh, ">", $outname;
open my $dath, "<", $dat;

print $outh <<HEADER;
data_$name

_pd_block_id    '$name'

loop_
      _pd_proc_point_id                           
      _pd_meas_2theta_scan
      _pd_meas_counts_total                     
      _pd_calc_intensity_total 
HEADER
while (<$dath>){
    printf $outh qq(  %-5d %4.5f %8g %7.5f\n), $., (split)[0..2]
}

print $outh "\n";

print "Dataset in CIF format output to $outname\n";
