#! /usr/bin/perl

use strict;
use warnings;

use Chandra::Tools::Common qw/ read_bintbl_cols /;
use Astro::FITS::CFITSIO;
use Data::Dumper;
use PDL;

@ARGV == 6 or die "Usage: $0 file row src1 src2 bg1 bg2\n";

my ($f, $row, $s1, $s2, $b1, $b2) = @ARGV;

my ($c, $bg1, $bg2, $bin_lo, $bin_hi) = read_bintbl_cols($f, qw/counts background_up background_down bin_lo bin_hi/, { extname => "spectrum" });
$_ = $_->slice(",($row)") for $c, $bg1, $bg2, $bin_lo, $bin_hi;

my $ii = which(($bin_lo>$b1) & ($bin_hi<$b2));
my $bg_sum = $bg1->index($ii)->sum + $bg2->index($ii)->sum;

$ii = which(($bin_lo>$s1) & ($bin_hi<$s2));
my $src_sum = $c->index($ii)->sum;

my $hdr = Astro::FITS::CFITSIO::fits_read_header($f."[spectrum]");

my $bg = $bg_sum/($hdr->{BACKSCDN} + $hdr->{BACKSCUP}) / ($b2-$b1) * ($s2-$s1); 
my $src = $src_sum;
my $net = $src - $bg;

print $src, " ", $bg, " ", $net, "\n";
