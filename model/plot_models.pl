#! /usr/bin/perl

use strict;
use warnings;

use PDL;
use PDL::Graphics::PGPLOT;
use PDL::Fit::Polynomial;

dev 'models.ps/vcps', 2, 3, { hardlw => 2, hardch => 1.2 };
for my $year (2010..2014) {
  for my $order (qw/neg pos/) {

    my $ratio_file = "ratios/${year}_${order}.ratio";
    my ($wav, $flux, $ratio, $err) = rcols($ratio_file);

    my $model_file = "${year}_${order}.model";
    my ($wavm, $ratiom) = rcols($model_file);


    points $wav, $ratio, { yrange => [0.85,1.1], title => "$year - $order", border => 1 };
    hold;
    errb $wav, $ratio, $err;
    line $wavm, $ratiom;

    my $o = 6;
    my (undef, $coeffs)  = fitpoly1d($wav, $ratio, $o, { Weights => 1/$err/$err });
    $wavm = sequence(101)/100 * ($wav->max-$wav->min) + $wav->min;
    $ratiom = sumover($coeffs * $wavm->dummy(0,$o)**sequence($o));
    line $wavm, $ratiom, { color => 2 };

    release;

  }
}
