#! /usr/bin/perl

use strict;
use warnings;

use PDL;
use Astro::FITS::CFITSIO;
use Carp;

my ($wav, $models, $T_eff, $log_g) = get_models();

my $template = 'template.txt';
my $fits = 'wd2.mod';

my ($elo, $ehi);
($elo, $ehi, $models) = restructure_models($wav, $models);

write_template($template, $elo, $ehi, $models, $T_eff, $log_g);


my $status = 0;
my $fptr = Astro::FITS::CFITSIO::create_file("!$fits($template)", $status);
check_status($status) or die;

$fptr->movabs_hdu(1, Astro::FITS::CFITSIO::BINARY_TBL(), $status);
$fptr->write_date($status);
$fptr->write_chksum($status);

write_parameters($fptr, $status, $T_eff, $log_g);
write_energies($fptr, $status, $elo, $ehi);
write_spectra($fptr, $status, $elo, $ehi, $models, $T_eff, $log_g);

$fptr->close_file($status);
check_status($status) or die;

sub restructure_models {
  my $wav = shift;
  my $models = $_[0]->copy;

  my $wav_bins = zeroes($wav->nelem + 1);
  (my $tmp = $wav_bins->slice('1:-2')) .= 0.5 * ($wav->slice('1:-1') + $wav->slice('0:-2'));
  $wav_bins->set(0, $wav->at(0) - 0.5 * ($wav->at(1)-$wav->at(0)));
  $wav_bins->set(-1, $wav->at(-1) + 0.5 * ($wav->at(-1)-$wav->at(-2)));

  my $wav_lo = $wav_bins->slice('0:-2')->copy;
  my $wav_hi = $wav_bins->slice('1:-1')->copy;

  my $e_bins = 12.3985 / $wav_bins->slice('-1:0:-1');
  my $e_lo = $e_bins->slice('0:-2')->copy;
  my $e_hi = $e_bins->slice('1:-1')->copy;

  my $ergs_per_kev = 1.6021773e-9;
  my $e = 0.5 * ($e_lo + $e_hi);

  for my $i (0..$models->getdim(1)-1) {
    my $s = ",($i)";
    (my $tmp = $models->slice($s)) .= ($models->slice($s) * ($wav_hi-$wav_lo))->slice('-1:0:-1') * $ergs_per_kev / $e; # photons / cm^2 / s in each bin
  }

  return $e_lo, $e_hi, $models;
}

sub write_template {
  my ($file, $elo, $ehi, $models, $T_eff, $log_g) = @_;

  my $T_eff_samp = $T_eff->uniq;
  my $log_g_samp = $log_g->uniq;

  my $numbvals = $T_eff_samp->nelem;
  $numbvals = $log_g_samp->nelem if $log_g_samp->nelem > $numbvals;

  open my $fh, '>', $file or die "could not open '$file' for writing: $!";
  print $fh <<EOP;
simple = T
bitpix = 16
naxis = 0
modlname = wd
modlunit = photons/cm^2/s
redshift = F
addmodel = T
hduclass = OGIP
hduclas1 = 'XSPEC TABLE MODEL'
hduvers = 1.0.0

xtension = bintable
extname = PARAMETERS
naxis2 = 2
nintparm = 2
naddparm = 0
ttype# = NAME / name of the parameter
tform# = 12a
ttype# = METHOD / interpolation method (0 if linear, 1 if logarithmic)
tform# = j
ttype# = INITIAL / initial value in the fit for the parameter
tform# = e
ttype# = DELTA / parameter delta used in fit (if negative paramter is frozen)
tform# = e
ttype# = MINIMUM / hard lower limit for parameter value
tform# = e
ttype# = BOTTOM / soft lower limit for parameter value
tform# = e
ttype# = TOP / soft upper limit for parameter value
tform# = e
ttype# = MAXIMUM / hard upper limit for parameter value
tform# = e
ttype# = NUMBVALS / the number of tabulated parameter values
tform# = j
ttype# = VALUE / the tabulated parameter values (n is NUMBVALS)
tform# = ${numbvals}e
hduclass = OGIP
hduclas1 = 'XSPEC TABLE MODEL'
hduclas2 = PARAMETERS
hduvers = 1.0.0

xtension = bintable
extname = ENERGIES
naxis2 = @{[$elo->nelem]}
nintparm = 2
naddparm = 0
ttype# = ENERG_LO
tform# = e
ttype# = ENERG_HI
tform# = e
hduclass = OGIP
hduclas1 = 'XSPEC TABLE MODEL'
hduclas2 = ENERGIES
hduvers = 1.0.0

xtension = bintable
extname = SPECTRA
naxis2 = @{[$T_eff_samp->nelem * $log_g_samp->nelem]}
ttype# = PARAMVAL
tform# = 2e
ttype# = INTPSPEC
tform# = @{[$elo->nelem]}e
tunit# = photons/cm^2/s
hduclass = OGIP
hduclas1 = 'XSPEC TABLE MODEL'
hduclas2 = 'MODEL SPECTRA'
hduvers = 1.0.0

EOP

  $fh->close();
}

sub write_parameters {
  my $fptr = shift;
  my ($T_eff, $log_g) = @_[1..2];

  my $T_eff_samp = $T_eff->uniq;
  my $log_g_samp = $log_g->uniq;

  my $nvals = $T_eff_samp->nelem;
  $nvals = $log_g_samp->nelem if $log_g_samp->nelem > $nvals;

  # what we're actually writing to the table
  my $value = zeroes($nvals, 2);
  my $tmp;
  ($tmp = $value->slice('0:'.($T_eff_samp->nelem-1).',(0)')) .= $T_eff_samp;
  ($tmp = $value->slice('0:'.($log_g_samp->nelem-1).',(1)')) .= $log_g_samp;

  my $name = ['T_eff', 'log_g'];
  my $method = [1, 0]; # 0 = linear, 1 = logarithmic
  my $initial = [$T_eff_samp->median, $log_g_samp->median];
  my $delta = [10, 0.01]; # ???
  my $minimum = [$T_eff_samp->min, $log_g_samp->min];
  my $bottom = $minimum;
  my $top = [$T_eff_samp->max, $log_g_samp->max];
  my $maximum = $top;
  my $numbvals = [$T_eff_samp->nelem, $log_g_samp->nelem];

  $fptr->movabs_hdu(2, Astro::FITS::CFITSIO::BINARY_TBL(), $_[0]);
  $fptr->write_col_str(1, 1, 1, 2, $name, $_[0]);
  $fptr->write_col_lng(2, 1, 1, 2, $method, $_[0]);
  $fptr->write_col_flt(3, 1, 1, 2, $initial, $_[0]);
  $fptr->write_col_flt(4, 1, 1, 2, $delta, $_[0]);
  $fptr->write_col_flt(5, 1, 1, 2, $minimum, $_[0]);
  $fptr->write_col_flt(6, 1, 1, 2, $bottom, $_[0]);
  $fptr->write_col_flt(7, 1, 1, 2, $top, $_[0]);
  $fptr->write_col_flt(8, 1, 1, 2, $maximum, $_[0]);
  $fptr->write_col_lng(9, 1, 1, 2, $numbvals, $_[0]);
  $fptr->write_col_flt(10, 1, 1, $value->nelem, $value->float->get_dataref, $_[0]);
  $fptr->write_chksum($_[0]);
}

sub write_energies {
  my $fptr = shift;
  my ($elo, $ehi) = @_[1..2];
  $fptr->movabs_hdu(3, Astro::FITS::CFITSIO::BINARY_TBL(), $_[0]);
  $fptr->write_col_flt(1, 1, 1, $elo->nelem, $elo->float->get_dataref, $_[0]);
  $fptr->write_col_flt(2, 1, 1, $ehi->nelem, $ehi->float->get_dataref, $_[0]);
  $fptr->write_chksum($_[0]);
}

sub write_spectra {
  my $fptr = shift;
  my ($elo, $ehi, $models, $T_eff, $log_g) = @_[1..5];

  my $T_eff_samp = $T_eff->uniq;
  my $log_g_samp = $log_g->uniq;

  $fptr->movabs_hdu(4, Astro::FITS::CFITSIO::BINARY_TBL(), $_[0]);

  for my $i (0..$T_eff_samp->nelem-1) {
    for my $j (0..$log_g_samp->nelem-1) {
      my $index = which( ($T_eff==$T_eff_samp->at($i)) &
			 ($log_g==$log_g_samp->at($j))
		       );
      $index->nelem == 1 or die;
      my $k = $index->at(0);
      $fptr->write_col_flt(1, $k+1, 1, 2, [$T_eff->at($k), $log_g->at($k)], $_[0]);
      $fptr->write_col_flt(2, $k+1, 1, $models->getdim(0), $models->slice(",($k)")->float->get_dataref, $_[0]);
    }
  }
  $fptr->write_chksum($_[0]);
}

sub get_models {

  my @files = glob('non_LTE/??.model');
  # do not include Nick's T=51000, log_g=7.9 model (15.model)`
  @files = @files[0..14];

  my $T_eff = zeroes(scalar @files);
  my $log_g = $T_eff->copy;

  my ($wav, $models);

  for my $i (0..$#files) {
    my $file = $files[$i];
    open my $fh, '<', $file or die "cannot open '$file': $!";
    while (defined( my $line = <$fh>)) {
      $line =~ /^#/ or last;
      $line =~ /T_eff=(\S+)/ and $T_eff->set($i, $1);
      $line =~ /log_g=(\S+)/ and $log_g->set($i, $1);
    }
    close $fh;

    my ($lambda, $flux) = rcols $file;

    if (! defined $models) {
      $models = zeroes( $flux->nelem, $T_eff->nelem );
      $wav = $lambda->copy;
    }

    which($wav != $lambda)->nelem and die $i;
    (my $tmp = $models->slice(",($i)")) .= $flux;
  }

  return $wav, $models, $T_eff, $log_g;
}

sub check_status {
  my $s = shift;
  if ($s != 0) {
    my $txt;
    Astro::FITS::CFITSIO::fits_get_errstatus($s,$txt);
    carp "CFITSIO error: $txt";
    return 0;
  }

  return 1;
}

