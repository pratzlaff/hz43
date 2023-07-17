#! /usr/bin/perl

use warnings;
use strict;

my $version = '0.1';

use FindBin;
use Config;
use Carp;
use Chandra::Tools::Common qw/ read_bintbl_cols /;
use PDL;
use PDL::Graphics::PGPLOT;

use Getopt::Long;
my %default_opts = (
		    device => '/xs',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'device=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

@ARGV == 7 or die "Usage: $0 file chipid chipx chipy wmin wmax n\n\tTry --help.";

my ($file, $chipid, $chipx, $chipy, $wmin, $wmax, $n) = @ARGV;
my $wav = sequence($n)*($wmax-$wmin)/($n-1)+$wmin;
my ($e, $qeu) = qeu_at($file."[$chipid]", $chipx, $chipy);
wcols "%.2f\t%.6g", 12.39854/$e, $qeu;
exit;

=begin comment

@ARGV == 8 or die "Usage: $0 f1 f2 chipid chipx chipy wmin wmax n\n\tTry --help.";
my ($f1, $f2, $chipid, $chipx, $chipy, $wmin, $wmax, $n) = @ARGV;
plot_ratios($f1, $f2, $chipid, $wav);
exit;

my ($e1, $q1) = qeu_at($f1."[$chipid]", $chipx, $chipy);
my ($e2, $q2) = qeu_at($f2."[$chipid]", $chipx, $chipy);

my $ratio = interpol($wav, 12.39854/$e2, $q2) / interpol($wav, 12.39854/$e1, $q1);
print $ratio;
line $wav, $ratio, { border => .1 };

=cut

sub plot_ratios {
  my ($f1, $f2, $chipid, $wav) = @_;
  my ($jnk, $e1, $q1, $e2, $q2);
  ($jnk, $jnk, $jnk, $jnk, $e1, $q1) = read_qeu($f1."[$chipid]");
  ($jnk, $jnk, $jnk, $jnk, $e2, $q2) = read_qeu($f2."[$chipid]");

  dev $opts{device}, 3, 3;

  for my $i (0..@{$q1}-1) {
    for my $xi (0..+($q1->[$i]->dims)[1]-1) {
      for my $yi (0..+($q1->[$i]->dims)[2]-1) {
	my $q1 = $q1->[$i]->slice(",($xi),($yi)");
	my $e1 = $e1->[$i];
	my $q2 = $q2->[$i]->slice(",($xi),($yi)");
	my $e2 = $e2->[$i];
	my $ratio = interpol($wav, 12.39854/$e2, $q2) / interpol($wav, 12.39854/$e1, $q1);
	#line $wav, $ratio, { border => .1, title => "$xi,$yi" };
	$ratio = $q2/$q1;
	print $ratio->min,"\t",$ratio->max,"\n";
	line 12.39854/$e1, $q2/$q1, { border => .1, title => "$xi,$yi" };
      }
    }
  }
  exit;
}

sub qeu_at {
  my ($f, $chipx, $chipy) = @_;

  my ($x0, $x1, $y0, $y1, $energy, $qeu) = read_qeu($f);

  for my $i (0..@{$x0}-1) {
    my $xi = which(($x0->[$i]<=$chipx) & ($x1->[$i]>$chipx));
    my $yi = which(($y0->[$i]<=$chipy) & ($y1->[$i]>$chipy));
    if ($xi->nelem and $yi->nelem) {
      my $xi = $xi->at(0);
      my $yi = $yi->at(0);
      #print "$xi, $yi\n";
      my $qeu = $qeu->[$i]->slice(",($xi),($yi)");
      my $energy = $energy->[$i];
      return $energy, $qeu;
    }
  }
}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}

sub read_qeu {
  my $fname = shift;

  my ($rid, $energy, $qeu, $tdim3,
      $_2crpx3, $_2cdlt3, $_2crvl3,
      $_3crpx3, $_3cdlt3, $_3crvl3)
    = read_bintbl_cols($fname,
		       qw( regionid energy qeu tdim3
			   2crpx3 2cdlt3 2crvl3
			   3crpx3 3cdlt3 3crvl3
			));

  # function outputs
  my (@x0, @x1, @y0, @y1, @energy, @qeu);

  for my $i (0..$rid->nelem-1) {
    my $tdim3 = $tdim3->[$i];
    my @dims = $tdim3 =~ /\((\d+),(\d+),(\d+)\)/ or die $tdim3;

    my $qeu = $qeu->slice(",($i)")->copy;
    $qeu->reshape(@dims);

    my $energy = $energy->slice(",($i)")->copy;
    my $wav = 12.39854 / $energy;

    # in WCS, the first pixel spans 0.5-1.5, with the center being 1
    die unless ($_2crpx3->at($i) == 1 and $_3crpx3->at($i) == 1);

    my $chipx = sequence($dims[1]) * $_2cdlt3->at($i) + $_2crvl3->at($i);
    my $chipy = sequence($dims[2]) * $_3cdlt3->at($i) + $_3crvl3->at($i);

    my $x0 = $chipx - $_2cdlt3->at($i)/2;
    my $x1 = $x0 + $_2cdlt3->at($i);

    my $y0 = $chipy - $_3cdlt3->at($i)/2;
    my $y1 = $y0 + $_3cdlt3->at($i);

    push @x0, $x0;
    push @x1, $x1;
    push @y0, $y0;
    push @y1, $y1;
    push @energy, $energy;
    push @qeu, $qeu;

  }

  return \(@x0, @x1, @y0, @y1, @energy, @qeu);

}

=head1 NAME

template - A template for Perl programs.

=head1 SYNOPSIS

cp template newprog

=head1 DESCRIPTION

blah blah blah

=head1 OPTIONS

=over 4

=item --help

Show help and exit.

=item --version

Show version and exit.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> May 2012

=head1 SEE ALSO

perl(1).

=cut

