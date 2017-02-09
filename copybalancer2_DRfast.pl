#! /usr/bin/perl -w
use strict;

die "Usage: $0 coverageFile k toleranceRatio>=1 splitMPprefix\n" unless scalar @ARGV >= 4;
my ($covfile, $k, $tolratio, $splitpre) = @ARGV;
my $invratio = 1/$tolratio;
die "Can't find splitMP files $splitpre.L.fa and $splitpre.R.fa\n" unless -f "$splitpre.L.fa" and -f "$splitpre.R.fa";
my (%cts, %mers, $mp, %revcomps, %mps);

#my %propers; open IN, "../../NM1/bridges/NM1both.main.NM1both.main.mask3.bridges"; while (<IN>) { next unless /^(\S+).*PROPER/; $propers{$1} ++;} close IN;

open IN, "$covfile";
while (<IN>) {
  chomp;
  next unless /^(\S+)\t(\S+)/;
  $mers{$1} = $2;
  $revcomps{$1} = $1;
  $revcomps{Revcomp($1)} = $1;
}
close IN;

for my $part ('L', 'R') {
  open IN, "$splitpre.$part.fa";
  while (<IN>) {
    if (/^>(\S+)/) {$mp = $1; next}
    next if $part eq 'R' and $mps{$mp}{reject};
    my $prop = ''; #my $prop = 'IMPRO'; $prop = 'PROPR' if $propers{$mp};
    chomp;
    my ($n, $sum, @copies);
    for ($_ =~ /(?=(.{$k}))/g) {
      $n ++;
      next unless $revcomps{$_};
      push @copies, $mers{$revcomps{$_}};
      $sum += $mers{$revcomps{$_}};
    }
    unless ($sum) {$mps{$mp}{reject} ++; $cts{"rejectNoMatch$prop"} ++; next}
    my $m = scalar @copies; 
    unless ($m/$n >= 0.9) {$mps{$mp}{reject} ++; $cts{"rejectPoorMatch$prop"} ++; next}
    my $avg = $sum/$m;
    for (@copies) {
      my $ratio = $_/$avg;
      if ($ratio < $invratio or $ratio > $tolratio) {$mps{$mp}{reject} ++; $cts{"rejectInconsistent$prop"} ++; last}
    }
    next if $mps{$mp}{reject};
    $mps{$mp}{$part}{copy} = $avg;
    $mps{$mp}{$part}{seq} = $_;
  }
  close IN;
}

my @cuts = qw/1.05 1.11 1.18 1.29 1.43 1.66 2.10 3.32/;

open L, ">$splitpre.Lbalance.fa";
open R, ">$splitpre.Rbalance.fa";
open CT, ">$splitpre.ctbalance";
for my $mp (keys %mps) {
  my $prop = ''; #'IMPRO'; $prop = 'PROPR' if $propers{$mp};
  $cts{"mps$prop"} ++;
  next if $mps{$mp}{reject};
  $cts{"mpsMatch$prop"} ++;
  my $ratio = $mps{$mp}{L}{copy}/$mps{$mp}{R}{copy};
  my $test=$ratio; $test=1/$ratio if $ratio<1; my $bin = 1; for my $i (0 .. $#cuts) {$bin = $cuts[$i] if $test>$cuts[$i]} $cts{"$bin$prop"} ++;
  if ($ratio < $invratio or $ratio > $tolratio) {$mps{$mp}{reject} ++; $cts{"rejectImbalance$prop"} ++; next}
  $cts{"yield$prop"} ++;
  print L ">$mp ", sprintf("%.1f", $mps{$mp}{L}{copy}), "\n$mps{$mp}{L}{seq}\n";
  print R ">$mp ", sprintf("%.1f", $mps{$mp}{R}{copy}), "\n$mps{$mp}{R}{seq}\n";
  print CT "$mp\t", sprintf("%.1f", ($mps{$mp}{L}{copy}+$mps{$mp}{R}{copy})/2), "\n";
}
close L;
close R;
close CT;

print "covfile=$covfile; k=$k; tolratio=$tolratio; splitMPprefix=$splitpre\n";
for (sort keys %cts) {print "$cts{$_}\t$_\n"}


sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}
