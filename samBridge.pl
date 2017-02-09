#! /usr/bin/perl -w
use strict;

unless (@ARGV >= 1) {die "Usage: $0 insam ctbalance-file [maskfile]\n"}
my ($insam, $ctfile, $maskfile) = @ARGV;

# SAMFLAGS
# BIT DEC MEANING
#  1    1 read paired
#  2    2 read mapped in proper pair
#  3    4 read unmapped
#  4    8 mate unmapped
#  5   16 read reverse strand
#  6   32 mate reverse strand
#  7   64 first in pair
#  8  128 second in pair
#  9  256 not primary alignment
# 10  512 read fails platform/vendor quality checks
# 11 1024 read is PCR or optical duplicate

my ($min, $max) = (501, 15000); # Allowable distance range
my %saminterp = (1 => 'paired', 2 => 'proper', 3 => 'unmapped', 4 => 'Munmapped', 5 => 'rev', 6 => 'Mrev', 7 => 'read1', 8 => 'read2', 9 => 'rehit', 10 => 'fail', 11 => 'dupe');
my (%pairs, %dnas, %cts, %bridges);

if ($maskfile) {
 open IN, $maskfile;
 while (<IN>) {
  next unless /^\S+\s+\S+\s+\S+\s+(\S+)/;
  $pairs{$1}{reject}{mask} ++;
 }
 close IN;
}
#die scalar(keys %pairs), " pairs masked\n";

open IN, $insam or die "Can't open sam file $insam\n";
while (<IN>) {
 if (/^\@SQ\s+SN:(\S+)\s+LN:(\d+)/) {
  # @SQ     SN:NODE_12_length_119100_cov_398.273_ID_23      LN:119100
  $dnas{$1} = $2;
 }
 next if /^\@/;
 #die scalar(keys %dnas), " dnas\n";
 #$cts{line} ++; last if $cts{line} > 1000;
 my @f = split "\t";
 my ($pair, $flag) = ($f[0], $f[1]);
 next if $pairs{$pair}{reject};
 if ($flag & 4 or $flag & 8) {$pairs{$pair}{reject}{unhit} ++; next}
 if ($flag & 64 and $flag & 128) {die "both forward and reverse? $_\n"}
 unless ($flag & 64 or $flag & 128) {die "neither forward nor reverse? $_\n"}
 my $mate = 1; $mate = 2 if $flag & 128;
 if ($pairs{$pair}{$mate}) {$pairs{$pair}{reject}{multi} ++; next}
 my $dir = 1; $dir = -1 if $flag & 16;
 my $len = 0;
 $f[5] =~ s/(\d+)[MIN]/$len+=$1/eg;
 %{$pairs{$pair}{$mate}} = (dir => $dir, dna => $f[2], pos => $f[3], len => $len);
 next;
}
close IN;
#print scalar(keys %dnas), " dnas\n";

open IN, $ctfile;
while (<IN>) {
 chomp;
 next unless /^(\S+)\t(\S+)/ and $pairs{$1};
 $pairs{$1}{copy} = $2;
}
close IN;

$insam =~ s/\.sam$//;
if ($maskfile) {$insam .= ".$maskfile"; $insam =~ s/\.bed//}
open OUT, ">$insam.bridges";
for my $pair (sort {$a <=> $b} keys %pairs) {
 if ($pairs{$pair}{reject} and $pairs{$pair}{reject}{mask}) {$cts{mask} ++; next}
 if ($pairs{$pair}{reject} and $pairs{$pair}{reject}{unhit}) {$cts{unhit} ++; next}
 if ($pairs{$pair}{reject} and $pairs{$pair}{reject}{multi}) {$cts{multi} ++; next}
 unless ($pairs{$pair}{1} and $pairs{$pair}{2}) {$cts{missing} ++; next} # Shouldn't occur if flags work properly
 my %ends;
 for (1 .. 2) { # find distances 3' to DNA end, for both reads
  if ($pairs{$pair}{$_}{dir} == 1) {%{$ends{$_}} = (end => 'R', dist => $dnas{$pairs{$pair}{$_}{dna}} - $pairs{$pair}{$_}{pos} + 1, pos => $pairs{$pair}{$_}{pos})}
  else {%{$ends{$_}} = (end => 'L', dist => $pairs{$pair}{$_}{pos} + $pairs{$pair}{$_}{len} - 1, pos => $pairs{$pair}{$_}{pos} + $pairs{$pair}{$_}{len} - 1)}
  if ($pairs{$pair}{$_}{dna} =~ /NODE_(\d+)_/) {
   $ends{$_}{dna} = $1;
   $ends{$_}{name} = $1 . $ends{$_}{end};
  } else {
   $ends{$_}{dna} = $pairs{$pair}{$_}{dna};
   $ends{$_}{name} = $pairs{$pair}{$_}{dna} . $ends{$_}{end};
  }
 }
 my $out = "$pair\t$ends{1}{dna}\t$ends{1}{pos}\t$pairs{$pair}{1}{dir}\t$ends{2}{pos}";
 if ($pairs{$pair}{1}{dna} eq $pairs{$pair}{2}{dna}) {
  if ($pairs{$pair}{1}{dir} eq $pairs{$pair}{2}{dir}) {$cts{samedir} ++; print OUT "$out\tSAMEDIR\n"; next}
  my $dist = $pairs{$pair}{2}{pos} - $pairs{$pair}{1}{pos};
  if ($pairs{$pair}{1}{dir} == -1) {
   $dist *= -1;
   $dist += $pairs{$pair}{1}{len} -1;
  } else {$dist += $pairs{$pair}{2}{len} -1}
  if ($dist >= $min and $dist <= $max) {$cts{proper} ++; print OUT "$out\tPROPER\n"; next}
  if ($dist > $max) {$cts{">max"} ++; print OUT "$out\t>MAX\n"; next}
  if ($dist > 0 and $dist < $min) {$cts{"0-min"} ++; print OUT "$out\t0-MIN\n"; next}
  # remainder should be negative distance, but could reveal circularity if within range of end
  if ($ends{1}{dist} + $ends{2}{dist} > $max) {$cts{neg} ++; print OUT "$out\tNEG\n"; next}
 }
 $out = "$pair\t$ends{1}{name}\t$ends{1}{dist}\t$ends{2}{name}\t$ends{2}{dist}";
 if ($ends{1}{dist} + $ends{2}{dist} > $max) {$cts{toofar} ++; print OUT "$out\tTOOFAR\n"; next}
 $cts{bridge} ++;
 print OUT "$out\tBRIDGE\n";
 my $bridge = join "\t", sort($ends{1}{name}, $ends{2}{name});
 push @{$bridges{$bridge}}, $pairs{$pair}{copy};
}
close OUT;

open OUT, ">$insam.bridgesum";
for (qw/mask unhit multi missing samedir proper >max 0-min neg toofar bridge/) {
 if ($cts{$_}) {print OUT "$_\t$cts{$_}\n"} elsif ($_ eq 'missing') {next} else {print OUT "$_\t0\n"}
}
for (sort {@{$bridges{$b}} <=> @{$bridges{$a}}} keys %bridges) { 
 my ($n, $mean, $sd) = NMeanSD(@{$bridges{$_}}); 
 for ($mean, $sd) {$_ = sprintf("%.1f", $_)}
 print OUT "$_\t$n\t$mean\t$sd\n";
}

sub NMeanSD { # returns mean and SD for list of values
 my $n = @_;
 if ($n < 1) {return 0, 0, 0}
 my $mean;
 for (@_) {$mean += $_}
 $mean /= $n;
 if ($n <= 1) {return $n, $mean, 0}
 my $sumsquares = 0;
 for (@_) {$sumsquares += ($_ - $mean) ** 2}
 return $n, $mean, sqrt ($sumsquares/($n - 1));
}

