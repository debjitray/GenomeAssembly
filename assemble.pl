#! /usr/bin/perl
use strict; use warnings;

die "Usage: $0 brijsumfile contigfile [gfffile]\n" unless @ARGV >= 2;
my ($brijsumfile, $contigfile, $gfffile) = @ARGV;
my (%ends, @joints, $collect, %ass, %contigs, %used, @ignored, %gffin, %gffout);
my ($serial, $circle, $linear) = (0, 0, 0);

open IN, $brijsumfile;
while (<IN>) {
 if (/^bridge/) {$collect ++; next}
 next unless $collect;
 chomp;
 my ($L, $R, $ct, $copy, $copysd) = split "\t";
 push @joints, {L => $L, R => $R, brij => $ct, copy => $copy, copysd => $copysd};
 for ($L, $R) {
  next if $ends{$_}; # only top bridge occurrence
  /(.*)[LR]$/; 
  %{$ends{$_}} = (serial => $serial, contig => $1, oppcontig => $_);
  $ends{$_}{oppcontig} =~ tr/LR/RL/; # can't delimit to last char only
 }
}
close IN;
print scalar(@joints), " joints; ", scalar(keys %ends), " ends\n";

# assemblies are arrays named for oppcontig of 
for my $joint (@joints) {
 #for (qw/L R brij copy/) {print "$_=$$joint{$_}; "} print "\n";
 my (@seen, @new);
 my $name = "$$joint{L}-$$joint{R}/$$joint{brij},$$joint{copy}";
 for (qw/L R/) {push @seen, $$joint{$_} if $ends{$$joint{$_}}{used} and not $ass{$$joint{$_}}}
 if (@seen) {push @ignored, $name; next}
 for (qw/L R/) {
  push @seen, $$joint{$_} if $ass{$$joint{$_}};
  $ends{$$joint{$_}}{used} ++;
 }
 my %mate = ($$joint{L} => $$joint{R}, $$joint{R} => $$joint{L});
 if (@seen == 0) { # neither end in assembly
  print "$name NEW\n";
  for ($$joint{L}, $$joint{R}) { # make fwd and rev version of new scaffold
   @{$ass{$ends{$_}{oppcontig}}} = ($mate{$_});
   push @new, $ends{$_}{oppcontig};
  }
 } elsif (@seen == 1) { # append to existing scaffold, in both directions   
  print "$name APPENDS\n";
  for ($$joint{L}, $$joint{R}) { 
   next unless $seen[0] eq $_;
   my ($comp, $newass) = ($ends{$ass{$_}[-1]}{oppcontig}, $ends{$mate{$_}}{oppcontig});
   push @{$ass{$comp}}, $mate{$_}; # extend complementary assem
   unshift @{$ass{$_}}, $_; # extend seen assem
   $ass{$newass} = delete $ass{$_}; # rename seen assem
   @new = ($comp, $newass);
  }
 } elsif (@seen == 2 and $ends{$$joint{L}}{oppcontig} eq $ass{$$joint{R}}[-1]) { # circle
  print "$name $ends{$ass{$$joint{L}}[-1]}{oppcontig} MULTI-CONTIG CIRCULARIZES\n";
  for ($$joint{L}, $$joint{R}) {push @{$ass{$_}}, $_; push @new, $_}
 } else { # merge two scaffolds
  print "$name $ends{$ass{$$joint{L}}[-1]}{oppcontig} MERGES\n";
  for ($$joint{L}, $$joint{R}) { # make merger in both directions
   my $comp = $ends{$ass{$_}[-1]}{oppcontig};
   #print "$_, $ass{$_}[-1], $comp, $ass{$comp}[-1], $mate{$_}, $ass{$mate{$_}}[-1]\n";
   push @{$ass{$comp}}, $mate{$_}, @{$ass{$mate{$_}}}; # extend complementary assem
   push @new, $comp; 
  }
  for ($$joint{L}, $$joint{R}) {delete $ass{$_}}
 }
 if ($new[0] eq $ass{$new[0]}[-1]) { # close out and singularize recent scaffold if now circular
  $circle ++;
  delete $ass{$new[1]}; # remove reverse scaffold
  #print scalar(keys %ass), " assems $_ circle $circle\n"; for (keys %ass) {print "  $_: \n"}
  $ass{"Circle$circle"} = delete $ass{$new[0]}; # change to inert hash key name
  print "CIRCLE\n"; 
 }
 #for (keys %ass) {print "  $_: ", join(',', @{$ass{$_}}), "\n"}
}
print "Ignored due to previous incorporation: ", join(', ', @ignored), "\n" if @ignored;

open IN, $contigfile or die "No contig file $contigfile\n";
while (<IN>) {
 if (/^>(NODE_(\d+)_length_(\d+)_cov_(\d+)\S+)/) {%{$contigs{$2}} = (len => $3, cov => $4, name => $1); $collect = $2; next}
 chomp;
 $contigs{$collect}{seq} .= $_;
}
close IN;

my (%seen, %scaffs, %lens);
for my $ass (keys %ass) { # remove reversed duplicates
 next unless $ass{$ass};
 my $comp = $ends{$ass{$ass}[-1]}{oppcontig};
 if ($ass{$ass}[-1]) {$comp = $ends{$ass{$ass}[-1]}{oppcontig}} else {next}
 if ($ass{$comp}) {
  delete $ass{$ass};
  $linear ++;
  unshift @{$ass{$comp}}, $comp;
  $ass{"Linear$linear"} = delete $ass{$comp};
  next;
 }
}

if ($gfffile) {
 open IN, $gfffile;
 while (<IN>) {
  chomp;
  next if /^#/;
  my @f = split "\t";
  next if /^#/ or @f != 9;
  $f[0] = $1 if $f[0] = /^NODE_(\d+)/;
  push @{$gffin{$f[0]}}, [];
  for (@f) {push @{$gffin{$f[0]}[-1]}, $_}
 }
}

for my $ass (keys %ass) { # assemble
 my ($scaff, $offset) = ('', 0);
 for my $end (@{$ass{$ass}}) {
  next unless $end =~ /(.*)[LR]$/;
  my $contig = $1;
  $seen{$1} ++;
  my ($len, $sign) = (length($contigs{$contig}{seq}), '+');
  if ($scaff or $ass =~ /^Circle/) {$scaff .= 'N' x 50; $offset += 50}
  if ($end =~ /L$/) {$scaff .= $contigs{$contig}{seq}}
  else {$scaff .= Revcomp($contigs{$contig}{seq}); $sign = '-'}
  $gffout{$ass}[0] = join "\t", 'SPAdes', 'contig', $offset+1, $offset+$len, $contigs{$contig}{cov}, $sign, '.', "ID=contig$contig";
  if ($gffin{$contig} and $sign eq '+') {Gff($ass, $contig, $offset)}
  elsif ($gffin{$contig}) {GffRev($ass, $contig, $offset+$len+1)}
  $offset += $len;
 }
 $scaffs{$ass} = $scaff;
 $lens{$ass} = length $scaff;
}

my $longest = (sort {$contigs{$b}{len} <=> $contigs{$a}{len}} keys %contigs)[0];
for (sort {$a <=> $b} keys %contigs) {
 $contigs{$_}{relcov} = $contigs{$_}{cov}/$contigs{$longest}{cov};
 print "$_ $contigs{$_}{relcov}\n";
}

$brijsumfile =~ s/\.bridgesum//; $brijsumfile .= ".scaffolds";
open OUT, ">$brijsumfile";
$brijsumfile =~ s/\.scaffolds//; $brijsumfile .= ".gff"; #my $n; while (-f $brijsumfile) {$n ++; $brijsumfile =~ s/\.gff\d+$/.gff$n/}
my ($contigct, $unseenct) = (scalar(keys %contigs), scalar(keys %contigs) - scalar(keys %seen));
open GFF, ">$brijsumfile";
#print scalar(keys %ass), " scaffolds\n"; for (keys %ass) {print "$_: ", join(',', @{$ass{$_}}), "\n"}
print "\nFinal Assembly: $circle circular scaffolds + $linear linear scaffolds + $unseenct unassembled contigs of $contigct original contigs\n"; 
for my $type ('Circle', 'Linear') {
 my $ct;
 for my $ass (sort {$lens{$b} <=> $lens{$a}} keys %lens) {
  next unless $ass =~ /^$type/;
  $ct ++;
  print "$type$ct, $lens{$ass} bp: ";
  my ($longestDna, $longestLen) = ('', 0);
  for (@{$ass{$ass}}) {
   unless (/(.*)([LR])$/) {print "$_,"; next}
   print "$2$1"; if ($2 eq 'L') {print 'R,'} else {print 'L,'}
   ($longestDna, $longestLen) = ($1, $contigs{$1}{len}) if $contigs{$1}{len} > $longestLen;
  }
  my $relcov = 0; $relcov = $contigs{$longestDna}{relcov} if $contigs{$longestDna};
  print "\n"; 
  print OUT ">$type$ct ", join(',', @{$ass{$ass}}), " relcov=$relcov\n$scaffs{$ass}\n";
  for (@{$gffout{$ass}}) {print GFF "$type$ct\t$_\n"}
 }
}

for (sort {$contigs{$b}{len} <=> $contigs{$a}{len}} keys %contigs) {
 next if $seen{$_};
 my $name = "Contig$_";
 my $relcov = 0; $relcov = $contigs{$_}{relcov} if $contigs{$_};
 print "$name, $contigs{$_}{len} bp\n";
 print OUT ">$name $contigs{$_}{name} relcov=$relcov\n$contigs{$_}{seq}\n";
 print GFF join("\t", $name, 'SPAdes', 'contig', 1, $contigs{$_}{len}, $contigs{$_}{cov}, '+', '.', "ID=contig$_"), "\n";
 if ($gffin{$_}) {for (@{$gffin{$_}}) {$$_[0] = $name; print GFF join("\t", @{$_})}}
}
close OUT;
close GFF;

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGT/TGCA/; return $ret}

sub Gff {
 my ($ass, $contig, $offset) = @_;
 my @out;
 for (@{$gffin{$contig}}) {
  $$_[3] += $offset;
  $$_[4] += $offset;
  shift @{$_};
  push @out, [@{$_}];
 }
 for (sort {$$a[2] <=> $$b[2] || $$b[3] <=> $$a[3]} @out) { # note: have removed usual $$_[0], so coords are 2 & 3 now
  push @{$gffout{$ass}}, join("\t", @{$_});
 } 
}

sub GffRev {
 my ($ass, $contig, $top) = @_;
 my @out;
 for (@{$gffin{$contig}}) {
  ($$_[3], $$_[4]) = ($top-$$_[4], $top-$$_[3]);
  $$_[6] =~ tr/-+/+-/;
  shift @{$_};
  push @out, [@{$_}];
 }
 for (sort {$$a[2] <=> $$b[2] || $$b[3] <=> $$a[3]} @out) { # note: have removed usual $$_[0], so coords are 2 & 3 now
  push @{$gffout{$ass}}, join("\t", @{$_});
 } 
}

