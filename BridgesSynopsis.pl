#!/usr/bin/perl -w
# perl /data1/users/debray/scripts/BridgesSynopsis.pl No_of_contigs_in_contigs.fa

use strict;

my $inputfile     = "MP.contigs.bridges";
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";

my $outputfile     = "MP_ProperImProper_Sumary.txt";
open(FDW,">$outputfile") or die "Can't open $outputfile: $!\n";

my %hashGMAX;
my %hash0MIN;
my %hashBRIDGE;
my %hashNEG;
my %hashTOOFAR;
my %hashSAMEDIR;
my %hashPROPER;



for (my $i=1; $i <=$ARGV[0]; $i++) {
  for (my $j=1; $j <=$ARGV[0] ; $j++) {
    $hashGMAX{$i}{$j}=0;
    $hash0MIN{$i}{$j}=0;
    $hashBRIDGE{$i}{$j}=0;
    $hashNEG{$i}{$j}=0;
    $hashTOOFAR{$i}{$j}=0;
    $hashSAMEDIR{$i}{$j}=0;
    $hashPROPER{$i}{$j}=0;
  }
}

my $count=0;

foreach(<FDR>){
  (my $line = $_) =~ s/\r*\n//;
  next if ($line =~ /^Read_ID/);
  my (@array) = split (/\t/,$line);
  $count++;
  if ($array[3] =~ /^1$/ or $array[3] =~ /^\-1$/ ) {
    if ($array[5] =~ /\>MAX/) {
      $hashGMAX{$array[1]}{$array[1]}++;
    }
    if ($array[5] =~ /0\-MIN/) {
      $hash0MIN{$array[1]}{$array[1]}++;
    }
    if ($array[5] =~ /BRIDGE/) {
      $hashBRIDGE{$array[1]}{$array[1]}++;
    }
    if ($array[5] =~ /NEG/) {
      $hashNEG{$array[1]}{$array[1]}++;
    }
    if ($array[5] =~ /TOOFAR/) {
      $hashTOOFAR{$array[1]}{$array[1]}++;
    }
    if ($array[5] =~ /SAMEDIR/) {
      $hashSAMEDIR{$array[1]}{$array[1]}++;
    }
    if ($array[5] =~ /PROPER/) {
      $hashPROPER{$array[1]}{$array[1]}++;
    }
    
  }
  else {
    my %temp;
    $array[1] =~ s/[LR]//;
    $array[3] =~ s/[LR]//;
    if ($array[5] =~ /\>MAX/) {
      $hashGMAX{$array[1]}{$array[3]}++;
      $hashGMAX{$array[3]}{$array[1]}++;
    }
    if ($array[5] =~ /0\-MIN/) {
      $hash0MIN{$array[1]}{$array[3]}++;
      $hash0MIN{$array[3]}{$array[1]}++;
    }
    if ($array[5] =~ /BRIDGE/) {
      if ($array[1] == $array[3]) {
	$hashBRIDGE{$array[1]}{$array[3]}++;
      }
      else {
	$hashBRIDGE{$array[1]}{$array[3]}++;
	$hashBRIDGE{$array[3]}{$array[1]}++;
      }
    }
    if ($array[5] =~ /NEG/) {
      $hashNEG{$array[1]}{$array[3]}++;
      $hashNEG{$array[3]}{$array[1]}++;
    }
    if ($array[5] =~ /TOOFAR/) {
      $hashTOOFAR{$array[1]}{$array[3]}++;
      $hashTOOFAR{$array[3]}{$array[1]}++;
    }
    if ($array[5] =~ /SAMEDIR/) {
      $hashSAMEDIR{$array[1]}{$array[3]}++;
      $hashSAMEDIR{$array[3]}{$array[1]}++;
    }
    if ($array[5] =~ /PROPER/) {
      $hashPROPER{$array[1]}{$array[3]}++;
      $hashPROPER{$array[3]}{$array[1]}++;
    }
  }
 
}

my @aks;
my $k=0;
while ( $k < $ARGV[2]) {
	$k++;
	push(@aks, join("","NODE_",$k))
}
#print join("\t",@aks)."\n";


print FDW "\n\n";
print FDW "#BRIDGE\n";
print FDW "\t".join("\t",@aks)."\n";
my $o=0;

foreach my $i(sort {$a<=>$b} keys %hashBRIDGE){
  print FDW $aks[$o]."\t";
  $o++;
  foreach my $j(sort {$a<=>$b} keys %{$hashBRIDGE{$i}}){
    print FDW $hashBRIDGE{$i}{$j}."\t";
   # print $i."\.".$j."\n";
  }
  print FDW "\n";
}

print FDW "\n\n";

print FDW "#MAX\n";
print FDW "\t".join("\t",@aks)."\n";
my $o1=0;
foreach my $i(sort {$a<=>$b} keys %hashGMAX){
    print FDW $aks[$o1]."\t";
    $o1++;
    foreach my $j(sort {$a<=>$b} keys %{$hashGMAX{$i}}){
      print FDW $hashGMAX{$i}{$j}."\t";
    #  print $i."\.".$j."\n";
    }
    print FDW "\n";
}

print FDW "\n\n";
print FDW "#0-MIN\n";
print FDW "\t".join("\t",@aks)."\n";
my $o2=0;
foreach my $i(sort {$a<=>$b} keys %hash0MIN){
  print FDW $aks[$o2]."\t";
  $o2++;
  foreach my $j(sort {$a<=>$b} keys %{$hash0MIN{$i}}){
    print FDW $hash0MIN{$i}{$j}."\t";
   # print $i."\.".$j."\n";
  }
  print FDW "\n";
}

print FDW "\n\n";
print FDW "#NEG\n";
print FDW "\t".join("\t",@aks)."\n";
my $o3=0;
foreach my $i(sort {$a<=>$b} keys %hashNEG){
  print FDW $aks[$o3]."\t";
  $o3++;
  foreach my $j(sort {$a<=>$b} keys %{$hashNEG{$i}}){
    print FDW $hashNEG{$i}{$j}."\t";
#    print $i."\.".$j."\n";
  }
  print FDW "\n";
}

print FDW "\n\n";
print FDW "#TOOFAR\n";
print FDW "\t".join("\t",@aks)."\n";
my $o4=0;
foreach my $i(sort {$a<=>$b} keys %hashTOOFAR){
  print FDW $aks[$o4]."\t";
  $o4++;
  foreach my $j(sort {$a<=>$b} keys %{$hashTOOFAR{$i}}){
    print FDW $hashTOOFAR{$i}{$j}."\t";
  #  print $i."\.".$j."\n";
  }
  print FDW "\n";
}

print FDW "\n\n";
print FDW "#SAMEDIR\n";
print FDW "\t".join("\t",@aks)."\n";
my $o5=0;
foreach my $i(sort {$a<=>$b} keys %hashSAMEDIR){
  print FDW $aks[$o5]."\t";
  $o5++;
  foreach my $j(sort {$a<=>$b} keys %{$hashSAMEDIR{$i}}){
    print FDW $hashSAMEDIR{$i}{$j}."\t";
   # print $i."\.".$j."\n";
  }
  print FDW "\n";
}

print FDW "\n\n";
print FDW "#PROPER\n";
print FDW "\t".join("\t",@aks)."\n";
my $o6=0;
foreach my $i(sort {$a<=>$b} keys %hashPROPER){
  print FDW $aks[$o6]."\t";
  $o6++;
  foreach my $j(sort {$a<=>$b} keys %{$hashPROPER{$i}}){
    print FDW $hashPROPER{$i}{$j}."\t";
    #print $i."\.".$j."\n";
  }
  print FDW "\n";
}



print $count."\n";
