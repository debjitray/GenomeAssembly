#!/usr/bin/perl -w
# Split a interlaced fastq file and creates two separeted file _1 and _2 for each reads
# USAGE: perl interlacedSplitter.pl S3_MP.mp.fastq
#use IO::Compress::Gzip qw(gzip $GzipError) ;

my $inputfile     = $ARGV[0];
open(FDR,"<$inputfile") or die "Can't open $inputfile: $!\n";
#open(FDR, "gzip -dc $inputfile |");  #

# OUTPUT FOR FORWARD READS
my ($TEMP1,$TEMP2)= split(/\./,$ARGV[0]);

my $outputfile1     = join ("_",$TEMP1,"1.fastq");
open(FDW1, ">$outputfile1") or die "Can't open $outputfile1: $!\n";

my $outputfile2     = join ("_",$TEMP1,"2.fastq");
open(FDW2, ">$outputfile2") or die "Can't open $outputfile2: $!\n";

my $count=0;
my $count1=0;

while(( my @lines = map $_ = <FDR>, 1 .. 8 )[0]) {
  chomp($lines[0]);
  chomp($lines[1]);
  chomp($lines[2]);
  chomp($lines[3]);
  chomp($lines[4]);
  chomp($lines[5]);
  chomp($lines[6]);
  chomp($lines[7]);

  print FDW1 $lines[0]."\n".$lines[1]."\n".$lines[2]."\n".$lines[3]."\n";
  print FDW2 $lines[4]."\n".$lines[5]."\n".$lines[6]."\n".$lines[7]."\n";
  
  $count++;
  #print $count."\n";
}

print "Total line ".$count."\n";
