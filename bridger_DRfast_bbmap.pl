#! /usr/bin/perl -w
use strict;
# nohup perl /data1/users/debray/bridger/bridger_DRfast.pl MP contigs 20 21 1.34 501 15000 PE_1.fq,PE_2.fq &
unless (@ARGV == 8) {
 die "Usage: perl $0 MPreadPrefix referencePrefix AdditionalKmerSource threads\n" . 
  " readPrefix: <path/prefix> for <>_1.fq and <>_2.fq MatePair PE file pair\n" . 
  " referencePrefix: <path/prefix> for either <>.fa or <> bowtie2 index of ref contigs\n" .
  " kmwerSize\n" .
  " toleranceratio\n" .
  " minins\n" .
  " maxins\n" .
  " kmerSource: <comma-sep'd path/filename> for additional read files (eg, nonMP) for kmer balancing\n";
}

my $bowtiepath = 'bowtie2';

my $perlpath = '/data1/users/debray/bridger/';

#my $perlpath = $0; $perlpath =~ s/[^\/]+$//;

my ($readpfx, $refpfx, $threads, $kmers, $tolratio, $minins, $maxins, $kmerSource) = (@ARGV);

for ($readpfx, $refpfx) {s/^.*\///}


# ASSEMBLE END PAIRS FOR MATEPAIR READS [USE IT AS SWITCH]
#unless (-f "$readpfx.assembled.fastq") {system "$pearpath --threads $threads -f $readpfx\_1.fq -r $readpfx\_2.fq -o $readpfx > $readpfx.pear.log"}
#unlink glob("$readpfx.unassembled*");
#unlink glob("$readpfx.discarded*");
# SPLIT ASSEMBLIES INTO PAIRED LEFT AND RIGHT MATEPAIR PARTS
#unless (-f "$readpfx.L.fa") {system "perl $perlpath" . "MPsplit.pl $readpfx";}

# CREATE THE kmers for the input files
#unless (-f "$readpfx.cov") {system "perl $perlpath" . "readMersSC.pl $kmerSource,$readpfx.assembled.fastq 21 $readpfx";}
# cat PE_1.fq |jellyfish count -m 21 -s 10000000 -t 32 -L 11 -C -o TEMP /dev/fd/0
unless (-f "$readpfx.cov") {system "cat PE_*.fq MP_*.fq|jellyfish count -m 21 -s 10000000 -t $threads -L 11 -C -o TEMP_0 /dev/fd/0";}
#unless (-f "$readpfx.cov") {system "cat PE_*.fq MP_*.fq|jellyfish count -m 21 -s 10000M -t $threads -L 11 -C -o TEMP /dev/fd/0";} # Original
system "jellyfish dump -c -t TEMP_0 -o $readpfx.cov";
system "jellyfish histo TEMP_0 -o $readpfx.profile";

# FILTER FOR BALANCED MP L & R COPY NUMBER
# perl copybalancer2_DR.pl MP.cov 21 1.34 MP
unless (-f "$readpfx.Lbalance.fa") {system "perl $perlpath" . "copybalancer2_DRfast.pl $readpfx.cov 21 1.34 $readpfx";}

# ALIGN MATEPAIR PARTS TO CONTIGS
# bowtie2-build contigs.fa MP > /dev/null 2>&1
# bowtie2 --threads 20 --sensitive -k 2 -x MP -1 MP.Rbalance.fa -2 MP.Lbalance.fa -f --minins 500 --maxins 15000 -S MP.sam
#unless (-f "$refpfx.rev.1.bt2[l]") {system $bowtiepath . "-build $refpfx.fa $refpfx > /dev/null 2>&1"}
#unless (-f "$readpfx.$refpfx.sam") {system "$bowtiepath --threads $threads --sensitive -k 2 -x $refpfx -1 $readpfx.Rbalance.fa -2 $readpfx.Lbalance.fa -f --minins $minins --maxins $maxins -S $readpfx.$refpfx.sam > $readpfx.$refpfx.bowtie.log 2>&1";}
unless (-f "$readpfx.$refpfx.sam") {system "/data1/users/debray/SOFTWARES/bbmap/bbmap.sh -Xmx20g maxsites=2 pairlen=$maxins ref=$refpfx.fa in=$readpfx.Rbalance.fa in2=$readpfx.Lbalance.fa out=$readpfx.$refpfx.sam > $readpfx.$refpfx.bowtie.log 2>&1";}


# CLASSIFY UNIQUELY MAPPED MATEPAIRS AS BRIDGES
# perl samBridge_DR.pl MP.sam MP.ctbalance
unless (-f "$readpfx.$refpfx.bridges") {system "perl $perlpath" . "samBridge.pl $readpfx.$refpfx.sam $readpfx.ctbalance";}

# ASSEMBLE
# perl assemble_DR.pl MP.bridgesum contigs.fa > MP.ass
unless (-f "$readpfx.$refpfx.ass") {system "perl $perlpath" . "assemble.pl $readpfx.$refpfx.bridgesum $refpfx.fa > $readpfx.$refpfx.ass";}
