#! /usr/bin/perl
use strict; use warnings;
use Getopt::Long;

my $verbose;
my $readpfx;
my $refpfx;
my $threads;
my $kmers;
my $tolratio;
my $minins;
my $maxins;
my $kmerSource;
my $scriptname = $0;
my $scriptpath = '.';
$scriptpath = $1 if $0 =~ /^(.+)\/[^\/]+$/;
my $VERSION = '0.1';

if (@ARGV < 1){
    print "\n Try '$scriptname --man' for full info\n\n";
    exit(0);
}
else{
	GetOptions('help' => sub {pod2usage(1);},
		'version' => sub {print STDOUT "\n $scriptname version $VERSION\n"; exit()},
		'man' => sub {pod2usage(-exitstatus => 0, -verbose => 2);},
		'readpfx=s' => \$readpfx,
		'refpfx=s' => \$refpfx,
		'threads=i' => \$threads,
		'kmers=i' => \$kmers,
		'tolratio=f' => \$tolratio,
		'minins=i' => \$minins,
		'maxins=i' => \$maxins,
		'kmerSource=s' => \$kmerSource,
	);
}

for ($readpfx, $refpfx) {s/^.*\///}

# CREATE THE kmers for the input files
unless (-f "$readpfx.cov") {system "cat PE_*.fq MP_*.fq|jellyfish count -m $kmers -s 10000000 -t $threads -L 11 -C -o TEMP_0 /dev/fd/0";}
system "jellyfish dump -c -t TEMP_0 -o $readpfx.cov";
system "jellyfish histo TEMP_0 -o $readpfx.profile";

# FILTER FOR BALANCED MP L & R COPY NUMBER
unless (-f "$readpfx.Lbalance.fa") {system "perl $scriptpath/copybalancer2_DRfast.pl $readpfx.cov $kmers $tolratio $readpfx";}

# ALIGN MATEPAIR PARTS TO CONTIGS
#unless (-f "$readpfx.$refpfx.sam") {system "/data1/users/debray/SOFTWARES/bbmap/bbmap.sh -Xmx20g maxsites=2 pairlen=$maxins ref=$refpfx.fa in=$readpfx.Rbalance.fa in2=$readpfx.Lbalance.fa out=$readpfx.$refpfx.sam > $readpfx.$refpfx.bowtie.log 2>&1";}
unless (-f "$refpfx.rev.1.bt2[l]") {system "bowtie2-build $refpfx.fa $refpfx > /dev/null 2>&1"}
unless (-f "$readpfx.$refpfx.sam") {system "bowtie2 --threads $threads --sensitive -k 2 -x $refpfx -1 $readpfx.Rbalance.fa -2 $readpfx.Lbalance.fa -f --minins $minins --maxins $maxins -S $readpfx.$refpfx.sam > $readpfx.$refpfx.bowtie.log 2>&1";}

# CLASSIFY UNIQUELY MAPPED MATEPAIRS AS BRIDGES
unless (-f "$readpfx.$refpfx.bridges") {system "perl $scriptpath/samBridge.pl $readpfx.$refpfx.sam $readpfx.ctbalance";}

# ASSEMBLE
unless (-f "$readpfx.$refpfx.ass") {system "perl $scriptpath/assemble.pl $readpfx.$refpfx.bridgesum $refpfx.fa > $readpfx.$refpfx.ass";}

# Bridger summary
my $contig_count= `grep ">" $refpfx.fa | wc -l`;
unless (-f "MP_ProperImProper_Sumary.txt") {system "perl $scriptpath/BridgesSynopsis.pl $contig_count";}

print "Recommended next step to compare assembly to reference <reference.fa>\n\$ sh $scriptpath/Scaffold_Arranger.sh <reference.fa>\n";

