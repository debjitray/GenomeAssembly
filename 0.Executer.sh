#!/usr/bin/sh

ln -s /data1/users/debray/Cdiff_FINISHED/BBDUK/All_PE/Cd$1_PE_clean1.fq PE_1.fq
ln -s /data1/users/debray/Cdiff_FINISHED/BBDUK/All_PE/Cd$1_PE_clean2.fq PE_2.fq

ln -s /data1/users/debray/Cdiff_FINISHED/BBDUK/All_MP/Cd$1_MP_clean1.fq MP_1.fq
ln -s /data1/users/debray/Cdiff_FINISHED/BBDUK/All_MP/Cd$1_MP_clean2.fq MP_2.fq

#ln -s /data1/users/debray/Cdiff_FINISHED/BBDUK/All_MP/Cd$1_MP_1.fastq Cd$1_MP_1.fastq
#ln -s /data1/users/debray/Cdiff_FINISHED/BBDUK/All_MP/Cd$1_MP_2.fastq Cd$1_MP_2.fastq

#Reusing the unknown mates
scp -r /data1/users/debray/Cdiff_FINISHED/BBDUK/All_MP/Cd$1_MP_1.fastq Cd$1_MP_1.fastq
scp -r /data1/users/debray/Cdiff_FINISHED/BBDUK/All_MP/Cd$1_MP_2.fastq Cd$1_MP_2.fastq
cat /data1/users/debray/Cdiff_FINISHED/BBDUK/All_MP/Cd$1_Unknown_MPExtracted_1.fq >> Cd$1_MP_1.fastq
cat /data1/users/debray/Cdiff_FINISHED/BBDUK/All_MP/Cd$1_Unknown_MPExtracted_2.fq >> Cd$1_MP_2.fastq

perl /data1/users/debray/scripts/changerFAtoFQ.pl Cd$1_MP_1.fastq Cd$1_MP_1_renamed.fasta
perl /data1/users/debray/scripts/changerFAtoFQ.pl Cd$1_MP_2.fastq Cd$1_MP_2_renamed.fasta

mkdir bridger
ln -s /data1/users/debray/CdBk/Assembly/S$1_NxTrimmed/TEST_SPADES_TrustedContigs/scaffolds_Selected.fasta bridger/contigs.fa

# Non SPAdes output formatter with contig length and dummy coverage
#perl /data1/users/debray/scripts/MakeBridgerContigs /data1/users/debray/CdBk/Assembly/S$1_NxTrimmed/TEST_SPADES_TrustedContigs/scaffolds_Selected.fasta bridger/contigs.fa

cat Cd$1_MP_1.fastq Cd$1_MP_2.fastq  > bridger/MP.assembled.fastq

fastx_reverse_complement -Q 33 -i Cd$1_MP_1_renamed.fasta -o bridger/MP.L.fa
fastx_reverse_complement -Q 33 -i Cd$1_MP_2_renamed.fasta -o bridger/MP.R.fa

cd bridger

ln -s ../MP_1.fq MP_1.fq 
ln -s ../MP_2.fq MP_2.fq 
ln -s ../PE_1.fq PE_1.fq 
ln -s ../PE_2.fq PE_2.fq

#nohup perl /data1/users/debray/bridger_Publication/bridger_DRfast.pl MP contigs 20 21 1.34 501 15000 PE_1.fq,PE_2.fq &
perl /data1/users/debray/bridger_Publication/bridger_DRfast.pl --readpfx MP --refpfx contigs --threads 20 --kmers 21 --tolratio 1.34 --minins 501 --maxins 15000 --kmerSource PE_1.fq,PE_2.fq

