#!/usr/bin/sh

reference='/data1/users/debray/Cdiff_FINISHED/ClostridiumDB/BI1.fasta'

/data1/users/debray/CdBk/X.BLASTN/MUMmer3.23/nucmer --maxmatch -c 500 -p X.MummerPlot $reference MP.contigs.scaffolds 
/data1/users/debray/CdBk/X.BLASTN/MUMmer3.23/mummerplot -postscript -p MP.contigs.scaffolds.ps X.MummerPlot.delta

unset DISPLAY

java -Xmx500m -cp /data1/users/debray/SOFTWARES/mauve/Mauve.jar org.gel.mauve.contigs.ContigOrderer -output results -ref $reference -draft MP.contigs.scaffolds

x="ls -t1 results/ |  head -n 1"
y=`eval $x`
echo $y
scp -r results/$y/MP.contigs.scaffolds.fas MP.contigs.scaffolds.fas_Arranged.fasta

rm -r results

/data1/users/debray/CdBk/X.BLASTN/MUMmer3.23/nucmer --maxmatch -c 500 -p X.MummerPlot $reference MP.contigs.scaffolds.fas_Arranged.fasta
/data1/users/debray/CdBk/X.BLASTN/MUMmer3.23/mummerplot -postscript -p MP.contigs.scaffolds.fas_Arranged.ps X.MummerPlot.delta
