#!/bin/bash

FILE=$1

# bigwigs.sh #######

bamDIR="../m6A/tophat_"${FILE}"/"
outDIR="../m6A/"
chromSizes="../m6A/GRCh38.sizes"
#cat GRCh38.fa.fai | awk '{print $1"\t" $2}' > GRCh38.sizes

#1-bam to wig
peakranger wig -d ${bamDIR}${FILE}.bam --format bam -l 0 -o ${outDIR}${FILE}.wig;
#2-wig to bw
wigToBigWig ${outDIR}${FILE}.wig $chromSizes ${outDIR}${FILE}.bw;
#3-bw to norm_bgr
cd ../RSeQC-2.6/scripts;
normalize_bigwig.py -i ${outDIR}${FILE}.bw -o ${outDIR}${FILE}.normalised.bgr -s $chromSizes -f bgr;
#4-norm_bgr to norm_bw
cd $outDIR
bedGraphToBigWig ${outDIR}${FILE}.normalised.bgr $chromSizes ${outDIR}${FILE}.normalised.bw;
#5-remove other files
rm ${outDIR}${FILE}.wig;
rm ${outDIR}${FILE}.bw;
rm ${outDIR}${FILE}.normalised.bgr;


