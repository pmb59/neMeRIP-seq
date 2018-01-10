#!/bin/bash -l

## Created on 15 February 2016
## @author: Igor Ruiz de los Mozos

### Script to merge bam files
#       input: merged_file file1 file2 file3  !!! all without the bam extension
#
#       output: merged
#
#

merged=$1
file1=$2
file2=$3
file3=$4

ml SAMtools/0.1.19-foss-2016b

echo "samtools merge"
samtools merge -f ${merged}.bam ${file1}.bam ${file2}.bam ${file3}.bam

echo "samtools sort"
samtools sort ${merged}.bam ${merged}.sort

echo "samtools sort ${merged}.bam ${merged}.sort.bam"

echo "samtools index"
samtools index ${merged}.sort.bam

# Clean
rm ${merged}.bam