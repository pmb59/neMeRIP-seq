#!/bin/bash -l

## Created on 27th of December 2017
## @author: Igor Ruiz de los Mozos

# Main pipeline to align m6a NeMeRIP-seq data

# Get fastq from Array Express
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/007/ERR1725677/ERR1725677_1.fastq.gz -O Data/m6a_IP_A1_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/007/ERR1725677/ERR1725677_2.fastq.gz -O Data/m6a_IP_A1_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/008/ERR1725678/ERR1725678_1.fastq.gz -O Data/m6a_IP_A2_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/008/ERR1725678/ERR1725678_2.fastq.gz -O Data/m6a_IP_A2_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/009/ERR1725679/ERR1725679_1.fastq.gz -O Data/m6a_IP_A3_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/009/ERR1725679/ERR1725679_2.fastq.gz -O Data/m6a_IP_A3_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/000/ERR1725680/ERR1725680_1.fastq.gz -O Data/m6a_input2_A1_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/000/ERR1725680/ERR1725680_2.fastq.gz -O Data/m6a_input2_A1_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/001/ERR1725681/ERR1725681_1.fastq.gz -O Data/m6a_input2_A2_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/001/ERR1725681/ERR1725681_2.fastq.gz -O Data/m6a_input2_A2_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/002/ERR1725682/ERR1725682_1.fastq.gz -O Data/m6a_input2_A3_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/002/ERR1725682/ERR1725682_2.fastq.gz -O Data/m6a_input2_A3_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/003/ERR1725683/ERR1725683_1.fastq.gz -O Data/m6a_input2_S1_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/003/ERR1725683/ERR1725683_2.fastq.gz -O Data/m6a_input2_S1_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/004/ERR1725684/ERR1725684_1.fastq.gz -O Data/m6a_input2_S2_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/004/ERR1725684/ERR1725684_2.fastq.gz -O Data/m6a_input2_S2_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/005/ERR1725685/ERR1725685_1.fastq.gz -O Data/m6a_input2_S3_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/005/ERR1725685/ERR1725685_2.fastq.gz -O Data/m6a_input2_S3_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/006/ERR1725686/ERR1725686_1.fastq.gz -O Data/m6a_IP_S1_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/006/ERR1725686/ERR1725686_2.fastq.gz -O Data/m6a_IP_S1_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/007/ERR1725687/ERR1725687_1.fastq.gz -O Data/m6a_IP_S2_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/007/ERR1725687/ERR1725687_2.fastq.gz -O Data/m6a_IP_S2_R2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/008/ERR1725688/ERR1725688_1.fastq.gz -O Data/m6a_IP_S3_R1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR172/008/ERR1725688/ERR1725688_2.fastq.gz -O Data/m6a_IP_S3_R2.fastq.gz

# Align sequencing data using STAR aligner
sbatch align_paired_STAR_m6A.sh

# Merge, sort and index bam files sorted by coordinates
bash merge_bam.sh m6a_IP_A_merged UCSC_GRCh37_ens59/m6a_IP_A1/STAR_Output/m6a_IP_A1Aligned.sortedByCoord.out UCSC_GRCh37_ens59/m6a_IP_A2/STAR_Output/m6a_IP_A2Aligned.sortedByCoord.out UCSC_GRCh37_ens59/m6a_IP_A3/STAR_Output/m6a_IP_A3Aligned.sortedByCoord.out
bash merge_bam.sh m6a_IP_S_merged UCSC_GRCh37_ens59/m6a_IP_S1/STAR_Output/m6a_IP_S1Aligned.sortedByCoord.out UCSC_GRCh37_ens59/m6a_IP_S2/STAR_Output/m6a_IP_S2Aligned.sortedByCoord.out UCSC_GRCh37_ens59/m6a_IP_S3/STAR_Output/m6a_IP_S3Aligned.sortedByCoord.out
bash merge_bam.sh m6a_input2_A_merged UCSC_GRCh37_ens59/m6a_input2_A1/STAR_Output/m6a_input2_A1Aligned.sortedByCoord.out UCSC_GRCh37_ens59/m6a_input2_A2/STAR_Output/m6a_input2_A2Aligned.sortedByCoord.out UCSC_GRCh37_ens59/m6a_input2_A3/STAR_Output/m6a_input2_A3Aligned.sortedByCoord.out
bash merge_bam.sh m6a_input2_S_merged UCSC_GRCh37_ens59/m6a_input2_S1/STAR_Output/m6a_input2_S1Aligned.sortedByCoord.out UCSC_GRCh37_ens59/m6a_input2_S2/STAR_Output/m6a_input2_S2Aligned.sortedByCoord.out UCSC_GRCh37_ens59/m6a_input2_S3/STAR_Output/m6a_input2_S3Aligned.sortedByCoord.out

# Explore on IGV representative candidates

