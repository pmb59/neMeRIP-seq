#!/bin/bash
# Simple SLURM sbatch array job
#
## Created on 27th of December 2017
## @author: Igor Ruiz de los Mozos
#
#SBATCH --job-name=Bertero
#SBATCH -n 12                               # number of cores
#SBATCH -t 1-00:00 # time (D-HH:MM)         # one day
#SBATCH --mem 42G                           # memory pool for all cores
#SBATCH --partition=compute
#
#SBATCH --array=1-12                        # Array of 12 jobs
#SBATCH -o S.%A_%a.out # STDOUT
#SBATCH -e S.%A_%a.err # STDERR
#SBATCH --mail-type=END,FAIL                # notifications for job done & fail
#SBATCH --mail-user=i.mozos {at} ucl.ac.uk  # send-to mail address

### Mapping with STAR aligner    version STAR_2.5.2a
### Usage:  align_pairedend_STAR.sh
###
### fastq.gz acepted

# Number of treads per job
tread="1"

date

# Inputs
real_path=`pwd`    ### path_to_directory

chr_GTF="/camp/lab/ulej/working/Igor/STAR/Genomes/hg19_ensembl59/hg59.gtf"    ## UCSC annotation hg19 ENSEMBL v. 59
genome_dir="/camp/lab/ulej/working/Igor/STAR/Genomes/hg19_ensembl59"
genome_fasta_seq="/camp/lab/ulej/working/Igor/STAR/Genomes/hg19_ensembl59/hg59.fa"
genome="hg19_ensembl59"
chr_size="/camp/lab/ulej/working/Igor/STAR/Genomes/hg19_ensembl59/chrNameLength.txt"

# Results folder
results_folder="UCSC_GRCh37_ens59"

# Index for array job
INDEX=$((SLURM_ARRAY_TASK_ID-1)) ## For array job
echo "$INDEX"

fastq_extensionR1="_R1.fastq.gz"
fastq_extensionR2="_R2.fastq.gz"


# fastq files names e.i: m6a_IP_A1_R1.fastq.gz & m6a_IP_A1_R2.fastq.gz
FILES=("m6a_IP_A1"
    "m6a_IP_A2"
    "m6a_IP_A3"
    "m6a_input2_A1"
    "m6a_input2_A2"
    "m6a_input2_A3"
    "m6a_input2_S1"
    "m6a_input2_S2"
    "m6a_input2_S3"
    "m6a_IP_S1"
    "m6a_IP_S2"
"m6a_IP_S3")

# Selected file on array job
f="${FILES[$INDEX]}"

# Print file name for records
echo "Align_paired_bertero.sh $SGE_TASK_ID $f"

# Modify file name
file=$(basename "$f")
Complete_filename="${file%.*}"
filename="${Complete_filename/.fastq}"
extension="${file##*.}"

# Quality check
ml FastQC/0.11.5-Java-1.7.0_80
fastqc ${real_path}/Data/${filename}${fastq_extensionR1}
fastqc ${real_path}/Data/${filename}${fastq_extensionR2}

# Align with STAR
ml STAR/2.5.2a-foss-2016b
mkdir -p ${real_path}/${results_folder}/${filename}/STAR_Output/
time STAR --runThreadN ${tread} --genomeDir $genome_dir --readFilesCommand zcat \
    --readFilesIn ${real_path}/Data/${filename}${fastq_extensionR1} ${real_path}/Data/${filename}${fastq_extensionR2}\
    --winAnchorMultimapNmax 100 --outFilterMultimapNmax 20 --outFileNamePrefix ${real_path}/${results_folder}/${filename}/STAR_Output/${filename} \
    --outSAMtype BAM SortedByCoordinate --outWigType wiggle --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic

# Transform wig files to TDF compressed and sort
ml IGVTools/2.3.75-Java-1.7.0_80
time igvtools toTDF ${real_path}/${results_folder}/${filename}/STAR_Output/${filename}Signal.UniqueMultiple.str1.out.wig ${real_path}/${results_folder}/${filename}/STAR_Output/${filename}Signal.UniqueMultiple.str1.out.tdf ${chr_size}
igvtools toTDF ${real_path}/${results_folder}/${filename}/STAR_Output/${filename}Signal.Unique.str1.out.wig ${real_path}/${results_folder}/${filename}/STAR_Output/${filename}Signal.Unique.str1.out.tdf ${chr_size}

# Prepare index and .sam files"
ml SAMtools/0.1.19-foss-2016b
samtools index ${real_path}/${results_folder}/${filename}/STAR_Output/${filename}Aligned.sortedByCoord.out.bam ${real_path}/${results_folder}/${filename}/STAR_Output/${filename}Aligned.sortedByCoord.out.bam.bai

echo "Finished ~~~~~~~~~~~~~~~~~"


date

exit

