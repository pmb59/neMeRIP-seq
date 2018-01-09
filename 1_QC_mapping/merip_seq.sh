
# m6A-seq data analysis pipeline
# Pedro Madrigal pmb59[ at ]cam.ac.uk

#--------------------------------------------------#
# Description:
# QC of FASTQ sequencing data using trimmomatic PE 
# Mapping to transcriptome using TopHat ( 1st transcriptome, then genome )
#--------------------------------------------------#
# GTF dowloaded from: ftp://ftp.ensembl.org/pub/release-83/gtf/homo_sapiens/
# FASTA downloaded from: ftp://ftp.ensembl.org/pub/release-83/fasta/homo_sapiens/dna/ >  Homo_sapiens.GRCh38.dna.primary_assembly.fa > GRCh38.fa

#threads;
T=8;

#Bowtie 2 version 2.1.0:
#bowtie2 -- OK
#build index for GRCh38 release 83 (execute just once)
#bowtie2-build GRCh38.fa GRCh38;


#Create tophat dir output
for j in m6a_IP_A1 m6a_IP_A2 m6a_IP_A3 m6a_IP_S1 m6a_IP_S2 m6a_IP_S3   m6a_input2_A1 m6a_input2_A2 m6a_input2_A3 m6a_input2_S1 m6a_input2_S2 m6a_input2_S3 ; do mkdir tophat_${j} ; 
done


RESDIR='../m6A';

GENOME_INDEX="../m6A/GRCh38";
GTF_FILE="../m6A/Homo_sapiens.GRCh38.83.gtf";
GENOME_FILE="../m6A/GRCh38.fa";



for sample in $1;
do

  #QC
  #Trimmomatic
  java -jar trimmomatic-0.35.jar PE -threads $T ${sample}_R1.fastq ${sample}_R2.fastq ${sample}_paired_R1.fastq ${sample}_unpaired_R1.fastq ${sample}_paired_R2.fastq ${sample}_unpaired_R2.fastq  LEADING:3 TRAILING:3 SLIDINGWINDOW:5:10 MINLEN:40 ;
  cat ${sample}_unpaired_R1.fastq ${sample}_unpaired_R2.fastq > ${sample}_unpaired.fastq


  #Tophat
  TOPHATDIR="/tophat_"${sample};

  # Next line only for the 1st run;
  #.tophat2 -p $T --GTF $GTF_FILE --transcriptome-index transcriptome $GENOME_INDEX;
   .tophat2 -p $T --transcriptome-index ${RESDIR}/transcriptome/Homo_sapiens.GRCh38.83  --library-type fr-firststrand -o ${RESDIR}${TOPHATDIR} $GENOME_INDEX ${RESDIR}/${sample}_paired_R1.fastq,${RESDIR}/${sample}_paired_R2.fastq,${RESDIR}/${sample}_unpaired.fastq ; 
  cd ${RESDIR}${TOPHATDIR};
  samtools index accepted_hits.bam;
  #Mapping Quality > 20;
  samtools view -b -h -q 20 -o ${sample}.bam accepted_hits.bam;
  samtools index ${sample}.bam;

done


