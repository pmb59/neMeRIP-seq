#--------------------------------------------------------------------
# METDIFF
#--------------------------------------------------------------------

# gitHub:
# https://github.com/compgenomics/MeTDiff

setwd("../m6A/MetDiff")

library(devtools)
library(MeTDiff)

# example data included in MetDiff package:

##Toy Example
## using the data included in the package
## in the real case, change the gtf to what you need
#gtf <- system.file('extdata','example.gtf',package='MeTDiff')
#
#ip1 <- system.file('extdata','IP1.bam',package='MeTDiff')
#ip2 <- system.file('extdata','IP2.bam',package='MeTDiff')
#ip3 <- system.file('extdata','IP3.bam',package='MeTDiff')
#input1 <- system.file('extdata','Input1.bam',package='MeTDiff')
#input2 <- system.file('extdata','Input2.bam',package='MeTDiff')
#input3 <- system.file('extdata','Input3.bam',package='MeTDiff')
#treated_ip <- system.file('extdata','treated_IP1.bam',package='MeTDiff')
#treated_input <- system.file('extdata','treated_Input1.bam',package='MeTDiff')
#
#IP_BAM <- c(ip1,ip2,ip3)
#INPUT_BAM <- c(input1,input2,input3)
#TREATED_IP_BAM <- c(treated_ip)
#TREATED_INPUT_BAM <- c(treated_input)
#
#metdiff(GENE_ANNO_GTF=gtf,IP_BAM = IP_BAM,INPUT_BAM = INPUT_BAM, TREATED_IP_BAM = TREATED_IP_BAM,TREATED_INPUT_BAM=TREATED_INPUT_BAM, EXPERIMENT_NAME="example")


#--------------------------------------------------------------------
# ACT vs SB
#--------------------------------------------------------------------
# in the real case, change the gtf to what you need
gtf <- "../m6A/Homo_sapiens.GRCh38.83.gtf"

ip1 <- "../m6A/tophat_m6a_IP_A1/m6a_IP_A1.bam"
ip2 <- "../m6A/tophat_m6a_IP_A2/m6a_IP_A2.bam"
ip3 <- "../m6A/tophat_m6a_IP_A3/m6a_IP_A3.bam"

input  <- "../m6A/tophat_m6a_input2_A/m6a_input2_A.bam"

treated_ip1 <- "../m6A/tophat_m6a_IP_S1/m6a_IP_S1.bam"
treated_ip2 <- "../m6A/tophat_m6a_IP_S2/m6a_IP_S2.bam"
treated_ip3 <- "../m6A/tophat_m6a_IP_S3/m6a_IP_S3.bam"

treated_input  <- "../m6A/tophat_m6a_input2_S/m6a_input2_S.bam"



IP_BAM <- c(ip1,ip2,ip3)
INPUT_BAM <- input  
TREATED_IP_BAM <- c(treated_ip1, treated_ip2, treated_ip3)
TREATED_INPUT_BAM <- treated_input  


EXPERIMENT_NAME <- "act_vs_sb_param1"
WINDOW_WIDTH <- 40 #an integer, which specifies the bin width of the sliding window, default: 50
SLIDING_STEP <- 20 #an integer, which specifies the step of the sliding window, usually set as bin width as HMM models the dependcy between continous bins, default: 50
FRAGMENT_LENGTH <- 250 #which specifies the fragment length in the library preparation, default: 100
PEAK_CUTOFF_PVALUE <- 1e-3  # a decimal number, which specifies the p-value cut-off in the peak detection algorithm, default: 1e-5
FOLD_ENRICHMENT <- 2 # a decimal number, which specifies the minimal fold enrichment in the peak calling process. default: 1 
DIFF_PEAK_ABS_FOLD_CHANGE <- 1.5
DIFF_PEAK_CUTOFF_FDR <- 0.1
MINIMAL_MAPQ <- 20

metdiff(GENE_ANNO_GTF=gtf,IP_BAM = IP_BAM,INPUT_BAM = INPUT_BAM, TREATED_IP_BAM = TREATED_IP_BAM,TREATED_INPUT_BAM=TREATED_INPUT_BAM, EXPERIMENT_NAME=EXPERIMENT_NAME, MINIMAL_MAPQ= MINIMAL_MAPQ, FRAGMENT_LENGTH=FRAGMENT_LENGTH, PEAK_CUTOFF_PVALUE=PEAK_CUTOFF_PVALUE,FOLD_ENRICHMENT=FOLD_ENRICHMENT, DIFF_PEAK_ABS_FOLD_CHANGE =DIFF_PEAK_ABS_FOLD_CHANGE, DIFF_PEAK_CUTOFF_FDR =DIFF_PEAK_CUTOFF_FDR, WINDOW_WIDTH=WINDOW_WIDTH, SLIDING_STEP=SLIDING_STEP )


print("...The End...")
