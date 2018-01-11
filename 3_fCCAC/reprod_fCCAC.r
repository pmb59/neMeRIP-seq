# 
# Reproducibility assessment of neMeRIP-seq with Bioconductor package fCCAC 
#
# Madrigal P (2017). “fCCAC: functional canonical correlation analysis to 
# evaluate covariance between nucleic acid sequencing datasets.” Bioinformatics
#
# https://bioconductor.org/packages/3.6/bioc/html/fCCAC.html
#
######################################################################
#Use 24,877 peaks #
library("fCCAC")
######################################################################
# m6A Ip
######################################################################
peaks <- "subset_fccac.bed"    
bigwigs <- c("m6a_IP_A1.normalised.bw", "m6a_IP_A2.normalised.bw", "m6a_IP_A3.normalised.bw" , "m6a_IP_S1.normalised.bw",  "m6a_IP_S2.normalised.bw",  "m6a_IP_S3.normalised.bw" )  
labels <- c( "m6A Activin", "m6A Activin","m6A Activin"  , "m6A SB", "m6A SB","m6A SB"  )  
x <- fCCAC(main="m6A MeRIP-seq peaks", peaks=peaks, bigwigs=bigwigs, labels=labels ,  splines=15, nbins=100, ncan=2)   #ncan < splines
#head(x)

#heatmap
pdf("reproducibility_F_values_heatmap_IP_24877peaks.pdf", height=5, width=5.5)
heatmapfCCAC(x)
dev.off()

