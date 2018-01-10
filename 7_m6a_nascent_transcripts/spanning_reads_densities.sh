#!/bin/bash -l

## Created on 10th of November 2017
## @author: Igor Ruiz de los Mozos

# Main script to plot densities and hetmaps of m6a NeMeRIP-seq peaks on Exons, Introns and spanning reads between both 

# Dependencies 
# Below programs must be instaled on the PATH. Check manuals for installation and usage
#
#	bedtools v2.26.0 			http://bedtools.readthedocs.io/en/latest/ 
#	samtools v0.1.19-44428cd	http://samtools.sourceforge.net/
#	deeptools 2.5.3 			https://deeptools.readthedocs.io/en/latest/index.html


########################################################################
# Data modification and 
########################################################################

# Merge bam files
bash merge_bam.sh m6a_IP_A_merge.TopHap.bam m6a_IP_A1.sorted.bam m6a_IP_A2.sorted.bam m6a_IP_A3.sorted.bam
bash merge_bam.sh m6a_IP_S_merge.TopHat.bam m6a_IP_S1.sorted.bam m6a_IP_S2.sorted.bam  m6a_IP_S3.sorted.bam 
bash merge_bam.sh m6a_input2_A_merge.TopHat.bam m6a_input2_A1.sorted.bam m6a_input2_A2.sorted.bam m6a_input2_A3.sorted.bam
bash merge_bam.sh m6a_input2_S_merge.TopHat.bam m6a_input2_S1.sorted.bam m6a_input2_S2.sorted.bam m6a_input2_S3.sorted.bam

# Calculate coverage of bam files to BigWig. 
bamCoverage --normalizeTo1x 2730870000 --bam /Volumes/lab-ulej/working/Igor/SMAD/m6a_input2_A_merge.TopHat.bam --outFileName m6a_input2_A_merge.sort.bigWig --outFileFormat bigwig
bamCoverage --normalizeTo1x 2730870000 --bam /Volumes/lab-ulej/working/Igor/SMAD/m6a_input2_S_merge.TopHat.bam --outFileName m6a_input2_S_merge.sort.bigWig --outFileFormat bigwig
bamCoverage --normalizeTo1x 2730870000 --bam /Volumes/lab-ulej/working/Igor/SMAD/m6a_IP_A_merge.TopHap.bam --outFileName m6a_IP_A_merge.sort.bigWig --outFileFormat bigwig
bamCoverage --normalizeTo1x 2730870000 --bam /Volumes/lab-ulej/working/Igor/SMAD/m6a_IP_S_merge.TopHat.bam --outFileName m6a_IP_S_merge.sort.bigWig --outFileFormat bigwig

# Insert ENSEML gene name ()
bedtools intersect -wao -a Subset3A.bed -b ENSEMBL_genes_hg38.bed > Subset3A_ann.bed
cat Subset3A_ann.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $16 "\t" $17}' > tmp.bed && mv tmp.bed Subset3A_ann.bed
sort -k1,1 -k2,2 -k3,3n -u Subset3A_ann.bed > tmp.bed && mv tmp.bed Subset3A_ann.bed
# Tranform to bed format for genome browser
cat Subset3A.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 }' > tmp.bed && mv tmp.bed Subset3A_Bed.bed

########################################################################
# Obtaining Exons with spanning reads 
########################################################################

# m6a peaks overlaping with Introns
# Remove chr string from chromosome column
sed 's/chr//g' Introns_hg38.bed > tmp.bed && mv tmp.bed Introns_hg38.bed
# Intersect with Introns (Introns genomic position downloaded from UCSC table browser GRCh38 GENCODE v24)
bedtools intersect -nonamecheck -wa -u -a Subset3A_ann.bed -b Introns_hg38.bed > Subset3A_ann_Introns.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Subset3A_ann_Introns.bed > tmp.bed && mv tmp.bed Subset3A_ann_Introns.bed
# Insert "INtrons" string on the table
cat Subset3A_ann_Introns.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" "Intron"}' > tmp.bed && mv tmp.bed Subset3A_ann_Introns.bed

# m6a peaks overlaping with Exons
# Remove chr string from chromosome column
sed 's/chr//g' Exons_hg38.bed > tmp.bed && mv tmp.bed Exons_hg38.bed
# Intersect with Exons (Exons genomic position downloaded from UCSC table browser GRCh38 GENCODE v24)
bedtools intersect -nonamecheck -v -wa -a Subset3A_ann_Introns.bed -b Exons_hg38.bed > Subset3A_ann_Introns_NoExon.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Subset3A_ann_Introns_NoExon.bed > tmp.bed && mv tmp.bed Subset3A_ann_Introns_NoExon.bed
cat Subset3A_ann_Introns_NoExon.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" "Intron"}' > tmp.bed && mv tmp.bed Subset3A_ann_Introns_NoExon.bed

# Spannig reads
bedtools intersect -nonamecheck -wa -a Subset3A_ann_Introns.bed -b Exons_hg38.bed > Subset3A_ann_spann_reads.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Subset3A_ann_spann_reads.bed > tmp.bed && mv tmp.bed Subset3A_ann_spann_reads.bed
cat Subset3A_ann_spann_reads.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" "SpannReads"}' > tmp.bed && mv tmp.bed Subset3A_ann_spann_reads.bed
# exons within spannig reads
bedtools intersect -nonamecheck -wa -a Exons_hg38.bed -b Subset3A_ann_spann_reads.bed > Subset3A_ann_spann_reads_Exons.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Subset3A_ann_spann_reads_Exons.bed > tmp.bed && mv tmp.bed Subset3A_ann_spann_reads_Exons.bed 

bedtools intersect -nonamecheck -wa -a Exons_hg38.bed -b Subset3A_ann_spann_reads_Down.bed > Subset3A_ann_spann_reads_Exons_m6Down.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Subset3A_ann_spann_reads_Exons_m6Down.bed > tmp.bed && mv tmp.bed Subset3A_ann_spann_reads_Exons_m6Down.bed


########################################################################
# Exploratory plots with Down and Up regulated m6a peaks
########################################################################

# Plot m6a peaks on Exons Introns Downregulated and Upregulated
computeMatrix scale-regions -S m6a_input2_A_merge.sort.bigWig m6a_IP_A_merge.sort.bigWig m6a_input2_S_merge.sort.bigWig m6a_IP_S_merge.sort.bigWig -R Subset3A_ann_m6a_Exons.bed Subset3A_ann_m6a_Introns.bed Downregulated.bed Upregulated.bed -out matrix_SMAD.tab.gz --missingDataAsZero --startLabel TSS --endLabel TES --skipZeros -a 5000 -b 5000 -bs 10

plotHeatmap -m matrix_matrix_SMAD_all.tab.gz -out SMAD.pdf --missingDataColor "#FFFFFF" --colorMap magma_r --samplesLabel Input_A m6a_A Input_S m6a_S --regionsLabel Exons Introns Downregulated_m6a Upregulated_m6a --startLabel start --endLabel end --plotTitle 'SMAD m6a' \
--heatmapHeight 15 --heatmapWidth 8 \
--sortUsing mean --sortUsingSamples 2 --sortRegions descend


## Only exons and introns with fold chance peaks
bedtools intersect -nonamecheck -wa -u -a Exons_hg38.bed -b Downregulated.bed > Downregulated_m6a_Exons.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Downregulated_m6a_Exons.bed > tmp.bed && mv tmp.bed Downregulated_m6a_Exons.bed

bedtools intersect -nonamecheck -wb -u -a Introns_hg38.bed -b Downregulated.bed > Downregulated_m6a_Introns.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Downregulated_m6a_Introns.bed > tmp.bed && mv tmp.bed Downregulated_m6a_Introns.bed

bedtools intersect -nonamecheck -wa -u -a Exons_hg38.bed -b Upregulated.bed > Upregulated_m6a_Exons.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Upregulated_m6a_Exons.bed > tmp.bed && mv tmp.bed Upregulated_m6a_Exons.bed

bedtools intersect -nonamecheck -wb -u -a Introns_hg38.bed -b Upregulated.bed > Upregulated_m6a_Introns.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Upregulated_m6a_Introns.bed > tmp.bed && mv tmp.bed Upregulated_m6a_Introns.bed


# Plot exons and introns with fold chance peaks Downregulated
computeMatrix scale-regions -S m6a_input2_A_merge.sort.bigWig m6a_IP_A_merge.sort.bigWig m6a_input2_S_merge.sort.bigWig m6a_IP_S_merge.sort.bigWig -R Downregulated_m6a_Exons.bed Downregulated_m6a_Introns.bed  -out matrix_SMAD_down_only.tab.gz --missingDataAsZero --startLabel Start --endLabel End --skipZeros -a 5000 -b 5000 -bs 10 

plotHeatmap -m matrix_SMAD_down_only.tab.gz -out SMAD_down_only.pdf --missingDataColor "#FFFFFF" --colorMap magma_r --samplesLabel Input_A m6a_A Input_S m6a_S --regionsLabel Down_m6a_Exons Down_m6a_Introns --startLabel Start --endLabel End --plotTitle 'SMAD m6a' \
--heatmapHeight 15 --heatmapWidth 8 \
--sortUsing mean --sortUsingSamples 2 --sortRegions descend


########################################################################
# Select m6a peaks FDR < 0.05 & FC < 2
########################################################################

cat Subset3A_ann.bed | awk -F "\t" '{if($5 <= -1.31 && $11 <= -0.58) {print}}' > DownregulatedHH2.bed
cat Subset3A_ann_spann_reads_Exons_m6Down.bed | awk -F "\t" '{if($5 <= -1.31 && $11 <= -0.58) {print}}' > Subset3A_ann_spann_reads_Down_HHFC2.bed

# Intersect high probabity peaks with Exons
bedtools intersect -nonamecheck -wa -u -a Exons_hg38.bed -b DownregulatedHH2.bed > DownregulatedHH_2_m6a_Exons.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u DownregulatedHH_2_m6a_Exons.bed > tmp.bed && mv tmp.bed DownregulatedHH_2_m6a_Exons.bed

# Intersect high probabity peaks with Exons
bedtools intersect -nonamecheck -wa -a Exons_hg38.bed -b Subset3A_ann_spann_reads_Down_HHFC2.bed > Exons_spanning_read_downHHFC_2.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Exons_spanning_read_downHHFC_2.bed > tmp.bed && mv tmp.bed Exons_spanning_read_downHHFC_2.bed

# Exclude Exons from the spanning reads Exons
bedtools intersect -nonamecheck -v -wa -a DownregulatedHH_2_m6a_Exons.bed -b Exons_spanning_read_downHHFC_2.bed > Substracted_exons_HHFC2.bed


# Plot substracted Exons and Spanning read Exons ## 
computeMatrix scale-regions -S m6a_input2_A_merge.sort.bigWig m6a_IP_A_merge.sort.bigWig m6a_input2_S_merge.sort.bigWig m6a_IP_S_merge.sort.bigWig -R Substracted_exons_HHFC2.bed Exons_spanning_read_downHHFC_2.bed -out matrix_SMAD_down_can_SpannHH2.tab.gz --missingDataAsZero --startLabel Start --endLabel End --skipZeros -a 500 -b 500 -bs 10 

plotHeatmap -m matrix_SMAD_down_can_SpannHH2.tab.gz -out SMAD_down_Can_Spann_DowHH2_500.pdf --missingDataColor "#FFFFFF" --colorMap magma_r --samplesLabel "Input Activin" "m6a Activin" "Input SB 2h" "m6a SB 2h" --regionsLabel "Downregulated Exons" "Downregulated Spanning reads Exons" --startLabel "3'ss" --endLabel "5'ss" --plotTitle 'SMAD m6a Spanning Read exon boundaries FDR < 0.05 FC < -2' \
--heatmapHeight 20 --heatmapWidth 8 --zMax 150 \
--sortUsing mean --sortUsingSamples 2 --sortRegions descend

## Join exons with amd without and spanning reads 
cat Substracted_exons_HHFC2.bed Exons_spanning_read_downHHFC_2.bed > Exons_merged.bed
sort -k1,1 -k2,2 -k3,3n -k5,5n -u Exons_merged.bed > tmp.bed && mv tmp.bed Exons_merged.bed

# Plot highly significant m6a peaks over Exons and Exons with spanning reads all together 
computeMatrix scale-regions -S m6a_input2_A_merge.sort.bigWig m6a_IP_A_merge.sort.bigWig m6a_input2_S_merge.sort.bigWig m6a_IP_S_merge.sort.bigWig -R Exons_merged.bed -out matrix_SMAD_down_can_SpannHH2_exon_merged.tab.gz --missingDataAsZero --startLabel Start --endLabel End --skipZeros -a 500 -b 500 -bs 10 

plotHeatmap -m matrix_SMAD_down_can_SpannHH2_exon_merged.tab.gz -out SMAD_down_Can_Spann_DowHH2_exon_merged_00.pdf --missingDataColor "#FFFFFF" --colorMap magma_r --samplesLabel "Input Activin" "m6a Activin" "Input SB 2h" "m6a SB 2h" --regionsLabel "All Downregulated Exons" --startLabel "3'ss" --endLabel "5'ss" --plotTitle 'SMAD m6a Spanning Read exon boundaries FDR < 0.05 FC < -2' \
--heatmapHeight 20 --heatmapWidth 8 --zMax 150 \
--sortUsing mean --sortUsingSamples 2 --sortRegions descend

# Final figure
mv SMAD_down_Can_Spann_DowHH2_exon_merged_00.pdf SMAD_2-3_figS6g_h.pdf
