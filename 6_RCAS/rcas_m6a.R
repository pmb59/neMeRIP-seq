# Pipeline RCAS


library(RCAS)
library(knitr)
library(rmarkdown)
library(htmlwidgets)

# set working directory
wd <- "../RCAS2017rev/data" 
setwd( wd )

REGIONS <- "m6A_new_up.bed" 

# for i in REGIONS

######################################################################
# 2. Data input
######################################################################
# 2.2 Importing custom data

queryRegions <- importBed(filePath = paste(wd, REGIONS, sep="/"), sampleN = 0, keepStandardChr = TRUE) 
queryRegions

gff <- importGtf(filePath = paste(wd, "Homo_sapiens.GRCh38.83.gtf", sep="/"), keepStandardChr = TRUE)

######################################################################
# 3. Summarizing Overlaps of Query Regions with Genomic Annotation Features 
######################################################################
# 3.1 Querying the annotation file
overlaps <- queryGff(queryRegions = queryRegions, gffData = gff)
#data.table is used to do quick summary operations
overlaps.dt <- data.table(as.data.frame(overlaps)) 

#3.1.1 Finding targeted gene types
#To find out the distribution of the query regions across gene types:

biotype_col <- grep('gene_biotype', colnames(overlaps.dt), value = T)

df <- overlaps.dt[,length(unique(overlappingQuery)), by = biotype_col]
colnames(df) <- c("feature", "count")
df$percent <- round(df$count / length(queryRegions) * 100, 1)
df <- df[order(count, decreasing = TRUE)]


#Recent versions of the R package include an export() function, which does image export locally, but requires the webshot package:
#if (!require("webshot")) install.packages("webshot")
#library(webshot)

p <- plot_ly(data = df, 
             type = "bar",
             x = df$feature,
             y = df$percent,
             text = paste("count:", df$count), color=df$feature, mode = "markers")  #pedro added: "mode = "markers"
layout(p = p, 
       margin = list(l=100, r=100, b=150), 
       xaxis = list(showticklabels = TRUE,  tickangle = 90), 
       yaxis = list(title = paste("percentage of query regions,", 
                                  "n =", length(queryRegions))))


#htmlwidgets::saveWidget(as.widget(p), "figure_3.1.1.html")
saveWidget(widget=as.widget(p), file="figure_3.1.1.html", selfcontained = FALSE  )  #, knitrOptions = list() 

rm(p)



#3.2 Extending the annotation feature space

txdbFeatures <- getTxdbFeaturesFromGRanges(gff)

# 3.2.1 Plotting overlap counts between query regions and gene features
summary <- summarizeQueryRegions(queryRegions = queryRegions, 
                                 txdbFeatures = txdbFeatures)

df <- data.frame(summary)
df$percent <- round((df$count / length(queryRegions)), 3) * 100



p <- plot_ly( data = df, 
              x = rownames(df), 
              y = df$percent, 
              type = 'bar',
              text = paste("count:", df$count), 
              color = rownames(df)
              )
layout(p = p, 
       xaxis = list(title = 'features'),
       yaxis = list(title = paste("percentage of query regions,", 
                                  "n =", length(queryRegions)
                                  )
                    ), 
       margin = list(b = 150, r = 50)
       )

htmlwidgets::saveWidget(as.widget(p), "figure_3.2.1.html", selfcontained = FALSE )



# 3.2.2 Obtaining a table of overlap counts between query regions and genes
dt <- getTargetedGenesTable(queryRegions = queryRegions, 
                           txdbFeatures = txdbFeatures)
dt <- dt[order(transcripts, decreasing = TRUE)]

#save table:
write.table(x=dt, file = "Table_3.2.2.txt", append = FALSE, quote = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)


#3.2.3 Profiling the coverage of query regions across gene features
#3.2.3.1 Coverage profile of query regions for 3’ UTRs
#cov <- calculateCoverageProfile(queryRegions = queryRegions, 
#                               targetRegions = txdbFeatures$threeUTRs, 
#                               sampleN = 0)
cov <- calculateCoverageProfile(queryRegions = queryRegions, 
                               targetRegions = txdbFeatures$threeUTRs, 
                               sampleN = 0)

rm(p)

p <- plot_ly(data = cov, x = ~bins, y = ~coverage, type = 'scatter', mode = 'lines')


htmlwidgets::saveWidget(as.widget(p), "figure_3.2.3.1.html", selfcontained = FALSE )

rm(p)

#3.2.3.2 Coverage profile of query regions for all gene features
covList <- calculateCoverageProfileList(queryRegions = queryRegions, 
                                       targetRegionsList = txdbFeatures, 
                                       sampleN = 0)

df <- do.call('cbind', covList)
df <- df[,!grepl(colnames(df), pattern = '*.bins')]
df$bins <- c(1:100)
colnames(df) <- gsub(pattern = ".coverage", replacement = "", x = colnames(df))
mdf <- reshape2::melt(df, id.vars = c('bins'))
colnames(mdf) <- c('bins', 'feature', 'coverage')
p = plot_ly(data = mdf, x = ~bins, y = ~coverage,  type = 'scatter', mode = 'lines', color = ~feature)
layout(p)

htmlwidgets::saveWidget(as.widget(p), "figure_3.2.3.2.html", selfcontained = FALSE )


#3.2.3.3 Coverage profile of query regions at/around Feature Boundaries
cvg <- getFeatureBoundaryCoverage(queryRegions = queryRegions,
                                  featureCoords = txdbFeatures$transcripts,
                                  flankSize = 2500, 
                                  sampleN = 0)

yLimit <- (as.integer(max(c(cvg$fivePrime, cvg$threePrime))/10)+1)*10

p <- subplot(
  plot_ly(data = cvg, x = ~bases, y = ~fivePrime, type = 'scatter', mode = 'lines'),
  plot_ly(data = cvg, x = ~bases, y = ~threePrime, type = 'scatter', mode = 'lines'),
  margin = 0.05
) %>% layout (xaxis = list(title = 'Distance (bp) to TSS'), 
              xaxis2 = list(title = 'Distance (bp) to TES'), 
              yaxis = list(title = 'coverage', range = c(0, yLimit)),
              yaxis2 = list(title = 'coverage', range = c(0, yLimit)),
                           showlegend = FALSE) 
layout(p) 
htmlwidgets::saveWidget(as.widget(p), "figure_3.2.3.3.html", selfcontained = FALSE )


#4 Motif Analysis using motifRG
#4.1 Calculating enriched motifs  (hg38 not available!!!)

#motifResults <- runMotifRG(queryRegions = queryRegions, 
#                           genomeVersion = 'hg19', 
#                           motifN = 2, nCores = 2)
#
#
#par(mfrow = c(1,2), mar = c(2,2,2,2))
#for (i in 1:length(motifResults$motifs)) {
#  motifPattern <- motifResults$motifs[[i]]@pattern
#  motifRG::plotMotif(match = motifResults$motifs[[i]]@match$pattern, 
#                     main = paste0('Motif-',i,': ',motifPattern),
#                     entropy = TRUE)
#}



#4.2 motif analysis: getting motif summary statistics

#hg38 not available!!!

#5 GO term analysis
#5.1 Biological processes enriched among targeted genes
#get all genes from the GTF data
backgroundGenes <- unique(gff$gene_id)
#get genes that overlap query regions
targetedGenes <- unique(overlaps$gene_id)

#run TopGO
goBP <- runTopGO(ontology = 'BP', 
                      species = 'human', 
                      backgroundGenes = backgroundGenes, 
                      targetedGenes = targetedGenes)

goBP <- goBP[order(goBP$foldEnrichment, decreasing = TRUE),]
rownames(goBP) <- goBP$GO.ID
goBP <- subset(goBP, select = -c(Annotated,classicFisher, bh, GO.ID))

#save table:
write.table(x=goBP, file = "Table_5.1.txt", append = FALSE, quote = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)




