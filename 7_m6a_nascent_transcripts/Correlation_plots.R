library("ggplot2")
library("smoother")
library("ggrepel")


####################################################################################
# Coorelation model function
####################################################################################

# Given a model, predict values of yvar from xvar
# This supports one predictor and one predicted variable.
# xrange: If NULL, determine the x range from the model object. If a vector with
#   two numbers, use those as the min and max of the prediction range.
# samples: Number of samples across the x range.
# ...: Further arguments to be passed to predict()
predictvals <- function(model, xvar, yvar, xrange=NULL, samples=100, ...) {
  
  # If xrange isn't passed in, determine xrange from the models.
  # Different ways of extracting the x range, depending on model type
  if (is.null(xrange)) {
    if (any(class(model) %in% c("lm", "glm")))
      xrange <- range(model$model[[xvar]])
    else if (any(class(model) %in% "loess"))
      xrange <- range(model$x)
  }
  
  newdata <- data.frame(x = seq(xrange[1], xrange[2], length.out = samples))
  names(newdata) <- xvar
  newdata[[yvar]] <- predict(model, newdata = newdata, ...)
  newdata
}


####################################################################################
# Features per ENSEMBL ID correlation min exon peak with mean of exon peaks 
####################################################################################

point1table <-read.table("./Subset3A_ann_ChrisPoint1.bed", sep="\t", header = TRUE, stringsAsFactors=TRUE)

head(point1table)
# Scatter plot
ggplot(point1table, aes(x=Mean, y=diff.log2.fc)) + geom_point()+
  geom_smooth(method=lm, se=FALSE)+
  labs(title="Correlation exonic peaks")

ggsave("figureS6_j.pdf", width = 20, height = 20, units = "cm")

model = lm(diff.log2.fc ~ Mean, point1table)
pred = predictvals(model, 'Mean', 'diff.log2.fc')

summary(model)



####################################################################################
# Features per ENSEMBL ID correlation min exon peak with mean of intron peaks 
####################################################################################

point2table <-read.table("./Subset3A_ann_ChrisPoint2.bed", sep="\t", header = TRUE, stringsAsFactors=TRUE)

head(point2table)
# Scatter plot
ggplot(point2table, aes(x=Mean, y=diff.log2.fc)) + geom_point()+
  geom_smooth(method=lm, se=FALSE)+
  labs(title="Correlation exonic peaks with intronic")


####################################################################################
# Features per ENSEMBL ID correlation min exon peak with mean everything else 
####################################################################################

point3table <-read.table("./Subset3A_ann_ChrisPoint3.bed", sep="\t", header = TRUE, stringsAsFactors=TRUE)

head(point3table)

model = lm(diff.log2.fc ~ Mean, point3table)
pred = predictvals(model, 'Mean', 'diff.log2.fc')

summary(model)


cr <- ggplot(point3table, aes(x=Mean, y=diff.log2.fc)) + geom_point()+
  theme_classic()+
  #geom_smooth(method=lm, se=FALSE, size=.5)+
  labs(title="Conservation of peaks variation within transcripts")+
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain"))

cr 


####################################################################################
# Features per ENSEMBL ID find position of min exon peak and distance to any other peak 
####################################################################################

# All point with exons and introns (RC_variations_all)

point4table <-read.table("./Subset3A_ann_ChrisPoint4-2_IG.txt", sep="\t", header = TRUE, stringsAsFactors=TRUE)
head(point4table)


sp_dist <- ggplot(point4table, aes(x=FC, y=diff.log2.fc, color=Functional_Location)) + geom_point(aes( fill=Functional_Location, size=Distance), alpha = 0.5)+
  theme_bw()+
  geom_smooth(method=lm, se=FALSE)+
  labs(title="Correlation exonic peaks with all peaks taking in acount distance to the most variable peak") +
  scale_size_continuous(range = c(0.1,2))+
  scale_color_manual(values=c('#999999','#E69F00','black', '#56B4E9'))+
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain"))+
  xlab("log2 m6A peak fold-cahnge")+
  ylab("Strongest -log2 m6A peak fold-change")+
  geom_text(aes(y=diff.log2.fc,x=FC-.2, label=Gene), size=2, vjust=0)

sp_dist



####################################################################################
# Features per ENSEMBL ID find position of min exon peak and distance to any other peak 
# Subset of transcripts that contain peaks falling on introns
####################################################################################

point5table <-read.table("./Subset3A_ann_ChrisPointOnlyIntons.txt", sep="\t", header = TRUE, stringsAsFactors=TRUE)

sp_introns <- ggplot(point5table, aes(x=FC, y=diff.log2.fc, color=Functional_Location)) + geom_point(aes( fill=Functional_Location, size=Distance), alpha = 0.5)+
  theme_bw()+
  labs(title="Correlation exonic peaks with intronic peaks taking in acount distance to the most variable peak") +
  scale_size_continuous(range = c(0.2,2))+
  scale_color_manual(values=c('white','gray48',"black", 'red'))+
  geom_text(aes(y=diff.log2.fc,x=FC-.2, label=Gene), size=2, vjust=0)

sp_introns


# Detail area
# This plot shows unequivocally how peaks on introns vary at the same level than the min exon peak and the rest inside exons. The distance from the min peak also account for the variation but even long distances (>10Kb) peaks on introns vary on the same direction and with the similar fold change.

sp_introns_detail <- ggplot(point5table, aes(x=FC, y=diff.log2.fc, color=Functional_Location)) + geom_point(aes( fill=Functional_Location, size=Distance), alpha = 0.5)+
  theme_classic()+
  labs(title="Correlation exonic peaks with intronic peaks taking in acount distance to the most variable peak") +
  scale_size_continuous(range = c(1,5))+
  scale_color_manual(values=c('white','gray48',"black", 'red'))+
  scale_y_continuous(limits = c(-1.6, -0.7)) +
  scale_x_continuous(limits = c(-2, 1.3)) +
  theme(text=element_text(size=6),axis.text=element_text(size=6), axis.title=element_text(size=6,face="plain"))+
  xlab("log2 m6A peak fold-cahnge")+
  ylab("Strongest -log2 m6A peak fold-change")+
  geom_text_repel(aes(y=diff.log2.fc,x=FC-0.3, label=Gene), size=3 , segment.size = 0, ) 

sp_introns_detail


