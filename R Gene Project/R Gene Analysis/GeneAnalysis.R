# Technical Assessment
# Author: Ryan Ellis
# Purpose: To demonstrate my technical skills and capabilities in the areas of data processing,
#          R programming and data visualization / reporting

library(nanostringr)
library(tidyverse)
library(gplots)
library(edgeR)
setwd("C:/Users/raemu/Documents/Technical Code")

# Import RCC files
samples<-read_rcc(path = "C:/Users/raemu/Documents/Technical Code/RCC")

# Rename 
names(samples$raw)[names(samples$raw) == "GSM2055823_01_4353_PD_mRNA"] <- "Subject 1 Base"
names(samples$raw)[names(samples$raw) == "GSM2055824_02_4355_PD_mRNA"] <- "Subject 1 Post"
names(samples$raw)[names(samples$raw) == "GSM2055825_03_3366_PD_mRNA"] <- "Subject 2 Base"
names(samples$raw)[names(samples$raw) == "GSM2055826_04_4078_PD_mRNA"] <- "Subject 2 Post"
names(samples$raw)[names(samples$raw) == "GSM2055827_05_4846_PD_mRNA"] <- "Subject 3 Base"
names(samples$raw)[names(samples$raw) == "GSM2055828_06_3746_PD_mRNA"] <- "Subject 3 Post"
names(samples$raw)[names(samples$raw) == "GSM2055829_07_3760_PD_mRNA"] <- "Subject 4 Base"
names(samples$raw)[names(samples$raw) == "GSM2055830_08_3790_PD_mRNA"] <- "Subject 4 Post"
names(samples$raw)[names(samples$raw) == "GSM2055831_09_4436_PD_mRNA"] <- "Subject 5 Base"
names(samples$raw)[names(samples$raw) == "GSM2055832_10_4050_PD_mRNA"] <- "Subject 5 Post"

# Get just the data
data<-samples$raw
# Set up matrix
ndata<-as.matrix(subset(data, select = -c(Name,Accession) )[,-1])
rownames(ndata)<-data$Name

# Normalize data
edger<-calcNormFactors(ndata, method = "TMM")
libSize<-effectiveLibSizes(ndata,log=FALSE)
nb<-libSize*edger
allNormData<-ndata/nb
allNormData<-t(allNormData)

# Filter Data
PNdata<-filter(data, Code.Class=="Positive" | Code.Class=="Negative")
NormPN<- subset(allNormData, select = PNdata$Name)
# Heat Map 
# Make heat map with Normalized values and clustering 
hr <- hclust(as.dist(1-cor(t(NormPN), method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(NormPN, method="kendall")), method="complete")
jpeg('heatmap.jpg',width=866, height=634, pointsize=10,
     type="windows", antialias="cleartype")
heatmap.2(NormPN, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), scale="col", density.info="none",margins=c(10,8), trace="none")
dev.off()

# Box Plot Analysis 
# Filter Genes 
Name<-c("MCL1","CXCL1")
BPdata<-subset(allNormData, select=colnames(allNormData)%in%Name)
MCL1<-subset(BPdata, select=colnames(BPdata)%in%c("MCL1"))
CXCL1<-subset(BPdata, select=colnames(BPdata)%in%c("CXCL1"))
MCL1<-t(MCL1)
CXCL1<-t(CXCL1)
MCL1<-as.data.frame(MCL1)
CXCL1<-as.data.frame(CXCL1)

# MCL1
MCL1Base<-MCL1%>%select(contains("Base"))
MCL1Post<-MCL1%>%select(contains("Post"))
BoxMCL1<-MCL1Base
BoxMCL1[nrow(BoxMCL1) + 1,]=MCL1Post
rownames(BoxMCL1)<-c("Base","Post")
colnames(BoxMCL1)<-gsub("Base","",colnames(BoxMCL1))
BoxMCL1<-t(BoxMCL1)
# CXCL1
CXCL1Base<-CXCL1%>%select(contains("Base"))
CXCL1Post<-CXCL1%>%select(contains("Post"))
BoxCXCL1<-CXCL1Base
BoxCXCL1[nrow(BoxCXCL1) + 1,]=CXCL1Post
rownames(BoxCXCL1)<-c("Base","Post")
colnames(BoxCXCL1)<-gsub("Base","",colnames(BoxCXCL1))
BoxCXCL1<-t(BoxCXCL1)

# Box Plots

#MCL1
jpeg('BoxPlot.jpg',width=733, height=686, pointsize=10,
     type="windows", antialias="cleartype")
par(mfrow=c(2,1))
MCL1Means<-t(as.data.frame(colMeans(BoxMCL1)))
x<-list('Base'=MCL1Base,'Post'=MCL1Post)
boxplot(BoxMCL1,
        horizontal = TRUE,
        main = "MCL1 Treatment Comparision",
        xlab = "Counts",
        ylab = "Treatment",
        xaxt='n')
points(x = MCL1Means,                             
       y = 1:2,
       col = "red",
       pch = 16,)
stripchart(x, vertical = FALSE, method = "jitter",
           pch = 19, add = TRUE,col=1)
axis(1, at = seq(0, 1, .005))


# CXCL1
# jpeg('CXCL1.jpg',width=866, height=634, pointsize=10,
#      type="windows", antialias="cleartype")
CXCL1Means<-t(as.data.frame(colMeans(BoxCXCL1)))
x<-list('CXCL1 Base'=CXCL1Base,'CXCL1 Post'=CXCL1Post)
boxplot(BoxCXCL1,
        horizontal = TRUE,
        main = "CXCL1 Treatment Comparision",
        xlab = "Counts",
        ylab = "Treatment",
        xaxt='n')
points(x = CXCL1Means,                             
       y = 1:2,
       col = "red",
       pch = 16)
stripchart(x, vertical = FALSE, method = "jitter",
           pch = 19, add = TRUE,col=1)
axis(1, at = seq(0, 1, .00025))
dev.off()

# CSV For Report
csvs<-list(data,ndata,allNormData,PNdata,NormPN,BPdata,MCL1,CXCL1,MCL1Base,MCL1Post,BoxMCL1,CXCL1Base,CXCL1Post,BoxCXCL1)
csvsNames<-c('data','ndata','allNormData','PNdata','NormPN','BPdata','MCL1','CXCL1','MCL1Base','MCL1Post','BoxMCL1','CXCL1Base','CXCL1Post','BoxCXCL1')
for (cv in 1:14) {
        write.csv(csvs[cv],paste(paste(as.character(cv),csvsNames[cv], sep='-'),'.csv',sep=''), row.names = TRUE)   
}
