library(AnnotationHub)
library(AnnotationDbi)
library(intervals)
library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)
library(devtools)
library(bedr)
library(scales)
source("asynt.R")
library(splitstackshape)
#install.packages("rstatix")
library(rstatix)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("multcomp")
library(multcomp)
#install.packages("dgof")
library(dgof)
#install.packages("Matching")
library(Matching)

############ Scurria scurra ################
#read in gtf file with pkg GenomicFeatures
#Scurra_txdb_parsed <- makeTxDbFromGFF(file="/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/01-SG_genome/assembly_annotation_v1/Scurria_scurra_annotation_v1_ch10_top10.gtf", format="gtf")
#saveDb(Scurra_txdb_parsed, "Scurria_scurra_parsed")
S.scurra.p <- loadDb("Scurria_scurra_parsed")
columns(S.scurra.p)

#create a table of gene and transcript IDs
txdf <- AnnotationDbi::select(S.scurra.p,
                              keys=keys(S.scurra.p, "GENEID"),
                              columns=c("GENEID", "TXID", "TXCHROM"),
                              keytype="GENEID")
head(txdf, 20)

#collect exons by transcript id
exons.list.per.transcript <- exonsBy(S.scurra.p, by="tx", use.names=TRUE)
mcols(exons.list.per.transcript) #get information about data collected

#coerce to dataframe
exons.list.per.transcript.df.ss <- as.data.frame(exons.list.per.transcript)
head(exons.list.per.transcript.df.ss, n=50)
colnames(exons.list.per.transcript.df.ss) <- c("group", "GENEID","CHRMID","start","end","width","strand","tx_id","tx_name")
head(exons.list.per.transcript.df.ss)

# create vector with species names
spID <- as.factor("sscurra")
exons.list.per.transcript.df.ss <- cbind(exons.list.per.transcript.df.ss, spID)

############ Scurria viridula ################
#read in gtf file with pkg GenomicFeatures
#Viridula_txdb_parsed <- makeTxDbFromGFF(file="/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/02-VG_genome/assembly_annotation_v1/Scurria_viridula_annotation_v1_ch10_top10.gtf", format="gtf")
#saveDb(Viridula_txdb_parsed, "Scurria_viridula_parsed")
S.viridula.p <- loadDb("Scurria_viridula_parsed")
columns(S.viridula.p)

#create a table of gene and transcript IDs
txdf <- AnnotationDbi::select(S.viridula.p,
                              keys=keys(S.viridula.p, "GENEID"),
                              columns=c("GENEID", "TXID", "TXCHROM"),
                              keytype="GENEID")
head(txdf, 20)

#collect exons by transcript id
exons.list.per.transcript <- exonsBy(S.viridula.p, by="tx")
mcols(exons.list.per.transcript) #get information about data collected

#coerce to dataframe
exons.list.per.transcript.df.sv <- as.data.frame(exons.list.per.transcript)
head(exons.list.per.transcript.df.sv, n=50)
colnames(exons.list.per.transcript.df.sv) <- c("group", "GENEID","CHRMID","start","end","width","strand","tx_id","tx_name")
head(exons.list.per.transcript.df.sv)

# create vector with species names
spID <- as.factor("sviridula")
exons.list.per.transcript.df.sv <- cbind(exons.list.per.transcript.df.sv, spID)

############ Scurria zebrina ################
#read in gtf file with pkg GenomicFeatures
#Zebrina_txdb_parsed <- makeTxDbFromGFF(file="/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/03-ZG_genome/assembly_annotation_v1/Scurria_zebrina_annotation_v1_ch10_top10.gtf", format="gtf")
#saveDb(Zebrina_txdb_parsed, "Scurria_zebrina_parsed")
S.zebrina.p <- loadDb("Scurria_zebrina_parsed")
columns(S.zebrina.p)

#create a table of gene and transcript IDs
txdf <- AnnotationDbi::select(S.zebrina.p,
                              keys=keys(S.zebrina.p, "GENEID"),
                              columns=c("GENEID", "TXID", "TXCHROM"),
                              keytype="GENEID")
head(txdf, 20)

#collect exons by transcript id
exons.list.per.transcript <- exonsBy(S.zebrina.p, by="tx")
mcols(exons.list.per.transcript) #get information about data collected

#coerce to dataframe
exons.list.per.transcript.df.sz <- as.data.frame(exons.list.per.transcript)
head(exons.list.per.transcript.df.sz, n=50)
colnames(exons.list.per.transcript.df.sz) <- c("group", "GENEID","CHRMID","start","end","width","strand","tx_id","tx_name")
head(exons.list.per.transcript.df.sz)

# create vector with species names
spID <- as.factor("szebrina")
exons.list.per.transcript.df.sz <- cbind(exons.list.per.transcript.df.sz, spID)

#######Merge dataframes ###
merged_df <- rbind(exons.list.per.transcript.df.ss, 
                   exons.list.per.transcript.df.sv, 
                   exons.list.per.transcript.df.sz)

#Add column with Chromosome number
merged_df$CHRMNUM <- substring(merged_df$CHRMID, 3)

#Add column with Orthologous group distinction based on MCSCAN collinearity plot
merged_df$ORTHGRP <- merged_df$CHRMNUM

merged_df$ORTHGRP <- ifelse(merged_df$CHRMID == "sv3", 4, 
                            ifelse(merged_df$CHRMID == "sv4", 3,
                                   ifelse(merged_df$CHRMID == "sz6", 7,
                                          ifelse(merged_df$CHRMID == "sz7", 6,
                                                 merged_df$CHRMNUM))) )

head(merged_df)
tail(merged_df, n=100)

str(merged_df)
summary(merged_df)

merged_df %>%
  group_by(ORTHGRP) %>%
  get_summary_stats(width, type = "common")

boxplot(width ~ ORTHGRP,
        data = merged_df)
#ANOVA
res_aov <- aov(length ~ spID,
               data = merged_df)
summary(res_aov)
#Post hoc pairwise tests
post_test <- glht(res_aov,
                  linfct = mcp(spID = "Tukey"))
summary(post_test)

#Kruskal Wallis tests
kruskal.test(width ~ spID,
             data = merged_df)

pairwise.wilcox.test(merged_df$width, merged_df$spID,
                     p.adjust.method = "fdr")

####PLOT ####
ggplot(merged_df) +
  aes(x = spID, y = width, color= spID) +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values=c("#FB61D7", "#A58AFF", "#00B6EB"))+
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="black", pch=16, size=.7) +
  stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="Exon Length (bp) ") + 
  theme(legend.title=element_blank(),
        axis.text.x=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))


#####Merge dataframes
merged_df <- rbind(exonic.transcript.sizes.ss, 
                   exonic.transcript.sizes.sv, 
                   exonic.transcript.sizes.sz)
head(merged_df)

#####Summary of data
str(merged_df)
summary(merged_df)

merged_df %>%
  group_by(spID) %>%
  get_summary_stats(length, type = "common")

boxplot(length ~ spID,
        data = merged_df)

######ANOVA
res_aov <- aov(length ~ spID,
               data = merged_df)
summary(res_aov)

post_test <- glht(res_aov,
                  linfct = mcp(spID = "Tukey"))
summary(post_test)

#######Kruskal Wallac Test and Pairwise Wilcox Test
kruskal.test(length ~ spID,
             data = merged_df)

pairwise.wilcox.test(merged_df$length, merged_df$spID,
                     p.adjust.method = "fdr")

########Kolmogorov-Smirnov Test
ks1<- ks.test(exonic.transcript.sizes.ss$length, 
        exonic.transcript.sizes.sv$length, alternative = c('two.sided', 'less', 'greater'), exact= NULL, tol= 1e-8, 
        simulate.p.value = FALSE, B=2000)

ks2<- ks.test(exonic.transcript.sizes.ss$length, 
              exonic.transcript.sizes.sv$length)

ks3<- ks.test(exonic.transcript.sizes.sv$length, 
              exonic.transcript.sizes.sz$length)

#plotting the results of the ks tests
plot(ecdf(exonic.transcript.sizes.ss$length),
     xlim = range(c(exonic.transcript.sizes.ss$length, exonic.transcript.sizes.sv$length)),
     col = "blue")
plot(ecdf(exonic.transcript.sizes.sv$length),
     add = TRUE,
     lty = "dashed",
     col = "red")

####PLOT of all lengths####
ggplot(merged_df) +
  aes(x = spID, y = length, color = spID) +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values=c("#FB61D7", "#A58AFF", "#00B6EB"))+
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="black", pch=16, size=.7) +
  stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="Exon Lengths (bp)") + 
  theme(legend.title=element_blank(),
        axis.text.x=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))
