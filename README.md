# Exon_Number

#### Script for generating a plot of number of exons per transcript for multiple genome comparisons
This script was generated for recovering exons from genomes of three different species of limpets. Parts of this script were retrieved from various parts of the internet. Reminder: What’s the distinction between genes and transcripts? A gene is generally defined as a region in the genome which is transcribed. However, the specific combination of exons is referred to as a transcript. The different transcripts that a gene produces can be referred to as isoforms of the gene. To give an example, suppose we have a gene with exons E1, E2, E3, E4. And there are two transcripts, or isoforms, of this gene: one which includes E3, and one which excludes E3. 

```source("asynt.R")```

```install.packages("rstatix")```

```install.packages("tidyverse")```

```install.packages("multcomp")```

```install.packages("dgof")```

```install.packages("Matching")```

```library(AnnotationHub)```

```library(AnnotationDbi)```

```library(intervals)```

```library(ggplot2)```

```library(dplyr)```

```library(data.table)```

```library(GenomicFeatures)```

```library(GenomicRanges)```

```library(rtracklayer)```

```library(devtools)```

```library(bedr)```

```library(scales)```

```library(splitstackshape)```

```library(rstatix)```

```library(tidyverse)```

```library(multcomp)```

```library(dgof)```

```library(Matching)```

### Scurria scurra 
read in gtf file with pkg GenomicFeatures

```Scurra_txdb_parsed <- makeTxDbFromGFF(file="{PWD}/Scurria_scurra_annotation_v1_ch10_top10.gtf", format="gtf")```

```saveDb(Scurra_txdb_parsed, "Scurria_scurra_parsed")```

```S.scurra.p <- loadDb("Scurria_scurra_parsed")```

```columns(S.scurra.p)```

create a table of gene and transcript IDs

```txdf <- AnnotationDbi::select(S.scurra.p,```
                              ```keys=keys(S.scurra.p, "GENEID"),```
                              ```columns=c("CDSNAME", "GENEID", "TXID", "TXCHROM","TXNAME","EXONID", "EXONNAME"),```
                              ```keytype="GENEID")```

```head(txdf, 20)```

collect exons by transcript id

```exons.list.per.transcript <- exonsBy(S.scurra.p, by="tx", use.names=TRUE)```

```mcols(exons.list.per.transcript) #get information about data collected```

```head(exons.list.per.transcript)```

coerce to dataframe

exons.list.per.transcript.df.ss.all <- as.data.frame(exons.list.per.transcript)
head(exons.list.per.transcript.df.ss.all, n=50)
colnames(exons.list.per.transcript.df.ss.all) <- c("GENEID", "TXID","CHRMID","start","end","width","strand","exon_id","exon_name", "exon_rank")
head(exons.list.per.transcript.df.ss.all)

#determine number of exons per transcript
mcols(exons.list.per.transcript)$num_exons <- lengths(exons.list.per.transcript)
exons_per_transcript.ss <- as.data.frame(mcols(exons.list.per.transcript))
exons_per_transcript.ss$TXID <- row.names(exons_per_transcript.ss)
head(exons_per_transcript.ss)
class(exons_per_transcript.ss)

#merge based on TXID
merged.ss <- dplyr::inner_join(exons_per_transcript.ss, exons.list.per.transcript.df.ss.all, by = "TXID")
head(merged.ss)

#subset only the columns that are needed
merged.ss.reduced <- subset(merged.ss, select=c(num_exons, TXID, GENEID, CHRMID))
head(merged.ss.reduced)

# get only distinct rows based on TXID
merged.ss.final <- merged.ss.reduced %>% 
  distinct(TXID, .keep_all = T)
head(merged.ss.final)

# create vector with species names
spID <- as.factor("sscurra")
merged.ss.final <- cbind(merged.ss.final, spID)

############ Scurria viridula ################
#read in gtf file with pkg GenomicFeatures
Viridula_txdb_parsed <- makeTxDbFromGFF(file="/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/02-VG_genome/assembly_annotation_v1/Scurria_viridula_annotation_v1_ch10_top10.gtf", format="gtf")
saveDb(Viridula_txdb_parsed, "Scurria_viridula_parsed")
S.viridula.p <- loadDb("Scurria_viridula_parsed")
columns(S.viridula.p)

#create a table of gene and transcript IDs
txdf <- AnnotationDbi::select(S.viridula.p,
                              keys=keys(S.viridula.p, "GENEID"),
                              columns=c("GENEID", "TXID", "TXCHROM","TXNAME","EXONID", "EXONNAME"),
                              keytype="GENEID")
head(txdf, 20)

#collect exons by transcript id
exons.list.per.transcript <- exonsBy(S.viridula.p, by="tx", use.names=TRUE)
mcols(exons.list.per.transcript) #get information about data collected
head(exons.list.per.transcript)

#coerce to dataframe
exons.list.per.transcript.df.sv.all <- as.data.frame(exons.list.per.transcript)
head(exons.list.per.transcript.df.sv.all, n=50)
colnames(exons.list.per.transcript.df.sv.all) <- c("GENEID", "TXID","CHRMID","start","end","width","strand","exon_id","exon_name", "exon_rank")
head(exons.list.per.transcript.df.sv.all)

#determine number of exons per transcript
mcols(exons.list.per.transcript)$num_exons <- lengths(exons.list.per.transcript)
exons_per_transcript.sv <- as.data.frame(mcols(exons.list.per.transcript))
exons_per_transcript.sv$TXID <- row.names(exons_per_transcript.sv)
head(exons_per_transcript.sv)
class(exons_per_transcript.sv)

#merge based on TXID
merged.sv <- dplyr::inner_join(exons_per_transcript.sv, exons.list.per.transcript.df.sv.all, by = "TXID")
head(merged.sv)

#subset only the columns that are needed
merged.sv.reduced <- subset(merged.sv, select=c(num_exons, TXID, GENEID, CHRMID))
head(merged.sv.reduced)

# get only distinct rows based on TXID
merged.sv.final <- merged.sv.reduced %>% 
  distinct(TXID, .keep_all = T)
head(merged.sv.final)

# create vector with species names
spID <- as.factor("sviridula")
merged.sv.final <- cbind(merged.sv.final, spID)

############ Scurria zebrina ################
#read in gtf file with pkg GenomicFeatures
Zebrina_txdb_parsed <- makeTxDbFromGFF(file="/Users/emily/Dropbox/School/Thesis/Genomics-Ch1/03-ZG_genome/assembly_annotation_v1/Scurria_zebrina_annotation_v1_ch10_top10.gtf", format="gtf")
saveDb(Zebrina_txdb_parsed, "Scurria_zebrina_parsed")
S.zebrina.p <- loadDb("Scurria_zebrina_parsed")
columns(S.zebrina.p)

#create a table of gene and transcript IDs
txdf <- AnnotationDbi::select(S.zebrina.p,
                              keys=keys(S.zebrina.p, "GENEID"),
                              columns=c("GENEID", "TXID", "TXCHROM","TXNAME","EXONID", "EXONNAME"),
                              keytype="GENEID")
head(txdf, 20)

#collect exons by transcript id
exons.list.per.transcript <- exonsBy(S.zebrina.p, by="tx", use.names=TRUE)
mcols(exons.list.per.transcript) #get information about data collected
head(exons.list.per.transcript)

#coerce to dataframe
exons.list.per.transcript.df.sz.all <- as.data.frame(exons.list.per.transcript)
head(exons.list.per.transcript.df.sz.all, n=50)
colnames(exons.list.per.transcript.df.sz.all) <- c("GENEID", "TXID","CHRMID","start","end","width","strand","exon_id","exon_name", "exon_rank")
head(exons.list.per.transcript.df.sz.all)

#determine number of exons per transcript
mcols(exons.list.per.transcript)$num_exons <- lengths(exons.list.per.transcript)
exons_per_transcript.sz <- as.data.frame(mcols(exons.list.per.transcript))
exons_per_transcript.sz$TXID <- row.names(exons_per_transcript.sz)
head(exons_per_transcript.sz)
class(exons_per_transcript.sz)

#merge based on TXID
merged.sz <- dplyr::inner_join(exons_per_transcript.sz, exons.list.per.transcript.df.sz.all, by = "TXID")
head(merged.sz)

#subset only the columns that are needed
merged.sz.reduced <- subset(merged.sz, select=c(num_exons, TXID, GENEID, CHRMID))
head(merged.sz.reduced)

# get only distinct rows based on TXID
merged.sz.final <- merged.sz.reduced %>% 
  distinct(TXID, .keep_all = T)
head(merged.sz.final)

# create vector with species names
spID <- as.factor("szebrina")
merged.sz.final <- cbind(merged.sz.final, spID)

#######Merge dataframes #########
merged_df <- rbind(merged.ss.final, 
                   merged.sv.final, 
                   merged.sz.final)

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

boxplot(num_exons ~ ORTHGRP,
        data = merged_df)

#Kruskal Wallis tests
kruskal.test(num_exons ~ spID,
             data = merged_df)

pairwise.wilcox.test(merged_df$num_exons, merged_df$spID,
                     p.adjust.method = "fdr")

####PLOT ####
ggplot(merged_df) +
  aes(x = spID, y = num_exons, color= spID) +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values=c("#FB61D7", "#A58AFF", "#00B6EB"))+
  stat_summary(fun.data=mean_sdl, 
               geom="pointrange", color="black", pch=16, size=.7) +
  stat_summary(fun.y=median, geom="point", size=3, color="black", pch=17) +
  scale_x_discrete(name=" ") +
  scale_y_continuous(name="Number of Exons/Transcript") + 
  theme(legend.title=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=10),
        panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10))
