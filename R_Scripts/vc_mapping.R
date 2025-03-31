library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)
library(readxl)
library(lubridate)
library(tidyr)
library(ggpubr)
library(rprojroot)
library(wesanderson)
library(plotly)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#produces a sample to reference mapping file (with following columns) which is used in variant calling pipeline. 
#column 1: SampleID of the sample
#column 2: Reference SampleID the sample is mapped to for variant calling. 

metadata<-read_excel('../p1707_NTZ_metadata.xlsx')
summary_stats <-  read.table("../Deduplicated/NoV.summaryStats.txt", sep = "\t", header = TRUE) %>%
                    mutate(SampLognum = as.character(sub(".*-", "", SampleID)))

v1_only <- metadata %>%
            filter(VisitID == 'V1')

mapping<-metadata %>%
  filter(VisitID != 'V1') %>%
  inner_join(v1_only, by=join_by(SubjectID==SubjectID))%>%
  mutate(SampLognum.x = as.character(SampLognum.x))%>%
  inner_join(summary_stats, by=join_by(SampLognum.x==SampLognum))%>%
  filter(Length > 7200 & Completness >= 90)%>%
 mutate(SampLognum.x = paste('p1707-', SampLognum.x, sep=''))%>%
  mutate(SampLognum.y = paste('p1707-', SampLognum.y, sep=''))%>%
select(SampLognum.x, SampLognum.y)

write.table(mapping, file='mapping.tsv', col.names=FALSE, quote=FALSE, sep='\t', row.names = FALSE)