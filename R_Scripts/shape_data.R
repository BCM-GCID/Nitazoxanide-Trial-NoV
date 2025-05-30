library(tidyverse)
library(ggplot2)
library(readxl)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#This script will combine variant calls from all samples into a single file per subject/patient. 

# https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
# https://pcingola.github.io/SnpEff/adds/VCFannotationformat_v1.0.pdf
# https://pcingola.github.io/SnpEff/se_inputoutput/

metadata<- read_xlsx('../p1707_NTZ_metadata.xlsx')
metadata$SampleID <- paste('p1707', metadata$SampLognum,sep='-')

filedir<- './output_by_subject'

filelist <- data.frame(list.files(path='../vc_output', full.names = TRUE, recursive = TRUE, pattern = '.snpEff.vcf'))
colnames(filelist) <- c('path')
filelist$SampleID <- sapply(filelist$path, function(x){unlist(strsplit(x, "/"))[3]})

metadata <- merge(metadata, filelist, by='SampleID')


for (patientID in unique(metadata$SubjectID)){
  print(patientID)
  patientDF<-metadata[metadata$SubjectID==patientID,]
  print(nrow(patientDF))
  
  finalDF<-NULL
  for(path in patientDF$path)
  {
    df <- try(read.table(path,  header=FALSE,  sep="\t"))
    if(inherits(df, "try-error"))
      next 
   
    colnames(df) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "DATA")
    sampleDF <-  df[,c( "POS", "REF", "ALT")]
    sampID <- patientDF[patientDF$path == path, ]$SampleID
    print(sampID)
    
    #AF=0.315436;GFF_FEATURE=NA;REF_CODON=NA;REF_AA=NA;ALT_CODON=NA;ALT_AA=NA;
    sampleDF$ALLEL_FREQ<- sapply(df$INFO, function(x){unlist(strsplit(unlist(strsplit(x, ";"))[1], "="))[2]})
    sampleDF$GFF_FEATURE<- sapply(df$INFO, function(x){unlist(strsplit(unlist(strsplit(x, ";"))[2], "="))[2]})
    sampleDF$REF_CODON<- sapply(df$INFO, function(x){unlist(strsplit(unlist(strsplit(x, ";"))[3], "="))[2]})
    sampleDF$REF_AA<- sapply(df$INFO, function(x){unlist(strsplit(unlist(strsplit(x, ";"))[4], "="))[2]})
    sampleDF$ALT_CODON<- sapply(df$INFO, function(x){unlist(strsplit(unlist(strsplit(x, ";"))[5], "="))[2]})
    sampleDF$ALT_AA<- sapply(df$INFO, function(x){unlist(strsplit(unlist(strsplit(x, ";"))[6], "="))[2]})
    
    # GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ
    # 1 :2     :0     :68      :9     :9     :37      :0.818182
    sampleDF$REF_DP <- sapply(df$DATA, function(x){unlist(strsplit(x, ":"))[2]})
    sampleDF$REF_QUAL <- sapply(df$DATA, function(x){unlist(strsplit(x, ":"))[4]})
    sampleDF$ALT_DP <- sapply(df$DATA, function(x){unlist(strsplit(x, ":"))[5]})
    sampleDF$ALT_QUAL <- sapply(df$DATA, function(x){unlist(strsplit(x, ":"))[7]})
    
    is.na(sampleDF) <- sampleDF == "NA"
    
    snpeff_annots<- sapply(df$INFO, function(x){unlist(strsplit(unlist(strsplit(x, ";"))[7], "="))[2]})
    snpeff_annots <- sapply(snpeff_annots, function(x){unlist(strsplit(x, ","))[1]})
    snpeff_annots <- paste(snpeff_annots, "|" ,sep="")
    
    #T|missense_variant|MODERATE|Gene_4_5373|Gene_4_5373|transcript|AAB50465.1|protein_coding|1/1|c.4052C>T|p.Ala1351Val|4052/5370|4052/5370|1351/1789||
    sampleDF$SNPEFF_ANNOT <- sapply(snpeff_annots, function(x){unlist(strsplit(x, "\\|"))})[2,]
    sampleDF$SNPEFF_IMPACT <- sapply(snpeff_annots, function(x){unlist(strsplit(x, "\\|"))})[3,]
    sampleDF$SNPEFF_GENE <- sapply(snpeff_annots, function(x){unlist(strsplit(x, "\\|"))})[4,]
    sampleDF$SNPEFF_FTYPE <- sapply(snpeff_annots, function(x){unlist(strsplit(x, "\\|"))})[6,]
    sampleDF$SNPEFF_FID <- sapply(snpeff_annots, function(x){unlist(strsplit(x, "\\|"))})[7,]
    sampleDF$SNPEFF_BIOTYPE <- sapply(snpeff_annots, function(x){unlist(strsplit(x, "\\|"))})[8,]
    sampleDF$SNPEFF_HGVSC <- sapply(snpeff_annots, function(x){unlist(strsplit(x, "\\|"))})[10,]
    sampleDF$SNPEFF_HGVSP <- sapply(snpeff_annots, function(x){unlist(strsplit(x, "\\|"))})[11,]
    
    colnames(sampleDF)[4:21] <- paste(sampID, colnames(sampleDF[,c(4:21)]),sep=".")
    
    if (is.null(finalDF)) {
      finalDF <- sampleDF
    }
    else {
      finalDF <- merge(finalDF, sampleDF, all=TRUE, by=c( "POS", "REF", "ALT"))
    }
    
    write.table(finalDF, file=paste(filedir, '/', patientID, '.tsv', sep=''), quote=FALSE, sep='\t', row.names = FALSE)
  }
  
}








