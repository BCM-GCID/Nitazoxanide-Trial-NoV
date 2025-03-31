library(tidyverse)
library(reshape2)
library(scales)
library(RColorBrewer)
library(patchwork)

## input information
# folder with per sample variant tables
filePath <- "output_by_subject"
pat <- "*\\d.tsv"

# data import
files <- list.files(path = filePath, pattern = pat, full.names = TRUE)
metadata<- read_xlsx('../p1707_NTZ_metadata.xlsx')
metadata$SampleID <- paste('p1707', metadata$SampLognum,sep='-')

# Create a single spreadsheet combined across all patient variant files - raw, unfiltered
new_data <- purrr::map_df(files, function(x) {
  mydata <- read_tsv(x, col_names = TRUE)
  mydata %>%
    select(POS | REF | ALT | ends_with("ALLEL_FREQ") | contains("SNPEFF") | ends_with("DP") | ends_with("DP") | ends_with("QUAL") | ends_with("AA") | ends_with("CODON")) %>% 
    mutate(across(ends_with("ALLEL_FREQ") | ends_with("DP") | ends_with("QUAL"), as.character)) %>% 
    pivot_longer(!c(POS, REF, ALT), names_to = "type", values_to = "values") %>% 
    filter(values != "NA") %>% 
    separate(type, into = c("SampleID", "type"), sep = "\\.") %>% 
    spread(type, values) %>% 
    type_convert() %>% 
    group_by(SampleID) %>% 
    mutate(Annotation = gsub("_", " ", SNPEFF_ANNOT),
           Annotation = str_replace(Annotation, "^\\w{1}", toupper),
           totalDP = ALT_DP + REF_DP) %>% 
    ungroup() %>% 
    type_convert()%>% 
    # if there is any metadata to add
    inner_join(., metadata, by = "SampleID")
})
write_delim(new_data, paste(filePath, "/allSampleVariantCalls-AA_metadata_UNFILTERED.tsv", sep=""), delim = "\t")


## QC figures
# new_data == UNFILTERED variant file
# filtered data
variants_filt <- 
  new_data %>% filter(ALLEL_FREQ >= 0.1 & totalDP >= 20)

## Figure 1: Detection of single nucleotide variants in study samples pre- and post-filtering
# Allele Frequency Detection Per Volunteer - QC boxplot of allele frequency and depth (top - unfiltered, bottom - filtered)
# final adjustments in Illustrator
# Panel A
fig1a <- new_data %>% select(POS, Volnum, SampleID, ALLEL_FREQ, totalDP) %>% 
  mutate(freqBin = case_when(ALLEL_FREQ < 0.1 ~ 'lowFreq',
                             ALLEL_FREQ >= 0.1 & ALLEL_FREQ <= 0.5 ~ 'medFreq',
                             ALLEL_FREQ > 0.5 & ALLEL_FREQ <= 0.9 ~ 'common',
                             ALLEL_FREQ > 0.9 ~ 'fixed')) %>% 
  ggplot(aes(x=factor(new_data$Volnum, levels=c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "Stool")), y=totalDP)) +
  geom_boxplot(outlier.shape = NA, color = "grey") +
  geom_jitter(aes(color = freqBin), alpha=0.8, width = 0.2, size = 2) +
  geom_hline(yintercept = 20, color = "black", linetype = "dashed") +
  scale_y_log10(labels = label_number_si()) +
  scale_color_manual(name = "Alternative Allele \nFrequency",
                     breaks = c("lowFreq", "medFreq", "common", "fixed"),
                     labels = c("lowFreq" = "Low frequency (< 10%)",
                                "medFreq" = "Medium frequency (10% - 50%)", "common" = "High frequency (50% - 90%)",
                                "fixed" = "Fixed allele ( > 90%)"),
                     values = c("orchid", "dodgerblue", "blue", "darkolivegreen")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_text(size=10),
    axis.text.y  = element_text(size=10),
    axis.title.y = element_text(size=11),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 11)) +
  labs(y="Total mapping depth (log)") 

ggsave(file=paste(filePath, "/Figure1a.pdf", sep=""), plot = fig1a, width = 11.69, height = 8.27, units = c("in"))

# Panel B
fig1b <- variants_filt %>% select(POS, Volnum, SampleID, ALLEL_FREQ, totalDP) %>% 
  mutate(freqBin = case_when(ALLEL_FREQ < 0.1 ~ 'lowFreq',
                             ALLEL_FREQ >= 0.1 & ALLEL_FREQ <= 0.5 ~ 'medFreq',
                             ALLEL_FREQ > 0.5 & ALLEL_FREQ <= 0.9 ~ 'common',
                             ALLEL_FREQ > 0.9 ~ 'fixed')) %>% 
  ggplot(aes(x=factor(variants_filt$Volnum, levels=c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10", "Stool")), y=totalDP)) +
  geom_boxplot(outlier.shape = NA, color = "grey") +
  geom_jitter(aes(color = freqBin), alpha=0.8, width = 0.2, size = 2) +
  geom_hline(yintercept = 20, color = "black", linetype = "dashed") +
  scale_y_log10(labels = label_number_si()) +
  scale_color_manual(name = "Alternative Allele \nFrequency",
                     breaks = c("lowFreq", "medFreq", "common", "fixed"),
                     labels = c("lowFreq" = "Low frequency (< 10%)",
                                "medFreq" = "Medium frequency (10% - 50%)", "common" = "High frequency (50% - 90%)",
                                "fixed" = "Fixed allele ( > 90%)"),
                     values = c("orchid", "dodgerblue", "blue", "darkolivegreen")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_text(size=10),
    axis.text.y  = element_text(size=10),
    axis.title.y = element_text(size=11),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 11)) +
  labs(y="Total mapping depth (log)")

ggsave(file=paste(filePath, "/Figure1b.pdf", sep=""), plot = fig1b, width = 11.69, height = 8.27, units = c("in"))



## Figure 6 (to become suppl): Single nucleotide variants in the GI.1 genome for each person in the study
# SNPs_all-DPC_DP-AF-filter
# Bubble plot - all samples, days post challenge, filtered
#variants_filt<-variants_filt[variants_filt$SampleID=='p1540-BCM16-16-AP',]
fig2 <- variants_filt %>% 
  # select(SampleID, Vist2, Annotation, POS, ALLEL_FREQ, Volnum) %>% 
  # mutate(Vist2 = factor(Vist2, levels = sort(unique(as.numeric(Vist2)),decreasing=TRUE))) %>% 
  ggplot() +
  geom_point(aes(x=POS, y=factor(Volnum, levels=c( "P10", "P9", "P8", "P7", "P6", "P5", "P4", "P3", "P2", "P1", "Stool")),
                 color = Annotation, size = ALLEL_FREQ), shape = 16) +
 # facet_wrap(~ Volnum, scales = "free_x", ncol = 5) +
  scale_x_continuous(limits = c(0,7654), breaks = seq(0, 7654, by = 1000)) +
  scale_color_manual(name = "Variant Effect",
                     breaks = c("Synonymous variant", "Missense variant",
                                "Downstream gene variant", "Stop gained",
                                "Frameshift variant", "Disruptive inframe deletion", "Stop lost&splice region variant", 
                                "Upstream gene variant", "Start lost"),
                     labels=c("Synonymous variant", "Non-Synonymous variant",
                              "Downstream gene variant", "Stop gained",
                              "Frameshift variant", "Disruptive inframe deletion", "Stop lost&splice region variant",
                              "Upstream gene variant", "Start lost"),
                     values = c("#D81B60", "#1E88E5", "#E2C878", "#61E2CC", "#905715", "#9C9A95", "#CCCBC9", "#373737", "#020202")) +
  scale_size(range = c(1, 4),
             limits = c(0.1, 1),
             breaks = seq(0.1, 1, by = 0.2),
             name="Frequency\n(bubble size)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x  = element_text(size=8),
    axis.text.y  = element_text(size=8),
    axis.title = element_blank()) +
  # edit
  labs(y="Days Post Challenge", x="Position in genome") 

ggsave(file=paste(filePath, "/Figure2.pdf", sep=""), plot = fig2, width = 11.69, height = 8.27, units = c("in"))

