library(gggenes)
#https://www.ncbi.nlm.nih.gov/nuccore/NC_001959.2
genes <- as_tibble(data.frame(molecule = c("GI.1","GI.1","GI.1","GI.1","GI.1"),
                              gene = c("5'-UTR","ORF1","ORF2","ORF3", "3'-UTR"), 
                              start = c(1,5,5358,6950,7589), 
                              end = c(4,5374,6950,7588,7654),
                              strand = c("forward","forward","forward","forward","forward"),
                              orientation = c("-1","-1","-1","-1","-1")))

proteins <- as_tibble(data.frame(molecule = c("GI.1","GI.1","GI.1","GI.1","GI.1","GI.1","GI.1","GI.1"),
                                    gene = c("ORF1", "ORF1","ORF1","ORF1","ORF1","ORF1", "ORF2","ORF3"), 
                                    start = c(5,5,5,5,5,5,5358,6950), 
                                    end = c(5374,5374,5374,5374,5374,5374,6950,7588),
                                    protein = c("p48", "NTPase", "p22", "VPg", "Pro", "RdRp", "VP1","VP2"), 
                                    from = c(5,1199,2288,2891,3305,3848,5358,6950), 
                                    to = c(1198,2287,2890,3304,3847,5371,6950,7588),
                                    strand = c("forward","forward","forward","forward","forward","forward","forward","forward"),
                                    orientation = c("-1","-1","-1","-1","-1","-1","-1","-1")))


geneMap <- ggplot(genes, aes(xmin = start, xmax = end, y = molecule)) +
  geom_gene_arrow(fill = "white") +
  geom_subgene_arrow(data = proteins,
                     aes(xmin = start, xmax = end, y = molecule, fill = protein,
                         xsubmin = from, xsubmax = to), color="black", alpha=.7) +
  geom_subgene_label(data = proteins,
                     aes(xsubmin = from, xsubmax = to, label = protein),
                     min.size = 0) +
  scale_fill_manual(name = "Gene Annotation",
                    breaks = c("p48", "NTPase", "p22", "VPg", "Pro", "RdRp", "VP1","VP2"),
                    #values = c("mediumpurple", "mediumorchid", "lightpink4","maroon4", "plum", "rosybrown3", "turquoise4", "darkgoldenrod2"),
                    values = c("#D8CDE5", "#D8CDE5", "#D8CDE5","#D8CDE5", "#D8CDE5", "#D8CDE5", "#D8CDE5", "#D8CDE5")) +
  scale_x_continuous(limits = c(0,7654), breaks = seq(0, 7654, by = 500)) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none"
  )
ggsave(geneMap, file="figures_pre-edit/GI.1-geneMap-oneColour.pdf", width = 11.69, height = 8.27, units = c("in"))

