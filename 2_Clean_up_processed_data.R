library(phyloseq)
library(tidyverse)

# Load data ####
ps <- readRDS("phyloseq_object_16S.RDS")

# quick peak at ps object

sample_names(ps)
rank_names(ps)
sample_variables(ps)

library("phyloseq"); packageVersion("phyloseq")

library("ggplot2"); packageVersion("ggplot2")

# Housekeeping ####
# Rename metadata columns
names(ps@sam_data) <- c("Multiplex ID","Library ID","SampleID","Location","Country","Species","CoralAge","CoralAgeBinned","Average_LE_mm","GPS","Control","Genotype")
# Add Lat/Lon columns
LAT <- unlist(map(str_split(ps@sam_data$GPS,pattern = " "),1))
LON <- unlist(map(str_split(ps@sam_data$GPS,pattern = " "),2))
LAT <- str_remove(LAT,"N")
LON <- str_remove(LON,"E")
ps@sam_data$LAT <- as.numeric(LAT)
ps@sam_data$LON <- as.numeric(LON)
meta = as.data.frame(sample_data(ps))


# Subset taxa to bacteria only (remove eukaryotes)
table(tax_table(ps)[,1])
ps <- subset_taxa(ps, Kingdom == "Bacteria")
colnames(tax_table(ps)) # No species-level assignments (see previous script)


# remove empty (and singleton) ESVs and Samples
summary(rowSums(otu_table(ps))) # no empty samples
otu_table(ps) <- otu_table(ps)[,colSums(otu_table(ps)) > 1]

# Save raw abundance phyloseq object
saveRDS(ps,"./output/phyloseq_object_16S_cleaned.RDS")

# Normalize with relative abundance
ps_ra <- transform_sample_counts(ps, function(x) x / sum(x) )

# Save relabund phyloseq object
saveRDS(ps_ra,"./output/phyloseq_cleaned_relabund.RDS")