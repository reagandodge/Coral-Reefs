library("phyloseq")
library("ggplot2")
library(tidyverse)
library(vegan)
theme_set(theme_bw())
data("GlobalPatterns")
data("esophagus")
data("enterotype")
data(soilrep)

?GlobalPatterns


example(enterotype, ask = FALSE)

?phyloseq
?otu_table
?sample_data
?tax_table
# Create a pretend OTU table 
otumat = matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol = 10)
otumat
# add row names to the OTU table 
rownames(otumat) <- paste0("OTU", 1:nrow(otumat))
colnames(otumat) <- paste0("Sample", 1:ncol(otumat))
otumat


# pretend taxonomic table 

taxmat = matrix(sample(letters, 70, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat

class(otumat)
class(taxmat)

# now we are going to combine into a phyloseq object 
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU

physeq = phyloseq(OTU, TAX )

physeq

plot_bar(physeq, fill = "Family")

# creating a random sample data set to add to this 
sampledata = sample_data(data.frame(
  Location = sample(LETTERS[1:4], size=nsamples(physeq), replace=TRUE),
  Depth = sample(50:1000, size=nsamples(physeq), replace=TRUE),
  row.names=sample_names(physeq),
  stringsAsFactors=FALSE
))
sampledata

physeq
physeq = phyloseq(otu_table(physeq),tax_table(physeq),sample_data(sampledata))

plot_bar(physeq, fill = "Phylum") + facet_wrap(~Location)

# creating a random phylogentic tree 

library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

# now merge new data with current phyloseq data 
physeq1 = merge_phyloseq(physeq, sampledata, random_tree)
physeq1


# are they identical?
identical(physeq1, physeq) # no 

# lets build a couple of trees 

plot_tree(physeq1, color="Location", label.tips="taxa_names", ladderize="left", plot.margin=0.3)

plot_tree(physeq1, color="Depth", shape="Location", label.tips="taxa_names", ladderize="right", plot.margin=0.3)


# heat maps 

plot_heatmap(physeq1)

plot_heatmap(physeq1, taxa.label="Phylum")


#### CHAGOS PRACTICE  ####

ps <- readRDS("phyloseq_object_16S.RDS")

sample_names(ps)
rank_names(ps)
sample_variables(ps)


head(otu_table(ps)) 

head(sample_data(ps))

head(tax_table(ps))

ps
ps_ra = transform_sample_counts(ps,function(x){x/sum(x)})

# quick ordination
ord1=ordinate(ps_ra,method="DCA")

plot_ordination(ps,ord1,color="Island")

vegan::adonis(otu_table(ps_ra) ~ ps_ra@sam_data$`Colony Colour`)

# shannon diversity
shannon = vegan::diversity(otu_table(ps_ra),index = "shannon")
sample_data(ps_ra)$Shannon = shannon
?specnumber
# richness
rich = specnumber(otu_table(ps_ra))
sample_data(ps_ra)$Richness = rich

meta = as.data.frame(sample_data(ps_ra))

meta$salinity = as.numeric(meta$salinity)

shannonmod = glm(Shannon ~ salinity,data=meta)
summary(shannonmod)



x = sample_data(ps_ra)
class(x)
ps_ra@sam_data$Island
sample_data(ps_ra)$Island

# look at it
plot_bar(ps, fill = "Phylum") + facet_wrap(~Island) + theme(legend.position="bottom")


# phylotree
library(phangorn)
library(msa)

#### Build Phylogenetic Tree ####
seqs <- rownames(tax_table(ps))
names(seqs) <- seqs # This propagates to the tip labels of the tree

# alignment (took 37 hours to complete!)
alignment <- msa(seqs,method = "Muscle", type = "dna",verbose = TRUE,order = "input",maxiters = 12)

# save progress 
saveRDS(alignment,"../Chagos_2/WWP_dna_alignment_muscle.RDS")


# re-load point
alignment <- readRDS("../Chagos_2/WWP_dna_alignment_muscle.RDS")
phang.align = as.phyDat(alignment, type = "DNA")

# distance max likelihood
dm <- dist.ml(phang.align)

#save
saveRDS(dm,"../Chagos_2/ML_Distance.RDS")

# neighbor-joining tree
treeNJ <- NJ(dm) # Note, tip order != sequence order
treeNJ$tip.label <- seqs

#save
saveRDS(treeNJ, "../Chagos_2/treeNJ.RDS")
treeNJ <- readRDS("../Chagos_2/treeNJ.RDS")


fit = pml(treeNJ, data=phang.align)

#save
saveRDS(fit,"../Chagos_2/fit_treeNJ.RDS")

## negative edges length changed to 0!

fitJC <- optim.pml(fit, TRUE)

# save
saveRDS(fitJC, "../Chagos_2/WWP_fitJC.RDS") # This is the new tree using optim.pml
write.tree(fitJC$tree, file = "../Chagos_2/WWP_bact_tree_JC.nwk")

# reload point
fitJC <- readRDS("../Chagos_2/WWP_fitJC.RDS")

(unique(names(seqs)))



# Bootstrap 
bs = bootstrap.pml(fitJC, bs=100, control = pml.control(trace = 0),multicore = TRUE, mc.cores = 1)

#?bootstrap.pml
# check tree with plot

plotBS(fitJC$tree,type="p")
plot_tree()

ps = phyloseq(otu_table(ps@otu_table),
              tax_table(ps@tax_table),
              sample_data(ps@sam_data),
              phy_tree(fitJC$tree))

#save this as rds
#load tidyverse and objects in a new script

tree<-plot_tree(ps)

saveRDS(tree,"../Chagos_2/PHYLO_tree.RDS")



#### redo with GTR model?
# fitGTR <- update(fit, k=4, inv=0.2)
# fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                               control = phangorn::pml.control(trace = 1L),rearrangement = "stochastic")


# saveRDS(fitGTR, "./output/WWP_fitGTR2.RDS") # This is the new tree using optim.pml
# fitGTR = readRDS("./output/WWP_fitGTR2.RDS") # loading new tree. does it work??
# write.tree(fitGTR$tree, file = "./output/WWP_bact_tree.nwk")

# bs = bootstrap.pml(fitGTR, bs=100, optNni=TRUE, multicore=TRUE)

detach("package:phangorn", unload=TRUE)












# Subset taxa to bacteria only (remove eukaryotes)
table(tax_table(ps)[,1]) # 21 eukaryotes, remove them 
ps <- subset_taxa(ps, Kingdom == "Bacteria")
colnames(tax_table(ps))

# combine... 
otu_table(ps) <- otu_table(ps)[,colSums(otu_table(ps)) > 1]



