# Alpha Diversity Tutorial 


data("GlobalPatterns")
library(phyloseq)
library(ggplot2)

GP <- prune_species(speciesSums(GlobalPatterns) > 0, GlobalPatterns)

#plot richness
plot_richness(GP)


# specify richness 
plot_richness(GP,measures = c("Choal1", "Shannon"))

# choose meaningful experimental value "SampleType"
plot_richness(GP, x="SampleType", measures=c("Chao1", "Shannon"))

# adding an external variable 
sampleData(GP)$human <- getVariable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")

# now plot it using that new variable
plot_richness(GP, x="human", color="SampleType", measures=c("Chao1", "Shannon"))

#merge samples
GPst = merge_samples(GP, "SampleType")

# repair them that were damaged during the merge 
sample_data(GPst)$SampleType <- factor(sample_names(GPst))
sample_data(GPst)$human <- as.logical(sample_data(GPst)$human)

# plot this version of merged data 
p = plot_richness(GPst, x="human", color="SampleType", measures=c("Chao1", "Shannon"))
p + geom_point(size=5, alpha=0.7)


#CHAGOS PRACTICE ALPHA DIVERSITY 

ps <- prune_species(speciesSums(ps) > 0, ps)

#plot richness

plot_richness(ps)

#specify richness 
plot_richness(ps,measures = c("Simpson", "Shannon"))

# choose some meaninful experimental values 
sample_variables(ps)

plot_richness(ps, x="site", measures=c("Simpson", "Shannon"))

plot_richness(ps, x="Island", measures=c("Simpson", "Shannon"))

plot_richness(ps, x="depth..m.", measures=c("Simpson", "Shannon"))

plot_richness(ps, x="Temp.Range..oC.", measures=c("Simpson", "Shannon"))



