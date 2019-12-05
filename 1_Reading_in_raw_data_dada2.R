#### 00_READING IN RAW DATA 



library(dada2); packageVersion("dada2")
library(vegan)
library(tidyverse)
library(phyloseq)
library(decontam)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")


# Find raw fastq files and prepare workspace ####
path <- "../../Chagos_Project/raw_data" 
list.files(path)

# Parse fwd and rev reads
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Get Sample Names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 3)


# Peek at quality profiles
plotQualityProfile(fnFs[c(1,2)]) # fwd reads
plotQualityProfile(fnRs[c(1,2)]) # rev reads

# Make filtered outfile names
length(fnFs) 
length(fnRs) 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

# Filter and trim ####
#go back and run this again and rerun at 250 and 160 
# cut fwd reads at 250 and rev reads at 160
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
#save this step as an rds object 
saveRDS(out,"Chagos_Project_filt250_160.RDS")
out<- readRDS("Chagos_Project_filt250_160.RDS")



head(out)

# reassign filts for lost samples, if any
filtFs <- list.files("raw_data/filtered", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("raw_data/filtered", pattern = "_R_filt", full.names = TRUE)


# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# plot error rates for sanity
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate ####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

# Sample Inferrence ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge fwd and rev reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
# Save progress
saveRDS(seqtab.nochim,"seqtab.nochim.RDS")

# re-load point for sequence table ####
seqtab.nochim <- readRDS("seqtab.nochim.RDS")

# Track reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
# Export results
write.csv(track, file = "tracked_reads.csv", quote = FALSE)
track <- read.csv("tracked_reads.csv")
plot(track$nonchim)
# remove all ASVs that don't have at least 2 hits ####
seqtab.nochim <- seqtab.nochim[,(colSums(seqtab.nochim) > 1)]

# Assign taxonomy - Silva v132 exact match / 80% bootstrap min ####
taxa <- assignTaxonomy(seqtab.nochim,"../../silva_nr_v132_train_set.fa.gz",multithread = TRUE)


# More Species addition steps for 16S
taxa <- addSpecies(taxa, "../../silva_species_assignment_v132.fa.gz")
saveRDS(taxa,"taxa.RDS")


# Inspect taxonomic assignments 
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# re-load point for taxonomy ####
taxa <- readRDS("taxa.RDS")


# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df) <- map(strsplit(row.names(seqtab.df), "_"),1)

# Prepare for PhyloSeq object ####

# read in metadata
library(readxl)
meta = read_excel("metadata/Chagos_PlateMpas_MetaData_1_Geoff.xlsx")

# subset meta and seqtab to match
in.meta = which(names(seqtab.nochim[,1]) %in% meta$`Pooled library ID (from GIS)`== TRUE)
seqtab.nochim = seqtab.nochim[in.meta,]

in.seqtab = which(meta$`Pooled library ID (from GIS)`%in% names(seqtab.nochim[,1]))
meta = meta[in.seqtab,]

# re-order
meta = meta[order(meta$`Pooled library ID (from GIS)`),]
row.names(meta) <- meta$`Pooled library ID (from GIS)`
# Check 
#identical(row.names(seqtab.nochim), meta$`Pooled library ID (from GIS)`)
row.names(seqtab.nochim)
meta$`Pooled library ID (from GIS)` <-  as.character(meta$`Pooled library ID (from GIS)`)

currentorder = order(as.numeric(names(seqtab.nochim[,1])))
seqtab.nochim=seqtab.nochim[currentorder,]
names(seqtab.nochim[,1])

identical(row.names(seqtab.nochim), meta$`Pooled library ID (from GIS)`)

meta$Negatives = (meta$salinity == "Blank")
# Remove contaminants ####
# find contaminants
contams = isContaminant(seqtab.nochim, neg = meta$Negatives, normalize = TRUE)
table(contams$contaminant)
write.csv(contams, file = "likely_contaminants.csv", row.names = TRUE)

# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$Negatives == FALSE,]
meta = meta[meta$Negatives == FALSE,]

# subset meta and seqtab to match
in.meta = which(names(seqtab.nochim[,1]) %in% meta$`Pooled library ID (from GIS)` == TRUE)
seqtab.nochim = seqtab.nochim[in.meta,]

in.seqtab = which(meta$`Pooled library ID (from GIS)` %in% names(seqtab.nochim[,1]))
meta = meta[in.seqtab,]

#re-order
meta = meta[order(meta$`Pooled library ID (from GIS)`),]
row.names(meta) <- meta$`Pooled library ID (from GIS)`
# Check
identical(row.names(seqtab.nochim), meta$`Pooled library ID (from GIS)`)


# make phyloseq object ####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))
# save it
saveRDS(ps, "phyloseq_object_16S.RDS")


