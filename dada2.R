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
path <- "raw_data" 
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
# cut fwd reads at 200 and rev reads at 200
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,200),
                     #maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     #compress=TRUE, multithread=TRUE)
#save this step as an rds object 
saveRDS(out,"Chagos_Project.RDS")
out<- readRDS("Chagos_Project.RDS")



head(out)

# reassign filts for lost samples, if any
filtFs <- list.files("raw_data/filtered/", pattern = "_F_filt", full.names = TRUE)
filtRs <- list.files("raw_data/filtered/", pattern = "_R_filt", full.names = TRUE)


# Learn error rates ####
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e8)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e8)

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

# Save progress
saveRDS(seqtab.nochim,"./output/seqtab.nochim.RDS")

# re-load point for sequence table ####
seqtab.nochim <- readRDS("./output/seqtab.nochim.RDS")

# Track reads through pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

# Export results
write.csv(track, file = "./output/tracked_reads.csv", quote = FALSE)


# remove all ASVs that don't have at least 2 hits ####
seqtab.nochim <- seqtab.nochim[,(colSums(seqtab.nochim) > 1)]

# Assign taxonomy - Silva v132 exact match / 80% bootstrap min ####
taxa <- assignTaxonomy(seqtab.nochim,"./taxonomy/silva_nr_v132_train_set.fa.gz", minBoot = 80,multithread = TRUE)
saveRDS(taxa,"./output/taxa.RDS")

# re-load point for taxonomy ####
taxa <- readRDS("./output/taxa.RDS")


# rename seqtab object samples
seqtab.df <- as.data.frame(seqtab.nochim)
row.names(seqtab.df) <- map(strsplit(row.names(seqtab.df), "_"),1)

# Prepare for PhyloSeq object ####

# read in metadata
meta = read_csv("./metadata_edited.csv")

# subset meta and seqtab to match
in.meta = which(names(seqtab.nochim[,1]) %in% meta$`Library ID` == TRUE)
seqtab.nochim = seqtab.nochim[in.meta,]

in.seqtab = which(meta$`Library ID` %in% names(seqtab.nochim[,1]))
meta = meta[in.seqtab,]

# re-order
meta = meta[order(meta$`Library ID`),]
row.names(meta) <- meta$`Library ID`
# Check 
identical(row.names(seqtab.nochim), meta$`Library ID`)


# Remove contaminants ####
# find contaminants
contams = isContaminant(seqtab.nochim, neg = meta$Control, normalize = TRUE)
table(contams$contaminant)
write.csv(contams, file = "./output/likely_contaminants.csv", row.names = TRUE)


# remove them
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$Control == FALSE,]
meta = meta[meta$Control == FALSE,]

# subset meta and seqtab to match
in.meta = which(names(seqtab.nochim[,1]) %in% meta$`Library ID` == TRUE)
seqtab.nochim = seqtab.nochim[in.meta,]

in.seqtab = which(meta$`Library ID` %in% names(seqtab.nochim[,1]))
meta = meta[in.seqtab,]

#re-order
meta = meta[order(meta$`Library ID`),]
row.names(meta) <- meta$`Library ID`
# Check
identical(row.names(seqtab.nochim), meta$`Library ID`)


# make phyloseq object ####
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxa))
# save it
saveRDS(ps, "./output/phyloseq_object_16S.RDS")


© 2019 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
About
