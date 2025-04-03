## Install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.19")

## Load packages
library("dada2"); packageVersion("dada2")
library("DECIPHER"); packageVersion("DECIPHER")
library("phyloseq"); packageVersion("phyloseq")
library("Biostrings"); packageVersion("Biostrings")
library("ggplot2"); packageVersion("ggplot2")
library("decontam"); packageVersion("decontam")
library("DESeq2"); packageVersion("DESeq2")
library("doParallel"); packageVersion("doParallel")
library("foreach"); packageVersion("foreach")
library("ape"); packageVersion("ape")
library("vegan"); packageVersion("vegan")
library("dplyr"); packageVersion("dplyr")
library("nortest"); packageVersion("nortest")
library("dunn.test"); packageVersion("dunn.test")
library("lawstat"); packageVersion("lawstat")
library("nlme"); packageVersion("nlme")
library("multcomp"); packageVersion("multcomp")
library("WeightedTreemaps"); packageVersion("WeightedTreemaps")
library("reshape2"); packageVersion("reshape2")
library("beeswarm"); packageVersion("beeswarm")
library("AICcmodavg"); packageVersion("AICcmodavg")
library("VennDiagram")
library("FUNGuildR")
library("MiscMetabar")
library("microeco")
library("file2meco")
library("rCSCS")
library("RCurl")

## Set working directory
setwd("~/UL/Recherche/Chimie_des_lichens/ADN/Aviti/")


#################### BACTERIA TAXONOMY ##################
# Set path to bacteria data
path <- "~/UL/Recherche/Chimie_des_lichens/ADN/Aviti/Metabolome/16S/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", 
                        full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", 
                        full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

## Inspect read quality profiles
# Forward - Great quality - no need to trim
plotQualityProfile(fnFs[1:2])
# Reverse - Great quality - no need to trim
plotQualityProfile(fnRs[1:2])

## Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", 
                    paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", 
                    paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(300,260), maxN=0, 
                     maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 
# On Windows set multithread=FALSE
out

## Learn the error rates
# Forward
errF <- learnErrors(filtFs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
# Reverse
errR <- learnErrors(filtRs, multithread=FALSE)
plotErrors(errR, nominalQ=TRUE)

## Sample inference
# Forward
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaFs[[1]]
# Reverse
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaRs[[1]]

## Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                      verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## remove sequences that did not correctly merged
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% seq(400,470)]
table(nchar(getSequences(seqtab)))

## Remove chimeras
seqtab.nochim_16S <- removeBimeraDenovo(seqtab,
                                        method="consensus",
                                        multithread=FALSE,
                                        verbose=TRUE)
dim(seqtab.nochim_16S)
sum(seqtab.nochim_16S)/sum(seqtab)

## Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim_16S))
colnames(track) <- c("input", "filtered", 
                     "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample.names
track

## IDTAXA with Decipher
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DECIPHER")
library(DECIPHER); packageVersion("DECIPHER")

# Assign taxonomy
# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim_16S)) 
load("Taxa/SILVA_SSU_r138_2019_DECIPHER.RData") 
ids <- IdTaxa(dna, trainingSet, strand="top", 
              processors=NULL, verbose=FALSE)
ranks <- c("domain", "phylum", "class", "order", 
           "family", "genus", "species")
# Convert the output object of class "Taxa" to 
# a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim_16S)

## Carry on using taxid
taxa_16S <- taxid


#################### FUNGI TAXONOMY ##################
# Set path to fungi data
path <- "~/UL/Recherche/Chimie_des_lichens/ADN/Aviti/Metabolome/ITS/"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", 
                        full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", 
                        full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Inspect read quality profiles
# Forward
plotQualityProfile(fnFs[1:2])
# Reverse
plotQualityProfile(fnRs[1:2])

## Filter and trim
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0,
                     maxEE = c(2, 2), truncQ = 2, 
                     minLen = 50, rm.phix = TRUE, 
                     compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE
out

## Learn error rates
# Forward
errF <- learnErrors(filtFs, multithread=FALSE)
plotErrors(errF, nominalQ=TRUE)
# Reverse
errR <- learnErrors(filtRs, multithread=FALSE)
plotErrors(errR, nominalQ=TRUE)

## Sample inference
# Forward
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaFs[[1]]
# Reverse
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
dadaRs[[1]]

## Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                      verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

## Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

## Remove chimeras
seqtab.nochim_ITS <- removeBimeraDenovo(seqtab,
                                        method="consensus",
                                        multithread=FALSE,
                                        verbose=TRUE)
dim(seqtab.nochim_ITS)
sum(seqtab.nochim_ITS)/sum(seqtab)

## Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim_ITS))
colnames(track) <- c("input", "filtered", 
                     "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample.names
track

## IDTAXA with Decipher
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DECIPHER")
library(DECIPHER); packageVersion("DECIPHER")

# Assign taxonomy
# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim_ITS)) 
load("Taxa/UNITE_v2023_July2023_DECIPHER.RData") 
ids <- IdTaxa(dna, trainingSet, strand="top", 
              processors=NULL, verbose=FALSE)
ranks <- c("domain", "phylum", "class", "order", 
           "family", "genus", "species")
# Convert the output object of class "Taxa" to 
# a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim_ITS)

## Carry on using taxid
taxa_ITS <- taxid



################### PHYLOSEQ OBJECTS ###############
# BiocManager::install("phyloseq", force = TRUE )
# BiocManager::install("decontam")
theme_set(theme_bw())

## BACTERIA
#Creation du df avec les donnees des echantillons
samples.out<- rownames(seqtab.nochim_16S)
samples.out
samples.out <- rownames(seqtab.nochim_16S)
echantillon <- sapply(strsplit(samples.out, "_"), `[`, 1)
type <- substr(echantillon,1,1)
species <- substr(echantillon,3,6)
substrate <- substr(echantillon,8,11)
subtype <- substr(echantillon,13,15)
number <- substr(echantillon,17,18)
samdf <- data.frame(echantillon=echantillon, 
                    type=type, 
                    species=species, 
                    substrate=substrate,
                    subtype=subtype,
                    number=number)
rownames(samdf) <- samples.out

ps_16S <- phyloseq(otu_table(seqtab.nochim_16S,
                             taxa_are_rows=FALSE),
                   sample_data(samdf),
                   tax_table(taxa_16S))
# Have a look on the data structure
ps_16S
head(sample_data(ps_16S), 10)

# Use short name for ASVs rather than the full DNA sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps_16S))
names(dna) <- taxa_names(ps_16S)
ps_16S <- merge_phyloseq(ps_16S, dna)
taxa_names(ps_16S) <- paste0("ASV", 
                             seq(ntaxa(ps_16S)))

# Remove ambiguous phylum, chloroplast, and mitochondria
ps_16S <- subset_taxa(ps_16S, !is.na(phylum)
                      & !order %in% c("Chloroplast")
                      & !family %in% c("Mitochondria"))
ps_16S

## Create table, number of features for each phyla
table(tax_table(ps_16S)[, "phylum"], exclude = NULL)

# Exporter la table ASV dans un fichier Excel 
# write.csv2(tax_table(ps_16S), "Taxa/ASV_table_16S.csv",
#            row.names = TRUE)

# Check the number of reads per sample
hist(sample_sums(ps_16S), 
     main="Histogram: Read Counts", 
     xlab="Total Reads", border="black", 
     col="darkgreen", las=1, breaks=15)

# Check the number of reads per taxa
head(taxa_sums(ps_16S), 50)


##FUNGI
#Creation du df avec les donnees des echantillons
samples.out<- rownames(seqtab.nochim_ITS)
samples.out
samples.out <- rownames(seqtab.nochim_ITS)
echantillon <- sapply(strsplit(samples.out, "_"), `[`, 1)
type <- substr(echantillon,1,1)
species <- substr(echantillon,3,6)
substrate <- substr(echantillon,8,11)
subtype <- substr(echantillon,13,15)
number <- substr(echantillon,17,18)
samdf <- data.frame(echantillon=echantillon, 
                    type=type, 
                    species=species, 
                    substrate=substrate,
                    subtype=subtype,
                    number=number)
rownames(samdf) <- samples.out

ps_ITS <- phyloseq(otu_table(seqtab.nochim_ITS,
                             taxa_are_rows=FALSE),
                   sample_data(samdf),
                   tax_table(taxa_ITS))
# Have a look on the data structure
ps_ITS
head(sample_data(ps_ITS), 10)

# Use short name for ASVs rather than the full DNA sequence
dna <- Biostrings::DNAStringSet(taxa_names(ps_ITS))
names(dna) <- taxa_names(ps_ITS)
ps_ITS <- merge_phyloseq(ps_ITS, dna)
taxa_names(ps_ITS) <- paste0("ASV", seq(ntaxa(ps_ITS)))

# Remove ambiguous phylum, chloroplast, and mitochondria
ps_ITS <- subset_taxa(ps_ITS, !is.na(phylum)
                      & !order %in% c("Chloroplast")
                      & !family %in% c("Mitochondria"))
ps_ITS

## Create table, number of features for each phyla
table(tax_table(ps_ITS)[, "family"], exclude = NULL)

# Exporter la table ASV dans un fichier Excel 
# write.csv2(tax_table(ps_ITS), "Taxa/ASV_table_ITS.csv",
#            row.names = TRUE)

# Check the number of reads per sample
hist(sample_sums(ps_ITS), main="Histogram: Read Counts", 
     xlab="Total Reads", border="black", col="darkgreen", 
     las=1, breaks=15)

# Check the number of reads per taxa
head(taxa_sums(ps_ITS), 50)


############### DECONTAMINATION ###################
## BACTERIA
# Ajouter a ps une variable pour distinguer les echantillons negatifs 
# (Control sample) des vrais echantillons (True sample)
sample_data(ps_16S)
sample_control <- c("True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample",
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample",
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample",
                    "True sample", "True sample",
                    "True sample", "True sample",
                    "True sample", "True sample",
                    "Control sample", "Control sample",
                    "Control sample", "Control sample", 
                    "Control sample", "Control sample")

sample_or_control <- factor(sample_control)
sample_data(ps_16S)$sample_or_control <- sample_or_control
sample_data(ps_16S) #Verification

# Inspect library sizes
df <- as.data.frame(sample_data(ps_16S)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_16S)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, 
                    color=sample_or_control)) + 
  geom_point()

# Ajouter la concentration ADN (ng/ul)
sample_data(ps_16S)
DNA_Concentration <- c(6.50, 6.30, 10.0, 17.1, 9.70,
                       8.50, 8.90, 18.7, 14.7, 6.30,
                       12.7, 8.00, 11.5, 17.6, 9.10, 
                       7.10, 18.3, 8.40, 6.80, 8.70, 
                       8.00, 8.50, 13.7, 19.9, 15.3, 
                       8.70, 8.60, 17.1, 22.1, 16.9, 
                       6.40, 8.80, 5.90, 10.6, 7.10, 
                       6.50)

# Specifier quel est le champ incluant la concentration ADN
sample_data(ps_16S)$DNA_Concentration <- DNA_Concentration
sample_data(ps_16S) #Verification

#Specifier quels sont les controles negatifs
sample_data(ps_16S)$is.neg <- sample_data(ps_16S)$sample_or_control == "Control sample"

## FREQUENCE AND PREVALENCE - Identification des contaminants
# Method = "either" permet de faire les deux methodes simultanement
set.seed(420)
contam <- isContaminant(ps_16S, method="either",
                        neg="is.neg",
                        threshold=0.25,
                        conc = "DNA_Concentration")

# Voir allure du resultat
head(contam, 10)
# Combien de contaminants
table(contam$contaminant)
# Voir les contaminants
which(contam$contaminant)

## Visualisation des contaminants
# Ne permet pas de retirer les contaminants!
# Frequence
plot_frequency(ps_16S, 
               taxa_names(ps_16S)[c(1:9)],
               conc="DNA_Concentration") + 
  xlab("DNA Concentration (Nanodrop)")

set.seed(420) #Inspection des AVSs "contaminants"
plot_frequency(ps_16S, 
               taxa_names(ps_16S)[sample(which(contam$contaminant),
                                         10)], 
               conc="DNA_Concentration") +
  xlab("DNA Concentration (Nanodrop)")

# Prevalence
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps_16S,
                                 function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_or_control == "Control sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_or_control == "True sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), 
                    pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contam$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, 
                       y=pa.pos, 
                       color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True samples)")

# RETIRER LES CONTAMINANTS
ps_16S
ps_16S <- prune_taxa(!contam$contaminant, ps_16S)
ps_16S


## FUNGI
# Ajouter a ps une variable pour distinguer les echantillons negatifs 
# (Control sample) des vrais echantillons (True sample)
sample_data(ps_ITS)
sample_control <- c("True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample",
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample",
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample", 
                    "True sample", "True sample",
                    "True sample", "True sample",
                    "True sample", "True sample",
                    "True sample", "True sample",
                    "Control sample", "Control sample")

sample_or_control <- factor(sample_control)
sample_data(ps_ITS)$sample_or_control <- sample_or_control
sample_data(ps_ITS) #Verification

# Inspect library sizes
df <- as.data.frame(sample_data(ps_ITS)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps_ITS)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, 
                    color=sample_or_control)) + 
  geom_point()

# Ajouter la concentration ADN (ng/ul)
sample_data(ps_ITS)
DNA_Concentration <- c(18.1, 19.9, 20.9, 26.4, 23.5,
                       22.8, 18.9, 9.30, 12.5, 19.2,
                       20.8, 16.8, 20.6, 20.5, 24.2, 
                       25.8, 25.3, 31.5, 24.1, 24.0,
                       25.8, 29.9, 26.3, 32.9, 24.8,
                       21.3, 25.9, 27.3, 21.9, 24.2,
                       6.80, 7.50)

# Specifier quel est le champ incluant la concentration ADN
sample_data(ps_ITS)$DNA_Concentration <- DNA_Concentration
sample_data(ps_ITS) #Verification

#Specifier quels sont les controles negatifs
sample_data(ps_ITS)$is.neg <- sample_data(ps_ITS)$sample_or_control == "Control sample"

## FREQUENCE AND PREVALENCE - Identification des contaminants
# Method = "either" permet de faire les deux methodes simultanement
set.seed(420)
contam <- isContaminant(ps_ITS, method="either",
                        neg="is.neg",
                        threshold=0.25,
                        conc = "DNA_Concentration")

# Voir allure du resultat
head(contam, 10)
# Combien de contaminants
table(contam$contaminant)
# Voir les contaminants
which(contam$contaminant)

## Visualisation des contaminants
# Ne permet pas de retirer les contaminants!
# Frequence
plot_frequency(ps_ITS, 
               taxa_names(ps_ITS)[c(1:9)],
               conc="DNA_Concentration") + 
  xlab("DNA Concentration (Nanodrop)")

set.seed(420) #Inspection des AVSs "contaminants"
plot_frequency(ps_ITS, 
               taxa_names(ps_ITS)[sample(which(contam$contaminant),
                                         9)], 
               conc="DNA_Concentration") +
  xlab("DNA Concentration (Nanodrop)")

# Prevalence
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps_ITS,
                                 function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_or_control == "Control sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_or_control == "True sample", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), 
                    pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contam$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, 
                       y=pa.pos, 
                       color=contaminant)) + 
  geom_point() +
  xlab("Prevalence (Negative Controls)") + 
  ylab("Prevalence (True samples)")

# RETIRER LES CONTAMINANTS
ps_ITS
ps_ITS <- prune_taxa(!contam$contaminant, ps_ITS)
ps_ITS

# write.csv2(otu_table(ps_ITS),"otu_table_ITS.csv")
# Inspect manually and remove leftover contaminants
# i.e. ASVs only present in negative controls
# ps_ITS <- subset_taxa(ps_ITS, "ASV613")


################# PREVALENCE ####################
## BACTERIA
# Subset excluding negative controls
ps2_16S <- subset_samples(ps_16S, subtype != "NAN")
ps2_16S
# Remove taxa not seen at least 5 times in two samples.
# This protects against an OTU with small mean &
# trivially large coefficient of variation
# It  was applied to the OTU counts prior to 
# creating the figures in the main phyloseq manuscript.
ps3_16S <- phyloseq::filter_taxa(ps2_16S, 
                                 function(x) sum(x > 25) > (0.02*length(x)),
                                 TRUE)
ps3_16S
# you can also subset_taxa to keep only taxonomic groups of interest
# GP.chl = subset_taxa(GlobalPatterns, Phylum=="Chlamydiae")
# Compute prevalence ONLY FUNGI SEEN AT LEAST 5 TIMES IN TWO SAMPLES 
prevdf = apply(X = otu_table(ps3_16S),
               MARGIN = ifelse(taxa_are_rows(ps3_16S), 
                               yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps3_16S),
                    tax_table(ps3_16S))
plyr::ddply(prevdf, "phylum", 
            function(df1){cbind(mean(df1$Prevalence),
                                sum(df1$Prevalence))})
# Subset to the remaining phyla
prevdf1 = subset(prevdf, 
                 phylum %in% get_taxa_unique(ps3_16S,
                                             "phylum"))
ggplot(prevdf1, aes(TotalAbundance, 
                    Prevalence / nsamples(ps3_16S),
                    color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, 
             linetype = 2) +  
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Reads") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")
# ONLY BACTERIA SEQUENCED AT LEAST 5 TIMES IN TWO SAMPLES

## FUNGI
# Subset excluding negative controls
ps2_ITS <- subset_samples(ps_ITS, subtype != "NAN")
ps2_ITS
# Remove taxa not seen at least 2 times in one samples.
# This protects against an OTU with small mean & trivially large 
# Coefficient of Variation
# It  was applied to the OTU counts prior to creating the figures 
# in the main phyloseq manuscript.
ps3_ITS <- phyloseq::filter_taxa(ps2_ITS, 
                             function(x) sum(x > 25) > (0.02*length(x)),
                             TRUE)
ps3_ITS

## Create table, number of features for each phyla
table(tax_table(ps3_ITS)[, "genus"], exclude = NULL)

# you can also subset_taxa to keep only taxonomic groups of interest
# GP.chl = subset_taxa(GlobalPatterns, Phylum=="Chlamydiae")
# Compute prevalence ONLY FUNGI SEEN AT LEAST 5 TIMES IN TWO SAMPLES 
prevdf = apply(X = otu_table(ps3_ITS),
               MARGIN = ifelse(taxa_are_rows(ps3_ITS), 
                               yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps3_ITS),
                    tax_table(ps3_ITS))
plyr::ddply(prevdf, "phylum", 
            function(df1){cbind(mean(df1$Prevalence),
                                sum(df1$Prevalence))})
# Subset to the remaining phyla
prevdf1 = subset(prevdf, 
                 phylum %in% get_taxa_unique(ps3_ITS,
                                             "phylum"))
ggplot(prevdf1, aes(TotalAbundance, 
                    Prevalence / nsamples(ps3_ITS),
                    color=phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, 
             linetype = 2) +  
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Reads") + 
  ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~phylum) + theme(legend.position="none")
# ONLY FUNGI SEQUENCED AT LEAST 5 TIMES IN TWO SAMPLES


############## NORMALISATION DESEQ2 #######################

## BACTERIA
# AS DESeq2 sends back negative value for raw count = 0
# I will replace manually all 0 by 1 in the tax_table
ps3_16S <- transform_sample_counts(ps3_16S,
                                   function(x) {x + 1})
otu_table(ps3_16S)
# Stabilisation en fonction du type echantillon
ps4_16S <- phyloseq_to_deseq2(ps3_16S, ~type)
ps4_16S <- estimateSizeFactors(ps4_16S)
ps4_16S <- estimateDispersions(ps4_16S, 
                               fitType = "mean")
ps5_16S <- DESeq2::getVarianceStabilizedData(ps4_16S)

# Save the untransformed data as a separate variable so you can go back to it
ps3.1_16S <- ps3_16S

# Formating the DESeq object back to be used in standard phyloseq extensions
# Replace the OTU table in ps3 by the variance stabilized counts of ps6
otu_table(ps3_16S) <- otu_table(ps5_16S, 
                                taxa_are_rows = TRUE)
head(otu_table(ps3_16S),10)
# AS DESeq2 sends back negative value for raw count = 0 (see above)
# I will replace manually all < 1 (= true zero in original) by 0
# Extract the OTU table
otu <- otu_table(ps3_16S)
# Replace values < 1 with 0
otu[otu < 1] <- 0
# Reassign the modified OTU table to the phyloseq object
otu_table(ps3_16S) <- otu_table(otu, 
                                taxa_are_rows = taxa_are_rows(otu))
head(otu_table(ps3_16S),10)

# Abondance relative de chaque ASV dans chaque echantillon !!!
# This modified ps3 variable now has the results of DESeq2 variance-stabilization 
# of counts instead of the original counts.
# Check the number of normalized reads per sample
sample_sums(ps3_16S)
hist(sample_sums(ps3_16S), 
     main="Histogram: DESeq2 Variance-stabilization", 
     xlab="Variance-stabilization of reads", 
     border="black", col="darkgreen", 
     las=1, breaks=15)

# Make sure that no ASV was lost in the variance-stabilization process
table(tax_table(ps3.1_16S)[, "phylum"], exclude = NULL)
table(tax_table(ps3_16S)[, "phylum"], exclude = NULL)
# Exact same number of ASV per phylum after normalisation of data
# Now you can use object ps3 for downstream analysis

## FUNGI
# AS DESeq2 sends back negative value for raw count = 0
# I will replace manually all 0 by 1 in the tax_table
ps3_ITS <- transform_sample_counts(ps3_ITS,
                                   function(x) {x + 1})
otu_table(ps3_ITS)
# Stabilisation en fonction du type echantillon
ps4_ITS <- phyloseq_to_deseq2(ps3_ITS, ~type)
ps4_ITS <- estimateSizeFactors(ps4_ITS)
ps4_ITS <- estimateDispersions(ps4_ITS, fitType = "mean")
ps5_ITS <- DESeq2::getVarianceStabilizedData(ps4_ITS)

# Save the untransformed data as a separate variable so you can go back to it
ps3.1_ITS <- ps3_ITS

# Formating the DESeq object back to be used in standard phyloseq extensions
# Replace the OTU table in ps3 by the variance stabilized counts of ps6
otu_table(ps3_ITS) <- otu_table(ps5_ITS, 
                                taxa_are_rows = TRUE)

head(otu_table(ps3_ITS),10)
# AS DESeq2 sends back negative value for raw count = 0 (see above)
# I will replace manually all < 1 (= true zero in original) by 0
# Extract the OTU table
otu <- otu_table(ps3_ITS)
# Replace values < 1 with 0
otu[otu < 1] <- 0
# Reassign the modified OTU table to the phyloseq object
otu_table(ps3_ITS) <- otu_table(otu, taxa_are_rows = taxa_are_rows(otu))
head(otu_table(ps3_ITS),10)
# Abondance relative de chaque ASV dans chaque echantillon !!!
# This modified ps3 variable now has the results of DESeq2 variance-stabilization 
# of counts instead of the original counts.
# Check the number of normalized reads per sample
sample_sums(ps3_ITS)
hist(sample_sums(ps3_ITS), 
     main="Histogram: DESeq2 Variance-stabilization", 
     xlab="Variance-stabilization of reads", 
     border="black", col="darkgreen", 
     las=1, breaks=15)

# Make sure that no ASV was lost in the variance-stabilization process
table(tax_table(ps3.1_ITS)[, "phylum"], exclude = NULL)
table(tax_table(ps3_ITS)[, "phylum"], exclude = NULL)
# Exact same number of ASV per phylum after normalisation of data
# Now you can use object ps3 for downstream analysis


################### FUNGUILD & FAPROTAX #######################
## FUNGUILD (fungi)
ps3_ITS <- add_funguild_info(ps3_ITS,
                             taxLevels = c("domain", "phylum",
                                           "class", "order",
                                           "family", "genus",
                                           "species"))

## FAPROTAX (bacteria)
# Necessary to change headers for it to work
tax <- tax_table(ps3_16S)
tax_matrix <- as.matrix(tax)
colnames(tax_matrix) <- c("Kingdom", "Phylum", 
                          "Class", "Order", "Family",
                          "Genus", "Species")
new_tax_table <- tax_table(tax_matrix)
tax_table(ps3_16S) <- new_tax_table
# Run Faprotax
meco_obj <- phyloseq2meco(ps3_16S)
func_obj <- trans_func$new(meco_obj)
func_obj$for_what <- "prok"
func_obj$cal_spe_func(prok_database = "FAPROTAX")
func_obj$cal_spe_func_perc(abundance_weighted = TRUE)
func_obj$trans_spe_func_perc()
# Export results
#write.csv(func_obj$res_spe_func_perc_trans, 
#          "faprotax_results.csv", row.names = FALSE)


################### EXPORT MATRICES #######################
# Export bacteria data
# write.csv2(tax_table(ps3_16S), "bacteria_taxa.csv")
# write.csv(otu_table(ps3_16S), "bacteria_abundance.csv")
# Export fungi data
# write.csv2(tax_table(ps3_ITS), "fungi_taxa.csv")
# write.csv2(otu_table(ps3_ITS), "fungi_abundance.csv")


################## METABOLOME IMPORT ################
# Import phyloseq object from metabolome dataset
## QUALITY CONTROL WITH EDITED DATASET SIRIUS_w_METADATA
# Load data
MET <- read.csv(file = "metabolome_phyloseq_import.csv",
                header = TRUE, row.names = 1)

## Creation du dataframe avec les donnees des echantillons
data.out <- rownames(MET)
#Creation du df avec les donnees des echantillons
echantillon <- sapply(strsplit(data.out, "_"), `[`, 1)
type <- substr(echantillon,1,1)
species <- substr(echantillon,3,6)
substrate <- substr(echantillon,8,11)
subtype <- substr(echantillon,13,15)
number <- substr(echantillon,17,18)
samdf <- data.frame(echantillon=echantillon, 
                    type=type, 
                    species=species, 
                    substrate=substrate,
                    subtype=subtype,
                    number=number)
rownames(samdf) <- data.out
# Verification
samdf # OK

## Creation objet phyloseq pour analyses de composition
ps_MET <- phyloseq(otu_table(MET, taxa_are_rows=FALSE),
                   sample_data(samdf))
ps_MET

# Add taxa table with annotations into the phyloseq object
import_tax_table_ps_MET <- read.csv("metabolome_taxa.csv", 
                                    row.names = 1)
all(taxa_names(ps_MET) == rownames(import_tax_table_ps_MET))
tax_table(ps_MET) <- as.matrix(import_tax_table_ps_MET)
ps_MET

# Must transpose OTU_table to have taxa as rows
ps_MET_transposed <- ps_MET %>%
  otu_table() %>%
  t() %>%
  otu_table(taxa_are_rows = FALSE) %>%
  phyloseq(tax_table(ps_MET), sample_data(ps_MET))

## DESEQ2
# AS DESeq2 sends back negative value for raw count = 0
# I will replace manually all 0 by 1 in the tax_table
ps2_MET <- transform_sample_counts(ps_MET_transposed,
                                   function(x) {x + 1})
otu_table(ps2_MET) <- as.matrix(otu_table(ps2_MET))

# Check for NA values
any(is.na(otu_table(ps2_MET))) # OK
# Ensure all values in the count matrix are integers. 
# If they're not, round them:
otu_table(ps2_MET) <- round(otu_table(ps2_MET))
# Check for non-integer values:
any(!is.integer(otu_table(ps2_MET)))
# Verify that no values exceed the maximum integer limit:
any(otu_table(ps2_MET) > .Machine$integer.max)
# If large values are present, consider scaling down your data:
otu_table(ps2_MET) <- floor(otu_table(ps2_MET) / 10)
# Verify that no values exceed the maximum integer limit:
any(otu_table(ps2_MET) > .Machine$integer.max)

# Stabilisation en fonction du type echantillon
ps3_MET <- phyloseq_to_deseq2(ps2_MET, ~type)
ps3_MET <- estimateSizeFactors(ps3_MET)
ps3_MET <- estimateDispersions(ps3_MET, fitType = "mean")
ps4_MET <- DESeq2::getVarianceStabilizedData(ps3_MET)

# Save the untransformed data as a separate variable so you can go back to it
ps2.1_MET <- ps2_MET

# Formating the DESeq object back to be used 
# in standard phyloseq extensions
# Replace the OTU table by the variance stabilized counts
otu_table(ps2_MET) <- otu_table(ps4_MET, taxa_are_rows = TRUE)

head(otu_table(ps2_MET),10)
# AS DESeq2 sends back negative value for raw count = 0 (see above)
# I will replace manually all < 1 (= true zero in original) by 0
# Extract the OTU table
otu <- otu_table(ps2_MET)
# Replace values < 1 with 0
otu[otu < 1] <- 0
# Reassign the modified OTU table to the phyloseq object
otu_table(ps2_MET) <- otu_table(otu, 
                                taxa_are_rows = taxa_are_rows(otu))
head(otu_table(ps2_MET),10)
# Abondance relative de chaque ASV dans chaque echantillon !!!
# This modified ps3 variable now has the results of DESeq2 variance-stabilization 
# of counts instead of the original counts.
# Check the number of normalized reads per sample
sample_sums(ps2_MET)
hist(sample_sums(ps2_MET), 
     main="Histogram: DESeq2 Variance-stabilization", 
     xlab="Variance-stabilization of reads", 
     border="black", col="darkgreen", 
     las=1, breaks=15)



################# ACCUMULATION CURVES ##################
# Import csv files formatted for R
fungi <- read.csv("fungi_presence-absence.csv")
bacteria <- read.csv("bacteria_presence-absence.csv")
metabolome <- read.csv("metabolome_presence-absence.csv")

# Make sure they are considered dataframes to avoid problems later
fungi <- as.data.frame(fungi)
bacteria <- as.data.frame(bacteria)
metabolome <- as.data.frame(metabolome)

### Fungi
acc.ITS <- specaccum(fungi, method = "random", 
                     permutations = 1000)
plot(acc.ITS, xlab = "Number of samples", 
     ylab = "Number of fungi ASVs",
     col = "black", ylim = c(0,400), xlim = c(0,30), 
     mgp=c(2,1,0))

### Bacteria
acc.16S <- specaccum(bacteria, method = "random", 
                     permutations = 1000)
plot(acc.16S, xlab = "Number of samples", 
     ylab = "Number of bacteria ASVs",
     col = "black", ylim = c(0,2500), xlim = c(0,30), 
     mgp=c(2,1,0))

### Metabolome
acc.MET <- specaccum(metabolome, method = "random", 
                     permutations = 1000)
plot(acc.MET, xlab = "Number of samples", 
     ylab = "Number of metabolites",
     col = "black", ylim = c(0,7000), xlim = c(0,30), 
     mgp=c(2,1,0))

# Construct and export figure
#jpeg("Figure_A1_Accumulation_curves.jpeg", units="in",
#     width=6, height=3, res=600)
# Define margin size
par(mfrow = c(1, 3))
# Fungi
par(mar = c(3.0, 3.0, 0.1, 0.5), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
plot(acc.ITS, xlab = "", 
     ylab = "",
     col = "black", ylim = c(0,400),
     xlim = c(0,30), cex.axis = 0.6,
     cex.lab = 0.8)
title(ylab = "Number of features", mgp = c(1.5, 1, 0))    # Add y-axis text
text(x = 15, y = 10, "Fungi ASVs", cex = 0.8, lwd = 1.5)

# Bacteria
par(mar = c(3.0, 1.75, 0.1, 1.75), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
plot(acc.16S, xlab = "", 
     ylab = "",
     col = "black", ylim = c(0,2500), 
     xlim = c(0,30), cex.axis = 0.6,
     cex.lab = 0.8)
title(xlab = "Number of samples", mgp = c(1.5, 1, 0))
text(x = 15, y = 62.5, "Bacteria ASVs", cex = 0.8, lwd = 1.5)

# Metabolites
par(mar = c(3.0, 0.5, 0.1, 3.0), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
plot(acc.MET, xlab = "", 
     ylab = "",
     col = "black", ylim = c(0,7000), 
     xlim = c(0,30), cex.axis = 0.6,
     cex.lab = 0.8)
text(x = 15, y = 175, "Metabolites", cex = 0.8, lwd = 1.5)
#dev.off()



##################### RICHNESS #########################
# Load data
richness <- read.csv("richness_calicioids.csv")
richness <- as.data.frame(richness)

# Variables as objects
TYPE <- richness$Type
SUBSTRATE <- richness$Substrate
SAMPLE <- richness$Sample
BACTERIA <- richness$Bacteria
FUNGI <- richness$Fungi
METABOLITES <- richness$Metabolites

## MAMAZ - Janvier 2025
# Effectivement, c'est mieux de considérer un seul facteur 
# qui désigne le traitement (culture-agar, fungus-wood, 
# lichen-bark, lichen-wood), puisque tu as un gros 
# débalancement et les deux facteurs ne sont pas
# parfaitement croisés. Si tu utilises lm( ) pour voir 
# les estimés du modèle, tu vas voir qu'il y a un paquet 
# de NA's comme valeur des estimations qui est une 
# indication qu'il y a un problème (débalancement et 
# certaines combinaisons de facteurs jamais observées). 
# Même si R arrive à le "gérer", c'est plus clean 
# et correct d'utiliser un seul facteur à 4 niveaux.

# Si tu appelle ce nouveau facteur "inter",
richness$inter <- paste(richness$Type, ".", 
                        richness$Substrate, sep = "")
INTER <- richness$inter
# Graph parameters
par(mfrow = c(1, 2))

### BACTERIA
## See data
hist(BACTERIA)
boxplot(BACTERIA~INTER,las=3)
# Mean and standard deviation
mean(BACTERIA[INTER == "Culture.Agar"])
sd(BACTERIA[INTER == "Culture.Agar"])
mean(BACTERIA[INTER == "Fungus.Wood"])
sd(BACTERIA[INTER == "Fungus.Wood"])
mean(BACTERIA[INTER == "Lichen.Bark"])
sd(BACTERIA[INTER == "Lichen.Bark"])
mean(BACTERIA[INTER == "Lichen.Wood"])
sd(BACTERIA[INTER == "Lichen.Wood"])
## Modèle linéaire
m1 <- lm(Bacteria ~ INTER, data = richness)
# Check assumptions with residual diagnostics
# Homogeneity of variances
plot(residuals(m1) ~ fitted(m1)) # Graph
levene.test(BACTERIA, INTER) # Levene
fligner.test(BACTERIA ~ INTER) # Fligner
bartlett.test(BACTERIA ~ INTER) # Bartlett
# Normality of residuals
qqnorm(residuals(m1, type = "pearson")) # Graph
qqline(residuals(m1, type = "pearson")) # Line
ad.test(residuals(m1)) # Anderson Darling
shapiro.test(residuals(m1)) # Shapiro Wilk
# Assumptions OK
# See results
summary(m1)
## Confidence intervals
# Data for predictions
predDat <- expand.grid(INTER = factor(c("Culture.Agar",
                                        "Fungus.Wood",
                                        "Lichen.Bark",
                                        "Lichen.Wood")))
# Predictions
predict(m1, newdata = predDat,
        interval = "confidence",
        level = 0.95)
predDat


### FUNGI
## See data
hist(FUNGI)
boxplot(FUNGI~INTER,las=3)
# Mean and standard deviation
mean(FUNGI[INTER == "Culture.Agar"])
sd(FUNGI[INTER == "Culture.Agar"])
mean(FUNGI[INTER == "Fungus.Wood"])
sd(FUNGI[INTER == "Fungus.Wood"])
mean(FUNGI[INTER == "Lichen.Bark"])
sd(FUNGI[INTER == "Lichen.Bark"])
mean(FUNGI[INTER == "Lichen.Wood"])
sd(FUNGI[INTER == "Lichen.Wood"])
## Linear model
m2 <- lm(FUNGI ~ INTER, data = richness)
# Check assumptions with residual diagnostics
# Homogeneity of variances
plot(residuals(m2) ~ fitted(m2)) # Graph
levene.test(FUNGI, INTER) # Levene
fligner.test(FUNGI ~ INTER) # Fligner
bartlett.test(FUNGI ~ INTER) # Bartlett
# Normality of residuals
qqnorm(residuals(m2, type = "pearson")) # Graph
qqline(residuals(m2, type = "pearson")) # Line
ad.test(residuals(m2)) # Anderson Darling
shapiro.test(residuals(m2)) # Shapiro Wilk
# Assumptions NOT OK: multivariances heterogeneous
# Transformation LOG
m2 <- lm(log(FUNGI) ~ INTER, data = richness)
# Check assumptions with residual diagnostics
# Homogeneity of variances
plot(residuals(m2) ~ fitted(m2)) # Graph
levene.test(log(FUNGI), INTER) # Levene
fligner.test(log(FUNGI) ~ INTER) # Fligner
bartlett.test(log(FUNGI) ~ INTER) # Bartlett
# Normality of residuals
qqnorm(residuals(m2, type = "pearson")) # Graph
qqline(residuals(m2, type = "pearson")) # Line
ad.test(residuals(m2)) # Anderson Darling
shapiro.test(residuals(m2)) # Shapiro Wilk
# Assumptions OK, but confidence intervals are weird
## GENERALIZED LEAST SQUARE MODEL FOR HETEROGENEOUS VARIANCES
# Model heterogeneity explicitly
# Add term to allow variances to differ among groups
glsFUNGI <- gls(FUNGI ~ INTER, data = richness,
                      weights = varIdent(form = ~ 1 | INTER))

# Check homogeneity of variances
plot(residuals(glsFUNGI, 
               type = "pearson") ~ fitted(glsFUNGI))
# MUCH BETTER
# Check normality of residuals
qqnorm(residuals(glsFUNGI, type = "pearson"))
qqline(residuals(glsFUNGI, type = "pearson"))
# OK
# Check results
summary(glsFUNGI) # Beta estimates
anova(glsFUNGI)  # ANOVA table
# Post-hoc multiple comparisons
# Necessary to convert to factors fo post-hoc
INTER2 <- as.factor(INTER)
# Rerun model heterogeneity explicitly with factors
glsFUNGI2 <- gls(FUNGI ~ INTER2,
                 weights = varIdent(form = ~ 1 | INTER2))
# Run the following lines. These introduce methods for 'gls' objects.
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)}
# The line below gives pairwise comparisons now.
confint(glht(glsFUNGI2, 
             linfct=mcp(INTER2="Tukey")))
## Confidence intervals
# Create dataframe
full <- data.frame(FUNGI = FUNGI,
                   INTER = INTER)
# Pour un GLS, tu peux aller le chercher 
# avec une fonction de mon package:
preds <- predictSE(glsFUNGI, newdata = predDat)
DF <- nrow(full) - attr(logLik(glsFUNGI), "df")
low95 <- preds$fit - qt(df = DF, p = 0.975) * preds$se.fit
upp95 <- preds$fit + qt(df = DF, p = 0.975) * preds$se.fit
preds
low95
upp95
predDat


### METABOLITES
## See data
hist(METABOLITES)
boxplot(METABOLITES~INTER,las=3)
# Mean and standard deviation
mean(METABOLITES[INTER == "Culture.Agar"])
sd(METABOLITES[INTER == "Culture.Agar"])
mean(METABOLITES[INTER == "Fungus.Wood"])
sd(METABOLITES[INTER == "Fungus.Wood"])
mean(METABOLITES[INTER == "Lichen.Bark"])
sd(METABOLITES[INTER == "Lichen.Bark"])
mean(METABOLITES[INTER == "Lichen.Wood"])
sd(METABOLITES[INTER == "Lichen.Wood"])
## Linear model
m3 <- lm(METABOLITES ~ INTER, data = richness)
# Check assumptions with residual diagnostics
# Homogeneity of variances
plot(residuals(m3) ~ fitted(m3)) # Graph
levene.test(METABOLITES, INTER) # Levene
fligner.test(METABOLITES ~ INTER) # Fligner
bartlett.test(METABOLITES ~ INTER) # Bartlett
# Normality of residuals
qqnorm(residuals(m3, type = "pearson")) # Graph
qqline(residuals(m3, type = "pearson")) # Line
ad.test(residuals(m3)) # Anderson Darling
shapiro.test(residuals(m3)) # Shapiro Wilk
# Assumptions NOT OK: multivariances heterogeneous
# Transformation LOG
m3 <- lm(log(METABOLITES) ~ INTER, data = richness)
# Check assumptions with residual diagnostics
# Homogeneity of variances
plot(residuals(m3) ~ fitted(m3)) # Graph
levene.test(log(METABOLITES), INTER) # Levene
fligner.test(log(METABOLITES) ~ INTER) # Fligner
bartlett.test(log(METABOLITES) ~ INTER) # Bartlett
# Normality of residuals
qqnorm(residuals(m3, type = "pearson")) # Graph
qqline(residuals(m3, type = "pearson")) # Line
ad.test(residuals(m3)) # Anderson Darling
shapiro.test(residuals(m3)) # Shapiro Wilk
# Assumptions NOT OK
## GENERALIZED LEAST SQUARE MODEL FOR HETEROGENEOUS VARIANCES
# Model heterogeneity explicitly
# Add term to allow variances to differ among groups
glsMETABOLITES <- gls(METABOLITES ~ INTER, data = richness,
                      weights = varIdent(form = ~ 1 | INTER))

# Check homogeneity of variances
plot(residuals(glsMETABOLITES, 
               type = "pearson") ~ fitted(glsMETABOLITES))
# MUCH BETTER
# Check normality of residuals
qqnorm(residuals(glsMETABOLITES, type = "pearson"))
qqline(residuals(glsMETABOLITES, type = "pearson"))
# OK
# Check results
summary(glsMETABOLITES) # Beta estimates
anova(glsMETABOLITES)  # ANOVA table
# Post-hoc multiple comparisons
# Necessary to convert to factors fo post-hoc
INTER2 <- as.factor(INTER)
# Rerun model heterogeneity explicitly with factors
glsMETABOLITES2 <- gls(METABOLITES ~ INTER2,
                 weights = varIdent(form = ~ 1 | INTER2))
# Run the following lines. These introduce methods for 'gls' objects.
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)}
# The line below gives pairwise comparisons now.
confint(glht(glsMETABOLITES2, 
             linfct=mcp(INTER2="Tukey")))
## Confidence intervals
# Create dataframe
full <- data.frame(METABOLITES = METABOLITES,
                   INTER = INTER)
# Pour un GLS, tu peux aller le chercher 
# avec une fonction de mon package:
preds <- predictSE(glsMETABOLITES, newdata = predDat)
DF <- nrow(full) - attr(logLik(glsMETABOLITES), "df")
low95 <- preds$fit - qt(df = DF, p = 0.975) * preds$se.fit
upp95 <- preds$fit + qt(df = DF, p = 0.975) * preds$se.fit
preds
low95
upp95
predDat


##### CREATE GRAPH FOR PUBLICATION
library("beeswarm")

# Importer le csv 
colnames(richness) <- c("Type", "Substrate", "Number",
                        "Bacteria", "Fungi", 
                        "Metabolites", "Treatment")
richness
# Format labels, add \n where appropriate.
my.labels1 <- c("Cultures\nAgar", "Fungi\nWood",
                "Lichens\nBark", "Lichens\nWood") 
# Create variable for SUBSTRATE
z <- as.numeric(factor(richness$Treatment), replace = TRUE)
# Create colour variable for Group
my.pwcol1 <- c("red","red","red","red","red",
               "red","red","red","red","red",
               "red","red","red","red",
               "blue","blue","blue",
               "black","black","black","black",
               "blue","black","blue",
               "black","black","black","black","black",
               "blue")
# Create point shape variable for Group
my.pwpch1 <- c(17,17,17,17,17,17,17,17,17,17,
               17,17,17,17,15,15,15,19,19,19,
               19,15,19,15,19,19,19,19,19,15)

# 95% Confidence intervals and means
# Bacteria
MEAN_Bacteria <- c(159.57, 171.00, 157.40, 214.00)
CI95_Bacteria_UP <- c(200.10, 258.57, 205.36, 301.57)
CI95_Bacteria_LOW <- c(119.03, 83.43, 109.44, 126.43)
# Fungi
MEAN_Fungi <- c(11.64, 38.33, 29.60, 44.67)
CI95_Fungi_UP <- c(16.07, 58.74, 37.54, 72.42)
CI95_Fungi_LOW <- c(7.22, 17.93, 21.66, 16.91)
# Metabolites
MEAN_Metabolites <- c(3500.79, 5225.67, 5468.60, 5271.33)
CI95_Metabolites_UP <- c(3727.89, 5429.91, 5525.92, 5344.92)
CI95_Metabolites_LOW <- c(3273.69, 5021.42, 5411.28, 5197.74)

# Construct and export figure
#jpeg("Figure_Richness.jpeg", units="in",
#      width=6, height=3, res=600)

# Define margin size
par(mfrow = c(1, 3))

# Plot Fungi
par(mar = c(3.0, 3.0, 0.1, 0.5), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
beeswarm(richness$Fungi ~ richness$Treatment, 
         xlab ="", ylab ="",
         cex = 0.6, ylim = c(0,80),
         pwpch = my.pwpch1,
         pwcol = my.pwcol1, cex.axis = 0.7, 
         cex.main = 0.8, labels = my.labels1)
title(ylab = "Richness", mgp = c(1.5, 1, 0))    # Add y-axis text
abline(v = 1.5, lty = 2, col = "gray35")
text(x = 2.5, y = 80, "Fungi", cex = 0.8, lwd = 1.5)
text(x = 1, y = 77.33, "A", cex = 0.8, lwd = 1.5)
text(x = 3, y = 77.33, "B", cex = 0.8, lwd = 1.5)
# Add standard error and mean
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Fungi_UP,
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Fungi,
       code=3, angle=90, length=0.04, 
       lwd = 1.0, col = "gray15")
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Fungi_LOW, 
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Fungi,
       code=3, angle=90, length=0.05, 
       lwd = 1.0, col = "gray15")

# Plot Bacteria
par(mar = c(3.0, 1.75, 0.1, 1.75), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
beeswarm(richness$Bacteria ~ richness$Treatment, 
         xlab ="", ylab ="",
         cex = 0.6, ylim = c(0,350),
         pwpch = my.pwpch1,
         pwcol = my.pwcol1, cex.axis = 0.7, 
         cex.main = 0.8, labels = my.labels1)
title(xlab = "Tissue origin", mgp = c(1.75, 3, 0)) 
text(x = 2.5, y = 350, "Bacteria", cex = 0.8, lwd = 1.5)
# Add standard error and mean
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Bacteria_UP,
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Bacteria,
       code=3, angle=90, length=0.04, 
       lwd = 1.0, col = "gray15")
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Bacteria_LOW, 
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Bacteria,
       code=3, angle=90, length=0.05, 
       lwd = 1.0, col = "gray15")

# Plot Metabolites
par(mar = c(3.0, 0.5, 0.1, 3.0), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
beeswarm(richness$Metabolites ~ richness$Treatment, 
         xlab ="", ylab ="",
         cex = 0.6, ylim = c(2000,6000),
         pwpch = my.pwpch1,
         pwcol = my.pwcol1, cex.axis = 0.7, 
         cex.main = 0.8, labels = my.labels1)
abline(v = 1.5, lty = 2, col = "gray35")
text(x = 2.5, y = 6000, "Metabolites", cex = 0.8, lwd = 1.5)
text(x = 1, y = 5800, "A", cex = 0.8, lwd = 1.5)
text(x = 3, y = 5800, "B", cex = 0.8, lwd = 1.5)
# Add standard error and mean
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Metabolites_UP,
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Metabolites,
       code=3, angle=90, length=0.04, 
       lwd = 1.0, col = "gray15")
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Metabolites_LOW, 
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Metabolites,
       code=3, angle=90, length=0.05, 
       lwd = 1.0, col = "gray15")


# Add legend
par(xpd = TRUE, mar = c(0, 0, 0, 0))
legend("right", inset=c(-0.25,0), title = "Substrate", 
       legend = c("Agar", "Bark", "Wood"),
       bty = "n", # remove box around legend 
       col = c("red","black","blue"), 
       pch = c(17,19,15), 
       cex = 0.7, horiz = FALSE)
#dev.off()


###################### VENN DIAGRAMS #######################
## Fungi
Venn_ITS <- read.csv("Venn_diagram_ITS.csv",
                     header = TRUE)
Culture_ITS <- (Venn_ITS$Culture[1:79])
Fungus_ITS <- (Venn_ITS$Fungus[1:92])
Lichen_ITS <- (Venn_ITS$Lichen[1:178])
# Graph
#venn.diagram(x = list(Culture_ITS, Fungus_ITS, Lichen_ITS),
#             filename = "Venn_diagram_ITS.jpg",
#             category.names = c("Cultures", "Fungi", "Lichens"),
#             height = 4000, width = 4000, resolution = 600,
#             lwd = 2, cex = 1.5,
#             cat.cex = 1.5, cat.default.pos = "outer",
#             cat.pos = c(-27, 27, 180),
#             cat.dist = c(0.05, 0.05, 0.03))

## Bacteria
Venn_16S <- read.csv("Venn_diagram_16S.csv",
                     header = TRUE)
Culture_16S <- (Venn_16S$Culture[1:1190])
Fungus_16S <- (Venn_16S$Fungus[1:437])
Lichen_16S <- (Venn_16S$Lichen[1:1330])
# Graph
#venn.diagram(x = list(Culture_16S, Fungus_16S, Lichen_16S),
#             filename = "Venn_diagram_16S.jpg",
#             category.names = c("Cultures", "Fungi", "Lichens"),
#             height = 4000, width = 4000, resolution = 600,
#             lwd = 2, cex = 1.5,
#             cat.cex = 1.5, cat.default.pos = "outer",
#             cat.pos = c(-27, 27, 180),
#             cat.dist = c(0.05, 0.05, 0.03))

## Metabolites
Venn_MET <- read.csv("Venn_diagram_MET.csv",
                     header = TRUE)
Culture_MET <- (Venn_MET$Culture[1:5704])
Fungus_MET <- (Venn_MET$Fungus[1:5700])
Lichen_MET <- (Venn_MET$Lichen[1:6065])
# Graph
#venn.diagram(x = list(Culture_MET, Fungus_MET, Lichen_MET),
#             filename = "Venn_diagram_MET.jpg",
#             category.names = c("Cultures", "Fungi", "Lichens"),
#             height = 4000, width = 4000, resolution = 600,
#             lwd = 2, cex = 1.5,
#             cat.cex = 1.5, cat.default.pos = "outer",
#             cat.pos = c(-27, 27, 180),
#             cat.dist = c(0.05, 0.05, 0.03))




################# DIVERSITY LINEAR REGRESSION #################
# Must transpose taxa as rows to taxa as columns in ps3
# FUNGI
ps3_ITS_transposed <- ps3_ITS %>%
  otu_table() %>%
  t() %>%
  otu_table(taxa_are_rows = FALSE) %>%
  phyloseq(tax_table(ps3_ITS), sample_data(ps3_ITS))
# Bacteria
ps3_16S_transposed <- ps3_16S %>%
  otu_table() %>%
  t() %>%
  otu_table(taxa_are_rows = FALSE) %>%
  phyloseq(tax_table(ps3_16S), sample_data(ps3_16S))

## Now compute diversity indexes
## Bacteria
H_16S <- diversity(otu_table(ps3_16S_transposed),
                   index = "shannon")
D_16S <- diversity(otu_table(ps3_16S_transposed),
                   index = "simpson")
I_16S <- diversity(otu_table(ps3_16S_transposed),
                   index = "invsimpson")
## FUNGI
H_ITS <- diversity(otu_table(ps3_ITS_transposed),
                   index = "shannon")
D_ITS <- diversity(otu_table(ps3_ITS_transposed),
                   index = "simpson")
I_ITS <- diversity(otu_table(ps3_ITS_transposed),
                   index = "invsimpson")
## METABOLITES
H_MET <- diversity(otu_table(ps_MET),
                   index = "shannon")
D_MET <- diversity(otu_table(ps_MET),
                   index = "simpson")
I_MET <- diversity(otu_table(ps_MET),
                   index = "invsimpson")

# Correlation coefficient
# Bacteria vs metabolites
cor1 <- cor.test(H_16S, H_MET, method = "pearson")
print(cor1)
cor2 <- cor.test(D_16S, D_MET, method = "pearson")
print(cor2)
cor3 <- cor.test(I_16S, I_MET, method = "pearson")
print(cor3)

# Fungi vs metabolites
cor4 <- cor.test(H_ITS, H_MET, method = "pearson")
print(cor4)
cor5 <- cor.test(D_ITS, D_MET, method = "pearson")
print(cor5)
cor6 <- cor.test(I_ITS, I_MET, method = "pearson")
print(cor6)

## Modeles lineaires
lm1 <- lm(H_MET ~ H_16S)
summary(lm1)
lm1
confint(lm1, level = 0.95)
lm2 <- lm(H_MET ~ H_ITS)
summary(lm2)
lm2
intervals(object = lm2)
confint(lm2, level = 0.95)
lm3 <- lm(H_16S ~ H_ITS)
summary(lm3)

## Representation graphique
Tissue <- c("Culture","Culture","Culture",
          "Culture","Culture","Culture",
          "Culture","Culture","Culture",
          "Culture","Culture","Culture",
          "Culture","Culture",
          "Fungus","Fungus","Fungus",
          "Lichen","Lichen","Lichen",
          "Lichen","Lichen","Lichen",
          "Lichen","Lichen","Lichen",
          "Lichen","Lichen","Lichen",
          "Lichen")
Substrate <- c("Agar","Agar","Agar",
               "Agar","Agar","Agar",
               "Agar","Agar","Agar",
               "Agar","Agar","Agar",
               "Agar","Agar",
               "Wood","Wood","Wood",
               "Bark","Bark","Bark",
               "Bark","Wood","Bark",
               "Wood","Bark","Bark",
               "Bark","Bark","Bark",
               "Wood")
## Graphs
# Fungi vs metabolites - Shannon
#jpeg("Figure_linear_regression_fungi_metabolome_diversity.jpeg", 
#     units="in", width=6, height=3, res=600)
plot_data <- data.frame(x = H_MET, 
                        y = H_ITS)
ggplot(plot_data, aes(x = x, y = y)) +
  geom_smooth(method = "lm", color = "red") +
  geom_point(aes(colour = Tissue, shape = Substrate),
             alpha = 1, size = 3) +
  geom_point(colour = "black", size = 1) +
  scale_shape_manual(values = c(15,16,17))+
  scale_color_manual(values = c("darkgreen", 
                                "firebrick", 
                                "mediumblue"))+
  labs(x = "Metabolome Diversity (H')", 
       y = "Fungi Diversity (H')",
       title = "Shannon diversity scatter plot",
       subtitle = paste("R-squared = 0.306, ", 
                        "p-value < 0.001")) +
  theme_minimal()
#dev.off()

# Bacteria vs metabolites - Shannon
#jpeg("Figure_linear_regression_bacteria_metabolome_diversity.jpeg", 
#     units="in", width=6, height=3, res=600)
plot_data <- data.frame(x = H_MET, 
                        y = H_16S)
ggplot(plot_data, aes(x = x, y = y)) +
  geom_smooth(method = "lm", color = "red") +
  geom_point(aes(colour = Tissue, shape = Substrate),
             alpha = 1, size = 3) +
  geom_point(colour = "black", size = 1) +
  scale_shape_manual(values = c(15,16,17))+
  scale_color_manual(values = c("darkgreen", 
                                "firebrick", 
                                "mediumblue"))+
  labs(x = "Metabolome Diversity (H')", 
       y = "Bacteria Diversity (H')",
       title = "Shannon diversity scatter plot",
       subtitle = paste("R-squared < 0.001, ", 
                        "p-value = 0.878")) +
  theme_minimal()
#dev.off()


################# DIVERSITY LINEAR MODELS #################
# Load data
shannon <- read.csv("shannon_calicioids.csv")
shannon <- as.data.frame(shannon)

# Variables as objects
TYPE <- shannon$Type
SUBSTRATE <- shannon$Substrate
SAMPLE <- shannon$Sample
BACTERIA <- shannon$Bacteria
FUNGI <- shannon$Fungi
METABOLITES <- shannon$Metabolites

# Si tu appelle ce nouveau facteur "inter",
shannon$inter <- paste(shannon$Type, ".",
                       shannon$Substrate, sep = "")
INTER <- shannon$inter
# Graph parameters
par(mfrow = c(1, 2))

### BACTERIA
## See data
hist(BACTERIA)
boxplot(BACTERIA~INTER,las=3)
# Mean and standard deviation
mean(BACTERIA[INTER == "Culture.Agar"])
sd(BACTERIA[INTER == "Culture.Agar"])
mean(BACTERIA[INTER == "Fungus.Wood"])
sd(BACTERIA[INTER == "Fungus.Wood"])
mean(BACTERIA[INTER == "Lichen.Bark"])
sd(BACTERIA[INTER == "Lichen.Bark"])
mean(BACTERIA[INTER == "Lichen.Wood"])
sd(BACTERIA[INTER == "Lichen.Wood"])
## Modèle linéaire
m1 <- lm(Bacteria ~ INTER, data = shannon)
# Check assumptions with residual diagnostics
# Homogeneity of variances
plot(residuals(m1) ~ fitted(m1)) # Graph
levene.test(BACTERIA, INTER) # Levene
fligner.test(BACTERIA ~ INTER) # Fligner
bartlett.test(BACTERIA ~ INTER) # Bartlett
# Normality of residuals
qqnorm(residuals(m1, type = "pearson")) # Graph
qqline(residuals(m1, type = "pearson")) # Line
ad.test(residuals(m1)) # Anderson Darling
shapiro.test(residuals(m1)) # Shapiro Wilk
# Assumptions NOT OK, residuals not normally distributed
## GENERALIZED LEAST SQUARE MODEL FOR HETEROGENEOUS VARIANCES
# Model heterogeneity explicitly
# Add term to allow variances to differ among groups
glsBACTERIA <- gls(BACTERIA ~ INTER, data = shannon,
                   weights = varIdent(form = ~ 1 | INTER))

# Check homogeneity of variances
plot(residuals(glsBACTERIA, 
               type = "pearson") ~ fitted(glsBACTERIA))
# MUCH BETTER
# Check normality of residuals
qqnorm(residuals(glsBACTERIA, type = "pearson"))
qqline(residuals(glsBACTERIA, type = "pearson"))
# OK
# Check results
summary(glsBACTERIA) # Beta estimates
anova(glsBACTERIA)  # ANOVA table
# Post-hoc multiple comparisons
# Necessary to convert to factors fo post-hoc
INTER2 <- as.factor(INTER)
# Rerun model heterogeneity explicitly with factors
glsBACTERIA2 <- gls(BACTERIA ~ INTER2,
                 weights = varIdent(form = ~ 1 | INTER2))
# Run the following lines. These introduce methods for 'gls' objects.
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)}
# The line below gives pairwise comparisons now.
confint(glht(glsBACTERIA2, 
             linfct=mcp(INTER2="Tukey")))
## Confidence intervals
# Create dataframe
full <- data.frame(BACTERIA = BACTERIA,
                   INTER = INTER)
# Pour un GLS, tu peux aller le chercher 
# avec une fonction de mon package:
preds <- predictSE(glsBACTERIA, newdata = predDat)
DF <- nrow(full) - attr(logLik(glsBACTERIA), "df")
low95 <- preds$fit - qt(df = DF, p = 0.975) * preds$se.fit
upp95 <- preds$fit + qt(df = DF, p = 0.975) * preds$se.fit
preds
low95
upp95
predDat

### FUNGI
## See data
hist(FUNGI)
boxplot(FUNGI~INTER,las=3)
# Mean and standard deviation
mean(FUNGI[INTER == "Culture.Agar"])
sd(FUNGI[INTER == "Culture.Agar"])
mean(FUNGI[INTER == "Fungus.Wood"])
sd(FUNGI[INTER == "Fungus.Wood"])
mean(FUNGI[INTER == "Lichen.Bark"])
sd(FUNGI[INTER == "Lichen.Bark"])
mean(FUNGI[INTER == "Lichen.Wood"])
sd(FUNGI[INTER == "Lichen.Wood"])
## Linear model
m2 <- lm(FUNGI ~ INTER, data = richness)
# Check assumptions with residual diagnostics
# Homogeneity of variances
plot(residuals(m2) ~ fitted(m2)) # Graph
levene.test(FUNGI, INTER) # Levene
fligner.test(FUNGI ~ INTER) # Fligner
bartlett.test(FUNGI ~ INTER) # Bartlett
# Normality of residuals
qqnorm(residuals(m2, type = "pearson")) # Graph
qqline(residuals(m2, type = "pearson")) # Line
ad.test(residuals(m2)) # Anderson Darling
shapiro.test(residuals(m2)) # Shapiro Wilk
# Assumptions OK
# See results
summary(m2)
## Confidence intervals
# Data for predictions
predDat <- expand.grid(INTER = factor(c("Culture.Agar",
                                        "Fungus.Wood",
                                        "Lichen.Bark",
                                        "Lichen.Wood")))
# Predictions
predict(m2, newdata = predDat,
        interval = "confidence",
        level = 0.95)
predDat

### METABOLITES
## See data
hist(METABOLITES)
boxplot(METABOLITES~INTER,las=3)
# Mean and standard deviation
mean(METABOLITES[INTER == "Culture.Agar"])
sd(METABOLITES[INTER == "Culture.Agar"])
mean(METABOLITES[INTER == "Fungus.Wood"])
sd(METABOLITES[INTER == "Fungus.Wood"])
mean(METABOLITES[INTER == "Lichen.Bark"])
sd(METABOLITES[INTER == "Lichen.Bark"])
mean(METABOLITES[INTER == "Lichen.Wood"])
sd(METABOLITES[INTER == "Lichen.Wood"])
## Linear model
m3 <- lm(METABOLITES ~ INTER, data = shannon)
# Check assumptions with residual diagnostics
# Homogeneity of variances
plot(residuals(m3) ~ fitted(m3)) # Graph
levene.test(METABOLITES, INTER) # Levene
fligner.test(METABOLITES ~ INTER) # Fligner
bartlett.test(METABOLITES ~ INTER) # Bartlett
# Normality of residuals
qqnorm(residuals(m3, type = "pearson")) # Graph
qqline(residuals(m3, type = "pearson")) # Line
ad.test(residuals(m3)) # Anderson Darling
shapiro.test(residuals(m3)) # Shapiro Wilk
# Assumptions OK
# See results
summary(m3)
## Confidence intervals
# Data for predictions
predDat <- expand.grid(INTER = factor(c("Culture.Agar",
                                        "Fungus.Wood",
                                        "Lichen.Bark",
                                        "Lichen.Wood")))
# Predictions
predict(m3, newdata = predDat,
        interval = "confidence",
        level = 0.95)
predDat

### CREATE GRAPH FOR PUBLICATION
library("beeswarm")

# Importer le csv 
colnames(shannon) <- c("Type", "Substrate", "Number",
                       "Bacteria", "Fungi",
                       "Metabolites", "Treatment")
shannon
# Format labels, add \n where appropriate.
my.labels1 <- c("Cultures\nAgar", "Fungi\nWood",
                "Lichens\nBark", "Lichens\nWood") 
# Create variable for SUBSTRATE
z <- as.numeric(factor(shannon$Treatment), replace = TRUE)
# Create colour variable for Group
my.pwcol1 <- c("red","red","red","red","red",
               "red","red","red","red","red",
               "red","red","red","red",
               "blue","blue","blue",
               "black","black","black","black",
               "blue","black","blue",
               "black","black","black","black","black",
               "blue")
# Create point shape variable for Group
my.pwpch1 <- c(17,17,17,17,17,17,17,17,17,17,
               17,17,17,17,15,15,15,19,19,19,
               19,15,19,15,19,19,19,19,19,15)

# 95% Confidence intervals and means
# Bacteria
MEAN_Bacteria <- c(4.88, 4.96, 4.91, 5.29)
CI95_Bacteria_UP <- c(5.27, 5.84, 5.29, 5.72)
CI95_Bacteria_LOW <- c(4.48, 4.09, 4.52, 4.86)
# Fungi
MEAN_Fungi <- c(2.26, 3.57, 3.27, 3.66)
CI95_Fungi_UP <- c(2.52, 4.14, 3.59, 4.23)
CI95_Fungi_LOW <- c(1.99, 3.00, 2.96, 3.09)
# Metabolites
MEAN_Metabolites <- c(4.81, 6.19, 6.17, 6.22)
CI95_Metabolites_UP <- c(5.01, 6.61, 6.41, 6.64)
CI95_Metabolites_LOW <- c(4.61, 5.76, 5.94, 5.79)

# Construct and export figure
#jpeg("Figure_Shannon.jpeg", units="in",
#      width=6, height=3, res=600)

# Define margin size
par(mfrow = c(1, 3))

# Plot Fungi
par(mar = c(3.0, 3.0, 0.1, 0.5), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
beeswarm(shannon$Fungi ~ shannon$Treatment, 
         xlab ="", ylab ="",
         cex = 0.6, ylim = c(0,8),
         pwpch = my.pwpch1,
         pwcol = my.pwcol1, cex.axis = 0.7, 
         cex.main = 0.8, labels = my.labels1)
title(ylab = "Shannon diversity index (H')",
      mgp = c(1.5, 1, 0))    # Add y-axis text
abline(v = 1.5, lty = 2, col = "gray35")
text(x = 2.5, y = 7.8, "Fungi", cex = 0.8, lwd = 1.5)
text(x = 1, y = 7.3, "A", cex = 0.8, lwd = 1.5)
text(x = 3, y = 7.3, "B", cex = 0.8, lwd = 1.5)
# Add standard error and mean
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Fungi_UP,
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Fungi,
       code=3, angle=90, length=0.04, 
       lwd = 1.0, col = "gray15")
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Fungi_LOW, 
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Fungi,
       code=3, angle=90, length=0.05, 
       lwd = 1.0, col = "gray15")

# Plot Bacteria
par(mar = c(3.0, 1.75, 0.1, 1.75), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
beeswarm(shannon$Bacteria ~ shannon$Treatment, 
         xlab ="", ylab ="",
         cex = 0.6, ylim = c(0,8),
         pwpch = my.pwpch1,
         pwcol = my.pwcol1, cex.axis = 0.7, 
         cex.main = 0.8, labels = my.labels1)
title(xlab = "Tissue origin", mgp = c(1.75, 3, 0)) 
text(x = 2.5, y = 7.8, "Bacteria", cex = 0.8, lwd = 1.5)
# Add standard error and mean
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Bacteria_UP,
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Bacteria,
       code=3, angle=90, length=0.04, 
       lwd = 1.0, col = "gray15")
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Bacteria_LOW, 
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Bacteria,
       code=3, angle=90, length=0.05, 
       lwd = 1.0, col = "gray15")

# Plot Metabolites
par(mar = c(3.0, 0.5, 0.1, 3.0), 
    mgp = c(1.5,0.5,0), xpd = FALSE)
beeswarm(shannon$Metabolites ~ shannon$Treatment, 
         xlab ="", ylab ="",
         cex = 0.6, ylim = c(0,8),
         pwpch = my.pwpch1,
         pwcol = my.pwcol1, cex.axis = 0.7, 
         cex.main = 0.8, labels = my.labels1)
abline(v = 1.5, lty = 2, col = "gray35")
text(x = 2.5, y = 7.8, "Metabolites", cex = 0.8, lwd = 1.5)
text(x = 1, y = 7.3, "A", cex = 0.8, lwd = 1.5)
text(x = 3, y = 7.3, "B", cex = 0.8, lwd = 1.5)
# Add standard error and mean
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Metabolites_UP,
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Metabolites,
       code=3, angle=90, length=0.04, 
       lwd = 1.0, col = "gray15")
arrows(x0=c(0.75,1.75,2.75,3.75), y0=CI95_Metabolites_LOW, 
       x1=c(0.75,1.75,2.75,3.75), y1=MEAN_Metabolites,
       code=3, angle=90, length=0.05, 
       lwd = 1.0, col = "gray15")


# Add legend
par(xpd = TRUE, mar = c(0, 0, 0, 0))
legend("right", inset=c(-0.25,0), title = "Substrate", 
       legend = c("Agar", "Bark", "Wood"),
       bty = "n", # remove box around legend 
       col = c("red","black","blue"), 
       pch = c(17,19,15), 
       cex = 0.7, horiz = FALSE)
#dev.off()



##################### VORONOI TREEMAPS ######################
## Subset sample types (culture, fungus, lichen)
# Fungi
ps3_ITS_C <- subset_samples(ps3_ITS, type == "C")
ps3_ITS_F <- subset_samples(ps3_ITS, type == "F")
ps3_ITS_L <- subset_samples(ps3_ITS, type == "L")
# Bacteria
ps3_16S_C <- subset_samples(ps3_16S, type == "C")
ps3_16S_F <- subset_samples(ps3_16S, type == "F")
ps3_16S_L <- subset_samples(ps3_16S, type == "L")
# Metabolites
ps2_MET_C <- subset_samples(ps2_MET, type == "C")
ps2_MET_F <- subset_samples(ps2_MET, type == "F")
ps2_MET_L <- subset_samples(ps2_MET, type == "L")

## Export taxa tables to manually change unknown genus to 
## closest higher known classification, when needed
# Export taxa tables
#write.csv2(data.frame(OTU = rownames(tax_table(ps3_ITS_C)), 
#                      tax_table(ps3_ITS_C)), 
#           "taxa_ITS_voronoi_export.csv", 
#           row.names = FALSE)
#write.csv(data.frame(OTU = rownames(tax_table(ps3_16S_C)), 
#                      tax_table(ps3_16S_C)), 
#           "taxa_16S_voronoi_export.csv", 
#           row.names = FALSE)

# Reimport taxa table
# Fungi culture
import_tax_table_ps3_ITS_C <- read.csv("taxa_ITS_voronoi.csv", 
                                       row.names = 1)
all(taxa_names(ps3_ITS_C) == rownames(import_tax_table_ps3_ITS_C))
tax_table(ps3_ITS_C) <- as.matrix(import_tax_table_ps3_ITS_C)
# Fungi fungi
import_tax_table_ps3_ITS_F <- read.csv("taxa_ITS_voronoi.csv", 
                                       row.names = 1)
all(taxa_names(ps3_ITS_F) == rownames(import_tax_table_ps3_ITS_F))
tax_table(ps3_ITS_F) <- as.matrix(import_tax_table_ps3_ITS_F)
# Fungi lichen
import_tax_table_ps3_ITS_L <- read.csv("taxa_ITS_voronoi.csv", 
                                       row.names = 1)
all(taxa_names(ps3_ITS_L) == rownames(import_tax_table_ps3_ITS_L))
tax_table(ps3_ITS_L) <- as.matrix(import_tax_table_ps3_ITS_L)

# Bacteria culture
import_tax_table_ps3_16S_C <- read.csv("taxa_16S_voronoi_C.csv", 
                                       row.names = 1)
all(taxa_names(ps3_16S_C) == rownames(import_tax_table_ps3_16S_C))
tax_table(ps3_16S_C) <- as.matrix(import_tax_table_ps3_16S_C)
# Bacteria fungi
import_tax_table_ps3_16S_F <- read.csv("taxa_16S_voronoi_F.csv", 
                                       row.names = 1)
all(taxa_names(ps3_16S_F) == rownames(import_tax_table_ps3_16S_F))
tax_table(ps3_16S_F) <- as.matrix(import_tax_table_ps3_16S_F)
# Bacteria lichen
import_tax_table_ps3_16S_L <- read.csv("taxa_16S_voronoi_L.csv", 
                                       row.names = 1)
all(taxa_names(ps3_16S_L) == rownames(import_tax_table_ps3_16S_L))
tax_table(ps3_16S_L) <- as.matrix(import_tax_table_ps3_16S_L)

## Convert phyloseq object to data frame
# Fungi culture
ps3_ITS_C_df <- as.data.frame(otu_table(ps3_ITS_C))
ps3_ITS_C_df$OTU <- rownames(ps3_ITS_C_df)
ps3_ITS_C_df <- melt(ps3_ITS_C_df, id.vars = "OTU",
                   variable.name = "Sample",
                   value.name = "Abundance")
# Fungi fungi
ps3_ITS_F_df <- as.data.frame(otu_table(ps3_ITS_F))
ps3_ITS_F_df$OTU <- rownames(ps3_ITS_F_df)
ps3_ITS_F_df <- melt(ps3_ITS_F_df, id.vars = "OTU",
                     variable.name = "Sample",
                     value.name = "Abundance")
# Fungi lichen
ps3_ITS_L_df <- as.data.frame(otu_table(ps3_ITS_L))
ps3_ITS_L_df$OTU <- rownames(ps3_ITS_L_df)
ps3_ITS_L_df <- melt(ps3_ITS_L_df, id.vars = "OTU",
                     variable.name = "Sample",
                     value.name = "Abundance")
# Bacteria culture
ps3_16S_C_df <- as.data.frame(otu_table(ps3_16S_C))
ps3_16S_C_df$OTU <- rownames(ps3_16S_C_df)
ps3_16S_C_df <- melt(ps3_16S_C_df, id.vars = "OTU",
                   variable.name = "Sample",
                   value.name = "Abundance")
# Bacteria fungi
ps3_16S_F_df <- as.data.frame(otu_table(ps3_16S_F))
ps3_16S_F_df$OTU <- rownames(ps3_16S_F_df)
ps3_16S_F_df <- melt(ps3_16S_F_df, id.vars = "OTU",
                     variable.name = "Sample",
                     value.name = "Abundance")
# Bacteria lichen
ps3_16S_L_df <- as.data.frame(otu_table(ps3_16S_L))
ps3_16S_L_df$OTU <- rownames(ps3_16S_L_df)
ps3_16S_L_df <- melt(ps3_16S_L_df, id.vars = "OTU",
                     variable.name = "Sample",
                     value.name = "Abundance")
# Metabolites culture
ps2_MET_C_df <- as.data.frame(otu_table(ps2_MET_C))
ps2_MET_C_df$OTU <- rownames(ps2_MET_C_df)
ps2_MET_C_df <- melt(ps2_MET_C_df, id.vars = "OTU",
                    variable.name = "Sample",
                    value.name = "Abundance")
# Metabolites fungi
ps2_MET_F_df <- as.data.frame(otu_table(ps2_MET_F))
ps2_MET_F_df$OTU <- rownames(ps2_MET_F_df)
ps2_MET_F_df <- melt(ps2_MET_F_df, id.vars = "OTU",
                    variable.name = "Sample",
                    value.name = "Abundance")
# Metabolites lichen
ps2_MET_L_df <- as.data.frame(otu_table(ps2_MET_L))
ps2_MET_L_df$OTU <- rownames(ps2_MET_L_df)
ps2_MET_L_df <- melt(ps2_MET_L_df, id.vars = "OTU",
                    variable.name = "Sample",
                    value.name = "Abundance")

## Add taxonomic information
# Fungi culture
tax_table_ITS_C <- as.data.frame(tax_table(ps3_ITS_C))
ps3_ITS_C_df <- merge(ps3_ITS_C_df, tax_table_ITS_C,
                    by.x = "OTU", by.y = "row.names")
# Fungi fungi
tax_table_ITS_F <- as.data.frame(tax_table(ps3_ITS_F))
ps3_ITS_F_df <- merge(ps3_ITS_F_df, tax_table_ITS_F,
                      by.x = "OTU", by.y = "row.names")
# Fungi lichen
tax_table_ITS_L <- as.data.frame(tax_table(ps3_ITS_L))
ps3_ITS_L_df <- merge(ps3_ITS_L_df, tax_table_ITS_L,
                      by.x = "OTU", by.y = "row.names")
# Bacteria culture
tax_table_16S_C <- as.data.frame(tax_table(ps3_16S_C))
ps3_16S_C_df <- merge(ps3_16S_C_df, tax_table_16S_C,
                    by.x = "OTU", by.y = "row.names")
# Bacteria fungi
tax_table_16S_F <- as.data.frame(tax_table(ps3_16S_F))
ps3_16S_F_df <- merge(ps3_16S_F_df, tax_table_16S_F,
                      by.x = "OTU", by.y = "row.names")
# Bacteria lichen
tax_table_16S_L <- as.data.frame(tax_table(ps3_16S_L))
ps3_16S_L_df <- merge(ps3_16S_L_df, tax_table_16S_L,
                      by.x = "OTU", by.y = "row.names")
# Metabolites culture
tax_table_MET_C <- as.data.frame(tax_table(ps2_MET_C))
ps2_MET_C_df <- merge(ps2_MET_C_df, tax_table_MET_C,
                     by.x = "OTU", by.y = "row.names")
# Metabolites fungi
tax_table_MET_F <- as.data.frame(tax_table(ps2_MET_F))
ps2_MET_F_df <- merge(ps2_MET_F_df, tax_table_MET_F,
                     by.x = "OTU", by.y = "row.names")
# Metabolites lichen
tax_table_MET_L <- as.data.frame(tax_table(ps2_MET_L))
ps2_MET_L_df <- merge(ps2_MET_L_df, tax_table_MET_L,
                     by.x = "OTU", by.y = "row.names")

## Replace non-positive values with a small positive number
# Fungi culture
small_positive <- min(ps3_ITS_C_df$Abundance[ps3_ITS_C_df$Abundance > 0]) / 100
ps3_ITS_C_df$Abundance[ps3_ITS_C_df$Abundance <= 0] <- small_positive
# Fungi fungi
small_positive <- min(ps3_ITS_F_df$Abundance[ps3_ITS_F_df$Abundance > 0]) / 100
ps3_ITS_F_df$Abundance[ps3_ITS_F_df$Abundance <= 0] <- small_positive
# Fungi lichen
small_positive <- min(ps3_ITS_L_df$Abundance[ps3_ITS_L_df$Abundance > 0]) / 100
ps3_ITS_L_df$Abundance[ps3_ITS_L_df$Abundance <= 0] <- small_positive
# Bacteria culture
small_positive <- min(ps3_16S_C_df$Abundance[ps3_16S_C_df$Abundance > 0]) / 100
ps3_16S_C_df$Abundance[ps3_16S_C_df$Abundance <= 0] <- small_positive
# Bacteria fungi
small_positive <- min(ps3_16S_F_df$Abundance[ps3_16S_F_df$Abundance > 0]) / 100
ps3_16S_F_df$Abundance[ps3_16S_F_df$Abundance <= 0] <- small_positive
# Bacteria lichen
small_positive <- min(ps3_16S_L_df$Abundance[ps3_16S_L_df$Abundance > 0]) / 100
ps3_16S_L_df$Abundance[ps3_16S_L_df$Abundance <= 0] <- small_positive
# Metabolites culture
small_positive <- min(ps2_MET_C_df$Abundance[ps2_MET_C_df$Abundance > 0]) / 100
ps2_MET_C_df$Abundance[ps2_MET_C_df$Abundance <= 0] <- small_positive
# Metabolites fungi
small_positive <- min(ps2_MET_F_df$Abundance[ps2_MET_F_df$Abundance > 0]) / 100
ps2_MET_F_df$Abundance[ps2_MET_F_df$Abundance <= 0] <- small_positive
# Metabolites lichen
small_positive <- min(ps2_MET_L_df$Abundance[ps2_MET_L_df$Abundance > 0]) / 100
ps2_MET_L_df$Abundance[ps2_MET_L_df$Abundance <= 0] <- small_positive

## Create Voronoi treemap
# Fungi culture
Voronoi_ITS_C <- voronoiTreemap(
  data = ps3_ITS_C_df, levels = c("phylum", "order", "genus"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.001,
  maxIteration = 100, seed = 123)
# Fungi fungi
Voronoi_ITS_F <- voronoiTreemap(
  data = ps3_ITS_F_df, levels = c("phylum", "order", "genus"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.0001,
  maxIteration = 100, seed = 123)
# Fungi lichen
Voronoi_ITS_L <- voronoiTreemap(
  data = ps3_ITS_L_df, levels = c("phylum", "order", "genus"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.0001,
  maxIteration = 100, seed = 123)
# Bacteria culture
Voronoi_16S_C <- voronoiTreemap(
  data = ps3_16S_C_df, levels = c("Phylum", "Order"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.0001,
  maxIteration = 100, seed = 120)
# Bacteria fungi
Voronoi_16S_F <- voronoiTreemap(
  data = ps3_16S_F_df, levels = c("Phylum", "Order"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.0005,
  maxIteration = 100, seed = 120)
# Bacteria lichen
Voronoi_16S_L <- voronoiTreemap(
  data = ps3_16S_L_df, levels = c("Phylum", "Order"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.0001,
  maxIteration = 100, seed = 120)
# Metabolites culture
Voronoi_MET_C <- voronoiTreemap(
  data = ps2_MET_C_df, levels = c("Pathway", "Class"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.000008,
  maxIteration = 250, seed = 420)
# Metabolites fungi
Voronoi_MET_F <- voronoiTreemap(
  data = ps2_MET_F_df, levels = c("Pathway", "Class"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.00001,
  maxIteration = 250, seed = 420)
# Metabolites lichen
Voronoi_MET_L <- voronoiTreemap(
  data = ps2_MET_L_df, levels = c("Pathway", "Class"),
  cell_size = "Abundance", shape = "circle",
  filter = 0.000004,
  maxIteration = 250, seed = 420)

## Draw Voronoi treemaps
# Fungi culture
#jpeg("Voronoi_Fungi_culture.jpeg", units="in",
#      width=6, height=3, res=1500)
drawTreemap(Voronoi_ITS_C, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato"))
#dev.off()
# Fungi fungi
#jpeg("Voronoi_Fungi_fungi.jpeg", units="in",
#     width=6, height=3, res=1500)
drawTreemap(Voronoi_ITS_F, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato",
                              "gray76"))
#dev.off()
# Fungi lichen
#jpeg("Voronoi_Fungi_lichen.jpeg", units="in",
#     width=6, height=3, res=1500)
drawTreemap(Voronoi_ITS_L, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato",
                              "navajowhite"))
#dev.off()
# Bacteria culture
#jpeg("Voronoi_Bacteria_culture.jpeg", units="in",
#     width=6, height=3, res=1500)
drawTreemap(Voronoi_16S_C, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato",
                              "slateblue",
                              "seagreen1", "plum",
                              "navajowhite","olivedrab"))
#dev.off()
# Bacteria fungi
#jpeg("Voronoi_Bacteria_fungi.jpeg", units="in",
#     width=6, height=3, res=1500)
drawTreemap(Voronoi_16S_F, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato",
                              "plum", "navajowhite",
                              "olivedrab", "goldenrod1"))
#dev.off()
# Bacteria lichen
#jpeg("Voronoi_Bacteria_lichen.jpeg", units="in",
#     width=6, height=3, res=1500)
drawTreemap(Voronoi_16S_L, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato",
                              "gray76","slateblue", 
                              "seagreen1", "goldenrod4",
                              "plum", "navajowhite"))
#dev.off()
# Metabolites culture
#jpeg("Voronoi_Metabolites_culture.jpeg", units="in",
#     width=6, height=3, res=1500)
drawTreemap(Voronoi_MET_C, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato",
                              "goldenrod1","slateblue", 
                              "seagreen1", "plum",
                              "navajowhite", "gray76"))
#dev.off()
# Metabolites fungi
#jpeg("Voronoi_Metabolites_fungi.jpeg", units="in",
#     width=6, height=3, res=1500)
drawTreemap(Voronoi_MET_F, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato",
                              "goldenrod1","slateblue", 
                              "seagreen1", "plum",
                              "navajowhite","gray76"))
#dev.off()
# Metabolites lichen
#jpeg("Voronoi_Metabolites_lichen.jpeg", units="in",
#     width=6, height=3, res=1500)
drawTreemap(Voronoi_MET_L, 
            label_size = 1.2, label_color = "black",
            border_size = 2, border_color = "black",
            width = 0.55, height = 0.9,
            color_palette = c("skyblue","tomato",
                              "goldenrod1","slateblue", 
                              "seagreen1", "plum",
                              "navajowhite","gray76"))
#dev.off()



##################### DISTANCE MATRICES #####################
# Define a default theme for ggplot graphics.
theme_set(theme_bw())
fontsize = 18L
theme_update(axis.title.x = element_text(size=fontsize))
theme_update(axis.title.y = element_text(size=fontsize))
theme_update(plot.title = element_text(size=fontsize+2))

## BACTERIA
set.seed(420)
# Need to create a phylogenetic tree 
tree <- rtree(ntaxa(ps3_16S), 
              rooted=TRUE, 
              tip.label=taxa_names(ps3_16S))
plot_tree(tree)

ps3_16S <- merge_phyloseq(ps3_16S, tree)
ps3_16S

UF_16S <- UniFrac(ps3_16S, weighted=TRUE, 
                  normalized=FALSE, parallel=FALSE)
# weighted = TRUE because we consider relative abundance, 
# not presence/absence
# normalized = FALSE because our data is already 
# normalized through DESeq2

## FUNGI
set.seed(420)
# Need to create a phylogenetic tree 
tree <- rtree(ntaxa(ps3_ITS), 
              rooted=TRUE, 
              tip.label=taxa_names(ps3_ITS))
plot_tree(tree)

ps3_ITS <- merge_phyloseq(ps3_ITS, tree)
ps3_ITS

UF_ITS <- UniFrac(ps3_ITS, weighted=TRUE, 
                  normalized=FALSE, parallel=FALSE)
# weighted = TRUE because we consider relative abundance, 
# not presence/absence
# normalized = FALSE because our data is already 
# normalized through DESeq2

# METABOLITES
# Chemical structural and compositional similarity (CSCS)
# install.packages("RCurl")
# devtools::install_github("askerdb/rCSCS")

# Import GNPS data
gnps <- prepare_GNPS("76ef55d181c04e8889089c34ee72413e")

# Built CSCS distance matrix
CSCS <- cscs(gnps$features, gnps$css,
             dissimilarity = T, cosine_threshold = 0.6,
             weighted = T, normalize = T)


################# PERMANOVA & DB-RDA ################

## BACTERIA
# make a data frame from the sample_data
ps3.2_16S <- data.frame(sample_data(ps3_16S))
# Create new composite variable with type and substrate
ps3.2_16S$inter <- c("Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Fungus.Wood",  "Fungus.Wood", 
                     "Fungus.Wood",  "Lichen.Bark",
                     "Lichen.Bark",  "Lichen.Bark", 
                     "Lichen.Bark",  "Lichen.Wood",
                     "Lichen.Bark",  "Lichen.Wood", 
                     "Lichen.Bark",  "Lichen.Bark",
                     "Lichen.Bark",  "Lichen.Bark", 
                     "Lichen.Bark",  "Lichen.Wood")
# Run PERMANOVA
permanova_16S <- adonis2(UF_16S ~ ps3.2_16S$inter,
                         ps3.2_16S, permutations = 10000)
permanova_16S

# Homogeneity of dispersion test
beta_16S <- betadisper(UF_16S, ps3.2_16S$inter)
permutest(beta_16S) # OK

# Same question with dbrda instead of permanova
dbrda_16S <- dbrda(UF_16S ~ ps3.2_16S$inter, 
                   ps3.2_16S, permutations = 10000)
dbrda_16S
anova.cca(dbrda_16S)
anova.cca(dbrda_16S, by = "axis")

# Create variables for clean legend in figure
ps3.2_16S$Substrate <- c("Agar","Agar","Agar","Agar",
                         "Agar","Agar","Agar","Agar",
                         "Agar","Agar","Agar","Agar",
                         "Agar","Agar","Wood","Wood",
                         "Wood","Bark","Bark","Bark",
                         "Bark","Wood","Bark","Wood",
                         "Bark","Bark","Bark","Bark",
                         "Bark","Wood")
sample_data(ps3_16S)$SampleType <- c("Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Fungus","Fungus",
                                     "Fungus","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen")

# Representation graphique
#jpeg("Bacteria_dbRDA.jpeg", units="in", 
#     width=6, height=3, res=600)
cap_plot <- plot_ordination(physeq = ps3_16S, 
                            ordination = dbrda_16S, 
                            color = "SampleType", 
                            axes = c(1,2)) + 
  aes(shape = Substrate) + 
  geom_point(aes(colour = SampleType), alpha = 1, size = 3) + 
  geom_point(colour = "black", size = 1) + 
  scale_color_manual(values = c("firebrick", 
                                "darkgreen", 
                                "mediumblue"))
cap_plot
#dev.off()

## FUNGI
# make a data frame from the sample_data
ps3.2_ITS <- data.frame(sample_data(ps3_ITS))
# Create new composite variable with type and substrate
ps3.2_ITS$inter <- c("Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Fungus.Wood",  "Fungus.Wood", 
                     "Fungus.Wood",  "Lichen.Bark",
                     "Lichen.Bark",  "Lichen.Bark", 
                     "Lichen.Bark",  "Lichen.Wood",
                     "Lichen.Bark",  "Lichen.Wood", 
                     "Lichen.Bark",  "Lichen.Bark",
                     "Lichen.Bark",  "Lichen.Bark", 
                     "Lichen.Bark",  "Lichen.Wood")
# Run PERMANOVA
permanova_ITS <- adonis2(UF_ITS ~ ps3.2_ITS$inter,
                         ps3.2_ITS, permutations = 10000)
permanova_ITS

# Homogeneity of dispersion test
beta_ITS <- betadisper(UF_ITS, ps3.2_ITS$inter)
permutest(beta_ITS)
# NOT OK

# Same question with dbrda instead of permanova
dbrda_ITS <- dbrda(UF_ITS ~ ps3.2_ITS$inter, 
               ps3.2_ITS, permutations = 10000)
dbrda_ITS
anova.cca(dbrda_ITS)
anova.cca(dbrda_ITS, by = "axis")

# Create variables for clean legend in figure
ps3.2_ITS$Substrate <- c("Agar","Agar","Agar","Agar",
                         "Agar","Agar","Agar","Agar",
                         "Agar","Agar","Agar","Agar",
                         "Agar","Agar","Wood","Wood",
                         "Wood","Bark","Bark","Bark",
                         "Bark","Wood","Bark","Wood",
                         "Bark","Bark","Bark","Bark",
                         "Bark","Wood")
sample_data(ps3_ITS)$SampleType <- c("Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Fungus","Fungus",
                                     "Fungus","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen")

# Representation graphique
#jpeg("Fungi_dbRDA.jpeg", units="in", 
#     width=6, height=3, res=600)
cap_plot <- plot_ordination(physeq = ps3_ITS, 
                            ordination = dbrda_ITS, 
                            color = "SampleType", 
                            axes = c(1,2)) + 
  aes(shape = Substrate) + 
  geom_point(aes(colour = SampleType), alpha = 1, size = 3) + 
  geom_point(colour = "black", size = 1) + 
  scale_color_manual(values = c("firebrick", 
                                "darkgreen", 
                                "mediumblue"))
cap_plot
#dev.off()

## METABOLITES
# make a data frame from the sample_data
ps2.2_MET <- data.frame(sample_data(ps2_MET))
# Create new composite variable with type and substrate
ps2.2_MET$inter <- c("Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Culture.Agar", "Culture.Agar",
                     "Fungus.Wood",  "Fungus.Wood", 
                     "Fungus.Wood",  "Lichen.Bark",
                     "Lichen.Bark",  "Lichen.Bark", 
                     "Lichen.Bark",  "Lichen.Wood",
                     "Lichen.Bark",  "Lichen.Wood", 
                     "Lichen.Bark",  "Lichen.Bark",
                     "Lichen.Bark",  "Lichen.Bark", 
                     "Lichen.Bark",  "Lichen.Wood")
# Run PERMANOVA
permanova_MET <- adonis2(CSCS ~ ps2.2_MET$inter,
                         ps2.2_MET, permutations = 10000)
permanova_MET

# Homogeneity of dispersion test
beta_MET <- betadisper(CSCS, ps2.2_MET$inter)
permutest(beta_MET)

# Same question with dbrda instead of permanova
dbrda_MET <- dbrda(CSCS ~ ps2.2_MET$inter, 
                   ps2.2_MET, permutations = 10000)
dbrda_MET
anova.cca(dbrda_MET)
anova.cca(dbrda_MET, by = "axis")

# Create variables for clean legend in figure
ps2.2_MET$Substrate <- c("Agar","Agar","Agar","Agar",
                         "Agar","Agar","Agar","Agar",
                         "Agar","Agar","Agar","Agar",
                         "Agar","Agar","Wood","Wood",
                         "Wood","Bark","Bark","Bark",
                         "Bark","Wood","Bark","Wood",
                         "Bark","Bark","Bark","Bark",
                         "Bark","Wood")
sample_data(ps2_MET)$SampleType <- c("Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Culture","Culture",
                                     "Fungus","Fungus",
                                     "Fungus","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen",
                                     "Lichen","Lichen")

# Representation graphique
#jpeg("Metabolites_dbRDA.jpeg", units="in",
#     width=6, height=3, res=600)
MET_plot <- plot_ordination(physeq = ps2_MET, 
                            ordination = dbrda_MET, 
                            color = "SampleType",
                            #label ="number",
                            axes = c(1,2)) + 
  aes(shape = Substrate) + 
  geom_point(aes(colour = SampleType), alpha = 1, size = 3) + 
  geom_point(colour = "black", size = 1) +
  scale_color_manual(values = c("firebrick", 
                                "darkgreen", 
                                "mediumblue"))
MET_plot
#dev.off()



################### MANTEL TESTS ###################
## BACTERIA VS FUNGI
mantel1 = mantel(UF_16S, UF_ITS, 
                 method = "pearson", 
                 permutations = 10000, 
                 na.rm = TRUE)
# Pearson is parametric, Spearman is non-parametric (based of ranks)
mantel1 # Non-significant no correlation

## BACTERIA VS FULL METABOLITES
mantel2 = mantel(CSCS, UF_16S, 
                 method = "pearson", 
                 permutations = 10000, 
                 na.rm = TRUE)
mantel2

## FUNGI VS FULL METABOLITES
mantel3 = mantel(CSCS, UF_ITS, 
                 method = "pearson", 
                 permutations = 10000, 
                 na.rm = TRUE)
mantel3


## SCATTER PLOTS MANTEL TEST
# Convert distance matrices to vectors
UF_16S_vec <- as.vector(as.matrix(UF_16S)[lower.tri(as.matrix(UF_16S))])
UF_ITS_vec <- as.vector(as.matrix(UF_ITS)[lower.tri(as.matrix(UF_ITS))])
CSCS_vec <- as.vector(as.matrix(CSCS)[lower.tri(as.matrix(CSCS))])

## 16S vs ITS
plot_data <- data.frame(x = UF_16S_vec, 
                        y = UF_ITS_vec)
ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  labs(x = "Bacteria distances", 
       y = "Fungi distances",
       title = "Mantel Test Scatter Plot",
       subtitle = paste("Mantel r =", 
                        round(mantel1$statistic, 3),
                        ", p-value =", 
                        round(mantel1$signif, 4))) +
  theme_minimal()

## 16S vs METABOLITES
set.seed(420)
plot_data <- data.frame(x = UF_16S_vec, 
                        y = CSCS_vec)
#jpeg("Mantel_Metabolites_vs_16S.jpeg", units="in",
#     width=6, height=3, res=600)
ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  labs(x = "Bacterial distance", 
       y = "Metabolite distance",
       title = "Mantel Test Scatter Plot",
       subtitle = paste("Mantel r =", 
                        round(mantel2$statistic, 3),
                        ", p-value =", 
                        round(mantel2$signif, 4))) +
  theme_minimal()
#dev.off()

## ITS vs METABOLITES
set.seed(420)
plot_data <- data.frame(x = UF_ITS_vec, 
                        y = CSCS_vec)
#jpeg("Mantel_Metabolites_vs_ITS.jpeg", units="in",
#     width=6, height=3, res=600)
ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red") +
  labs(x = "Fungal distance", 
       y = "Metabolite distance",
       title = "Mantel Test Scatter Plot",
       subtitle = paste("Mantel r =", 
                        round(mantel3$statistic, 3),
                        ", p-value =", 
                        round(mantel3$signif, 4))) +
  theme_minimal()
#dev.off()


############### HEATMAP LICHEN SUBSTANCES #################
# Original intensity
library(RColorBrewer)

## Heatmap for intensity classes
## Intensity classes : Absence = 0 ; < 1M = 1 ; 
## 1M to 5M = 2 ; 5M to 20M = 3 ; 20M to 100M = 4 ; 
## 100M to 500M = 5 ; 500M to 2B = 6 ; > 2B = 7

# Import data
substances1 <- read.csv(file = "lichens_substances_codes.csv",
                        header = TRUE, row.names = 1)
# Transform as matrix
substances1 <- as.matrix(substances1)

# Color coding by intensity class
greens <- brewer.pal(8, "Greens")

# Export to jpeg for building figure in Inkscape
jpeg("Heatmap_lichen_substances_intensity.jpeg", units="in",
     width=10, height=5, res=600)
# Build heatmap ; scale = none is VERY important
heatmap(substances1, Rowv = NA, Colv = NA,
        col = greens, scale = "none")
# Add legend
legend("left", legend = c("1", "2", "3", "4", 
                          "5", "6", "7", "8"), 
       fill = greens)
# Export to build figure in Inkscape
dev.off()

## Heatmap for presence-absence 
## (percent of samples in which a substance is present)
## Percent : Absence = 0 ; 1-14% = 1 ; 
## 15-29% = 2 ; 30-44% = 3 ; 45-59% = 4 ; 
## 60-74% = 5 ; 75-89% = 6 ; 90-100% = 7

# Import data
substances2 <- read.csv(file = "lichens_substances_presence.csv",
                        header = TRUE, row.names = 1)
# Transform as matrix
substances2 <- as.matrix(substances2)

# Color coding by intensity class
reds <- brewer.pal(8, "YlOrRd")

# Export to jpeg for building figure in Inkscape
jpeg("Heatmap_lichen_substances_presence.jpeg", units="in",
     width=10, height=5, res=600)
# Build heatmap ; scale = none is VERY important
heatmap(substances2, Rowv = NA, Colv = NA,
        col = reds, scale = "none")
# Add legend
legend("left", legend = c("1", "2", "3", "4", 
                          "5", "6", "7", "8"), 
       fill = reds)
# Export to build figure in Inkscape
dev.off()


