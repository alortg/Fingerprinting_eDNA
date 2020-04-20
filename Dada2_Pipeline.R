##############################################################################
##############################################################################
##############################################################################
                       DADA2 Installation
  
Refer to https://benjjneb.github.io/dada2/dada-installation.html

source("https://bioconductor.org/biocLite.R") #Install binaries from Bioconductor. Refer to https://bioconductor.org/install/
biocLite("dada2")

Good luck! If it does not work the third time, 
take a break, meditate, go for a run. Take it easy. 
Be kind to yourself. 

Done? Congratualations!!! Reward yourself with a treat :)

DADA2 Tutorial: https://benjjneb.github.io/dada2/tutorial.html
##############################################################################
##############################################################################
##############################################################################

#Calling packages
library(ape)
library(taxize)
library(ShortRead)
library(ggplot2)
library(seqinr)
library(phyloseq)
library(dplyr)
library(dada2)

#Getting Ready
setwd("~/Barcoding/eDNARSS/") # CHANGE ME
dir<-("~/Barcoding/") # CHANGE ME

path <- "~/MiSeq" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2]) 

#Filter and trim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Lunch time! No worries, the next lines will run for a while.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(110,110),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE) #WAIT

errR <- learnErrors(filtRs, multithread=TRUE) #WAIT

plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample Inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE) #WAIT

dadaRs <- dada(derepRs, err=errR, multithread=TRUE) #WAIT

dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Write OTU table
write.table(seqtab.nochim, paste(dir,"~/OTU_Table.csv",sep=""),row.names=T, sep=",")

# Assign taxonomy using the curated Reference Library (fasta file)
path<-paste(dir,"/reference library/", sep="")
setwd(path)

taxa <- assignTaxonomy(seqtab.nochim, "silva_euka02_bothprim_dada_names_DNA.fasta", multithread=TRUE)
taxa <- addSpecies(taxa, "silva_euka02_bothprim_dada_names_DNA.fasta") 

# Write Taxa table
write.table(taxa, paste(dir,"Taxa.csv",sep=""),row.names=T, sep=",")


#Optional: Track reads throuhg the pipeline. As a final check of our progress, we’ll look at the number of reads that made it through each step in the pipeline:
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
hist(log(colSums(seqtab.nochim))) 



##################################################################################
##################################################################################
##################################################################################

END of code. 


Now, open the OTU and the TAXA tables, check that both tables have the same number
of sequences (silly but conforting, like checking in the calculator that 1+1=2).

-In one table, the name of the samples is in the wrong position (first sequence),
move the title row one cell, until the data match.
-Then, copy the OTU table into a new doc, then paste transposed the Taxa table 
into the OTU table (above). 

##################################################################################
##################################################################################
##################################################################################