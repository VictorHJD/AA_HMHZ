## Please uncomment the first time you run this and re-install packages

#require(devtools)
#devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)

library(ggplot2)
library(MultiAmplicon)
library(reshape)
library(phyloseq)
library(data.table)
library(taxonomizr)
library(taxize)
library(parallel)

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- TRUE

doMultiAmp <- TRUE

doTax <- TRUE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run Microbiome#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

#path <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/2018_22_HMHZ_1_1/"
#path <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/2018_22_HMHZ_1_2/" ## change according to where you downloaded
#path <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/2018_22_HMHZ_2_1/"
path <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/2018_22_HMHZ_2_2/"

fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 


samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
#samples<- gsub("s\\d+-", "\\1", basename(samples)) ##For Pool 1
samples<- gsub("S\\d+-", "\\1", basename(samples)) ##For Pool 2
samples<- gsub("-", "_", basename(samples))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[1]])
## plotQualityProfile(fastqF[[380]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[380]])

#Creation of a folder for filtrated reads 

#filt_path <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/filtered1_1/"
#filt_path <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/filtered_2_1/"
filt_path <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/filtered_2_2/"


#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))
names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering
if(doFilter){
  lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(250,250), minLen=c(250,250), 
                  maxN=0, maxEE=2, truncQ=2, 
                  compress=TRUE, verbose=TRUE)
  })
}


names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/primer_list.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "Seq_F"]
primerR <- ptable[, "Seq_R"]
names(primerF) <- as.character(ptable[, "Name_F"])
names(primerR) <- as.character(ptable[, "Name_R"])
primer <- PrimerPairsSet(primerF, primerR)


##Multi amplicon pipeline
if(doMultiAmp){
  MA <- MultiAmplicon(primer, files)
  #filedir <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/stratified_files_1_1"
  #filedir <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/stratified_files_2_1"
  filedir <- "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/stratified_files_2_2"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
  
  errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                       verbose=0, multithread = 12)
  errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                      verbose=0, multithread = 12)
  
  MA <- derepMulti(MA, mc.cores=12) 
  
  MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                  verbose=0, mc.cores=12)
  
  MA <- mergeMulti(MA, mc.cores=12) 
  
  propMerged <- MultiAmplicon::calcPropMerged(MA)
  
  MA <- mergeMulti(MA, justConcatenate=propMerged<0.8, mc.cores=12) 
  
  MA <- makeSequenceTableMulti(MA, mc.cores=12) 
  
  MA <- removeChimeraMulti(MA, mc.cores=12)
  
  #saveRDS(MA, "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA.RDS") ##Pool HMHZ 1 preliminary run
  #saveRDS(MA, "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1.RDS") ##Pool HMHZ 1 full run
  #saveRDS(MA, "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA2_1.RDS") ##Pool HMHZ 2 full run (1st batch of data)
  saveRDS(MA, "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA2_2.RDS") ##Pool HMHZ 2 full run (2nd batch of data)
} else{
  #MA <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA.RDS") ###START from here now!
  #MA <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1.RDS") 
  #MA <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA2_1.RDS")
  MA <- readRDS("/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA2_2.RDS")
}

trackingF <- getPipelineSummary(MA) 
plotPipelineSummary(trackingF) 
PipSum <- plotPipelineSummary(trackingF) + scale_y_log10()
#ggsave("Sequencing_summary_HMHZ_2.pdf", PipSum, path = "~/AA_HMHZ/", height = 15, width = 15)


Heatmap <- plotAmpliconNumbers(MA)
#ggsave("Sequencing_reads_HMHZ_2.pdf", Heatmap, path = "~/AA_HMHZ/", height = 15, width = 15)

#pdf("~/AA_HMHZ/Sequencing_reads.pdf", plotAmpliconNumbers(MA), height = 15, width = 15)
#dev.off()

###New taxonomic assignment

#Sys.setenv("BLASTDB" = "/SAN/db/blastdb/") #To make the annotation work, boss will fix this in the package
#library("vctrs", lib.loc="/usr/local/lib/R/site-library")
#MA <- blastTaxAnnot(MA,  dataBaseDir = Sys.getenv("BLASTDB"), negative_gilist = "/SAN/db/blastdb/uncultured.gi", num_threads = 20)


MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    #infasta = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/HMHZ1_1_in.fasta",
                    #infasta = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/HMHZ1_2_in.fasta",
                    #infasta = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/HMHZ2_1_in.fasta",
                    infasta = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/HMHZ2_2_in.fasta",
                    #outblast = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/blast1_1_out.fasta",
                    #outblast = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/blast1_2_out.fasta",
                    #outblast = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/blast2_1_out.fasta",
                    outblast = "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/blast2_2_out.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 20)

#saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1_1Tax.Rds") ##Just Test run HMHZ 1
#saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1_2Tax.Rds") ## Just Full run HMHZ 1
#saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA2_1Tax.Rds") ## 1st batch Full run HMHZ 2
#saveRDS(MA, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA2_2Tax.Rds") ## 2nd batch Full run HMHZ 2


##Start from here after the taxonomic annotation
#MA<- readRDS(file= "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1_1Tax.Rds") ###Test run
#MA<- readRDS(file= "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA1_2Tax.Rds") ###Full run
#MA<- readRDS(file= "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA2_1Tax.Rds") ###Full run (1st batch)
#MA<- readRDS(file= "/SAN/Victors_playground/Metabarcoding/AA_HMHZ/MA2_2Tax.Rds") ###Full run (2nd batch)

###Load sample information

sample.data <- read.csv("~/AA_HMHZ/Sample_selection_Metabarcoding_Complete.csv",
                        dec=",", stringsAsFactors=FALSE)

##Tiny adjustments to metadata
require(dplyr)

sample.data$X.1<- NULL
sample.data$X<- NULL
sample.data$Run_I_ID <- NULL
sample.data$Run_II_ID <- NULL
sample.data$Year <- as.factor(sample.data$Year)
sample.data$HI <- as.numeric(sample.data$HI)
sample.data$HI<-round(sample.data$HI, 2)
sample.data$OPG_Eimeria <- as.numeric(sample.data$OPG_Eimeria)
sample.data <- sample.data[!is.na(sample.data$HI),] ##Eliminate 3 samples without HI

sample.data$Concentration <- as.numeric(sample.data$Concentration)
sample.data$Chip_number <- as.factor(sample.data$Chip_number)
sample.data$Body_weight <- as.numeric(sample.data$Body_weight)
sample.data$Body_length <- as.numeric(sample.data$Body_length)
sample.data$Sex <- as.factor(sample.data$Sex)

sample.data%>%
  mutate(BMI = Body_weight/((Body_length)^2)) -> sample.data

sample.data%>%
  mutate(Seq_Run = case_when(Chip_number%in% c(1:8) ~ "Pool_1",
                             Chip_number%in% c(9:15) ~ "Pool_2")) -> sample.data

sample.data$Seq_Run <- as.factor(sample.data$Seq_Run)

sample.data$Longitude<- as.numeric(sample.data$Longitude)
sample.data$Latitude<- as.numeric(sample.data$Latitude)
sample.data$Longitude<-round(sample.data$Longitude, 4)
sample.data$Latitude<-round(sample.data$Latitude, 4)
sample.data$Locality <- paste(sample.data$Latitude, sample.data$Longitude)

sample.data$OPG_Eimeria <- as.numeric(sample.data$OPG_Eimeria)
sample.data$delta_ct_MminusE <- as.numeric(sample.data$delta_ct_MminusE)

sample.data$Transect <- gsub(" ", "", sample.data$Transect)
sample.data$Transect <- as.factor(sample.data$Transect)

sample.data %>%
  mutate(Genotype = case_when(HI >= 0.95 ~ "Mmm", 
                              HI <= 0.05 ~ "Mmd", HI < 0.95 | HI > 0.05 ~ "Hybrid")) -> sample.data

sample.data$Genotype <- as.factor(sample.data$Genotype)

sample.data %>%
  mutate(Pinworms = sample.data$Aspiculuris_tetraptera+sample.data$Syphacia_obvelata) -> sample.data ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

sample.data %>%
  mutate(Eimeria = case_when(Species_tissue== "E_ferrisi"| Species_tissue== "E_falciformis"| Species_tissue== "E_vermiformis"| Species_tissue== "Other" ~ "Positive", 
                             Species_tissue== "Negative" ~ "Negative",
                             Flot== "FALSE" & Ap5== "FALSE" ~ "Negative")) -> sample.data

sample.data$Eimeria <- as.factor(sample.data$Eimeria)

rownames(sample.data) <- make.unique(sample.data$Mouse_ID) ##Works when MA contains single run data 

##Adding sample data 

MAsample <- addSampleData(MA, sample.data)

#plotAmpliconNumbers(MAsample[, which(colnames(MAsample)%in%
#                                       sample.data$Mouse_ID)])

####Raw sums
rawcounts_HMHZ <- data.frame(colSums(getRawCounts(MAsample)))
rawcounts_HMHZ[,2] <- rownames(rawcounts_HMHZ)
colnames(rawcounts_HMHZ) <- c("Raw_counts", "Mouse_ID")
rownames(rawcounts_HMHZ) <- c(1:nrow(rawcounts_HMHZ))
rawcounts_HMHZ <- data.frame(Mouse_ID = rawcounts_HMHZ$Mouse_ID, Raw_counts = rawcounts_HMHZ$Raw_counts) 

hist(rawcounts_HMHZ$Raw_counts)
summary(rawcounts_HMHZ$Raw_counts)
boxplot(rawcounts_HMHZ$Raw_counts)

sum(rawcounts_HMHZ$Raw_counts) ###Total raw counts reads 

##To phyloseq
##Sample data
PS <- toPhyloseq(MAsample, colnames(MAsample)) ##Now it work
##saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_1.Rds") ###Results from preliminary analysis (Sample data)
##saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_2.Rds") ###Results from full run pool 1(Sample data)
##saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_2_1.Rds") ###Results from first batch run pool 2 (Sample data)
##saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_2_2.Rds") ###Results from first batch run pool 2 (Sample data)

sum(otu_table(PS)) ##Total denoised reads

##Primer data
PS.l <- toPhyloseq(MAsample, colnames(MAsample),  multi2Single=FALSE) ##It work
###For primer analysis (Victor)
##saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_1_1.Rds") ###Test run pool 1
##saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_1_2.Rds") ###Full run Pool 1
##saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_2_1.Rds") ###First batch run Pool 2
##saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_2_2.Rds") ###Second batch run Pool 2


###
lapply(getTaxonTable(MAsample), function (x) table(as.vector(x[, "phylum"])))
lapply(getTaxonTable(MA), function (x) table(as.vector(x[, "genus"])))
lapply(getTaxonTable(MAsample), function (x) table(as.vector(x[, "species"])))

####Merge phyloseq objects 
if(HMHZ_1){
PS1 <-  readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_1.Rds")
PS2 <-  readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_2.Rds")

PS <- merge_phyloseq(PS1, PS2) ###Works! 
saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_All.Rds") ###Results from full + test run pool 1 
}


if(HMHZ_2){
  PS1 <-  readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_2_1.Rds")
  PS2 <-  readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_2_2.Rds")
  
  PS <- merge_phyloseq(PS1, PS2) ###Works! 
  saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_2_All.Rds") ###Results from full + test run pool 1 
}

###Merging HMHZ_1 and HMHZ_2 

##Pool 1 Merged
PS1 <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_All.Rds")

##Pool 2 Merged
PS2 <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_2_All.Rds")

PS <- merge_phyloseq(PS1, PS2) ####WORKED!!!!!!!!! :D :D
saveRDS(PS, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_All.Rds") ###Results from full + test run pool 1 

rm(PS1)
rm(PS2)
rm(MA)
rm(MAsample)
rm(PS.l)
rm(PS)

##Merge PS.list

PS.l.1.1<- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_1_1.Rds")
PS.l.1.2<- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_1_2.Rds")
PS.l.2.1<- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_2_1.Rds")
PS.l.2.2<- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_2_2.Rds")

###Eliminate the empty primers first and then merge the phyloseq lists... empty list make the next function bug
PS.l.1.1[sapply(PS.l.1.1, is.null)]<- NULL
PS.l.1.2[sapply(PS.l.1.2, is.null)]<- NULL
PS.l.2.1[sapply(PS.l.2.1, is.null)]<- NULL
PS.l.2.2[sapply(PS.l.2.2, is.null)]<- NULL

##It is necessary to eliminate those that are empty in both 
along<- names(PS.l.1.1)
PS.l.1 <- lapply(along, function(i) merge_phyloseq(PS.l.1.1[[i]], PS.l.1.2[[i]])) ##Merge all the information from both experiments
names(PS.l.1) <- names(PS.l.1.1)

along<- names(PS.l.2.1)
PS.l.2 <- lapply(along, function(i) merge_phyloseq(PS.l.2.1[[i]], PS.l.2.2[[i]])) ##Merge all the information from both experiments
names(PS.l.2) <- names(PS.l.2.1)

along<- names(PS.l.2)
PS.l <- lapply(along, function(i) merge_phyloseq(PS.l.1[[i]], PS.l.2[[i]])) ##Merge all the information from both experiments
names(PS.l) <- names(PS.l.2)

x<- names(PS.l)
x[6]<-"BGf_132_F.BGr_132_R"
x[22]<-"LSU_Fwd_2_3Mod_55_F.LSU_Rev_4_54_R"
x[28]<-"NLF184cw_74_F.NL818cw_74_R" 
names(PS.l)<- x

saveRDS(PS.l, file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_All.Rds") 

rm(PS.l.1)
rm(PS.l.1.1)
rm(PS.l.1.2)
rm(PS.l.2)
rm(PS.l.2.1)
rm(PS.l.2.2)
rm(along)

along<- names(PS.l)
###Raw counts 
###Raw counts by primer pair 

rawcounts.1.1 <- rowSums(getRawCounts(MA11))
rawcounts.1.1 <- data.frame(rawcounts.1.1)
rawcounts.1.1[,2] <- rownames(rawcounts.1.1)
colnames(rawcounts.1.1) <- c("Raw_counts_1.1", "Primer_name")
rownames(rawcounts.1.1) <- c(1:nrow(rawcounts.1.1))
rawcounts.1.1 <- data.frame(Primer_name = rawcounts.1.1$Primer_name, Raw_counts_1.1 = rawcounts.1.1$Raw_counts_1.1) ###change the order of the columns
rawcounts.1.1$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcounts.1.1$Primer_name)
rawcounts.1.1$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcounts.1.1$Primer_name)

rawcounts.1.2 <- rowSums(getRawCounts(MA12))
rawcounts.1.2 <- data.frame(rawcounts.1.2)
rawcounts.1.2[,2] <- rownames(rawcounts.1.2)
colnames(rawcounts.1.2) <- c("Raw_counts_1.2", "Primer_name")
rownames(rawcounts.1.2) <- c(1:nrow(rawcounts.1.2))
rawcounts.1.2 <- data.frame(Primer_name = rawcounts.1.2$Primer_name, Raw_counts_1.2 = rawcounts.1.2$Raw_counts_1.2) ###change the order of the columns
rawcounts.1.2$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcounts.1.2$Primer_name)
rawcounts.1.2$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcounts.1.2$Primer_name)

rawcounts.2.1 <- rowSums(getRawCounts(MA21))
rawcounts.2.1 <- data.frame(rawcounts.2.1)
rawcounts.2.1[,2] <- rownames(rawcounts.2.1)
colnames(rawcounts.2.1) <- c("Raw_counts_2.1", "Primer_name")
rownames(rawcounts.2.1) <- c(1:nrow(rawcounts.2.1))
rawcounts.2.1 <- data.frame(Primer_name = rawcounts.2.1$Primer_name, Raw_counts_2.1 = rawcounts.2.1$Raw_counts_2.1) ###change the order of the columns
rawcounts.2.1$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcounts.2.1$Primer_name)
rawcounts.2.1$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcounts.2.1$Primer_name)

rawcounts.2.2 <- rowSums(getRawCounts(MA22))
rawcounts.2.2 <- data.frame(rawcounts.2.2)
rawcounts.2.2[,2] <- rownames(rawcounts.2.2)
colnames(rawcounts.2.2) <- c("Raw_counts_2.2", "Primer_name")
rownames(rawcounts.2.2) <- c(1:nrow(rawcounts.2.2))
rawcounts.2.2 <- data.frame(Primer_name = rawcounts.2.2$Primer_name, Raw_counts_2.2 = rawcounts.2.2$Raw_counts_2.2) ###change the order of the columns
rawcounts.2.2$Primer_name <- gsub(pattern = " ", replacement = "", x = rawcounts.2.2$Primer_name)
rawcounts.2.2$Primer_name <- gsub(pattern = "-", replacement = "_", x = rawcounts.2.2$Primer_name)

###Merged rawcounts
rawcounts <- plyr::join(rawcounts.1.1, rawcounts.1.2, by= "Primer_name")
rawcounts <- plyr::join(rawcounts, rawcounts.2.1, by= "Primer_name")
rawcounts <- plyr::join(rawcounts, rawcounts.2.2, by= "Primer_name")

rawcounts %>%
  rowwise() %>%
  mutate(Reads_Pool_1 = sum(Raw_counts_1.1, Raw_counts_1.2)) %>%
  mutate(Reads_Pool_2 = sum(Raw_counts_2.1, Raw_counts_2.2)) %>%
  mutate(Total_reads = sum(Reads_Pool_1, Reads_Pool_2))-> rawcounts



