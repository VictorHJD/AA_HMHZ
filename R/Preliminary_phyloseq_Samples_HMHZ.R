###Phyloseq pipeline for diversity analysis
library(ggplot2)
library(reshape)
library(phyloseq)
library(data.table)
library(parallel)
library(microbiome)
library("pheatmap")
library(plyr)
library(dplyr)
library(gridExtra)
library(grid)
library(vegan)
library(tidyr)

if(!exists("PS")){
  source("MA_General_HMHZ.R")
}

##Functions
nameOtuByTax <- function(PS, taxon="genus"){
  otable <- otu_table(PS)
  rownames(otable) <- make.unique(as.character(tax_table(PS)[, "genus"]))
  otable
}

##Use the functions from joey711
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}

summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  summarydt = mdt[, list(meanRA = mean(RelativeAbundance),
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

plot_taxa_summary = function(physeq, Rank, GroupBy = NULL){
  # Get taxa summary table 
  dt1 = summarize_taxa(physeq, Rank = Rank, GroupBy = GroupBy)
  # Set factor appropriately for plotting
  RankCol = which(colnames(dt1) == Rank)
  setorder(dt1, -meanRA)
  dt1[, RankFac := factor(dt1[[Rank]], 
                          levels = rev(dt1[[Rank]]))]
  dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
  dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
  # Set zeroes to one-tenth the smallest value
  ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
  ebarMinFloor <- ebarMinFloor / 10
  dt1[(ebarMin == 0), ebarMin := ebarMinFloor]
  
  pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
    scale_x_log10() +
    xlab("Mean Relative Abundance") +
    ylab(Rank) +
    theme_bw()
  if(!is.null(GroupBy)){
    # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
    pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                size = 5)
  } else {
    # Don't include error bars for faceted version
    pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                        xmin = ebarMin))
  }
  return(pRank)
}

##Load data

Rarefied <- TRUE

Normalized <- FALSE

TSS <- FALSE

Relative_abund <- TRUE

##Preliminary data
#PS <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_1.Rds")
#PS <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_2.Rds")

##Pool 1 Merged
#PS <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_1_All.Rds")

##Pool 2 Merged
#PS <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_2_All.Rds")

##Pool 1 + 2 Merged
PS <- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqCombi_HMHZ_All.Rds")

summarize_phyloseq(PS) ##Non-Normalized

##Check it!
nsamples(PS)
ntaxa(PS)
sample_variables(PS)
rank_names(PS)
#tax_table(PS)
##Some samples have NO reads!

###Let's eliminate samples with 0 read counts
PS <- prune_samples(sample_sums(PS)>0, PS)

############# General overview ##################
###Check taxa
table(tax_table(PS)[, "superkingdom"], exclude = NULL) ## 476 ASV's are not assigned as Eukaryotes or Bacteria

##Eliminate Unassigned to superkingdom level 
PS <- subset_taxa(PS, !is.na(superkingdom) & !superkingdom %in% c("", "uncharacterized"))

###Check how many reads have every superkingdom
###Raw counts 
rawcounts_HMHZ <- data.frame(rowSums(otu_table(PS)))
rawcounts_HMHZ[,2] <- rownames(rawcounts_HMHZ)
###Bacterial and Eukaryotic counts
rawcounts_HMHZ[,3] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, superkingdom%in%"Bacteria"))))
rawcounts_HMHZ[,4] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, superkingdom%in%"Eukaryota"))))
#rawcounts_HMHZ[,5] <- as.data.frame(rowSums(otu_table(subset_taxa(PS, superkingdom%in%NA))))
colnames(rawcounts_HMHZ) <- c("Raw_counts", "Mouse_ID", "Bacteria_reads", "Eukaryota_reads")
rownames(rawcounts_HMHZ) <- c(1:nrow(rawcounts_HMHZ))
rawcounts_HMHZ <- data.frame(Mouse_ID = rawcounts_HMHZ$Mouse_ID, 
                              Raw_counts = rawcounts_HMHZ$Raw_counts,
                              Bacteria_reads= rawcounts_HMHZ$Bacteria_reads,
                              Eukaryota_reads= rawcounts_HMHZ$Eukaryota_reads)#,
                             #Unassigned_reads= rawcounts_HMHZ$Unassigned_reads) 

hist(rawcounts_HMHZ$Raw_counts)
summary(rawcounts_HMHZ$Raw_counts)
sum(rawcounts_HMHZ$Raw_counts)

## Bacteria read numbes 
sum(otu_table(subset_taxa(PS, superkingdom %in% "Bacteria")))

sum(otu_table(subset_taxa(PS, superkingdom %in% "Bacteria")))/sum(otu_table(PS))
## Eukaryote read numbers 
sum(otu_table(subset_taxa(PS, superkingdom %in% "Eukaryota")))

sum(otu_table(subset_taxa(PS, superkingdom %in% "Eukaryota")))/sum(otu_table(PS))
## Host read numbers
hist(rowSums(otu_table(subset_taxa(PS, genus%in%"Mus"))))

sum(otu_table(subset_taxa(PS, genus%in%"Mus")))/sum(otu_table(PS))
###A lot of Mus :(

###Eliminate reads assigned as "Mus"
PS <- subset_taxa(PS, genus!= "Mus") ##Eliminate samples :S 

####Taxa detected
as.data.frame(table(tax_table(PS)[, "phylum"]))
as.data.frame(table(tax_table(PS)[, "genus"]))

###Rarefaction curve 
rarcurv <- vegan::rarecurve(otu_table(PS),
                            label = F)

## Eliminate samples with low counts
##At 6000 reads samples have reached the species saturation (from rarecurve2!)
PSHigh <- prune_samples(sample_sums(PS)>=6000, PS)
summarize_phyloseq(PSHigh)
nsamples(PSHigh)
ntaxa(PSHigh)
sum(otu_table(PSHigh))
sum(otu_table(subset_taxa(PSHigh, superkingdom %in% "Bacteria")))/sum(otu_table(PSHigh))
sum(otu_table(subset_taxa(PSHigh, superkingdom %in% "Eukaryota")))/sum(otu_table(PSHigh))
#sum(otu_table(subset_taxa(PSHigh, genus%in%"Mus")))/sum(otu_table(PSHigh))

table(tax_table(PSHigh)[, "superkingdom"], exclude = NULL)

rarcurv2 <- vegan::rarecurve(otu_table(PSHigh),
                             label = F)

###Rarefy### 
###To the minumim reads
if(Rarefied){
  PSRare <- rarefy_even_depth(PSHigh, rngseed=1, sample.size= min(sample_sums(PSHigh)), replace=F)
  
  summarize_phyloseq(PSRare)
  nsamples(PSRare)
  ntaxa(PSRare)
  sum(otu_table(PSRare))
  sum(otu_table(subset_taxa(PSRare, superkingdom %in% "Bacteria")))/sum(otu_table(PSRare))
  sum(otu_table(subset_taxa(PSRare, superkingdom %in% "Eukaryota")))/sum(otu_table(PSRare))
  #sum(otu_table(subset_taxa(PSRare, genus%in%"Mus")))/sum(otu_table(PSRare))
  table(tax_table(PSRare)[, "superkingdom"], exclude = NULL)
  
  rarcurv4 <- vegan::rarecurve(otu_table(PSRare),
                               label = F)
  
  ###To a fixed treshold (10,000 reads)
  #PSRare <- rarefy_even_depth(PSHigh, rngseed=1, sample.size=10000, replace=F)
  
  #summarize_phyloseq(PSRare2)
  #nsamples(PSRare2)
  #ntaxa(PSRare2)
  #sum(otu_table(PSRare2))
  #sum(otu_table(subset_taxa(PSRare2, superkingdom %in% "Bacteria")))/sum(otu_table(PSRare2))
  #sum(otu_table(subset_taxa(PSRare2, superkingdom %in% "Eukaryota")))/sum(otu_table(PSRare2))
  #sum(otu_table(subset_taxa(PSRare2, genus%in%"Mus")))/sum(otu_table(PSRare2))
  #table(tax_table(PSRare2)[, "superkingdom"], exclude = NULL)
}

###Normalize
if(Normalized){
PSNor <- transform_sample_counts(PSHigh, function(x) x/sum(x)) ##Preprocessed (compositional) 
summarize_phyloseq(PSNor)

###Keep ASVs with more than 10^⁻5 
PSNor <- filter_taxa(PSNor, function(x) mean(x) > 1e-5, TRUE)

ntaxa(PSNor)
sum(otu_table(PSNor))
sum(otu_table(subset_taxa(PSNor, superkingdom %in% "Bacteria")))/sum(otu_table(PSNor))
sum(otu_table(subset_taxa(PSNor, superkingdom %in% "Eukaryota")))/sum(otu_table(PSNor))
sum(otu_table(subset_taxa(PSNor, genus%in%"Mus")))/sum(otu_table(PSNor))
table(tax_table(PSNor)[, "superkingdom"], exclude = NULL)
}
#### Total-sum scaling (TSS) normalization
if(TSS){
  PS.TSS <- metagMisc::phyloseq_standardize_otu_abundance(PSHigh, method = "total") #strange things :S 
}
## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case phylum and genus)
PS.Gen <-  tax_glom(PSHigh, "genus", NArm = F)
summarize_phyloseq(PS.Gen)

PS.Phy <-  tax_glom(PSHigh, "phylum", NArm = TRUE)
summarize_phyloseq(PS.Phy)

###Summarize sequencing depths 
sdt <- data.table(as(sample_data(PSHigh), "data.frame"),
                  TotalReads= sample_sums(PSHigh), keep.rownames = T)

sample_data(PSHigh)$TotalReads <- as.numeric(sdt$TotalReads)

if(Relative_abund){
  sdt <- data.table(as(sample_data(PSRare), "data.frame"),
                    TotalReads= sample_sums(PSRare), keep.rownames = T)
  
  sample_data(PSRare)$TotalReads <- as.numeric(sdt$TotalReads)
}

setnames(sdt, "rn", "Sample_ID")

###Read distribution by chip
SeqDep <- ggplot(sdt, aes(TotalReads)) +
  geom_histogram() + 
  facet_wrap(~Chip_number)+
  labs(tag = "A)")

###Sequencing depth by year 
SeqYear <- ggplot(sdt, aes(Year, TotalReads, color= as.factor(Chip_number)))+
  geom_point(size= 5)+
  geom_smooth(method = lm)+
  #scale_y_log10() + 
  facet_grid(Chip_number~.)+
  labs(tag = "B)")

pdf(file = "~/AA_HMHZ/Seq_char.pdf", width = 15, height = 15)
grid.arrange(SeqDep, SeqYear, ncol= 1, nrow= 2)
dev.off()

###Year of the sample does not impact the sequencing depth 

###Total counts taxa
tdt <- data.table(tax_table(PSHigh),
                 TotalCounts = taxa_sums(PSHigh),
                 OTU = taxa_names(PSHigh))

tdt <- data.table(tax_table(PSRare),
                  TotalCounts = taxa_sums(PSRare),
                  OTU = taxa_names(PSRare))

ggplot(tdt, aes(TotalCounts)) + 
  geom_histogram() 

###A lot of taxa with low counts (negative binomial distributed)
# taxa cumulative sum
taxcumsum <- tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]

# Define the plot
ggplot(taxcumsum, aes(TotalCounts, CumSum)) + 
  geom_point() +
  xlab("Filtering Threshold, Minimum Total Counts") +
  ylab("OTUs Filtered") +
  ggtitle("OTUs that would be filtered vs. the minimum count threshold")#+
  #xlim(0,100)

##Prevalence estimations 
#Prevalence = number of times an ASV is observed at least once
source("~/AA_HMHZ/taxa_summary.R", local = TRUE)
mdt <- fast_melt(PSRare)
prevd <- mdt[, list(Prevalence = sum(count > 0), 
                  TotalCounts = sum(count)),
           by = TaxaID]
ggplot(prevd, aes(Prevalence)) + 
  geom_histogram() 

ggplot(prevd, aes(Prevalence, TotalCounts)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()

##Top 9 Phyla prevalence 
addPhylum <- unique(copy(mdt[, list(TaxaID, phylum)]))

# Join by TaxaID
setkey(prevd, TaxaID)
setkey(addPhylum, TaxaID)
prevd <- addPhylum[prevd]
showPhyla <- prevd[, sum(TotalCounts), by = phylum][order(-V1)][1:9]$phylum
setkey(prevd, phylum)
yy <- ggplot(prevd[showPhyla], 
       mapping = aes(Prevalence, TotalCounts, color = phylum)) + 
  geom_point(size = 4, alpha = 0.75) + 
  scale_y_log10()+
  xlim(0,80)+
  labs(tag = "A)")

rm(tdt)
rm(taxcumsum)
rm(mdt)
rm(prevd)
rm(showPhyla)
#####################Alpha diversity###################### 
###Richness and diveristy
alphaDiv <- as.data.frame(estimate_richness(PSRare, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
alphaDiv$Sample_ID <- rownames(alphaDiv)
sdt <- plyr::join(sdt, alphaDiv, by= "Sample_ID")

###Evenness
eveAll <- evenness(PSRare, c("pielou", "simpson", "evar", "bulla"))
eveAll$Sample_ID <- rownames(eveAll)
sdt <- plyr::join(sdt, eveAll, by= "Sample_ID")

rm(eveAll)
rm(alphaDiv)

#sdt$Sex <- as.factor(sdt$Sex)

##Tiny adjustments 
#alphaDiv$HI<- as.numeric(alphaDiv$HI)
#alphaDiv$HI<-round(alphaDiv$HI, 2)

#alphaDiv%>%
#  mutate(Seq_Run = case_when(Chip_number%in% c(1:8) ~ "Pool_1",
#                             Chip_number%in% c(9:15) ~ "Pool_2")) -> alphaDiv

#sample_data(PSHigh)$Seq_Run <- as.factor(alphaDiv$Seq_Run)
#sample_data(PSRare)$Seq_Run <- as.factor(alphaDiv$Seq_Run)

#alphaDiv$Longitude<- as.numeric(alphaDiv$Longitude)
#alphaDiv$Latitude<- as.numeric(alphaDiv$Latitude)
#alphaDiv$Longitude<-round(alphaDiv$Longitude, 4)
#alphaDiv$Latitude<-round(alphaDiv$Latitude, 4)
#alphaDiv$Locality <- paste(alphaDiv$Latitude, alphaDiv$Longitude)

#sdt%>%
#  mutate(Seq_Run = case_when(Chip_number%in% c(1:8) ~ "Pool_1",
#                             Chip_number%in% c(9:15) ~ "Pool_2")) -> sdt

#sdt %>%
#  mutate(Genotype = case_when(HI >= 0.95 ~ "Mmm", 
#                              HI <= 0.05 ~ "Mmd", HI < 0.95 | HI > 0.05 ~ "Hybrid")) -> sdt

#sdt$Genotype <- as.factor(sdt$Genotype)

#sample_data(PSHigh)$Genotype <- sdt$Genotype
#sample_data(PSRare)$Genotype <- sdt$Genotype


#sdt$Transect <- gsub(" ", "", sdt$Transect)
#sdt$Transect <- as.factor(sdt$Transect)
#sample_data(PSRare)$Transect <- sdt$Transect 

#sdt$Body_length <- as.numeric(sdt$Body_length)
#sdt%>%
#  mutate(BMI = Body_weight/((Body_length)^2)) -> sdt

##Alpha diversity by year
a<- plot_richness(PSRare, x= "Year", color = "Year" , measures = c("Chao1", "Shannon")) +
  geom_boxplot()+
  #geom_violin()+
  theme_bw()+
  labs(tag= "A)")

##Wilcoxon rank-sum test (Mann-Whitney) Non-parametric test
pairwise.wilcox.test(sdt$Chao1, sample_data(PSRare)$Year)
pairwise.wilcox.test(sdt$Shannon, sample_data(PSRare)$Year)

##Alpha diversity by sex
a2<- plot_richness(PSRare, x= "Sex", color = "Sex", measures = c("Chao1", "Shannon")) +
  geom_boxplot()+
  theme_bw()+
  labs(tag= "B)")

##Wilcoxon rank-sum test (Mann-Whitney) Non-parametric test
pairwise.wilcox.test(sdt$Chao1, sample_data(PSRare)$Sex)
pairwise.wilcox.test(sdt$Shannon, sample_data(PSRare)$Sex)

##Alpha diversity by year and HI
a3<- plot_richness(PSRare, x= "HI", color = "HI" , measures = c("Observed","Chao1", "Shannon")) +
  geom_point(size= 2, alpha= 0.005)+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")+
  theme_bw()+
  labs(tag= "C)")

pdf(file = "~/AA_HMHZ/Alpha_div_Sex_Year_Transect_Rare.pdf", width = 8, height = 10)
grid.arrange(a, a2, a4, ncol= 1, nrow= 3)
dev.off()

a4<- plot_richness(PSRare, x= "Transect", color = "Transect", measures = c("Chao1", "Shannon")) +
  geom_boxplot()+
  theme_bw()+
  labs(tag= "C)")

##Wilcoxon rank-sum test (Mann-Whitney) Non-parametric test
pairwise.wilcox.test(sdt$Chao1, sample_data(PSRare)$Transect)
pairwise.wilcox.test(sdt$Shannon, sample_data(PSRare)$Transect)

a5<-  ggplot(sdt, aes(Year, pielou, fill= Year))+
  geom_boxplot()+
  theme_bw()+
  ylab("Evenness (Pielou's Index)")+
  labs(tag= "A)")

a6<-  ggplot(sdt, aes(Sex, pielou, color= Sex))+
  geom_boxplot()+
  theme_bw()+
  ylab("Evenness (Pielou's Index)")+
  labs(tag= "B)")

a7<-  ggplot(sdt, aes(Transect, pielou, color= Transect))+
  geom_boxplot()+
  theme_bw()+
  ylab("Evenness (Pielou's Index)")+
  labs(tag= "C)")

pdf(file = "~/AA_HMHZ/Evenn_Sex_Year__Transect_Rare.pdf", width = 7, height = 7)
grid.arrange(a5, a6, a7, ncol= 1, nrow= 3)
dev.off()

###And HI?
##By year
b<- ggplot(sdt, aes(HI, Chao1, color= Year))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+ 
  facet_wrap(~Year)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Alpha diversity (Chao1)")+
  theme_bw() +
  labs(tag= "A)")

##All togheter 
b2 <- ggplot(sdt, aes(HI, Chao1, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Richness (Chao1 Index)")+
  labs(tag= "A)")+
  theme_bw() +
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

b3 <- ggplot(sdt, aes(HI, Shannon, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Diversity (Shannon Index)")+
  labs(tag= "B)")+
  theme_bw() +
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

b4 <- ggplot(sdt, aes(HI, pielou, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Evenness (Pielou's Index)")+
  labs(tag= "C)")+
  theme_bw() +
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")


pdf(file = "~/AA_HMHZ/Alpha_div_year_HI_Rare.pdf", width = 10, height = 15)
grid.arrange(b, b2, ncol= 1, nrow= 2)
dev.off()

pdf(file = "~/AA_HMHZ/Alpha_div_HI_Rare.pdf", width = 10, height = 15)
grid.arrange(b2, b3, b4, ncol= 1, nrow= 3)
dev.off()

require("sjPlot")
##Biological
glmRic <- glm(Chao1~HI+Genotype+Sex+Year+BMI+Longitude+Latitude+Transect, data = sdt)
plot_model(glmRic, vline.color = "grey", sort.est = TRUE, show.values = TRUE, value.offset = .3)+
  theme_bw()
summary(glmRic)

##Technical
#Seq_Run+Chip_number+TotalReads

glmDiv <- glm(Shannon~HI+Genotype+Sex+Year+BMI+Longitude+Latitude+Transect, data = sdt)
plot_model(glmDiv, vline.color = "grey", sort.est = TRUE, show.values = TRUE, value.offset = .3)+
  theme_bw()
summary(glmDiv)

glmEve <- glm(pielou~HI+Genotype+Sex+Year+BMI+Longitude+Latitude+Transect, data = sdt)
plot_model(glmEve, vline.color = "grey", sort.est = TRUE, show.values = TRUE, value.offset = .3)+
  theme_bw()
summary(glmEve)

#plot_composition(PSRare)
##Tiny adjustments to the metadata 
#sample_data(PS)$HI <- alphaDiv$HI
#alphaDiv$Year <- as.factor(alphaDiv$Year)
#sample_data(PS)$Year <- alphaDiv$Year
#sample_data(PSHigh.rare)$Chip_number <- as.factor(sample_data(PSHigh.rare)$Chip_number)
#sample_data(PSRare)$Concentration <- as.numeric(sample_data(PSRare)$Concentration)
#sample_data(PSRare)$Chip_number <- as.factor(sample_data(PSRare)$Chip_number)
#sample_data(PSRare)$Body_weight <- as.numeric(sample_data(PSRare)$Body_weight)
#sample_data(PSRare)$Transect <- sdt$Transect

###Beta diversity (just with bray-curtis dissimilarity)

ordi <- ordinate(PSRare, method="PCoA", distance="bray")

ordi2<- ordinate(PSRare, method="NMDS", distance="bray")

dis <- phyloseq::distance(PSRare, method="bray")

c<- plot_ordination(PSRare, ordi, color="HI")+
    scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
    theme_bw()+
    geom_point(size=5, alpha= 0.75)+
    labs(tag= "A)")

c2 <- plot_ordination(PSRare, ordi, color="Sex")+
      theme_bw()+
      geom_point(size=5, alpha= 0.75)+
      labs(tag= "B)")

c3 <- plot_ordination(PSRare, ordi, color="Concentration")+
      theme_bw() +
      geom_point(size=5, alpha= 0.75)+
      scale_color_gradient(low="green", high="red")+
      labs(tag= "C)")

c4 <- plot_ordination(PSRare, ordi, color="Chip_number")+
      theme_bw() +
      geom_point(size=5, alpha= 0.75)+
      labs(tag= "D)")

c5 <- plot_ordination(PSRare, ordi, color="Year")+
      theme_bw() +
      geom_point(size=5, alpha= 0.75)+
      labs(tag= "E)")

c6 <- plot_ordination(PSRare, ordi, color="Seq_Run")+
      theme_bw() +
      geom_point(size=5, alpha= 0.75)+
      labs(tag= "F)")

##Extra way of plotting with package "Microbiome"
c7 <- plot_landscape(PSRare, "MDS", "bray", col = "HI")+
      scale_color_gradient(low="blue", high="red", name= "Hybrid Index)")+
      labs(tag= "B)")

c8<- plot_ordination(PSRare, ordi, type="taxa", color="phylum") + 
     facet_wrap(~phylum, 9) +
     theme_bw()

c9<- plot_ordination(PSRare, ordi, color="TotalReads")+
  theme_bw()+
  scale_color_gradient(low="#F89441FF", high="#0D0887FF")+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "G)")

pdf(file = "~/AA_HMHZ/MDS_HI_Rare.pdf", width = 8, height = 10)
grid.arrange(c, ncol= 1, nrow= 1)
dev.off()

pdf(file = "~/AA_HMHZ/Beta_diversity_Phylum_Rare.pdf", width = 8, height = 10)
grid.arrange(c8, ncol= 1, nrow= 1)
dev.off()

pdf(file = "~/AA_HMHZ/Beta_diversity_HMHZ_Rare.pdf", width = 20, height = 15)
grid.arrange(c, c2, c3, c4, c5, c6, ncol= 2, nrow= 3)
dev.off()

pdf(file = "~/AA_HMHZ/Beta_diversity_HMHZ_Technical.pdf", width = 8, height = 10)
grid.arrange( c4, c6, c9, ncol= 1, nrow= 3)
dev.off()

rm(a, a2, a3, a4, a5, a6, a7, b, b2, b3, b4, c, c2, c3, c4, c5, c6, c7, c8, c9)
####Permanova
adonis(dis ~ HI+Longitude+Latitude+Year+Chip_number+Concentration+Seq_Run, data = sdt)
adonis(dis ~ HI+Genotype+Year+Longitude+Latitude+Transect, data = sdt)

####Mantel test

geo.dist<- dist(cbind(sdt$Longitude, sdt$Latitude))

HI.dist<- dist(sdt$HI)

pr_dis2simil

## Sum asvCount by taxa 
HMHZPhy<- summarize_taxa(PS.Phy, "phylum", "Mouse_ID")
phylnames<- unique(HMHZPhy$phylum)

##Create matrix for PCA
HMHZPhy%>%
  select_(1,2,3) %>%
  group_by(Mouse_ID)%>%
  spread(phylum, meanRA)-> foo

foo[is.na(foo)]<- 0
foo<-as.data.frame(foo)
rownames(foo)<-foo[,1]
foo.samples<-join(foo, sdt, by= "Mouse_ID")
rownames(foo.samples)<-foo.samples[,1]
foo[,1]<- NULL
foo.samples[,1]<- NULL ##foo.samples is the file required for the Maximum likelihood analysis

##Create data frame for composition plot 
HMHZPhy%>%
  select_(1,2,3) %>%
  group_by(Mouse_ID)%>%
  mutate(Main_taxa= meanRA>= 0.10)%>%
  group_by(Mouse_ID)%>%
  mutate(phylum= case_when(Main_taxa== FALSE ~ "Taxa less represented", TRUE ~ as.character(phylum))) %>%
  arrange(Mouse_ID, desc(phylum))-> HMHZPhy

### Let's plot 
ggplot(data=HMHZPhy, aes(x= Mouse_ID, y= meanRA, fill= phylum)) +
  scale_fill_manual(values = c( "#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99",
                                "#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99",
                                "#B15928","#1B9E77","#D95F02","#7570B3","#E7298A", "#666666", "#66A61E","#E6AB02")) +
  geom_bar(aes(), stat="identity", position="fill") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(x = "Mouse ID", y= "Relative abundance", tag = "E)")+
  guides(fill= guide_legend(nrow = 8))
  
require("factoextra")
require("FactoMineR")
foo.pca<- PCA(foo, graph = T)
foo.eig<- get_eigenvalue(foo.pca)
fviz_eig(foo.pca, addlabels = TRUE, ylim = c(0, 25))

fviz_pca_ind(foo.pca,
             geom.ind = "point", # show points only (but not "text")
             col.ind = foo.samples$HI, # color by groups
             gradient.cols = c("blue", "red"),
             alpha.ind = 0.7,
             addEllipses = F, # Concentration ellipses
             legend.title = "Hybrid \nIndex")+
  labs(tag = "A)")
  
fviz_pca_var(foo.pca, col.var = "cos2", select.var =  list(contrib = 10), ##top ten taxa contributing
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE )+ # Avoid text overlapping
  labs(tag = "B)")

###Beta diversity by Kingdom 
###Bacteria 

ordi.Bac <- ordinate(PS.bacteria, method="PCoA", distance="bray")
dis.Bac <- phyloseq::distance(PS.bacteria, method="bray")
plot_ordination(PS.bacteria, ordi.Bac, color="HI")+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  theme_bw()+
  facet_wrap(~Chip_number) +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

plot_ordination(PS.bacteria, ordi.Bac, color="Genotype")+
  theme_bw() +
  geom_point(size=5, alpha= 0.75)+
  #geom_path()+
  labs(tag= "D)")

plot_ordination(PS.bacteria, ordi.Bac, type="taxa", color="phylum") + 
  facet_wrap(~phylum, 9) +
  theme_bw()

plot_ordination(PS.bacteria, ordi.Bac, color="Chip_number")+
  theme_bw() +
  geom_point(size=5, alpha= 0.75)+
  #geom_path()+
  facet_wrap(~Chip_number) +
  labs(tag= "D)")

plot_scree(ordination = ordi.Bac)

###Eukaryota
ordi.Euk <- ordinate(PS.eukaryota, method="PCoA", distance="bray")
dis.Euk <- phyloseq::distance(PS.eukaryota, method="bray")
plot_ordination(PS.eukaryota, ordi.Euk, color="HI")+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  theme_bw()+
  facet_wrap(~Chip_number) +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

plot_ordination(PS.eukaryota, ordi.Euk, type="taxa", color="phylum") + 
  facet_wrap(~phylum, 9) +
  theme_bw()

plot_ordination(PS.eukaryota, ordi.Euk, color="Chip_number")+
  theme_bw() +
  geom_point(size=5, alpha= 0.75)+
  #geom_path()+
  facet_wrap(~Chip_number) +
  labs(tag= "D)")

require(ade4)
r1<- mantel.rtest(geo.dist, dis, nrepet = 9999)
r2<- mantel.rtest(HI.dist, dis, nrepet = 9999)

plot(r2, main= "Mantel's test")

####GLM to model Mean-variance relationship
#install.packages("mvabund", repos="http://R-Forge.R-project.org") ##Not possible to install :S 

###
require(DESeq2)
DSHigh<- phyloseq_to_deseq2(PSRare, ~Year)
DSHigh<- DESeq(DSHigh)
DSes <- results(DSHigh, alpha=0.01)

#adonis(dis ~ Chip_number+Year+Concentration, data = sdt)
#adonis(dis.TSS ~ Chip_number+Year+Concentration, data = sdt)

###HI is not working because of NAs :S 

rm(dis)
rm(dis2)
rm(ordi)
rm(ordi3)
##PCA by Genus 
library(FactoMineR)
library(factoextra)

df<- as.data.frame(mat)
foo <- PCA(df)
foo$eig
fviz_pca_ind(res.pca)
plot(foo,  choix = "var")
head(foo$ind$coord)

###Network
#net<- make_network(PSHigh, type = "samples", distance = "bray", max.dist = 0.85)
#plot_network(net, PSHigh, color = "Year", shape = "Sex", line_weight = 0.4, 
#             label = NULL)

#plot_bar(PSHigh, fill="phylum") #+ facet_wrap(~dpi, scales= "free_x", nrow=1)

#### Data by group
PS.bacteria <- subset_taxa(PSHigh, superkingdom%in%"Bacteria")
PS.bacteria <- subset_taxa(PSRare, superkingdom%in%"Bacteria")

PS.eukaryota <- subset_taxa(PSHigh, superkingdom%in%"Eukaryota") 
PS.eukaryota <- subset_taxa(PSRare, superkingdom%in%"Eukaryota") 


PS.para <- subset_taxa(PSHigh, phylum%in%c("Nematoda",
                                       "Apicomplexa",
                                       "Platyhelminthes"))

PS.nema <- subset_taxa(PSHigh, phylum%in%"Nematoda")

PS.api <- subset_taxa(PSHigh, phylum%in%"Apicomplexa")

PS.plat <- subset_taxa(PSHigh, phylum%in%"Platyhelminthes")

PS.cili <- subset_taxa(PSHigh, phylum%in%"Ciliophora")

#PS.para2 <- subset_taxa(PSHigh, order%in%c("Diplomonadida", "Trichomonadida"))

#PS.para3 <- merge_phyloseq(PS.para, PS.para2)

PS.Trich <- subset_taxa(PSHigh, genus%in%c("Trichomonas",
                                            "Tritrichomonas",
                                            "Tetratrichomonas"))

PS.Forni <- subset_taxa(PSHigh, genus%in%c("Spironucleus",
                                           "Giardia",
                                           "Octomitus", 
                                           "Chilomastix"))

PS.eimeria <- subset_taxa(PSHigh, genus%in%"Eimeria")

PS.eimcyc <-subset_taxa(PSHigh, genus%in%c("Eimeria", "Cyclospora"))

PS.Crypto <- subset_taxa(PSHigh, genus%in%"Cryptosporidium")

PS.pinworms <- subset_taxa(PSHigh, genus%in%c("Syphacia", "Aspiculuris"))

PS.syp <- subset_taxa(PSHigh, genus%in%"Syphacia")

PS.asp <- subset_taxa(PSHigh, genus%in%"Aspiculuris")

PS.fungi <-subset_taxa(PSHigh, phylum%in%c("Ascomycota",
                           "Basidiomycota",
                           "Zygomycota", "Mucoromycota", "Zoopagomycota", "Blastocladiomycota", "Cryptomycota", "Microsporidia"))

PGHigh <- tax_glom(PSHigh, "genus", NArm = TRUE)
PG.para <- tax_glom(PS.para, "species", NArm = TRUE)
PG.para <- tax_glom(PS.para3, "genus", NArm = TRUE)
PG.fungi <- tax_glom(PS.fungi, "genus", NArm = TRUE)
PG.bacteria <- tax_glom(PS.bacteria, "phylum", NArm = TRUE)
#PG.bacteria <- tax_glom(PS.bacteria, "family", NArm = TRUE)
PG.Eimeria <- tax_glom(PS.eimeria, "species", NArm = TRUE)
PG.Tricho <- tax_glom(PS.Trich, "genus", NArm = TRUE)
PG.Forni <- tax_glom(PS.Forni, "genus", NArm = TRUE)
PG.nema <- tax_glom(PS.nema, "genus", NArm = TRUE)
PG.api <- tax_glom(PS.api, "genus", NArm = TRUE)
PG.plat <- tax_glom(PS.plat, "genus", NArm = TRUE)
PG.cili <- tax_glom(PS.cili, "species", NArm = TRUE)
PG.eimcyc <- tax_glom(PS.eimcyc, "species", NArm = TRUE)
PG.Proteo <- tax_glom(PS.Proteo, "genus", NArm = TRUE)
PG.Asco <- tax_glom(PS.Asco, "genus", NArm = TRUE)
PG.Bas <- tax_glom(PS.Bas, "genus", NArm = TRUE)
PG.Muc<- tax_glom(PS.Muc, "genus", NArm = TRUE)

table(tax_table(PG.Asco)[,5])
table(tax_table(PG.Proteo)[,5])
table(tax_table(PG.eimcyc)[,6])
table(tax_table(PG.para)[,5])
table(tax_table(PG.fungi)[,5])
table(tax_table(PG.bacteria)[,2])
table(tax_table(PG.Eimeria)[,6])
table(tax_table(PG.Tricho)[,6])
table(tax_table(PG.Forni)[,6])
table(tax_table(PG.nema)[,5])
table(tax_table(PG.api)[,5])
table(tax_table(PG.plat)[,5])
table(tax_table(PG.cili)[,6])
table(tax_table(PGHigh)[,5])

mat <- as.matrix(t(otu_table(PG.Asco)))
mat <- as.matrix(t(otu_table(PG.Bas)))
mat <- as.matrix(t(otu_table(PG.Muc)))
mat <- as.matrix(t(otu_table(PG.Proteo)))
mat <- as.matrix(t(otu_table(PG.para)))
mat <- as.matrix(t(otu_table(PG.fungi)))
mat<-as.matrix(t(otu_table(PG.Eimeria)))
mat<-as.matrix(t(otu_table(PG.bacteria)))
mat<-as.matrix(t(otu_table(PG.Tricho)))
mat<-as.matrix(t(otu_table(PG.Forni)))
mat<-as.matrix(t(otu_table(PG.cili)))
mat<- as.matrix(otu_table(PGHigh))
mat<-as.matrix(t(otu_table(PG.eimcyc)))
mat<-as.matrix(t(otu_table(PG.plat)))

rownames(mat) <- make.names(tax_table(PG.plat)[, 5])
rownames(mat) <- make.names(tax_table(PG.Asco)[, 5])
rownames(mat) <- make.names(tax_table(PG.Bas)[, 5])
rownames(mat) <- make.names(tax_table(PG.Muc)[, 5])
rownames(mat) <- make.names(tax_table(PG.Proteo)[, 5])
rownames(mat) <- make.names(tax_table(PG.eimcyc)[, 6])
rownames(mat) <- make.names(tax_table(PG.para)[, 5])
rownames(mat) <- make.names(tax_table(PG.fungi)[, 5])
rownames(mat) <- make.names(tax_table(PG.Eimeria)[, 6])
rownames(mat) <- make.names(tax_table(PG.bacteria)[, 2])
rownames(mat) <- make.names(tax_table(PG.Tricho)[, 5])
rownames(mat) <- make.names(tax_table(PG.Forni)[, 5])
rownames(mat) <- make.names(tax_table(PG.cili)[, 6])
colnames(mat)<- make.names(tax_table(PGHigh)[,5])
#pdf("figures/parasite_genera_heat.pdf", width=8, height=14)

pheatmap(log10(mat+1))

######Alpha diversity (Bacteria, Parasites, Fungi)
alphaBac <- as.data.frame(estimate_richness(PS.bacteria, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
alphaBac$Sample_ID <- rownames(alphaBac)
alphaBac <- plyr::join(sdt, alphaBac, by= "Sample_ID")

alphaPar <- as.data.frame(estimate_richness(PS.para3, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
alphaPar$Sample_ID <- rownames(alphaPar)
alphaPar <- plyr::join(sdt, alphaPar, by= "Sample_ID")

alphaFun <- as.data.frame(estimate_richness(PS.fungi, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
alphaFun$Sample_ID <- rownames(alphaFun)
alphaFun <- plyr::join(sdt, alphaFun, by= "Sample_ID")

pdf(file = "~/AA_HMHZ/AlphaDiv_HI_Groups.pdf", width = 35, height = 20)
d <- ggplot(alphaBac, aes(HI, Chao1, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Bacteria diversity (Shannon Index)")+
  labs(tag= "A)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

d2 <- ggplot(alphaPar, aes(HI, Chao1, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Parasites diversity (Shannon Index)")+
  labs(tag= "B)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

d3 <- ggplot(alphaFun, aes(HI, Shannon, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Fungi diversity (Shannon Index)")+
  labs(tag= "C)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

grid.arrange(d, d2, d3, ncol= 3, nrow= 1)
dev.off()

### DNA vs alpha diveristy bacteria 
alphaBac$Concentration<- as.numeric(alphaBac$Concentration) #Adjust
ggplot(alphaBac, aes(Concentration, Chao1, color= Year))+
  geom_point(size= 3)+
  #geom_smooth(method = lm)+
  xlab("DNA concentration (ng/µL)")+
  ylab("Richness (Chao1 Index)")+
  theme_bw()+
  labs(tag= "A)")


#Add total amount of bacterial reads 
###Total Bacterial
sdtBac <- data.table(as(sample_data(PS.bacteria), "data.frame"),
                     ReadsBac= sample_sums(PS.bacteria), keep.rownames = T)

sdtBac <- dplyr::select(sdtBac, 4,66)

sdt <- plyr::join(sdt, sdtBac, by= "Mouse_ID")

rm(sdtBac)

sdt %>%
  mutate(Bac_abundance = sdt$ReadsBac/sdt$TotalReads) -> sdt ### Add a new variable 

require(ggpubr)
ggplot(sdt, aes(Concentration, ReadsBac))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("DNA concentration (ng/µL)")+
  ylab("Bacteria read count")+
  theme_bw()+
  labs(tag= "A)")+
  stat_cor(label.x = 300, label.y = 10000) +
  stat_regline_equation(label.x = 300, label.y = 11000)

cor.test(sdt$ReadsBac, sdt$Concentration,
         use="pairwise.complete.obs", method="spearman")

###Main Phylum
##Firmicutes
PS.Firmi <- subset_taxa(PS.bacteria, phylum%in%"Firmicutes")

sdtFirmi <- data.table(as(sample_data(PS.Firmi), "data.frame"),
                       ReadsFirm= sample_sums(PS.Firmi), keep.rownames = T)


sdtFirmi <- dplyr::select(sdtFirmi, 4,66)

sdt <- plyr::join(sdt, sdtFirmi, by= "Mouse_ID")

rm(sdtFirmi)

sdt %>%
  mutate(Firmicutes_abundance = sdt$ReadsFirm/sdt$ReadsBac) -> sdt

##Proteobacteria
PS.Proteo <- subset_taxa(PS.bacteria, phylum%in%"Proteobacteria")

sdtProteo <- data.table(as(sample_data(PS.Proteo), "data.frame"),
                        ReadsProteo= sample_sums(PS.Proteo), keep.rownames = T)

sdtProteo <- dplyr::select(sdtProteo, 4,66)

sdt <- plyr::join(sdt, sdtProteo, by= "Mouse_ID")

rm(sdtProteo)

sdt %>%
  mutate(Proteobacteria_abundance = sdt$ReadsProteo/sdt$TotalReads) -> sdt

##Helicobacter reads 
PS.Heli<- subset_taxa(PS.Proteo, genus%in%"Helicobacter")

sdtHeli <- dplyr::select((data.table(as(sample_data(PS.Heli), "data.frame"),
            ReadsHeli= sample_sums(PS.Heli), keep.rownames = T)),4,65)

sdt <- plyr::join(sdt, sdtHeli, by= "Mouse_ID")

rm(PS.Heli)
rm(sdtHeli)

##Desulfovibrio reads 
PS.Des<- subset_taxa(PS.Proteo, genus%in%"Desulfovibrio")


sdtDes <- dplyr::select((data.table(as(sample_data(PS.Des), "data.frame"),
                                    ReadsDes= sample_sums(PS.Des), keep.rownames = T)),4,65)

sdt <- plyr::join(sdt, sdtDes, by= "Mouse_ID")


rm(PS.Des)
rm(sdtDes)

##Bacteroidetes
PS.Badts <- subset_taxa(PS.bacteria, phylum%in%"Bacteroidetes")

sdtBacdts <- data.table(as(sample_data(PS.Badts), "data.frame"),
                        ReadsBacdts= sample_sums(PS.Badts), keep.rownames = T)


sdtBacdts <- dplyr::select(sdtBacdts, 4,66)

sdt <- plyr::join(sdt, sdtBacdts, by= "Mouse_ID")

rm(sdtBacdts)

sdt %>%
  mutate(Bacteroidetes_abundance = sdt$ReadsBacdts/sdt$ReadsBac) -> sdt

##Deferribacteres
PS.Def <- subset_taxa(PS.bacteria, phylum%in%"Deferribacteres")

sdtDef <- data.table(as(sample_data(PS.Def), "data.frame"),
                     ReadsDef= sample_sums(PS.Def), keep.rownames = T)

sdtDef <- dplyr::select(sdtDef, 4,66)

sdt <- plyr::join(sdt, sdtDef, by= "Mouse_ID")

rm(sdtDef)

sdt %>%
  mutate(Deferribacteres_abundance = sdt$ReadsDef/sdt$ReadsBac) -> sdt

##Tenericutes
PS.Ten <- subset_taxa(PS.bacteria, phylum%in%"Tenericutes")

sdtTen <- data.table(as(sample_data(PS.Ten), "data.frame"),
                     ReadsTen= sample_sums(PS.Ten), keep.rownames = T)

sdtTen <- dplyr::select(sdtTen, 4,66)

sdt <- plyr::join(sdt, sdtTen, by= "Mouse_ID")

rm(sdtTen)

sdt %>%
  mutate(Tenericutes_abundance = sdt$ReadsTen/sdt$ReadsBac) -> sdt

##Actinobacteria
PS.Act <- subset_taxa(PS.bacteria, phylum%in%"Actinobacteria")

sdtAct <- data.table(as(sample_data(PS.Act), "data.frame"),
                     ReadsAct= sample_sums(PS.Act), keep.rownames = T)

sdtAct <- dplyr::select(sdtAct, 4,66)

sdt <- plyr::join(sdt, sdtAct, by= "Mouse_ID")

rm(sdtAct)

sdt %>%
  mutate(Actinobacteria_abundance = sdt$ReadsAct/sdt$ReadsBac) -> sdt

###Plot

e<- ggplot(sdt, aes(HI, Firmicutes_abundance, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Relative abundance of Firmicutes")+
  labs(tag= "A)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")


e2<- ggplot(sdt, aes(HI, Bacteroidetes_abundance, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Relative abundance of Bacteroidetes")+
  labs(tag= "B)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

e3<- ggplot(sdt, aes(HI, Proteobacteria_abundance, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Relative abundance of Proteobacteria")+
  labs(tag= "C)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

e4<- ggplot(sdt, aes(HI, Tenericutes_abundance, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Relative abundance of Tenericutes")+
  labs(tag= "D)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

e5<- ggplot(sdt, aes(HI, Deferribacteres_abundance, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Relative abundance of Deferribacteres")+
  labs(tag= "E)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

e6<- ggplot(sdt, aes(HI, Actinobacteria_abundance, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Relative abundance of Actinobacteria")+
  labs(tag= "F)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

pdf(file = "~/AA_HMHZ/Bacteria_phyla_HI.pdf", width = 15, height = 15)
grid.arrange(e, e2, e3, e4, e5, e6, ncol= 2, nrow= 3)
dev.off()

#########Total reads pathogens #########
###Eimeria 
sdtEim <- data.table(as(sample_data(PS.eimeria), "data.frame"),
                     ReadsEim= sample_sums(PS.eimeria), keep.rownames = T)


sdtEim <- dplyr::select(sdtEim, 4,66)

sdt <- plyr::join(sdt, sdtEim, by= "Mouse_ID")

rm(sdtEim)

sdt %>%
  mutate(Eimeria_abundance = sdt$ReadsEim/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair
###Eimeria Cryptospora 
sdtEimC <- data.table(as(sample_data(PS.eimcyc), "data.frame"),
                     ReadsEimC= sample_sums(PS.eimcyc), keep.rownames = T)


sdtEimC <- dplyr::select(sdtEimC, 4,66)

sdt <- plyr::join(sdt, sdtEimC, by= "Mouse_ID")

rm(sdtEimC)

sdt %>%
  mutate(Eimeria_abundance2 = sdt$ReadsEimC/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair


###Crypto 
sdtCry <- data.table(as(sample_data(PS.Crypto), "data.frame"),
                     ReadsCry= sample_sums(PS.Crypto), keep.rownames = T)


sdtCry <- dplyr::select(sdtCry, 4,66)

sdt <- plyr::join(sdt, sdtCry, by= "Mouse_ID")

rm(sdtCry)

sdt %>%
  mutate(Crypto_abundance = sdt$ReadsCry/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

###Pinworms
sdtPin <- data.table(as(sample_data(PS.pinworms), "data.frame"),
                     ReadsPin= sample_sums(PS.pinworms), keep.rownames = T)


sdtPin <- dplyr::select(sdtPin, 4,65)

sdt <- plyr::join(sdt, sdtPin, by= "Mouse_ID")

rm(sdtPin)

sdt %>%
  mutate(Pin_abundance = sdt$ReadsPin/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

###Syp
sdtSyp <- data.table(as(sample_data(PS.syp), "data.frame"),
                     ReadsSyp= sample_sums(PS.syp), keep.rownames = T)


sdtSyp <- dplyr::select(sdtSyp, 4,65)

sdt <- plyr::join(sdt, sdtSyp, by= "Mouse_ID")

rm(sdtSyp)

sdt %>%
  mutate(Syp_abundance = sdt$ReadsSyp/sdt$TotalReads) -> sdt

###Asp
sdtAsp <- data.table(as(sample_data(PS.asp), "data.frame"),
                     ReadsAsp= sample_sums(PS.asp), keep.rownames = T)


sdtAsp <- dplyr::select(sdtAsp, 4,65)

sdt <- plyr::join(sdt, sdtAsp, by= "Mouse_ID")

rm(sdtAsp)

sdt %>%
  mutate(Asp_abundance = sdt$ReadsAsp/sdt$TotalReads) -> sdt ### Add a new variable 

###Trichomonadida
sdtTricho <- data.table(as(sample_data(PS.Trich), "data.frame"),
                        ReadsTricho= sample_sums(PS.Trich), keep.rownames = T)


sdtTricho <- dplyr::select(sdtTricho, 4,66)

sdt <- plyr::join(sdt, sdtTricho, by= "Mouse_ID")

rm(sdtTricho)

sdt %>%
  mutate(Tricho_abundance = sdt$ReadsTricho/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

###Fornicata
sdtForni <- data.table(as(sample_data(PS.Forni), "data.frame"),
                       ReadsForni= sample_sums(PS.Forni), keep.rownames = T)


sdtForni <- dplyr::select(sdtForni, 4,66)

sdt <- plyr::join(sdt, sdtForni, by= "Mouse_ID")

rm(sdtForni)

sdt %>%
  mutate(Forni_abundance = sdt$ReadsForni/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

###Ciliophora
sdtcili <- data.table(as(sample_data(PS.cili), "data.frame"),
                       Readscili= sample_sums(PS.cili), keep.rownames = T)


sdtcili <- dplyr::select(sdtcili, 4,66)

sdt <- plyr::join(sdt, sdtcili, by= "Mouse_ID")

rm(sdtcili)

sdt %>%
  mutate(Cili_abundance = sdt$Readscili/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair


###Trichuris 

PS.Trich<- subset_taxa(PS.nema, genus%in%"Trichuris")

sdtTrich <- dplyr::select((data.table(as(sample_data(PS.Trich), "data.frame"),
                                     ReadsTrich= sample_sums(PS.Trich), keep.rownames = T)),4,65)

sdt <- plyr::join(sdt, sdtTrich, by= "Mouse_ID")

rm(PS.Trich)
rm(sdtTrich)

table(count=sdt$Trichuris_muris>0, seq=sdt$ReadsTrich>0)

###Parasites (Nematoda, Apicomplexa, PLatihelminths)
sdtpara <- data.table(as(sample_data(PS.para), "data.frame"),
                       Readspara= sample_sums(PS.para), keep.rownames = T)


sdtpara <- dplyr::select(sdtpara, 4,66)

sdt <- plyr::join(sdt, sdtpara, by= "Mouse_ID")

rm(sdtpara)

sdt %>%
  mutate(Para_abundance = sdt$Readspara/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

##Sum all parasites (Apicomplexa, Nematoda, Platihelminths, Ciliophora, Trichomonas, Fornicata)
sdt %>%
  mutate(Para_reads_all = Readspara+ Readscili + ReadsTricho + ReadsForni) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

sdt %>%
  mutate(Para_all_abundance = sdt$Para_reads_all/sdt$TotalReads) -> sdt
###Nematoda
sdtnema <- data.table(as(sample_data(PS.nema), "data.frame"),
                      Readsnema= sample_sums(PS.nema), keep.rownames = T)

sdtnema <- dplyr::select(sdtnema, 4,66)

sdt <- plyr::join(sdt, sdtnema, by= "Mouse_ID")

rm(sdtnema)

sdt %>%
  mutate(Nema_abundance = sdt$Readsnema/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

##Apicomplexa
sdtapi <- data.table(as(sample_data(PS.api), "data.frame"),
                      Readsapi= sample_sums(PS.api), keep.rownames = T)


sdtapi <- dplyr::select(sdtapi, 4,66)

sdt <- plyr::join(sdt, sdtapi, by= "Mouse_ID")

rm(sdtapi)

sdt %>%
  mutate(Api_abundance = sdt$Readsapi/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

##Platyhelminthes
sdtplat <- data.table(as(sample_data(PS.plat), "data.frame"),
                     Readsplat= sample_sums(PS.plat), keep.rownames = T)


sdtplat <- dplyr::select(sdtplat, 4,66)

sdt <- plyr::join(sdt, sdtplat, by= "Mouse_ID")

rm(sdtplat)

sdt %>%
  mutate(Plat_abundance = sdt$Readsplat/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

###Fungi
sdtfungi <- data.table(as(sample_data(PS.fungi), "data.frame"),
                     Readsfungi= sample_sums(PS.fungi), keep.rownames = T)


sdtfungi <- dplyr::select(sdtfungi, 4,66)

sdt <- plyr::join(sdt, sdtfungi, by= "Mouse_ID")

rm(sdtfungi)

sdt %>%
  mutate(Fungi_abundance = sdt$Readsfungi/sdt$TotalReads) -> sdt ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

##Ascomycota
PS.Asco <- subset_taxa(PS.fungi, phylum%in%"Ascomycota")

sdtAsco <- data.table(as(sample_data(PS.Asco), "data.frame"),
                     ReadsAsco= sample_sums(PS.Asco), keep.rownames = T)

sdtAsco <- dplyr::select(sdtAsco, 4,66)

sdt <- plyr::join(sdt, sdtAsco, by= "Mouse_ID")

rm(sdtAsco)

sdt %>%
  mutate(Ascomycota_abundance = sdt$ReadsAsco/sdt$TotalReads) -> sdt

##Basidiomycota
PS.Bas <- subset_taxa(PS.fungi, phylum%in%"Basidiomycota")

sdtBas <- data.table(as(sample_data(PS.Bas), "data.frame"),
                      ReadsBas= sample_sums(PS.Bas), keep.rownames = T)

sdtBas <- dplyr::select(sdtBas, 4,66)

sdt <- plyr::join(sdt, sdtBas, by= "Mouse_ID")

rm(sdtBas)

sdt %>%
  mutate(Basidiomycota_abundance = sdt$ReadsBas/sdt$TotalReads) -> sdt

##Mucoromycota
PS.Muc <- subset_taxa(PS.fungi, phylum%in%"Mucoromycota")

sdtMuc <- data.table(as(sample_data(PS.Muc), "data.frame"),
                     ReadsMuc= sample_sums(PS.Muc), keep.rownames = T)

sdtMuc <- dplyr::select(sdtMuc, 4,66)

sdt <- plyr::join(sdt, sdtMuc, by= "Mouse_ID")

rm(sdtMuc)

sdt %>%
  mutate(Mucoromycota_abundance = sdt$ReadsMuc/sdt$TotalReads) -> sdt


###let's do a psudo bananas
pdf(file = "~/AA_HMHZ/Pseudobananas_Eimeria_Pinworms.pdf", width = 15, height = 15)
f <- ggplot(sdt, aes(HI, Eimeria_abundance, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Relative abundance of Eimeria")+
  labs(tag= "A)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

f2 <- ggplot(sdt, aes(HI, Pin_abundance, color= HI))+
  geom_point(size= 3)+
  geom_smooth(method = loess)+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Relative abundance of Syphacia/Aspiculuris")+
  labs(tag= "B)")+
  theme_bw()+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")

grid.arrange(f, f2, ncol= 1, nrow= 2)
dev.off()

########Correlation reads vs counts######
require(ggpubr)
g <- ggplot(sdt, aes(ReadsPin+1, Pinworms+1))+
  #geom_point()+
  geom_jitter(height = 0.5)+
  geom_smooth(method = lm)+
  scale_x_log10("Sequence reads (Syphacia/Aspiculuris)")+
  scale_y_log10("Counts of Syphacia/Aspiculuris")+
  labs(tag= "A)")+
  theme_bw()+
  stat_cor(label.x = 1000, label.y = 100, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 1000, label.y = 150)

table(count=sdt$Pinworms>0, seq=sdt$ReadsPin>0)

g2 <- ggplot(sdt, aes(ReadsEim, delta_ct_MminusE))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Sequence abundance (Eimeria)")+
  ylab("Relative intensity of Eimeria in tissue")+
  labs(tag= "B)")+
  theme_bw()+
  stat_cor(label.x = 5000, label.y = 3, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 5000, label.y = 4)

g3 <-ggplot(sdt, aes(ReadsEim+1, OPG_Eimeria+1))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  scale_x_log10("Sequence reads count (Eimeria)")+
  scale_y_log10("Oocyst of Eimeria per gram of feces")+
  labs(tag= "B)")+
  theme_bw()+
  stat_cor(label.x = 4500, label.y = 400000, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  stat_regline_equation(label.x = 4500, label.y = 450000)

ggplot(sdt, aes(ReadsEimC, OPG_Eimeria))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Sequence reads count (Eimeria)")+
  ylab("Oocyst of Eimeria per gram of feces")+
  labs(tag= "B)")+
  theme_bw()+
  stat_cor(label.x = 4500, label.y = 400000, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~")))+
  stat_regline_equation(label.x = 4500, label.y = 450000)

g4 <-ggplot(sdt, aes(Crypto_Oocyst_calculated, ReadsCry))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Oocyst of Cryptosporidium per gram of feces")+
  ylab("Sequence reads count (Cryptosporidium)")+
  labs(tag= "B)")+
  theme_bw()+
  stat_cor(label.x = 5000000000, label.y = 500) +
  stat_regline_equation(label.x = 5000000000, label.y = 550)

require(ggpubr)
ggplot(sdt, aes(Syp_abundance, Asp_abundance))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Sequence abundance (Syphacia)")+
  ylab("Sequence abundance (Aspiculuris)")+
  labs(tag= "A)")+
  theme_bw()+
  stat_cor(label.x = 0.15, label.y = 0.3, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 0.15, label.y = 0.32)
####

###By Eimeria species infection

pdf(file = "~/AA_HMHZ/Eimeria_alpha_div.pdf", width = 8, height = 10)
ggplot(sdt, aes(Species_tissue, Chao1))+
  geom_boxplot(fill= c("pink", "#FF3300","#66CC00", "#FFCC00", "gray", "blue", "gray"), colour = "black",
               alpha = 0.8, outlier.colour = "black", outlier.shape = 20)+
  geom_jitter(pch= 21)+
  xlab("Eimeria species identified")+
  ylab("Richness (Chao1 Index)")+
  theme(axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold"))+
  theme_bw()#+
  #geom_violin(trim = F, alpha = 0.1)
dev.off()

pdf(file = "~/AA_HMHZ/Eimeria_abundance_sp.pdf", width = 8, height = 10)
ggplot(sdt, aes(Species_tissue, Eimeria_abundance))+
  geom_boxplot(fill= c("pink", "#FF3300","#66CC00", "#FFCC00", "gray", "blue", "gray"), colour = "black",
               alpha = 0.8, outlier.colour = "black", outlier.shape = 20)+
  geom_jitter(pch= 21)+
  theme(axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold"))+
  theme_bw()#+
dev.off()

ggplot(sdt, aes(Eimeria, Chao1))+
  geom_boxplot(fill= c("red", "darkgreen", "gray"), colour = "black",
               alpha = 0.8, outlier.colour = "black", outlier.shape = 20)+
  geom_jitter(pch= 21)+
  xlab("Eimeria Infection")+
  ylab("Richness (Chao1 Index)")+
  theme(axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold"))+
  theme_bw()#


ggplot(sdt, aes(Eimeria, Shannon))+
  geom_boxplot(fill= c("red", "darkgreen", "gray"), colour = "black",
               alpha = 0.8, outlier.colour = "black", outlier.shape = 20)+
  geom_jitter(pch= 21)+
  xlab("Eimeria Infection")+
  ylab("Diversity (Shannon Index)")+
  theme(axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold"))+
  theme_bw()#

ggplot(sdt, aes(Eimeria, ReadsEim))+
  geom_boxplot()

alphaBac %>%
  mutate(Eimeria = case_when(Species_tissue== "E_ferrisi"| Species_tissue== "E_falciformis"| Species_tissue== "E_vermiformis"| Species_tissue== "Other" ~ "Positive", 
                             Species_tissue== "Negative" ~ "Negative",
                             Flot== "FALSE" & Ap5== "FALSE" ~ "Negative")) -> alphaBac

alphaDiv %>%
  mutate(Eimeria = case_when(Species_tissue== "E_ferrisi"| Species_tissue== "E_falciformis"| Species_tissue== "E_vermiformis"| Species_tissue== "Other" ~ "Positive", 
                             Species_tissue== "Negative" ~ "Negative",
                             Flot== "FALSE" & Ap5== "FALSE" ~ "Negative")) -> alphaDiv


sdt %>%
  mutate(Eimeria = case_when(Species_tissue== "E_ferrisi"| Species_tissue== "E_falciformis"| Species_tissue== "E_vermiformis"| Species_tissue== "Other" ~ "Positive", 
                             Species_tissue== "Negative" ~ "Negative",
                             Flot== "FALSE" & Ap5== "FALSE" ~ "Negative")) -> sdt


ggplot(alphaBac, aes(Eimeria, Shannon))+
  geom_boxplot(fill= c("red", "darkgreen"), colour = "black",
               alpha = 0.8, outlier.colour = "black", outlier.shape = 20)+
  #geom_violin()+
  geom_jitter(pch= 21)+
  xlab("Eimeria Infection")+
  ylab("Diversity (Shannon Index)")+
  theme(axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold"))+
  theme_bw()#


alphaBac %>%
  mutate(H_diminuta = case_when(Hymenolepis_diminuta >= 1 ~ "Positive", 
                                Hymenolepis_diminuta == 0 ~ "Negative")) -> alphaBac

ggplot(alphaBac, aes(H_diminuta, InvSimpson))+
  geom_boxplot(fill= c("red", "gray"), colour = "black",
               alpha = 0.8, outlier.colour = "black", outlier.shape = 20)+
  geom_jitter(pch= 21)+
  xlab("Hymenolepis diminuta Infection")+
  ylab("Diversity (Inverse Simpson Index)")+
  theme(axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold"))+
  theme_bw()#

alphaBac %>%
  mutate(H_micostoma = case_when(Hymenolepis_microstoma >= 1 ~ "Positive", 
                                Hymenolepis_microstoma == 0 ~ "Negative")) -> alphaBac

ggplot(alphaBac, aes(H_micostoma, InvSimpson))+
  geom_boxplot(fill= c("red", "darkgreen", "gray"), colour = "black",
               alpha = 0.8, outlier.colour = "black", outlier.shape = 20)+
  geom_jitter(pch= 21)+
  xlab("Hymenolepis microstoma presence")+
  ylab("Diversity (Inverse Simpson Index)")+
  theme(axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold"))+
  theme_bw()#

alphaBac %>%
  mutate(T_muris = case_when(Trichuris_muris >= 1 ~ "Positive", 
                                Trichuris_muris == 0 ~ "Negative")) -> alphaBac

ggplot(alphaBac, aes(T_muris, InvSimpson))+
  geom_boxplot(fill= c("red", "darkgreen", "gray"), colour = "black",
               alpha = 0.8, outlier.colour = "black", outlier.shape = 20)+
  geom_jitter(pch= 21)+
  xlab("Trichuris muris presence")+
  ylab("Diversity (Inverse Simpson Index)")+
  theme(axis.text= element_text(size = 18), axis.title=element_text(size=18,face="bold"))+
  theme_bw()#


alphaBac %>%
  mutate(Asp = case_when(Aspiculuris_tetraptera >= 1 ~ "Positive", 
                             Aspiculuris_tetraptera == 0 ~ "Negative")) -> alphaBac


alphaBac %>%
  mutate(Syp = case_when(Syphacia_obvelata >= 1 ~ "Positive", 
                         Syphacia_obvelata == 0 ~ "Negative")) -> alphaBac

##Stats
wilcox.exact(Chao1 ~ Eimeria, data = alphaBac, exact = TRUE)
wilcox.exact(Shannon ~ Eimeria, data = alphaBac, exact = TRUE)


###########Composition#########
library("microbiome")

PS.Core.Bac <- core(PSHigh, detection= 0, prevalence = 50/100)

plot_composition(PS.bacteria, average_by = "Genotype", transform = "compositional")

PS.Comp <- PS.Core %>%
  aggregate_taxa(level= "phylum") %>%
  microbiome::transform(transform = "compositional")

plot_composition(microbiome::transform(PSHigh, "compositional"),
                 plot.type = "heatmap",
                 sample.sort = "neatmap",
                 otu.sort = "neatmap")

microbiome::dominance(PSHigh)

ggplot(sdt, aes(Flot, ReadsEim, color= Species_tissue))+
  geom_jitter(size= 3, aes(alpha= 0.5))+
  xlab("Flotation status")+
  ylab("Reads assigned as Eimeria")+
  labs(tag= "A)")+
  theme_bw()

ggplot(sdt, aes(Ap5, ReadsEim, color= Species_tissue))+
  geom_jitter(size= 3, aes(alpha= 0.5))+
  xlab("Ap5 status")+
  ylab("Reads assigned as Eimeria")+
  labs(tag= "A)")+
  theme_bw()

h<- ggplot(sdt, aes(BMI, ReadsBac))+
  geom_jitter(size= 3, aes(alpha= 0.5))+
  xlab("Body mass index")+
  ylab("Reads count Bacteria")+
  labs(tag= "A)")+
  theme_bw()

h2<- ggplot(sdt, aes(BMI, Bac_abundance, color= Sex))+
  geom_jitter(size= 3, aes(alpha= 0.5))+
  xlab("Body mass index")+
  ylab("Relative abundance Bacteria")+
  labs(tag= "B)")+
  theme_bw()

#sdt %>%
#  mutate(Firm_Bac = sdt$ReadsFirm/sdt$ReadsBacdts) -> sdt

h3 <- ggplot(sdt, aes(BMI, Firm_Bac, color= Sex))+
  geom_jitter(size= 3, aes(alpha= 0.5))+
  xlab("Body mass index")+
  ylab("Firmicutes/Bacteroidetes ratio")+
  labs(tag= "C)")+
  theme_bw()

h4<- ggplot(sdt, aes(Latitude, Firm_Bac, color= Sex))+
  geom_jitter(size= 3, aes(alpha= 0.5))+
  xlab("Latitude")+
  ylab("Firmicutes/Bacteroidetes ratio")+
  labs(tag= "D)")+
  theme_bw()

pdf(file = "~/AA_HMHZ/BMI.pdf", width = 8, height = 10)
grid.arrange(h,h2, h3, h4, ncol= 2, nrow= 2)
dev.off()

i<- ggplot(sdt, aes(Latitude, Proteobacteria_abundance, color= HI))+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  geom_jitter(size= 3, aes(alpha= 0.5))+
  geom_smooth(method = lm)+
  xlab("Latitude")+
  ylab("Relative abundance Proteobacteria")+
  facet_wrap(~Transect, 2) +
  labs(tag= "A)")+
  theme_bw()

i2<- ggplot(sdt, aes(Longitude, HI, color= HI))+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  geom_jitter(size= 3, aes(alpha= 0.5))+
  xlab("Longitude")+
  ylab("Hybrid Index")+
  facet_wrap(~Transect, 2) +
  labs(tag= "B)")+
  theme_bw()

pdf(file = "~/AA_HMHZ/Transect.pdf", width = 8, height = 10)
grid.arrange(i,i2, ncol= 1, nrow= 2)
dev.off()

###Replicating the analysis strategy from Wang 2015 with our wild samples

###Beta diversity (just with bray-curtis dissimilarity)

ordi <- ordinate(PSRare, method="PCoA", distance="bray")
dis <- phyloseq::distance(PSRare, method="bray")

c1<- plot_ordination(PSRare, ordi, color="Genotype")+
  scale_color_manual(values = c("purple","blue","red"))+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "B)")

adonis(dis ~ HI+Genotype+Longitude+Latitude+Transect+Year+Chip_number+Concentration+Seq_Run, data = sdt)

pdf(file = "~/AA_HMHZ/Beta_Div_HI_Genotype.pdf", width = 8, height = 10)
grid.arrange(c,c1, ncol= 1, nrow= 2)
dev.off()

c1.1 <- plot_ordination(PSRare, ordi, color="Transect")+
        theme_bw()+
        geom_point(size=5, alpha= 0.75)+
        labs(tag= "C)")

pdf(file = "~/AA_HMHZ/Beta_Div_HI_Gen_Tra.pdf", width = 8, height = 10)
grid.arrange(c,c1, c1.1, ncol= 1, nrow= 3)
dev.off()

j<- ggplot(sdt, aes(ReadsSyp, ReadsHeli))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Sequence reads count (Syphacia)")+
  ylab("Sequence reads count (Helicobacter)")+
  labs(tag= "A)")+
  theme_bw()+
  stat_cor(label.x = 3500, label.y = 2000, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 3500, label.y = 2100)

j2<- ggplot(sdt, aes(ReadsSyp, ReadsEim))+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Sequence reads count (Syphacia)")+
  ylab("Sequence reads count (Eimeria)")+
  labs(tag= "B)")+
  theme_bw()+
  stat_cor(label.x = 3500, label.y = 6000, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 3500, label.y = 6300)

pdf(file = "~/AA_HMHZ/Syp-Heli-Eim_correlation.pdf", width = 16, height = 8)
grid.arrange(j,j2, ncol= 2, nrow= 1)
dev.off()

###Principal component regression analysis 
require("pls")


