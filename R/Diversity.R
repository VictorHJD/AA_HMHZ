###Phyloseq pipeline for alpha and beta diversity analysis
###Preparation of Phyloseq object with final data for analysis 
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

##Functions
source("GitProjects/AA_HMHZ/R/Functions.R")

###Data 
source("GitProjects/AA_HMHZ/R/PS_Samples_HMHZ.R")

##############Rarefy##############
###To the minumim reads
PSRare <- rarefy_even_depth(PSHigh, rngseed=1, sample.size= min(sample_sums(PSHigh)), replace=F)

summarize_phyloseq(PSRare)
nsamples(PSRare)
ntaxa(PSRare)
sum(otu_table(PSRare))
sum(otu_table(subset_taxa(PSRare, superkingdom %in% "Bacteria")))/sum(otu_table(PSRare))
sum(otu_table(subset_taxa(PSRare, superkingdom %in% "Eukaryota")))/sum(otu_table(PSRare))

table(tax_table(PSRare)[, "superkingdom"], exclude = NULL)
rarcurv4 <- vegan::rarecurve(otu_table(PSRare),
                             label = F)

rarecounts_HMHZ <- data.frame(rowSums(otu_table(PSRare)))
rarecounts_HMHZ[,2] <- rownames(rarecounts_HMHZ)

###Bacterial and Eukaryotic counts
rarecounts_HMHZ[,3] <- as.data.frame(rowSums(otu_table(subset_taxa(PSRare, superkingdom%in%"Bacteria"))))
rarecounts_HMHZ[,4] <- as.data.frame(rowSums(otu_table(subset_taxa(PSRare, superkingdom%in%"Eukaryota"))))
colnames(rarecounts_HMHZ) <- c("Rarefaction_counts", "Mouse_ID", "Bacteria_reads", "Eukaryota_reads")
rownames(rarecounts_HMHZ) <- c(1:nrow(rarecounts_HMHZ))
rarecounts_HMHZ <- data.frame(Mouse_ID = rarecounts_HMHZ$Mouse_ID, 
                              Rarefaction_counts = rarecounts_HMHZ$Rarefaction_counts,
                              Bacteria_reads= rarecounts_HMHZ$Bacteria_reads,
                              Eukaryota_reads= rarecounts_HMHZ$Eukaryota_reads)

#####################Separate Eukaryotic and Bacterial###########
PS.bacteria.H <- subset_taxa(PSHigh, superkingdom%in%"Bacteria") ###For Beta diversity
PS.bacteria.R <- subset_taxa(PSRare, superkingdom%in%"Bacteria") ##For Alpha diversity

PS.eukaryota.H <- subset_taxa(PSHigh, superkingdom%in%"Eukaryota") 
PS.eukaryota.R <- subset_taxa(PSRare, superkingdom%in%"Eukaryota") 

######################Separate by transect################
PS.HZ_BR.H<- subset_samples(PSHigh, Transect%in%"HZ_BR")
PS.HZ_BR.R<- subset_samples(PSRare, Transect%in%"HZ_BR")

PS.HZ_BAV.H<- subset_samples(PSHigh, Transect%in%"HZ_BR")
PS.HZ_BAV.R<- subset_samples(PSRare, Transect%in%"HZ_BR")

##################Extract sample information##############

sdt.bacteria.H <- data.table(as(sample_data(PS.bacteria.H), "data.frame"),
                  TotalReadsH= sample_sums(PS.bacteria.H), keep.rownames = T)

sample_data(PS.bacteria.H)$TotalReadsH <- as.numeric(sdt.bacteria.H$TotalReadsH)
setnames(sdt.bacteria.H, "rn", "Sample_ID")
colnames(sdt.bacteria.H)[16]<- "Mastaphorus_fdx"
colnames(sdt.bacteria.H)[64]<- "Eimeria_dx"

sdt.bacteria.R <- data.table(as(sample_data(PS.bacteria.R), "data.frame"),
                             TotalReadsR= sample_sums(PS.bacteria.R), keep.rownames = T)

sample_data(PS.bacteria.R)$TotalReadsR <- as.numeric(sdt.bacteria.R$TotalReadsR)
setnames(sdt.bacteria.R, "rn", "Sample_ID")
colnames(sdt.bacteria.R)[16]<- "Mastaphorus_fdx"
colnames(sdt.bacteria.R)[64]<- "Eimeria_dx"

sdt.eukaryota.H <- data.table(as(sample_data(PS.eukaryota.H), "data.frame"),
                             TotalReadsH= sample_sums(PS.eukaryota.H), keep.rownames = T)

sample_data(PS.eukaryota.H)$TotalReadsH <- as.numeric(sdt.eukaryota.H$TotalReadsH)
setnames(sdt.eukaryota.H, "rn", "Sample_ID")
colnames(sdt.eukaryota.H)[16]<- "Mastaphorus_fdx"
colnames(sdt.eukaryota.H)[64]<- "Eimeria_dx"

sdt.eukaryota.R <- data.table(as(sample_data(PS.eukaryota.R), "data.frame"),
                             TotalReadsR= sample_sums(PS.eukaryota.R), keep.rownames = T)

sample_data(PS.eukaryota.R)$TotalReadsR <- as.numeric(sdt.eukaryota.R$TotalReadsR)
setnames(sdt.eukaryota.R, "rn", "Sample_ID")
colnames(sdt.eukaryota.R)[16]<- "Mastaphorus_fdx"
colnames(sdt.eukaryota.R)[64]<- "Eimeria_dx"

#####################Alpha diversity###################### 
###Richness and diveristy
alphaDiv <- as.data.frame(estimate_richness(PSRare, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
alphaDiv$Sample_ID <- rownames(alphaDiv)
sdt <- plyr::join(sdt, alphaDiv, by= "Sample_ID")

alphaDiv.bacteria <- as.data.frame(estimate_richness(PS.bacteria.R, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
alphaDiv.bacteria$Sample_ID <- rownames(alphaDiv.bacteria)
sdt.bacteria.R <- plyr::join(sdt.bacteria.R, alphaDiv.bacteria, by= "Sample_ID")

alphaDiv.eukaryota <- as.data.frame(estimate_richness(PS.eukaryota.R, measures=c("Observed", "InvSimpson", "Shannon", "Chao1")))
alphaDiv.eukaryota$Sample_ID <- rownames(alphaDiv.eukaryota)
sdt.eukaryota.R <- plyr::join(sdt.eukaryota.R, alphaDiv.eukaryota, by= "Sample_ID")

###Evenness
eveAll <- evenness(PSRare, c("pielou", "simpson", "evar", "bulla"))
eveAll$Sample_ID <- rownames(eveAll)
sdt <- plyr::join(sdt, eveAll, by= "Sample_ID")

eveAll.bacteria <- evenness(PS.bacteria.R, c("pielou", "simpson", "evar", "bulla"))
eveAll.bacteria$Sample_ID <- rownames(eveAll.bacteria)
sdt.bacteria.R <- plyr::join(sdt.bacteria.R, eveAll.bacteria, by= "Sample_ID")

eveAll.eukaryota <- evenness(PS.eukaryota.R, c("pielou", "simpson", "evar", "bulla"))
eveAll.eukaryota$Sample_ID <- rownames(eveAll.eukaryota)
sdt.eukaryota.R <- plyr::join(sdt.eukaryota.R, eveAll.eukaryota, by= "Sample_ID")

rm(eveAll, alphaDiv)
rm(eveAll.bacteria, eveAll.eukaryota, alphaDiv.bacteria, alphaDiv.eukaryota)

#######################Beta diversity#######################
##Ordination Bray-Curtis
ordi <- ordinate(PSHigh, method="PCoA", distance="bray")
ordi2 <- ordinate(PSHigh, method="NMDS", distance="bray")
PS.bacteria.H<- prune_samples(sample_sums(PS.bacteria.H)>0, PS.bacteria.H) ##eliminate 2 samples with 0 bacterial taxa
ordi.bacteria <- ordinate(PS.bacteria.H, method="PCoA", distance="bray")
ordi.eukaryota <- ordinate(PS.eukaryota.H, method="PCoA", distance="bray")

##Distance matrix
dis <- phyloseq::distance(PSHigh, method="bray")
dis.bacteria <- phyloseq::distance(PS.bacteria.H, method="bray")
dis.eukaryota <- phyloseq::distance(PS.eukaryota.H, method="bray")

##PCoA plot
beta.plot <- plot_ordination(PSHigh, ordi, color="HI")+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybrid Index)")+
  theme_classic(base_size = 15, base_family = "Helvetica")+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

beta.plot.chip <- plot_ordination(PSHigh, ordi, color="HI")+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  theme_classic(base_size = 15, base_family = "Helvetica")+
  facet_wrap(~Chip_number) +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

ggsave(file = "~/AA_HMHZ/Beta_div/Beta_diversity_HMHZ_High_Chip.pdf", plot = beta.plot.chip, width = 10, height = 8)

plot_ordination(PS.eukaryota.H, ordi.eukaryota, color="HI")+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  theme_classic(base_size = 15, base_family = "Helvetica")+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

plot_ordination(PS.bacteria.H, ordi.bacteria, color="HI")+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  theme_classic(base_size = 15, base_family = "Helvetica")+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

ggsave(file = "~/AA_HMHZ/Beta_div/Beta_diversity_HMHZ_High.pdf", plot = beta.plot, width = 10, height = 8)

beta.plot2<-  plot_ordination(PSHigh, ordi, color="Genotype")+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)+
  scale_color_manual(values =  c("purple", "blue", "red"))+
  labs(tag= "B)")

ggsave(file = "~/AA_HMHZ/Beta_div/Beta_diversity_HMHZ_High_cate.pdf", plot = beta.plot2, width = 10, height = 8)

plot_ordination(PSHigh, ordi, color="Sex")+
  theme_bw()+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "B)")

plot_ordination(PSHigh, ordi, color="Concentration")+
  theme_bw() +
  geom_point(size=5, alpha= 0.75)+
  scale_color_gradient(low="green", high="red")+
  labs(tag= "C)")

plot_ordination(PSHigh, ordi, color="Chip_number")+
  theme_bw() +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "D)")

plot_ordination(PSHigh, ordi, color="Year")+
  theme_bw() +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "E)")

plot_ordination(PSHigh, ordi, color="Seq_Run")+
  theme_bw() +
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "F)")

plot_ordination(PSHigh, ordi2, color="HI")+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  theme_classic(base_size = 15, base_family = "Helvetica")+
  geom_point(size=5, alpha= 0.75)+
  labs(tag= "A)")

##Adonis 
####Permanova
adonis(dis ~ HI+Longitude+Latitude+Year+Chip_number+Concentration+Seq_Run, data = sdt, permutations = 9999)
adonis(dis ~ HI+Genotype+Year+Longitude+Latitude+Transect, data = sdt, permutations = 9999)

adonis(dis.bacteria ~ HI+Longitude+Latitude+Year+Chip_number+Concentration+Seq_Run, data = sdt.bacteria.H, permutations = 9999)
adonis(dis.bacteria ~ HI+Genotype+Year+Longitude+Latitude+Transect, data = sdt.bacteria.H, permutations = 9999)

adonis(dis.eukaryota ~ HI+Longitude+Latitude+Year+Chip_number+Concentration+Seq_Run, data = sdt.eukaryota.H,  permutations = 9999)
adonis(dis.eukaryota ~ HI+Genotype+Year+Longitude+Latitude+Transect, data = sdt.eukaryota.H, permutations = 9999)

##Multivariate extension of GLMs
install.packages("mvabund")
library("mvabund")
