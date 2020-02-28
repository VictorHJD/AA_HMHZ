###Phyloseq pipeline for diversity analysis
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
library(tidyverse)
library(geosphere)

##Functions
source("GitProjects/AA_HMHZ/R/Functions.R")

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
summarize_phyloseq(PS)
####Taxa detected
as.data.frame(table(tax_table(PS)[, "phylum"]))
as.data.frame(table(tax_table(PS)[, "genus"]))

###Rarefaction curve 
#rarcurv <- vegan::rarecurve(otu_table(PS),
#                            label = F)

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

#rarcurv2 <- vegan::rarecurve(otu_table(PSHigh),
#                             label = F)

## Merge ASVs that have the same taxonomy at a certain taxonomic rank (in this case phylum and genus)
PS.Gen <-  tax_glom(PSHigh, "genus", NArm = F)
summarize_phyloseq(PS.Gen)

PS.Fam<-  tax_glom(PSHigh, "family", NArm = F)
summarize_phyloseq(PS.Fam)

PS.Ord <-  tax_glom(PSHigh, "order", NArm = F)
summarize_phyloseq(PS.Ord)

PS.Phy <-  tax_glom(PSHigh, "phylum", NArm = TRUE)
summarize_phyloseq(PS.Phy)

###Summarize sequencing depths 
sdt <- data.table(as(sample_data(PSHigh), "data.frame"),
                  TotalReads= sample_sums(PSHigh), keep.rownames = T)

sample_data(PSHigh)$TotalReads <- as.numeric(sdt$TotalReads)
setnames(sdt, "rn", "Sample_ID")
colnames(sdt)[16]<- "Mastaphorus_fdx"
colnames(sdt)[64]<- "Eimeria_dx"

## Sum asvCount by taxa 
HMHZPhy<- summarize_taxa(PS.Phy, "phylum", "Mouse_ID")
phylnames<- unique(HMHZPhy$phylum)

HMHZOrd<- summarize_taxa(PS.Ord, "order", "Mouse_ID")
ordnames<- unique(HMHZOrd$order)

HMHZFam<- summarize_taxa(PS.Fam, "family", "Mouse_ID")
famnames<- unique(HMHZFam$family)

HMHZGen<- summarize_taxa(PS.Gen, "genus", "Mouse_ID")
genunames<- unique(HMHZGen$genus)

###Plot relative abundances by taxa level
theme_set(theme_classic(base_size = 15, base_family = "Helvetica"))

## Phylum
HMHZPhy_sorted <- HMHZPhy %>%
  mutate(phylum = fct_reorder(phylum, -meanRA))

HMHZPhy_sorted%>%
  group_by(phylum)%>%
  mutate(phylum_avg= mean(meanRA))%>%
  mutate(Main_phylum= phylum_avg>=0.01)->HMHZPhy_sorted 

HMHZPhy_sorted%>%
  summarise(no_rows = length(phylum))%>%
  mutate(Prev_phylum= (no_rows/438)*100)%>%
  mutate(High_prev= Prev_phylum>=30)-> HMHZPhy_prev

HMHZPhy_sorted<- merge(HMHZPhy_sorted, HMHZPhy_prev, by= "phylum")

phy_levels<- c("Proteobacteria", "Bacteroidetes", "Firmicutes",
              "Apicomplexa", "Nematoda", "Basidiomycota", "Ascomycota")

HMHZPhy_sorted$phylum<- factor(HMHZPhy_sorted$phylum, levels = phy_levels)

HMHZPhy_sorted%>%
filter(phylum %in% c("Firmicutes", "Bacteroidetes", "Proteobacteria", "Nematoda", "Apicomplexa",
                      "Ascomycota", "Basidiomycota")) ->HMHZPhy_sorted

sample_avg <- HMHZPhy %>%
  summarize(avg = mean(meanRA, na.rm = T)) %>%
  pull(avg)

g <- ggplot(subset(HMHZPhy_sorted, High_prev==T&Main_phylum==T), aes(x = phylum, y = meanRA, color= phylum)) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
  labs(x = NULL, y = "Relative abundance") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 10),
        panel.grid = element_blank())

set.seed(123)
a<- g + 
  #geom_hline(aes(yintercept = sample_avg), color = "gray70", size = 0.6) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.2)+
  stat_summary(fun.y = mean, geom = "point", size = 2, color= "black")+
  labs(tag = "A)")

ggsave(file = "~/AA_HMHZ/ML_Analysis/Main_taxa_A.pdf", plot = a, width = 10, height = 8)

##Order
HMHZOrd_sorted <- HMHZOrd %>%
  mutate(order = fct_reorder(order, -meanRA))

HMHZOrd_sorted%>%
  group_by(order)%>%
  mutate(order_avg= mean(meanRA))%>%
  mutate(Main_order= order_avg>=0.05)->HMHZOrd_sorted ###Select Orders with average higher relative abundace 

sample_avg <- HMHZOrd %>%
  summarize(avg = mean(meanRA, na.rm = T)) %>%
  pull(avg)

g <- ggplot(subset(HMHZOrd_sorted, Main_order==T), aes(x = order, y = meanRA, color= order)) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
  labs(x = NULL, y = "Relative abundance") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 10),
        panel.grid = element_blank())

set.seed(123)
b<- g + 
  geom_hline(aes(yintercept = sample_avg), color = "gray70", size = 0.6) +
  geom_jitter(size = 2, alpha = 0.05, width = 0.2)+
  stat_summary(fun.y = mean, geom = "point", size = 2, color= "black")+
  labs(tag = "B)")

##Family
HMHZFam_sorted <- HMHZFam %>%
  mutate(family = fct_reorder(family, -meanRA))

HMHZFam_sorted%>%
  group_by(family)%>%
  mutate(family_avg= mean(meanRA))%>%
  mutate(Main_family= family_avg>=0.05)->HMHZFam_sorted

sample_avg <- HMHZFam %>%
  summarize(avg = mean(meanRA, na.rm = T)) %>%
  pull(avg)

g <- ggplot(subset(HMHZFam_sorted, Main_family==T), aes(x = family, y = meanRA, color= family)) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
  labs(x = NULL, y = "Relative abundance") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 10),
        panel.grid = element_blank())

set.seed(123)
c<- g + 
  geom_hline(aes(yintercept = sample_avg), color = "gray70", size = 0.6) +
  geom_jitter(size = 2, alpha = 0.05, width = 0.2)+
  stat_summary(fun.y = mean, geom = "point", size = 2, color= "black")+
  labs(tag = "C)")

##Genus
HMHZGen_sorted <- HMHZGen %>%
  mutate(genus = fct_reorder(genus, -meanRA))

HMHZGen_sorted%>%
  group_by(genus)%>%
  mutate(genus_avg= mean(meanRA))%>%
  mutate(Main_genus= genus_avg>=0.01)->HMHZGen_sorted 

HMHZGen_sorted%>%
  summarise(no_rows = length(genus))%>%
  mutate(Prev_genus= (no_rows/438)*100)%>%
  mutate(High_prev= Prev_genus>=30)-> HMHZGen_prev

HMHZGen_sorted<- merge(HMHZGen_sorted, HMHZGen_prev, by= "genus")

HMHZGen_sorted%>%
  filter(genus %in% c("Kazachstania", "Aspiculuris", "Eimeria", "Syphacia", "Tritrichomonas", "Bacteroides", "Helicobacter", "Blautia"))->HMHZGen_sorted 

gen_levels<- c("Blautia", "Helicobacter", "Bacteroides", "Tritrichomonas", 
                "Syphacia","Eimeria","Aspiculuris","Kazachstania")

HMHZGen_sorted$genus<- factor(HMHZGen_sorted$genus, levels = gen_levels)

sample_avg <- HMHZGen %>%
  summarize(avg = mean(meanRA, na.rm = T)) %>%
  pull(avg)

g <- ggplot(HMHZGen_sorted, aes(x = genus, y = meanRA, color= genus)) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), expand = c(0.005, 0.005)) +
  labs(x = NULL, y = "Relative abundance") +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text.x = element_text(family = "Helvetica", size = 10),
        panel.grid = element_blank())

set.seed(123)
d<- g + 
  #geom_hline(aes(yintercept = sample_avg), color = "gray70", size = 0.6) +
  geom_jitter(size = 2, alpha = 0.5, width = 0.2)+
  stat_summary(fun.y = mean, geom = "point", size = 2, color= "black")+
  labs(tag = "D)")

#ggsave(file = "~/AA_HMHZ/ML_Analysis/Main_taxa_B.pdf", plot = d, width = 10, height = 8)

e<- grid.arrange(a, b, c, d, ncol= 2, nrow= 2)

ggsave(file = "~/AA_HMHZ/ML_Analysis/Main_taxa.pdf", plot = e, width = 15, height = 10)

####Analysis### (Move later to a different R script)
###Square root transformation of the data to assess Bray-Curtis dissimilarity
##Phylum
HMHZPhy%>%
  mutate(transRA= sqrt(meanRA))%>%
  select(1,2,7) %>%
  group_by(Mouse_ID)%>%
  spread(phylum, transRA)-> Phylum.tRA 

Phylum.tRA[is.na(Phylum.tRA)]<- as.numeric(0.000000000)
Phylum.tRA<- column_to_rownames(Phylum.tRA, var= "Mouse_ID")

##Genus
HMHZGen%>%
  mutate(transRA= sqrt(meanRA))%>%
  select(1,2,7) %>%
  group_by(Mouse_ID)%>%
  spread(genus, transRA)-> Genus.tRA 

Genus.tRA[is.na(Genus.tRA)]<- as.numeric(0.000000000)
Genus.tRA<- column_to_rownames(Genus.tRA, var= "Mouse_ID")

##Check for correlation between Bray Curtis dissimilarity in composition and geographical distance
##Create a matrix with the Longitude and latitude information
sdt%>%
  select_(1,7,8)-> geo.dist
geo.dist<- column_to_rownames(geo.dist, var= "Sample_ID")
geo.dist<- as.matrix(geo.dist)
geo.dist<- distm(geo.dist, fun = distGeo)

##Plot BC dissimiliraty and geographic distance 
combi.dist<- data.frame(Geo= as.vector(geo.dist)/1000,
                        Bray= as.vector(Genus.bctRA))

geo.plot<- ggplot(combi.dist, aes(Geo, Bray)) +
  geom_jitter(alpha=0.03, width=0.3, height=0, color= "lightblue") + 
  stat_smooth(se=T, method="lm") +
  scale_x_continuous("Geographical distance (km)") +
  scale_y_continuous("Bray-Curtis dissimilarity between samples") +
  labs(tag = "B)")+
  theme_classic(base_size = 15, base_family = "Helvetica")

#ggsave(file = "~/AA_HMHZ/Beta_div/BC_Geo_dist_HMHZ.pdf", plot = geo.plot, width = 10, height = 8)

##For genetic distance based on HI differences
sdt%>%
  select_(1,41)-> gen.dist

gen.dist<- column_to_rownames(gen.dist, var= "Sample_ID")
gen.dist<- as.matrix(gen.dist)
gen.dist<- dist(gen.dist)

combi.dist2<- data.frame(Bray= as.vector(dis),
                        Gen= as.vector(gen.dist))

gen.plot<- ggplot(combi.dist2, aes(Gen, Bray)) +
  geom_jitter(alpha=0.03, width=0.3, height=0, color= "pink") + 
  stat_smooth(se=T, method="lm") +
  scale_x_continuous("Genetic distance (Hybrid Index)") +
  scale_y_continuous("Bray-Curtis dissimilarity between samples") +
  labs(tag = "C)")+
  theme_classic(base_size = 15, base_family = "Helvetica")

#ggsave(file = "~/AA_HMHZ/Beta_div/BC_Gen_dist_HMHZ.pdf", plot = gen.plot, width = 10, height = 8)


##Create a matrix for BC dissimilarity 
##Phylum
Phylum.matRA<- as.matrix(Phylum.tRA)
Phylum.bctRA<- vegdist(Phylum.matRA, method = "bray")
Phylum.bctRA<- as.matrix(Phylum.bctRA)

##Genus
Genus.matRA<- as.matrix(Genus.tRA)
Genus.bctRA<- vegdist(Genus.matRA, method = "bray")
Genus.bctRA<- as.matrix(Genus.bctRA)

adonis(Genus.bctRA ~ HI+Genotype+Year+Longitude+Latitude+Transect, data = Genus.samples)

###Mantel test
##Phylum
bc.geo<- mantel(Phylum.bctRA, geo.dist, method = "spearman", permutations = 9999, na.rm = T)

##Genus
bc.geo<- mantel(Genus.bctRA, geo.dist, method = "spearman", permutations = 9999, na.rm = T)

##Surprise surprise Bray-Curtis dissimilarity has not a significant relationship with the geographical separation of the samples
## The more separated the mice are, their bacterial/eukaryotic communities (at Phylum level) do not become more dissimilar 

##Non-metric Multidimensional Scaling (NMDS)
NMDS.scree(Phylum.bctRA)
set.seed(2)
# Here, we perform the final analysis and check the result
NMDS1 <- metaMDS(Phylum.bctRA, k = 5, trymax = 1000, trace = F)
stressplot(NMDS1)
orditorp(NMDS1, display = "sites")

##Create matrix for PCA and NMDS
HMHZPhy%>%
  select_(1,2,3) %>%
  group_by(Mouse_ID)%>%
  spread(phylum, meanRA)-> Phylum.RA ##RA= relative abundances

Phylum.RA[is.na(Phylum.RA)]<- as.numeric(0)
Phylum.RA<-as.data.frame(Phylum.RA)
rownames(Phylum.RA)<-Phylum.RA[,1]
Phylum.samples<-join(Phylum.RA, sdt, by= "Mouse_ID")
rownames(Phylum.samples)<-Phylum.samples[,1]
Phylum.RA[,1]<- NULL
Phylum.samples[,1]<- NULL ##Phylum.samples is the file required for the Maximum likelihood analysis

NMDS2<- metaMDS(Phylum.RA, k=5, trymax = 1000, autotransform = FALSE, distance="bray")
stressplot(NMDS2)

plot(NMDS2, display = "species")
points(NMDS2, display = "species", col= "green", cex= 1.5, pch= 20)
text(NMDS2, display = "species", pos= 1, cex= 0.75)

plot(NMDS2, display = "sites")
points(NMDS2, display = "sites", col= "purple", cex= 1.5, pch= 20)
text(NMDS2, display = "species", pos= 1, cex= 0.75)

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
  guides(fill= guide_legend(nrow = 11))

require("factoextra")
require("FactoMineR")
Phylum.pca<- PCA(Phylum.RA, graph = T)
Phylum.eig<- get_eigenvalue(Phylum.pca)
fviz_eig(Phylum.pca, addlabels = TRUE, ylim = c(0, 25))

fviz_pca_ind(Phylum.pca,
             geom.ind = "point", # show points only (but not "text")
             col.ind = Phylum.samples$HI, # color by groups
             gradient.cols = c("blue", "red"),
             alpha.ind = 0.7,
             addEllipses = F, # Concentration ellipses
             legend.title = "Hybrid \nIndex")+
  labs(tag = "A)")

fviz_pca_var(Phylum.pca, col.var = "cos2", select.var =  list(contrib = 10), ##top ten taxa contributing
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE )+ # Avoid text overlapping
  labs(tag = "B)")

###GLM
#Can the relative abundance at Phylum level predict the genotype of host?

Phylum.samples -> Phylum.model
Phylum.model<- subset(Phylum.model, select = -c(31:69))
Phylum.model<- subset(Phylum.model, select = c(1:31))

modelPhy <- glm(HI~ ., data = Phylum.model)
summary(modelPhy)

##Genus all
HMHZGen%>%
  select(1,2,3) %>%
  group_by(Mouse_ID)%>%
  spread(genus, meanRA)-> Genus.RA ##RA= relative abundances

Genus.RA[is.na(Genus.RA)]<- as.numeric(0)
Genus.RA<-as.data.frame(Genus.RA)
rownames(Genus.RA)<-Genus.RA[,1]
Genus.samples<-join(Genus.RA, sdt, by= "Mouse_ID")
rownames(Genus.samples)<-Genus.samples[,1]
Genus.RA[,1]<- NULL
Genus.samples[,1]<- NULL

Genus.samples -> Genus.model
Genus.model<- subset(Genus.model, select = -c(644:682))
Genus.model<- subset(Genus.model, select = c(1:644))

modelGen <- glm(HI~ ., data = Genus.model)
summary(modelGen)
plot_model(modelPhy)

###Correlation 
require(ggpubr)
ggplot(Phylum.samples, aes(Nematoda, Proteobacteria, color= HI))+
  scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
  geom_point(size= 3)+
  geom_smooth(method = lm)+
  xlab("Relative abundance (Nematoda)")+
  ylab("Relative abundance (Proteobacteria)")+
  #labs(tag= "A)")+
  theme_bw()+
  stat_cor(label.x = 0.5, label.y = 0.25, aes(label= paste(..rr.label.., ..p.label.., sep= "~`,`~"))) +
  stat_regline_equation(label.x = 0.5, label.y = 0.27)

###Prevalence estimation
###Kazachstania
(sum(Genus.samples$Kazachstania>0)/sum(Genus.samples$Kazachstania>=0))*100

###Syphacia
(sum(Genus.samples$Syphacia>0)/sum(Genus.samples$Syphacia>=0))*100

##Aspiculuris
(sum(Genus.samples$Aspiculuris>0)/sum(Genus.samples$Aspiculuris>=0))*100

###Eimeria
(sum(Genus.samples$Eimeria>0)/sum(Genus.samples$Eimeria>=0))*100

###Bacteroides
(sum(Genus.samples$Bacteroides>0)/sum(Genus.samples$Bacteroides>=0))*100

###Helicobacter
(sum(Genus.samples$Helicobacter>0)/sum(Genus.samples$Helicobacter>=0))*100

###Blautia
(sum(Genus.samples$Blautia>0)/sum(Genus.samples$Blautia>=0))*100

##Cryptosporidium
(sum(Genus.samples$Cryptosporidium>0)/sum(Genus.samples$Cryptosporidium>=0))*100

##Tritrichomonas
(sum(Genus.samples$Tritrichomonas>0)/sum(Genus.samples$Tritrichomonas>=0))*100

##Hymenolepis
(sum(Genus.samples$Hymenolepis>0)/sum(Genus.samples$Hymenolepis>=0))*100
