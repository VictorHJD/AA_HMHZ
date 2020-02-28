###Phyloseq pipeline for diversity analysis

library(ggplot2)
library(reshape)
library(phyloseq)
library(data.table)
library(parallel)
library(microbiome)
library("pheatmap")
library(dplyr)
library(gridExtra)
library(grid)
library(vegan)
library("TmCalculator")

##Functions
sumSeqByTax <- function(PS.l, rank) {
  lapply(PS.l, function(x){ 
    counts <- data.frame(cbind(asvCount=colSums(otu_table(x)), tax_table(x)))
    counts$asvCount <- as.numeric(as.character(counts$asvCount))
    tapply(counts$asvCount, counts[, rank], sum)})
}

readDis <- function(df){
  ggplot(df, aes(TotalReads)) +
    geom_histogram() + 
    facet_wrap(~Chip_number)#+
  #labs(tag = x)
}

corPlot<- function(df){
  ggplot(df, aes(ReadsPin+1, Pinworms+1))+
    geom_jitter(height = 0.5)+
    geom_smooth(method = lm)+
    scale_x_log10("Sequence reads (Syphacia/Aspiculuris)")+
    scale_y_log10("Counts of Syphacia/Aspiculuris")+
    labs(tag= x)+
    theme_bw()
}

##Load data
PS.l<- readRDS(file="/SAN/Victors_playground/Metabarcoding/AA_HMHZ/PhyloSeqList_HMHZ_All.Rds")

##curatedata 

PS.l <- lapply(PS.l, function(x){ 
  prune_samples(sample_sums(x)>0, x)
})

####Primer information
primerInput <- read.csv("/localstorage/victor/AA_Primer_evaluation/primerInputUnique.csv") 
primerInput$Primer_name <- gsub(pattern = " ", replacement = "", x = primerInput$Primer_name)
primerInput$X <- NULL
primerInput$Primer_name <- gsub(pattern = " ", replacement = "", x = primerInput$Primer_name)
##Primer characteristics (Tm and GC)
primerInput$Seq_F <- as.vector(primerInput$Seq_F)
primerInput$Seq_R <- as.vector(primerInput$Seq_R)

primerInput <- primerInput[-c(20,21,82),] ###Eliminating primers with inosine ... To make the loop work :S 

output <- data.frame(matrix(nrow = nrow(primerInput), ncol = 5))
colnames(output) <- c("Primer_name","Tm_F", "Tm_R", "GC_F", "GC_R")
output[,1]<-primerInput$Primer_name

for (row in 1: nrow(primerInput))
{
  Fwd <- primerInput[row, "Seq_F"]
  Rev <- primerInput[row, "Seq_R"]
  
  TmF <- Tm_NN(Fwd, ambiguous = T, imm_table = "DNA_IMM1")
  TmR <- Tm_NN(Rev, ambiguous = T, imm_table = "DNA_IMM1")
  GCF <- GC(Fwd, ambiguous = T)
  GCR <- GC(Rev, ambiguous = T)
  
  output[row, 2] <- TmF
  output[row, 3] <- TmR
  output[row, 4] <- GCF
  output[row, 5] <- GCR
  
}

primerInput <- read.csv("~/AA_Primer_evaluation/primerInputUnique.csv")
primerInput$Primer_name <- gsub(pattern = " ", replacement = "", x = primerInput$Primer_name)
primerInput <- merge(primerInput, output, by= "Primer_name", all= T)
primerInput$X <- NULL
primerInput$Primer_name <- gsub(pattern = " ", replacement = "", x = primerInput$Primer_name)

rm(GCF)
rm(GCR)
rm(Fwd)
rm(Rev)
rm(row)
rm(TmF)
rm(TmR)
rm(output)

###Raw counts by amplicon
rawcounts <- as.data.frame(unlist(lapply(PS.l, function(x){
  data.frame(cbind(asvCount=colSums(data.frame(rowSums(otu_table(x))))))
})))

rawcounts[,2] <- rownames(rawcounts)
rownames(rawcounts) <- c(1:nrow(rawcounts))
colnames(rawcounts) <- c("Reads","Primer_name")
rawcounts <- data.frame(Primer_name = rawcounts$Primer_name, Reads = rawcounts$Reads) 
rawcounts$Primer_name <- gsub(".asvCount", "\\1", rawcounts$Primer_name)


###Reads by taxa
readNumByPhylum <- sumSeqByTax(PS.l = PS.l, rank = "phylum")
readNumByGenus <- sumSeqByTax(PS.l = PS.l, rank = "genus")
readNumByFamily <- sumSeqByTax(PS.l = PS.l, rank = "family")
readNumByOrder <- sumSeqByTax(PS.l = PS.l, rank = "order")
readNumBySpecies <- sumSeqByTax(PS.l = PS.l, rank = "species")

###
AbPhy <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByPhylum)) ### Start a loop: fro every element in the list ...
{ 
  phyla <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByPhylum[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    phyla[1,1] <- 0    ### Add a zero in the first column
    phyla[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    phyla <- as.data.frame((readNumByPhylum[[i]]))  ###Make a data frame with the data included in each element of the list 
    phyla[,2] <- rownames(phyla) ### And use the rownames as information of the second column 
  }
  
  phyla[,3] <- names(readNumByPhylum)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(phyla) <- c("ASV", "Phyla", "Primer_name") ### change the names for the columns 
  AbPhy <- rbind(AbPhy, phyla) ### Join all the "individual" data frames into the final data frame 
  
}   ### close loop

rownames(AbPhy) <- c(1:nrow(AbPhy)) ### change the rownames to consecutive numbers 
AbPhy <- data.frame(Primer_name = AbPhy$Primer_name, Phyla = AbPhy$Phyla, ASV = AbPhy$ASV) ###change the order of the columns
AbPhy$ASV <- as.numeric(AbPhy$ASV)

AbPhy %>%
  group_by(Primer_name) %>% 
  mutate(Total_ASV = sum(ASV)) -> AbPhy ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

Relative_abundance = AbPhy$ASV/AbPhy$Total_ASV ### create a vector with the result of the operation 

AbPhy[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbPhy)[5] <- "Relative_abundance" ### Change the name of the column 

AbPhy$Primer_name <- gsub(pattern = " ", replacement = "", x = AbPhy$Primer_name) ### check if the primer names have extra spaces

AbPhy$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbPhy$Primer_name)

AbPhy <- merge(AbPhy, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

AbPhy <- plyr::join(AbPhy, rawcounts, by= "Primer_name")


###What about genus level?
AbGen <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByGenus)) ### Start a loop: for every element in the list ...
{ 
  genus <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByGenus[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    genus[1,1] <- 0    ### Add a zero in the first column
    genus[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    genus <- as.data.frame((readNumByGenus[[i]]))  ###Make a data frame with the data included in each element of the list 
    genus[,2] <- rownames(genus) ### And use the rownames as information of the second column 
  }
  
  genus[,3] <- names(readNumByGenus)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(genus) <- c("ASV", "Genus", "Primer_name") ### change the names for the columns 
  AbGen <- rbind(AbGen, genus) ### Join all the "individual" data frames into the final data frame 
  
}  ### close loop

rownames(AbGen) <- c(1:nrow(AbGen)) ### change the rownames to consecutive numbers 
AbGen <- data.frame(Primer_name = AbGen$Primer_name, Genus = AbGen$Genus, ASV = AbGen$ASV) ###change the order of the columns

AbGen %>%
  group_by(Primer_name) %>% 
  mutate(Total_ASV = sum(ASV)) -> AbGen ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

Relative_abundance = AbGen$ASV/AbGen$Total_ASV ### create a vector with the result of the operation 

AbGen[,5] <- Relative_abundance ### And put it in the same order in the colum 5

colnames(AbGen)[5] <- "Relative_abundance" ### Change the name of the column 

AbGen$Primer_name <- gsub(pattern = " ", replacement = "", x = AbGen$Primer_name) ### check if the primer names have extra spaces

AbGen$Primer_name <- gsub(pattern = "-", replacement = "_", x = AbGen$Primer_name)


AbGen <- merge(AbGen, primerInput, by= "Primer_name") ###merge the selected information with the origial data frame created 

AbGen <- plyr::join(AbGen, rawcounts, by= "Primer_name")

###Rarefaction curves by amplicon 
###Summary by amplicon
s.list <- lapply(PS.l, function (x) {
  summary(sample_sums(x))
})

###Distribution of reads by amplicon

hist.list <- lapply(PS.l, function (x) {
  histogram(rowSums(otu_table(x)))
})

pdf("~/AA_HMHZ/Primers/Distribution_by_amplicon.pdf", 
    width=8, height=10, onefile=T)
hist.list 
dev.off()

###Rarefaction curves 
rarecurv.l <- lapply(PS.l, function (x) {
              vegan::rarecurve(otu_table(x),
              label = F)
})

#pdf("~/AA_HMHZ/Primers/RareCurv_by_amplicon.pdf", 
#    width=8, height=10, onefile=T)
#rarecurv.l
#dev.off()

###Rarefied data
PS.l.rar <- lapply(PS.l, function (x) {
  rarefy_even_depth(x, rngseed=1, sample.size=0.3*mean(sample_sums(x)), replace=F)
})

##Alpha diversity 

alphaDiv <- data.frame()

for (i in 1:length(PS.l.rar)){

  a <- data.frame()

  b <- data.frame(estimate_richness(PS.l.rar[[i]], measures = c("Observed","Chao1", "Shannon")))

  a<-colMeans(b)

  alphaDiv<- rbind(alphaDiv, a)
}

alphaDiv[,5] <- names(PS.l.rar)

colnames(alphaDiv) <- c("Observed","Chao1","se_Chao1", "Shannon", "Primer_name")


##Beta diversity
ord.l <-lapply(PS.l.rar, function (x) {
  ordinate(x, method="PCoA", distance="bray")
})


pca.list <- for (i in 1:length(PS.l.rar)) {
  x<- names(PS.l.rar)
  pca <- plot_ordination(PS.l.rar[[i]], ord.l[[i]], color="HI")+
    scale_color_gradient(low="blue", high="red", name= "Mouse genotype\n(Hybird Index)")+
    theme_bw()+
    geom_point(size=5, alpha= 0.75)+
    labs(title = x[i])
}

#####Heatmap primer HMHZ
num.taxa <-sapply(c("species","genus", "family", "order", "phylum", "superkingdom"), function (rank){
  lapply(PS.l, function (x) 
    length(get_taxa_unique(x, taxonomic.rank = rank)))
})

num.reads <- unlist(lapply(PS.l, function (x) 
  sum(otu_table(x))))

PrimTax <- as.data.frame(cbind(num.taxa, num.reads))

PrimTax <- as.data.frame(apply(PrimTax, 2, unlist))
PrimTax[,8] <- row.names(PrimTax)  
colnames(PrimTax)[8] <- "Primer_name"
PrimTax<- plyr::join(PrimTax, primerInput, by= "Primer_name")

ps.numord.l <- PS.l[order(PrimTax$num.reads, decreasing=TRUE)]

###Select markers for Eukaryotes 
Primtax.comb <- subset(PrimTax, !Gen=="16S")
ps.numord.l.Euk <- PS.l[c(Primtax.comb$Primer_name)][order(Primtax.comb$num.reads, decreasing=TRUE)]

PM.Euk <- lapply(ps.numord.l.Euk, function (x) {
  subset_taxa(x, superkingdom%in%"Eukaryota")
})

PG.Euk <- lapply(PM.Euk, function (x) {
  tax_glom(x, "genus", NArm = TRUE)
})

PS.Euk <- lapply(PG.Euk, function (y) {
  make.names(tax_table(y)[, 5])
}) 

P.intersection <- sapply(seq_len(length(PS.Euk)), function(x)
  sapply(seq_len(length(PS.Euk)), function(y) length(intersect(unlist(PS.Euk[x]), unlist(PS.Euk[y]))))
)

rownames(P.intersection)<- names(PS.Euk)
colnames(P.intersection)<- names(PS.Euk)

P.clust <- hclust(dist(P.intersection), method = "complete") ##Dendogram
require(dendextend)
as.dendrogram(P.clust) %>%
  plot(horiz = TRUE)

P.col <- cutree(tree = P.clust, k = 2)
P.col  <- data.frame(cluster = ifelse(test = P.col  == 1, yes = "cluster 1", no = "cluster 2"))
P.col$Primer_name <- rownames(P.col)

P.col <- merge(P.col, PrimTax, by="Primer_name", sort= F)

col_groups <- P.col %>%
  select("Primer_name", "Gen") ##Here It is possible to add the "Region"

row.names(col_groups)<- col_groups$Primer_name

col_groups$Primer_name<- NULL

colour_groups <- list( Gen= c("18S"= "#440154FF", 
                              "28S"= "#21908CFF", 
                              "COI"= "#FDE725FF", 
                              #"16S"= "pink",
                              #"12S"= "#E3DAC9",
                              "ITS"= "#C46210",
                              #"rbcL" = "#D0FF14",
                              "GDH"= "#0BF2DD",
                              "HSP70"= "#F20BD3",
                              "tRNA"= "#2AF20B",
                              "BG"= "#613A08"))
require("viridis")
paraheatmap <- pheatmap(P.intersection, 
                        color = plasma(100),
                        border_color = NA,
                        annotation_col = col_groups, 
                        #annotation_row = col_groups,
                        annotation_colors = colour_groups,
                        #cutree_rows = 2,
                        #cutree_cols = 2,
                        show_rownames = F,
                        show_colnames = F,
                        main= "Redundant genus of Eukaryotes amplified")


TopPhylum <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByPhylum)) ### Start a loop: for every element in the list ...
{ 
  top <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByPhylum[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    top[1,1] <- 0    ### Add a zero in the first column
    top[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    top <- as.data.frame(readNumByPhylum[[i]])  ###Make a data frame with the data included in each element of the list 
    colnames(top)<- "j"
    top <-subset(top, top$j == max(top))
    top[,2] <- rownames(top) ### And use the rownames as information of the second column 
  }
  
  top[,3] <- names(readNumByPhylum)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(top) <- c("ASV_Phylum", "Phylum", "Primer_name") ### change the names for the columns 
  TopPhylum <- rbind(TopPhylum, top) ### Join all the "individual" data frames into the final data frame 
  
}


TopOrder <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByOrder)) ### Start a loop: for every element in the list ...
{ 
  top <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByPhylum[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    top[1,1] <- 0    ### Add a zero in the first column
    top[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    top <- as.data.frame(readNumByOrder[[i]])  ###Make a data frame with the data included in each element of the list 
    colnames(top)<- "j"
    top <-subset(top, top$j == max(top))
    top[,2] <- rownames(top) ### And use the rownames as information of the second column 
  }
  
  top[,3] <- names(readNumByOrder)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(top) <- c("ASV_Order", "Order", "Primer_name") ### change the names for the columns 
  TopOrder <- rbind(TopOrder, top) ### Join all the "individual" data frames into the final data frame 
  
}

TopFam <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByFamily)) ### Start a loop: for every element in the list ...
{ 
  top <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByFamily[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    top[1,1] <- 0    ### Add a zero in the first column
    top[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    top <- as.data.frame(readNumByFamily[[i]])  ###Make a data frame with the data included in each element of the list 
    colnames(top)<- "j"
    top <-subset(top, top$j == max(top))
    top[,2] <- rownames(top) ### And use the rownames as information of the second column 
  }
  
  top[,3] <- names(readNumByFamily)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(top) <- c("ASV_Family", "Family", "Primer_name") ### change the names for the columns 
  TopFam <- rbind(TopFam, top) ### Join all the "individual" data frames into the final data frame 
  
}

TopGen <- data.frame() ###Create the data frame 

for (i in 1: length(readNumByGenus)) ### Start a loop: for every element in the list ...
{ 
  top <- data.frame() #### make an individual data frame ...
  
  if(nrow(as.data.frame(readNumByGenus[[i]])) == 0) ### A tiny condition if the longitude of the list is = to 0 
  {
    top[1,1] <- 0    ### Add a zero in the first column
    top[1,2] <- "NA" ### And a NA in the second column.
  }else               ### For the rest of the elements: 
  {
    top <- as.data.frame(readNumByGenus[[i]])  ###Make a data frame with the data included in each element of the list 
    colnames(top)<- "j"
    top <-subset(top, top$j == max(top))
    top[,2] <- rownames(top) ### And use the rownames as information of the second column 
  }
  
  top[,3] <- names(readNumByGenus)[i] ### Take the names of every list and use them to fill column 3 as many times the logitude of the column 2
  colnames(top) <- c("ASV_Genus", "Genus", "Primer_name") ### change the names for the columns 
  TopGen <- rbind(TopGen, top) ### Join all the "individual" data frames into the final data frame 
  
}

require(plyr)
TopByRank <- join(TopPhylum, TopOrder, by= "Primer_name")
TopByRank <- join(TopFam, TopByRank, by= "Primer_name")
TopByRank <- join(TopByRank, TopGen, by= "Primer_name")

TopByRank<- TopByRank[c("Primer_name", "Phylum", "ASV_Phylum", "Order", "ASV_Order", "Family", "ASV_Family", "Genus", "ASV_Genus")]


###Composition plots by primer pair 
ggplot(data=AbGen, aes(x= Primer_name, y= Relative_abundance, fill= Genus)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
  #facet_wrap(~ AbyPp$Target_gene) +
  theme(legend.position="none") + guides(fill=guide_legend(nrow=50)) 

###Correlation pinworms

##sequence per amplicon
spa<- lapply(PS.l, function (x) {
  data.table(as(sample_data(x), "data.frame"),
  TotalReads= sample_sums(x), keep.rownames = T)
}) 

##Read destribution by amplicon by chip
pdf("~/AA_HMHZ/Primers/Reads_per_chip.pdf")
lapply(spa, readDis)
dev.off()

###Reads for Syp and Asp by amplicon
syphacia<- subset(AbGen, Genus%in%"Syphacia")
aspiculuris<- subset(AbGen, Genus%in%"Aspiculuris")

syphacia<-syphacia[,1:3]
aspiculuris<-aspiculuris[,1:3]
require(plyr)
pinworms<- merge(syphacia, aspiculuris, by= "Primer_name")
colnames(pinworms)<- c("Primer_name", "Syphacia", "Syp_reads", "Aspiculuris", "Asp_reads")

pinworms<- pinworms%>%select(1,3,5)

pinworms %>%
  mutate(Total_pin = Syp_reads+Asp_reads) -> pinworms ### Add a new variable that will contain the sum of all the sequencing reads by primer pair

###subset PS amplicons that amplified for Syp and/or Asp 
along<- pinworms$Primer_name

PS.SyAs <- PS.l[c(along)]

spa2<- lapply(PS.SyAs, function (x) {
  x<- data.table(as(sample_data(x), "data.frame"),
             TotalReads= sample_sums(x), keep.rownames = T)
  
  setnames(x, "rn", "Sample_ID")
}) 

###Get read number for Syp and Asp per sample 
PS.l.pinworms <- lapply(PS.SyAs, function (x) {
  subset_taxa(x, genus%in%c("Syphacia", "Aspiculuris"))
})

spa3 <- lapply(PS.l.pinworms, function (x) {
  x<- data.table(as(sample_data(x), "data.frame"),
                 ReadsPin= sample_sums(x), keep.rownames = T)
}) 


pdf("~/AA_HMHZ/Primers/Pinworms_correlations.pdf")
lapply(spa3, corPlot)
dev.off()

names(spa3)
