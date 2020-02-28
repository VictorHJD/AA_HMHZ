# load packages
library("VennDiagram")
library(grid)
library(gridExtra)
library(ggmap)
library(ggrepel)
library(RCurl)
library(ggmap)
library(pegas)
library(msa)
library(RColorBrewer)
library(Biostrings)
library(IRanges)
library(XVector)
library("ggplot2")
library("BarcodingR")
library(dplyr)
library("ggpubr")
library(multcomp)
library("Publish") 
library("AICcmodavg")
library(Rmisc)
library("ggpubr")
library("vegan")
library("picante")
library(phangorn) 
library(plotly)


# load HMHZ functions
source("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/R/functions/HMHZ_Functions.R")
source("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Databasing/master/R/functions/makeMiceTable.R")


#Get all the data 
###Detection data from flotation and Ap5
dx.data <- read.csv("~/Dokumente/Git_projects/Mouse_Eimeria_Databasing/data/Eimeria_detection/Inventory_contents_all.csv")
dx.data <- dx.data[-c(5:11,14:23)]
colnames(dx.data) <- c("Year", "Transect", "Code", "Mouse_ID", "Flot", "Ap5")
dx.data$Mouse_ID <- gsub(pattern = " ", replacement = "", x = dx.data$Mouse_ID)
dx.data <- dx.data[c(4:6)]
#dx.data <- na.omit(dx.data[c(1,2,3)])


raw.data1417 <- read.csv("~/Dokumente/Git_projects/Mouse_Eimeria_Databasing/data/MiceTable_fullEimeriaInfos_2014to2017.csv")
#colnames(raw.data1417)
raw.data2018 <- read.csv("~/Dokumente/Git_projects/Mouse_Eimeria_Databasing/data/Field_data/HZ18_September_Mice_Dissection.csv")

###Coordinates data
loc.data1417 <- raw.data1417[names(raw.data1417) %in% c("Mouse_ID", "Longitude", "Latitude")]
loc.data1417$Mouse_ID <- gsub(pattern = " ", replacement = "", x = loc.data$Mouse_ID)
#loc.data <- na.omit(loc.data[c(1,2,3)])
loc.data18 <- raw.data2018[names(raw.data2018) %in% c("Mouse_ID", "Longitude", "Latitude")]
loc.data <- na.omit(rbind(loc.data1417, loc.data18))


###Worms data
#colnames(raw.data1417)
worms.data1417 <- raw.data1417[names(raw.data1417) %in% c("Mouse_ID", "Hymenolepis_diminiuta", "Taenia_martis", "Heterakis_spumosa",
                                                          "Hymenolepis_microstoma", "Catenotaenia_pusilla", "Aspiculuris_tetraptera", "Syphacia_obvelata", "Trichuris_muris", 
                                                          "Mastophorus")]

worms.data1417 <- data.frame(Mouse_ID = worms.data1417$Mouse_ID, Aspiculuris_tetraptera = worms.data1417$Aspiculuris_tetraptera, 
                             Syphacia_obvelata = worms.data1417$Syphacia_obvelata, Heterakis_spumosa= worms.data1417$Heterakis_spumosa, 
                             Taenia_martis= worms.data1417$Taenia_martis, Catenotaenia_pusilla= worms.data1417$Catenotaenia_pusilla,
                             Hymenolepis_diminuta= worms.data1417$Hymenolepis_diminiuta, Hymenolepis_microstoma= worms.data1417$Hymenolepis_microstoma,
                             Mastaphorus= worms.data1417$Mastophorus, Trichuris_muris= worms.data1417$Trichuris_muris) ###change the order of the columns


#colnames(raw.data2018)

worms.data2018 <- raw.data2018[names(raw.data2018) %in% c("Mouse_ID", "ASP", "SYP", "HET", "MART", "CP", "HD", "HM", "MM", "TM")]
colnames(worms.data2018) <- c("Mouse_ID", "Aspiculuris_tetraptera", "Syphacia_obvelata", "Heterakis_spumosa", "Taenia_martis", "Catenotaenia_pusilla", 
                              "Hymenolepis_diminuta", "Hymenolepis_microstoma", "Mastophorus",  "Trichuris_muris")


colnames(worms.data1417) <- colnames(worms.data2018)
identical(colnames(worms.data1417), colnames(worms.data2018))
worms.data <- rbind(worms.data1417, worms.data2018)

###Start to merge 

all.data <- merge(dx.data, loc.data, all = T, by= "Mouse_ID")
all.data <- merge(all.data, worms.data, all = T, by= "Mouse_ID")


###Identification data 
### For COI (Cocci_primers + EimCOI_Primers)
phylogroup<- read.csv("~/Dokumente/Sequences/Manuscript/Species_assignment_14_18.csv")
phylogroup$Mouse_ID <- gsub(pattern = " ", replacement = "", x = phylogroup$Mouse_ID)
phylogroup$Species <- gsub(pattern = " ", replacement = "", x = phylogroup$Species)
phylogroup$E_ferrisi <- gsub(pattern = " ", replacement = "", x = phylogroup$E_ferrisi)
phylogroup$E_vermiformis <- gsub(pattern = " ", replacement = "", x = phylogroup$E_vermiformis)
table(phylogroup$Species) 

#phylogroup<- phylogroup[phylogroup$Species%in%c("E_falciformis", "E_ferrisi", "E_vermiformis"),]
#table(phylogroup$Species) ### Single infections


###All qPCR data 
allqpcr <-  read.csv("~/Dokumente/Git_projects/Mouse_Eimeria_Databasing/data/Eimeria_detection/FINALqpcrData_2016_2017_threshold5.csv")
allqpcr$Mouse_ID <- gsub(pattern = " ", replacement = "", x = allqpcr$Mouse_ID)


##DNA information

dnaConc <- read.csv("/home/victor/Dokumente/Metabarcoding/HMHZ_Metabarcoding/DNA_Concentration.csv")
dnaConc$Mouse_ID <- gsub(pattern = " ", replacement = "", x = dnaConc$Mouse_ID)
dnaConc$Transect <- gsub(pattern = " ", replacement = "", x = dnaConc$Transect)
dnaConc$Year <- gsub(pattern = " ", replacement = "", x = dnaConc$Year)
colnames(dnaConc) <- c("Year", "Transect","Mouse_ID", "ColContDNA", "Ap5", "Concentration", "p260_280", "p260_230")
dnaConc <- dnaConc[names(dnaConc) %in% c("Mouse_ID", "Concentration", "p260_280", "p260_230")]


####FINAL DATA TABLE###

###All detection information COI primer pair 1
all.data <- merge(all.data, allqpcr, all=T, by="Mouse_ID")

all.data <- merge(all.data, phylogroup, all=T, by="Mouse_ID")

all.data <- merge(all.data, dnaConc, by= "Mouse_ID")

all.data$year <- NULL

all.data$COI_Seq_1 <- ifelse(all.data$COI_Seq_1%in%c("positive","Tissue"), TRUE, FALSE)
all.data$COI_Seq_2 <- ifelse(all.data$COI_Seq_2%in%c("positive","Tissue"), TRUE, FALSE)
all.data$n18S_Seq <- ifelse(all.data$n18S_Seq%in%c("positive","Tissue"), TRUE, FALSE)
all.data$Ap5 <- ifelse(all.data$Ap5%in%"positive", TRUE, FALSE)
all.data$Flot <- ifelse(all.data$Flot%in%"positive", TRUE, FALSE)
all.data$qPCRstatus <- ifelse(all.data$qPCRstatus%in%"positive", TRUE, FALSE)
all.data$Group_18S <- ifelse(all.data$Group_18S%in%""|is.na(all.data$Group_18S), "negative", as.character(all.data$Group_18S))
all.data$Group_COI_1 <- ifelse(all.data$Group_COI_1%in%""|is.na(all.data$Group_COI_1), "negative", as.character(all.data$Group_COI_1))
all.data$Group_COI_2<- ifelse(all.data$Group_COI_2%in%""|is.na(all.data$Group_COI_2), "negative", as.character(all.data$Group_COI_2))


###Selection of samples

MetabarDNA <- subset(all.data, all.data$Concentration >= 30)
colnames(MetabarDNA)

#write.csv(MetabarDNA, "~/Dokumente/Metabarcoding/HMHZ_Metabarcoding/Sample_selection_Metabarcoding.csv")

## Plot samples with HI 

a1 <- ggplot(MetabarDNA, aes(x= HI, y= Concentration,  color= as.character(Year))) +
  geom_jitter(shape=16, position=position_jitter(0.0), alpha=0.8) + geom_hline(yintercept=30, linetype="dashed", color = "black") + 
  #geom_vline(xintercept = -6, linetype="dashed", color = "blue") + 
  labs(x="HI (Hybrid Index)", y = "DNA concentration ng/µL") +
  theme_classic()

b1 <- ggplot(MetabarDNA, aes(x= p260_280, y= Concentration,  color= as.character(Year))) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5) + geom_hline(yintercept=30, linetype="dashed", color = "black") + 
  #geom_vline(xintercept = -6, linetype="dashed", color = "blue") + 
  labs(x="Purity 260/280", y = "DNA concentration ng/µL") +
  theme_classic()

c1 <- ggplot(MetabarDNA, aes(x= p260_230, y= Concentration,  color= as.character(Year))) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5) + geom_hline(yintercept=30, linetype="dashed", color = "black") + 
  #geom_vline(xintercept = -6, linetype="dashed", color = "blue") + 
  labs(x="Purity 260/230", y = "DNA concentration ng/µL") +
  theme_classic()

d1 <- ggplot(MetabarDNA, aes(x= p260_280, y= p260_230,  color= as.character(Year))) +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha=0.5) + geom_hline(yintercept=1.0, linetype="dashed", color = "black") + 
  geom_vline(xintercept = 1.4, linetype="dashed", color = "black") + 
  labs(x="Purity 260/280", y = "Purity 260/230") +
  theme_classic()

grid.arrange(a1, b1, c1, d1, nrow= 2)

ggsave("DNA_characteristics_selec_samples.pdf", plot = grid.arrange(a1, b1, c1, d1, nrow= 2), device = "pdf", path = "~/Dokumente/Metabarcoding/HMHZ_Metabarcoding/",
       scale = 1, width = 25, height = 25, units ="cm",
       dpi = 450)

plot_ly(data= MetabarDNA, x = MetabarDNA$Concentration, y = MetabarDNA$p260_280, z = MetabarDNA$p260_230,
        marker = list(color = as.character (MetabarDNA$Year), showscale = TRUE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'DNA concentration ng/µl'),
                      yaxis = list(title = 'Purity 260/280'),
                      zaxis = list(title = 'Purity 260/230')),
         annotations = list(
           x = 1.13,
           y = 1.05,
           text = 'Hybrid Index',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))

### DNA worms 

at <- ggplot(MetabarDNA, aes(x= Concentration, y= Aspiculuris_tetraptera,  color= as.character(Year))) +
  geom_jitter(shape=16, position=position_jitter(0.0), alpha=0.8) + #geom_hline(yintercept=30, linetype="dashed", color = "black") + 
  #geom_vline(xintercept = -6, linetype="dashed", color = "blue") + 
  labs(x="DNA concentration ng/µL", y = "Aspiculuris count") +
  #geom_smooth()+
  theme_classic()

so <- ggplot(MetabarDNA, aes(x= Concentration, y= Syphacia_obvelata,  color= as.character(Year))) +
  geom_jitter(shape=16, position=position_jitter(0.0), alpha=0.8) + #geom_hline(yintercept=30, linetype="dashed", color = "black") + 
  #geom_vline(xintercept = -6, linetype="dashed", color = "blue") + 
  labs(x="DNA concentration ng/µL", y = "Syphacia count") +
  #geom_smooth()+
  theme_classic()

grid.arrange(at, so)

ggsave("DNA_select_worms.pdf", plot = grid.arrange(at, so), device = "pdf", path = "~/Dokumente/Metabarcoding/HMHZ_Metabarcoding/",
       scale = 1, width = 25, height = 25, units ="cm",
       dpi = 450)



####Map
areamus <- get_map(location =
                     c(min(MetabarDNA$Longitude - 0.5),
                       min(MetabarDNA$Latitude - 0.5),
                       max(MetabarDNA$Longitude + 0.5),
                       max(MetabarDNA$Latitude + 0.5)),
                   source = "stamen", maptype= "toner-lite", zoom = 8)



ggmap(areamus)+
  geom_point(data = MetabarDNA, shape = 21, size = 5, aes(Longitude, Latitude, fill= as.factor(Year)), alpha = 0.5)


ggsave("Localization_samples.pdf", plot = last_plot(), device = "pdf", path = "~/Dokumente/Metabarcoding/HMHZ_Metabarcoding/",
       scale = 1, width = 10, height = 15, units ="cm",
       dpi = 450)

####Randomization

Chip <- as.data.frame(split(MetabarDNA$Mouse_ID, sample(15))) 
colnames(Chip)<- c("Chip_1", "Chip_2", "Chip_3", "Chip_4", "Chip_5", "Chip_6", "Chip_7", "Chip_8", "Chip_9", "Chip_10", "Chip_11", "Chip_12", "Chip_13", "Chip_14", "Chip_15")
Chip[46:48,]<- c("N1_1", "N1_2", "N1_3")

for (i in 1:ncol(Chip)) {
  Chip[,i]<- sample(Chip[,i])
}

#write.csv(Chip, "~/Dokumente/Metabarcoding/HMHZ_Metabarcoding/Chips.csv")

#Eimeria_OPG <- read.csv("~/AA_HMHZ/Oocyst_counts.csv")
#Crypto_Oocyst <- read.csv("~/AA_HMHZ/Rawdata_Crypto_qPCR_results.csv")
#Crypto_Oocyst <- dplyr::select(Crypto_Oocyst, 1,9)

#all.data <- plyr::join(sample.data, Eimeria_OPG, by= "Mouse_ID")
#all.data <- plyr::join(sample.data, Crypto_Oocyst, by= "Mouse_ID")
#unique(all.data$Mouse_ID)

#Eimeria_OPG <- dplyr::select(all.data, 1,54) ### Select just the info for those samples used in the metabarcoding analysis
#Crypto_Oocyst<- dplyr::select(all.data, 1,54) ### Select just the info for those samples used in the metabarcoding analysis

#Extra.data <- plyr::join(Crypto_Oocyst, Eimeria_OPG, by= "Mouse_ID")
#sample.data <- plyr::join(sample.data, Extra.data, by= "Mouse_ID")
#write.csv(sample.data, "~/AA_HMHZ/Sample_selection_Metabarcoding_Complete.csv")

#rm(all.data)



### Sample selection for E. ferrisi genotyping 

all.data <- merge(dnaConc, phylogroup, by="Mouse_ID")

nrow(subset(all.data, all.data$Concentration >= 30 & all.data$Species_tissue == "E_ferrisi")) 
