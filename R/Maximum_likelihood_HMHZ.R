###Packages

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
library("fitdistrplus")
library("optimx")
library(FSA)

devtools::install_github("alicebalard/parasiteLoad@v2.0")
library(parasiteLoad)

# to remove CI in case bad optimisation
# source(whatever function script)

### Functions
source("GitProjects/AA_HMHZ/R/TestDistributions.R")
source("GitProjects/AA_HMHZ/R/bananaplotNoCI.R")

##Load data 
if(!exists("Phylum.samples")){
  source("PS_Samples_HMHZ.R")
}

if(!exists("Genus.samples")){
  source("PS_Samples_HMHZ.R")
}

if(!exists("sdt")){
  source("Diversity.R")
}
####Transform data to fit parasiteLoad
##Relative abundances Phylum level
Phylum.samples$Ascomycota.trans<- round(log10(Phylum.samples$Ascomycota+1)*10000)
Phylum.samples$Basidiomycota.trans<- round(log10(Phylum.samples$Basidiomycota+1)*10000)
Phylum.samples$Nematoda.trans<- round(log10(Phylum.samples$Nematoda+1)*10000)
Phylum.samples$Apicomplexa.trans<- round(log10(Phylum.samples$Apicomplexa+1)*10000)
Phylum.samples$Platyhelminthes.trans<- round(log10(Phylum.samples$Platyhelminthes+1)*10000)
Phylum.samples$Ciliophora.trans<- round(log10(Phylum.samples$Ciliophora+1)*10000)
Phylum.samples$Firmicutes.trans<- round(log10(Phylum.samples$Firmicutes+1)*10000)
Phylum.samples$Bacteroidetes.trans<- round(log10(Phylum.samples$Bacteroidetes+1)*10000)
Phylum.samples$Proteobacteria.trans<- round(log10(Phylum.samples$Proteobacteria+1)*10000)
Phylum.samples$Fusobacteria.trans<- round(log10(Phylum.samples$Fusobacteria+1)*10000)


##Relative abundances Genus level
Genus.samples$Eimeria.trans<- round(log10(Genus.samples$Eimeria+1)*10000)
Genus.samples$Kazachstania.trans<- round(log10(Genus.samples$Kazachstania+1)*10000)
Genus.samples$Aspiculuris.trans<- round(log10(Genus.samples$Aspiculuris+1)*10000)
Genus.samples$Syphacia.trans<- round(log10(Genus.samples$Syphacia+1)*10000)
Genus.samples$Tritrichomonas.trans<- round(log10(Genus.samples$Tritrichomonas+1)*10000)
Genus.samples$Cryptosporidium.trans<- round(log10(Genus.samples$Cryptosporidium+1)*10000)
Genus.samples$Hymenolepis.trans<- round(log10(Genus.samples$Hymenolepis+1)*10000)
Genus.samples$Bacteroides.trans<- round(log10(Genus.samples$Bacteroides+1)*10000)

Genus.samples$Helicobacter.trans<- round(log10(Genus.samples$Helicobacter+1)*10000)
Genus.samples$Blautia.trans<- round(log10(Genus.samples$Blautia+1)*10000)
Genus.samples$Lactobacillus.trans<- round(log10(Genus.samples$Lactobacillus+1)*10000)
Genus.samples$Clostridium.trans<- round(log10(Genus.samples$Clostridium+1)*10000)
Genus.samples$Fusimonas.trans<- round(log10(Genus.samples$Fusimonas+1)*10000)
Genus.samples$Muribaculum.trans<- round(log10(Genus.samples$Muribaculum+1)*10000)


##Diversity index
##Small adjustment for a sample with NA
sdt$pielou[is.na(sdt$pielou)]<- as.numeric(1)

####Choose distribution
##Total reads
x<- Phylum.samples$TotalReads

##Alpha diversity measurments 
x<- sdt$Chao1 ##Normal
x<- sdt$Shannon ##Normal
x<- sdt$pielou

##Phylum level
x<- Phylum.samples$Ascomycota.trans
x<- Phylum.samples$Basidiomycota.trans
x<- Phylum.samples$Nematoda.trans
x<- Phylum.samples$Apicomplexa.trans
x<- Phylum.samples$Platyhelminthes.trans
x<- Phylum.samples$Ciliophora.trans
x<- Phylum.samples$Firmicutes.trans
x<- Phylum.samples$Bacteroidetes.trans
x<- Phylum.samples$Proteobacteria.trans


##Genus level
x<- Genus.samples$Kazachstania.trans
x<- Genus.samples$Aspiculuris.trans
x<- Genus.samples$Syphacia.trans
x<- Genus.samples$Eimeria.trans
x<- Genus.samples$Tritrichomonas.trans
x<- Genus.samples$Cryptosporidium.trans
x<- Genus.samples$Hymenolepis.trans
x<- Genus.samples$Bacteroides.trans
x<- Genus.samples$Helicobacter.trans
x<- Genus.samples$Blautia.trans

x <- x[x>0]
hist(x, breaks = 100)
descdist(x)
findGoodDist(x, distribs = c("normal", "negative binomial", "poisson"), 
             distribs2 = c("norm", "nbinom", "pois"))

###ML anaylsis
##Total reads
fitReads <- parasiteLoad::analyse(data = Phylum.samples, response = "TotalReads", model = "weibull", group = "Seq_Run")

Reads.plot <- parasiteLoad::bananaPlot(mod = fitReads$H1,
                                     data = sdt,
                                     response = "TotalReads",
                                     islog10 = T, group = "Seq_Run",
                                     cols = c("#006A4E", "#006A4E"))+ labs(tag = "A)")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log10 Total read counts")+
  theme(legend.position ="none")

##Alpha diversity measurments
Chao1<- parasiteLoad::analyse(data = sdt, response = "Chao1", model = "normal", group = "Seq_Run")
Chao1.plot <- parasiteLoad::bananaPlot(mod = Chao1$H1,
                                            data = sdt,
                                            response = "Chao1",
                                            islog10 = F, group = "Seq_Run",
                                            cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "A)")+
  #ggtitle("Ascomycota")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Richness (Chao1 Index)")+
  theme(legend.position ="none")

Shannon<- parasiteLoad::analyse(data = sdt, response = "Shannon", model = "normal", group = "Seq_Run")
Shannon.plot <- parasiteLoad::bananaPlot(mod = Shannon$H1,
                                       data = sdt,
                                       response = "Shannon",
                                       islog10 = F, group = "Seq_Run",
                                       cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "B)")+
  #ggtitle("Ascomycota")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Diversity (Shannon Index)")+
  theme(legend.position ="none")

Pielou<- parasiteLoad::analyse(data = sdt, response = "pielou", model = "normal", group = "Seq_Run")
Pielou.plot <- parasiteLoad::bananaPlot(mod = Pielou$H1,
                                         data = sdt,
                                         response = "pielou",
                                         islog10 = F, group = "Seq_Run",
                                         cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "C)")+
  #ggtitle("Ascomycota")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Evenness (Pielou's Index)")+
  theme(legend.position ="none")

##Just bacteria
Chao1.bacteria<- parasiteLoad::analyse(data = sdt.bacteria.R, response = "Chao1", model = "normal", group = "Seq_Run")
Chao1.bacteria.plot <- parasiteLoad::bananaPlot(mod = Chao1.bacteria$H1,
                                       data = sdt.bacteria.R,
                                       response = "Chao1",
                                       islog10 = F, group = "Seq_Run",
                                       cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "A)")+
  ggtitle("Bacteria")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Richness (Chao1 Index)")+
  theme(legend.position ="none")

Shannon.bacteria<- parasiteLoad::analyse(data = sdt.bacteria.R, response = "Shannon", model = "normal", group = "Seq_Run")
Shannon.bacteria.plot <- parasiteLoad::bananaPlot(mod = Shannon.bacteria$H1,
                                         data = sdt.bacteria.R,
                                         response = "Shannon",
                                         islog10 = F, group = "Seq_Run",
                                         cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "B)")+
  ggtitle("Bacteria")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Diversity (Shannon Index)")+
  theme(legend.position ="none")

Pielou.bacteria<- parasiteLoad::analyse(data = sdt.bacteria.R, response = "pielou", model = "normal", group = "Seq_Run")
Pielou.bacteria.plot <- parasiteLoad::bananaPlot(mod = Pielou.bacteria$H1,
                                        data = sdt.bacteria.R,
                                        response = "pielou",
                                        islog10 = F, group = "Seq_Run",
                                        cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "C)")+
  ggtitle("Bacteria")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Evenness (Pielou's Index)")+
  theme(legend.position ="none")
##Just eukaryota
Chao1.eukaryota<- parasiteLoad::analyse(data = sdt.eukaryota.R, response = "Chao1", model = "normal", group = "Seq_Run")
Chao1.eukaryota.plot <- parasiteLoad::bananaPlot(mod = Chao1.eukaryota$H1,
                                       data = sdt.eukaryota.R,
                                       response = "Chao1",
                                       islog10 = F, group = "Seq_Run",
                                       cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "A)")+
  ggtitle("Eukaryota")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Richness (Chao1 Index)")+
  theme(legend.position ="none")

Shannon.eukaryota<- parasiteLoad::analyse(data = sdt.eukaryota.R, response = "Shannon", model = "normal", group = "Seq_Run")
Shannon.eukaryota.plot <- parasiteLoad::bananaPlot(mod = Shannon.eukaryota$H1,
                                         data = sdt.eukaryota.R,
                                         response = "Shannon",
                                         islog10 = F, group = "Seq_Run",
                                         cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "B)")+
  ggtitle("Eukaryota")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Diversity (Shannon Index)")+
  theme(legend.position ="none")

Pielou.eukaryota<- parasiteLoad::analyse(data = sdt.eukaryota.R, response = "pielou", model = "normal", group = "Seq_Run")
Pielou.eukaryota.plot <- parasiteLoad::bananaPlot(mod = Pielou.eukaryota$H1,
                                        data = sdt.eukaryota.R,
                                        response = "pielou",
                                        islog10 = F, group = "Seq_Run",
                                        cols = c("darkgoldenrod1", "darkgoldenrod1"))+ labs(tag = "C)")+
  ggtitle("Eukaryota")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("Evenness (Pielou's Index)")+
  theme(legend.position ="none")

##Phylum level
Ascomycota<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Ascomycota.trans>0,], response = "Ascomycota.trans", model = "weibull", group = "Seq_Run")
Ascomycota.plot <- parasiteLoad::bananaPlot(mod = Ascomycota$H1,
                         data = Phylum.samples[Phylum.samples$Ascomycota.trans>0,],
                         response = "Ascomycota.trans",
                         islog10 = T, group = "Seq_Run",
                         cols = c("#006A4E", "#006A4E"))+ labs(tag = "A)")+
  ggtitle("Ascomycota")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Basidiomycota<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Basidiomycota.trans>0,], response = "Basidiomycota.trans", model = "negbin", group = "Seq_Run")
Basidiomycota.plot <- parasiteLoad::bananaPlot(mod = Basidiomycota$H1,
                                            data = Phylum.samples[Phylum.samples$Basidiomycota.trans>0,],
                                            response = "Basidiomycota.trans",
                                            islog10 = T, group = "Seq_Run",
                                            cols = c("#006A4E", "#006A4E"))+ labs(tag = "B)")+
  ggtitle("Basidiomycota")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Nematoda<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Nematoda.trans>0,], response = "Nematoda.trans", model = "negbin", group = "Seq_Run")
Nematoda.plot <- parasiteLoad::bananaPlot(mod = Nematoda$H1,
                                               data = Phylum.samples[Phylum.samples$Nematoda.trans>0,],
                                               response = "Nematoda.trans",
                                               islog10 = T, group = "Seq_Run",
                                               cols = c("#006A4E", "#006A4E"))+ labs(tag = "A)")+
  ggtitle("Nematoda")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Apicomplexa<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Apicomplexa.trans>0,], response = "Apicomplexa.trans", model = "negbin", group = "Seq_Run")
Apicomplexa.plot <- parasiteLoad::bananaPlot(mod = Apicomplexa$H1,
                                          data = Phylum.samples[Phylum.samples$Apicomplexa.trans>0,],
                                          response = "Apicomplexa.trans",
                                          islog10 = T, group = "Seq_Run",
                                          cols = c("#006A4E", "#006A4E"))+ labs(tag = "B)")+
  ggtitle("Apicomplexa")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Platyhelminthes<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Platyhelminthes.trans>0,], response = "Platyhelminthes.trans", model = "negbin", group = "Seq_Run")
Platyhelminthes.plot <-  bananaplotNoCI(mod = Platyhelminthes$H1,
                                             data = Phylum.samples[Phylum.samples$Platyhelminthes.trans>0,],
                                             response = "Platyhelminthes.trans",
                                             islog10 = T, group = "Seq_Run",
                                             cols = c("#006A4E", "#006A4E"))+ labs(tag = "C)")+
  ggtitle("Platyhelminthes")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Ciliophora<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Ciliophora.trans>0,], response = "Ciliophora.trans", model = "negbin", group = "Seq_Run")
Ciliophora.plot <-  bananaplotNoCI(mod = Ciliophora$H1,
                                                 data = Phylum.samples[Phylum.samples$Ciliophora.trans>0,],
                                                 response = "Ciliophora.trans",
                                                 islog10 = T, group = "Seq_Run",
                                                 cols = c("#006A4E", "#006A4E"))+ labs(tag = "D)")+
  ggtitle("Ciliophora")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")


Firmicutes<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Firmicutes.trans>0,], response = "Firmicutes.trans", model = "negbin", group = "Seq_Run")
Firmicutes.plot <- parasiteLoad::bananaPlot(mod = Firmicutes$H1,
                                            data = Phylum.samples[Phylum.samples$Firmicutes.trans>0,],
                                            response = "Firmicutes.trans",
                                            islog10 = T, group = "Seq_Run",
                                            cols = c("#006A4E", "#006A4E"))+ labs(tag = "A)")+
  ggtitle("Firmicutes")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Bacteroidetes<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Bacteroidetes.trans>0,], response = "Bacteroidetes.trans", model = "negbin", group = "Seq_Run")
Bacteroidetes.plot <- parasiteLoad::bananaPlot(mod = Bacteroidetes$H1,
                                            data = Phylum.samples[Phylum.samples$Bacteroidetes.trans>0,],
                                            response = "Bacteroidetes.trans",
                                            islog10 = T, group = "Seq_Run",
                                            cols = c("#006A4E", "#006A4E"))+ labs(tag = "B)")+
  ggtitle("Bacteroidetes")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Proteobacteria<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Proteobacteria.trans>0,], response = "Proteobacteria.trans", model = "negbin", group = "Seq_Run")
Proteobacteria.plot <- parasiteLoad::bananaPlot(mod = Proteobacteria$H1,
                                               data = Phylum.samples[Phylum.samples$Proteobacteria.trans>0,],
                                               response = "Proteobacteria.trans",
                                               islog10 = T, group = "Seq_Run",
                                               cols = c("#006A4E", "#006A4E"))+ labs(tag = "C)")+
  ggtitle("Proteobacteria")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Fusobacteria<- parasiteLoad::analyse(data = Phylum.samples[Phylum.samples$Fusobacteria.trans>0,], response = "Fusobacteria.trans", model = "negbin", group = "Seq_Run")
Fusobacteria.plot <- bananaplotNoCI(mod = Fusobacteria$H1,
                                                data = Phylum.samples[Phylum.samples$Fusobacteria.trans>0,],
                                                response = "Fusobacteria.trans",
                                                islog10 = T, group = "Seq_Run",
                                                cols = c("#006A4E", "#006A4E"))+ labs(tag = "D)")+
  ggtitle("Fusobacteria")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

##Genus
Kazachstania<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Kazachstania.trans>0,], response = "Kazachstania.trans", model = "negbin", group = "Seq_Run")
Kazachstania.plot <- parasiteLoad::bananaPlot(mod = Kazachstania$H1,
                                            data = Genus.samples[Genus.samples$Kazachstania.trans>0,],
                                            response = "Kazachstania.trans",
                                            islog10 = T, group = "Seq_Run",
                                            cols = c("#006A4E", "#006A4E"))+ labs(tag = "A)")+
  ggtitle("Kazachstania")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Aspiculuris<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Aspiculuris.trans>0,], response = "Aspiculuris.trans", model = "negbin", group = "Seq_Run")
Aspiculuris.plot <- parasiteLoad::bananaPlot(mod = Aspiculuris$H1,
                                              data = Genus.samples[Genus.samples$Aspiculuris.trans>0,],
                                              response = "Aspiculuris.trans",
                                              islog10 = T, group = "Seq_Run",
                                              cols = c("#006A4E", "#006A4E"))+ labs(tag = "A)")+
  ggtitle("Apiculuris")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Syphacia<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Syphacia.trans>0,], response = "Syphacia.trans", model = "negbin", group = "Seq_Run")
Syphacia.plot <- parasiteLoad::bananaPlot(mod = Syphacia$H1,
                                             data = Genus.samples[Genus.samples$Syphacia.trans>0,],
                                             response = "Syphacia.trans",
                                             islog10 = T, group = "Seq_Run",
                                             cols = c("#006A4E", "#006A4E"))+ labs(tag = "B)")+
  ggtitle("Syphacia")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Eimeria<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Eimeria.trans>0,], response = "Eimeria.trans", model = "negbin", group = "Seq_Run")
Eimeria.plot <- parasiteLoad::bananaPlot(mod = Eimeria$H1,
                                          data = Genus.samples[Genus.samples$Eimeria.trans>0,],
                                          response = "Eimeria.trans",
                                          islog10 = T, group = "Seq_Run",
                                          cols = c("#006A4E", "#006A4E"))+ labs(tag = "C)")+
  ggtitle("Eimeria")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Tritrichomonas<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Tritrichomonas.trans>0,], response = "Tritrichomonas.trans", model = "negbin", group = "Seq_Run")
Tritrichomonas.plot <- parasiteLoad::bananaPlot(mod = Tritrichomonas$H1,
                                         data = Genus.samples[Genus.samples$Tritrichomonas.trans>0,],
                                         response = "Tritrichomonas.trans",
                                         islog10 = T, group = "Seq_Run",
                                         cols = c("#006A4E", "#006A4E"))+ labs(tag = "D)")+
  ggtitle("Tritrichomonas")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Cryptosporidium<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Cryptosporidium.trans>0,], response = "Cryptosporidium.trans", model = "negbin", group = "Seq_Run")
Cryptosporidium.plot <-  bananaplotNoCI(mod = Cryptosporidium$H1,
                                                data = Genus.samples[Genus.samples$Cryptosporidium.trans>0,],
                                                response = "Cryptosporidium.trans",
                                                islog10 = T, group = "Seq_Run",
                                                cols = c("#006A4E", "#006A4E"))+ labs(tag = "E)")+
  ggtitle("Cryptosporidium")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Hymenolepis<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Hymenolepis.trans>0,], response = "Hymenolepis.trans", model = "negbin", group = "Seq_Run" )
Hymenolepis.plot <- bananaplotNoCI(mod = Hymenolepis$H1,
                                                 data = Genus.samples[Genus.samples$Hymenolepis.trans>0,],
                                                 response = "Hymenolepis.trans",
                                                 islog10 = T, group = "Seq_Run",
                                                 cols = c("#006A4E", "#006A4E"))+ labs(tag = "F)")+
  ggtitle("Hymenolepis")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Bacteroides<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Bacteroides.trans>0,], response = "Bacteroides.trans", model = "negbin", group = "Seq_Run")
Bacteroides.plot <- parasiteLoad::bananaPlot(mod = Bacteroides$H1,
                                             data = Genus.samples[Genus.samples$Bacteroides.trans>0,],
                                             response = "Bacteroides.trans",
                                             islog10 = T, group = "Seq_Run",
                                             cols = c("#006A4E", "#006A4E"))+ labs(tag = "A)")+
  ggtitle("Bacteroides")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Helicobacter<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Helicobacter.trans>0,], response = "Helicobacter.trans", model = "negbin", group = "Seq_Run")
Helicobacter.plot <- parasiteLoad::bananaPlot(mod = Helicobacter$H1,
                                             data = Genus.samples[Genus.samples$Helicobacter.trans>0,],
                                             response = "Helicobacter.trans",
                                             islog10 = T, group = "Seq_Run",
                                             cols = c("#006A4E", "#006A4E"))+ labs(tag = "B)")+
  ggtitle("Helicobacter")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Blautia<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Blautia.trans>0,], response = "Blautia.trans", model = "negbin", group = "Seq_Run")
Blautia.plot <- parasiteLoad::bananaPlot(mod = Blautia$H1,
                                              data = Genus.samples[Genus.samples$Blautia.trans>0,],
                                              response = "Blautia.trans",
                                              islog10 = T, group = "Seq_Run",
                                              cols = c("#006A4E", "#006A4E"))+ labs(tag = "C)")+
  ggtitle("Blautia")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Lactobacillus<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Lactobacillus.trans>0,], response = "Lactobacillus.trans", model = "negbin", group = "Seq_Run")
Lactobacillus.plot <- parasiteLoad::bananaPlot(mod = Lactobacillus$H1,
                                         data = Genus.samples[Genus.samples$Lactobacillus.trans>0,],
                                         response = "Lactobacillus.trans",
                                         islog10 = T, group = "Seq_Run",
                                         cols = c("#006A4E", "#006A4E"))+ labs(tag = "C)")+
  ggtitle("Lactobacillus")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Clostridium<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Clostridium.trans>0,], response = "Clostridium.trans", model = "negbin", group = "Seq_Run")
Clostridium.plot <- parasiteLoad::bananaPlot(mod = Clostridium$H1,
                                               data = Genus.samples[Genus.samples$Clostridium.trans>0,],
                                               response = "Clostridium.trans",
                                               islog10 = T, group = "Seq_Run",
                                               cols = c("#006A4E", "#006A4E"))+ labs(tag = "C)")+
  ggtitle("Clostridium")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Fusimonas<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Fusimonas.trans>0,], response = "Fusimonas.trans", model = "negbin", group = "Seq_Run")
Fusimonas.plot <- parasiteLoad::bananaPlot(mod = Fusimonas$H1,
                                               data = Genus.samples[Genus.samples$Fusimonas.trans>0,],
                                               response = "Fusimonas.trans",
                                               islog10 = T, group = "Seq_Run",
                                               cols = c("#006A4E", "#006A4E"))+ labs(tag = "C)")+
  ggtitle("Fusimonas")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

Muribaculum<- parasiteLoad::analyse(data = Genus.samples[Genus.samples$Muribaculum.trans>0,], response = "Muribaculum.trans", model = "negbin", group = "Seq_Run")
Muribaculum.plot <- parasiteLoad::bananaPlot(mod = Muribaculum$H1,
                                               data = Genus.samples[Genus.samples$Muribaculum.trans>0,],
                                               response = "Muribaculum.trans",
                                               islog10 = T, group = "Seq_Run",
                                               cols = c("#006A4E", "#006A4E"))+ labs(tag = "C)")+
  ggtitle("Muribaculum")+
  xlab("Mouse genotype (Hybrid Index)")+
  ylab("log 10 transformed rel. abundance")+
  theme(legend.position ="none")

##Multiple comparisons correction
Input<- read.csv("~/AA_HMHZ/ML_Analysis/Maximum_likelihood_results.csv")
Input%>%
  slice(1:7) -> Data.phylum
Input%>%
  slice(8:15) -> Data.genus

Data<- Input[order(Input$Raw_p),]
Data.phylum<- Data.phylum[order(Data.phylum$Raw_p),]
Data.genus<- Data.genus[order(Data.genus$Raw_p),]

Data$Bonferroni <- p.adjust(Data$Raw_p, method = "bonferroni")
Data.phylum$Bonferroni <- p.adjust(Data.phylum$Raw_p, method = "bonferroni")
Data.genus$Bonferroni <- p.adjust(Data.genus$Raw_p, method = "bonferroni")

###Compile figure for: 
##Phylum
pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Phylum.pdf", width = 20, height = 50)
grid.arrange(Ascomycota.plot, Basidiomycota.plot, Nematoda.plot,
             Apicomplexa.plot, Platyhelminthes.plot, Ciliophora.plot,
             Firmicutes.plot, Bacteroidetes.plot, Proteobacteria.plot, Fusobacteria.plot, ncol= 2, nrow= 5)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Phylum_Fungi.pdf", width = 20, height = 10)
grid.arrange(Ascomycota.plot, Basidiomycota.plot, ncol= 2, nrow= 1)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Phylum_Basidiomycota.pdf", width = 10, height = 8)
Basidiomycota.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Phylum_Parasites.pdf", width = 20, height = 20)
grid.arrange(Nematoda.plot, Apicomplexa.plot, 
             Platyhelminthes.plot, Ciliophora.plot, ncol= 2, nrow= 2)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Phylum_Nematoda.pdf", width = 10, height = 8)
Nematoda.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Phylum_Bacteria.pdf", width = 20, height = 20)
grid.arrange(Firmicutes.plot, Bacteroidetes.plot,
             Proteobacteria.plot, Fusobacteria.plot, ncol= 2, nrow= 2)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Phylum_Proteobacteria.pdf", width = 10, height = 8)
Proteobacteria.plot
dev.off()

##Genus
pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus.pdf", width = 40, height = 30)
grid.arrange(Kazachstania.plot, Aspiculuris.plot, Syphacia.plot,
             Eimeria.plot, Tritrichomonas.plot, Cryptosporidium.plot, Hymenolepis.plot,
             Bacteroides.plot, Helicobacter.plot, Blautia.plot,ncol= 4, nrow= 3)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Parasites.pdf", width = 20, height = 30)
grid.arrange(Aspiculuris.plot, Syphacia.plot,
             Eimeria.plot, Tritrichomonas.plot, 
             Cryptosporidium.plot, Hymenolepis.plot,
             ncol= 2, nrow= 3)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Bacteria.pdf", width = 30, height = 10)
grid.arrange(Bacteroides.plot, Helicobacter.plot, Blautia.plot,ncol= 3, nrow= 1)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Kazachstania.pdf", width = 10, height = 10)
Kazachstania.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Syp.pdf", width = 10, height = 8)
Syphacia.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Eimeria.pdf", width = 10, height = 8)
Eimeria.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Tritrichomonas.pdf", width = 10, height = 8)
Tritrichomonas.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Helicobacter.pdf", width = 10, height = 8)
Helicobacter.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Blautia.pdf", width = 10, height = 8)
Blautia.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Genus_Bacteroides.pdf", width = 10, height = 8)
Bacteroides.plot
dev.off()


##Alpha diversity indexes
pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Alpha_div.pdf", width = 30, height = 10)
grid.arrange(Chao1.plot, Shannon.plot, Pielou.plot, ncol= 3, nrow= 1)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Chao1.pdf", width = 10, height = 10)
Chao1.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Shannon.pdf", width = 10, height = 10)
Shannon.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Pielou.pdf", width = 10, height = 10)
Pielou.plot
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Alpha_div_Bac.pdf", width = 30, height = 10)
grid.arrange(Chao1.bacteria.plot, Shannon.bacteria.plot, Pielou.bacteria.plot, ncol= 3, nrow= 1)
dev.off()

pdf(file = "~/AA_HMHZ/ML_Analysis/ML_Alpha_div_Euk.pdf", width = 30, height = 10)
grid.arrange(Chao1.eukaryota.plot, Shannon.eukaryota.plot, Pielou.eukaryota.plot, ncol= 3, nrow= 1)
dev.off()

####Raw reads parasites (Not useful any more)

#parasiteLoad::getParamBounds(model = "normal",
#              data = Phylum.samples,
#              response = "Ascomycota.trans")

#fitPara <- parasiteLoad::analyse(data = sdt, response = "Readspara", model = "negbin", group = "Seq_Run")
#fitPara <- parasiteLoad::analyse(data = sdt, response = "Para_abundance", model = "negbin", group = "Seq_Run", myparamBounds = speparam)
#fitPara <- parasiteLoad::analyse(data = sdt[sdt$Readspara>0,], response = "Readspara", model = "negbin", group = "Seq_Run")


#Paraban <- parasiteLoad::bananaPlot(mod = fitPara$H1,
#                                    data = sdt,
#                                    response = "Readspara",
#                                    islog10 = FALSE, group = "Seq_Run",
#                                    cols = c("#006A4E", "#006A4E"))+ labs(tag = "B)")+
#  xlab("Mouse genotype (Hybrid Index)")+
#  ylab("Read counts Parasites (Nem+Api+Platy)")+
#  theme(legend.position ="none" )

###Other parasites

#fitPara2 <- parasiteLoad::analyse(data = sdt, response = "Para_reads_all", model = "negbin", group = "Seq_Run")
#fitPara2 <- parasiteLoad::analyse(data = sdt[sdt$Para_reads_all>0,], response = "Para_reads_all", model = "negbin", group = "Seq_Run")


#Paraban2 <- parasiteLoad::bananaPlot(mod = fitPara2$H1,
#                                    data = sdt,
#                                    response = "Para_reads_all",
#                                    islog10 = FALSE, group = "Seq_Run",
#                                    cols = c("#006A4E", "#006A4E"))+ labs(tag = "A)")+
#  xlab("Mouse genotype (Hybrid Index)")+
#  ylab("Read counts Protozoa + Helminths")+
#  theme(legend.position ="none" )

#speparam <- c(L1start = 0,
#              L1LB = -20,
#              L1UB = 4000,
#              L2start = 0,
#              L2LB = -20,
#              L2UB = 4000,
#              alphaStart = 0, alphaLB = -100, alphaUB = 100,
#              A1start = 1, A1LB = 1.e-9, A1UB = 1000,
#              A2start = 1, A2LB = 1.e-9, A2UB = 1000,
#              Zstart = 0, ZLB = -100, ZUB = 100)

#fitcili <- parasiteLoad::analyse(data = sdt, response = "Readscili", model = "negbin", group = "Seq_Run") #, myparamBounds = speparam

#fitcili_pos <- parasiteLoad::analyse(data = sdt[sdt$Readscili>0,], response = "Readscili", model = "negbin", group = "Seq_Run", myparamBounds = speparam)

# improve paramBounds
#speparam <- c(L1start = 0,
#              L1LB = -20,
#              L1UB = 20,
#              L2start = 0,
#              L2LB = -20,
#              L2UB = 20,
#              alphaStart = 0, alphaLB = -100, alphaUB = 100,
#              A1start = 1, A1LB = 1.e-9, A1UB = 1000,
#              A2start = 1, A2LB = 1.e-9, A2UB = 1000,
#              Zstart = 0, ZLB = -2, ZUB = 2)


##All
#fitResiduals_eimeria <-
#  parasiteLoad::analyse(data = body_data_eimeria,
#                        response = "residuals",
#                        model = "normal",
#                        group = "presence_eimeria_tissues",
#                        myparamBounds = speparam)


#parasiteLoad::getParamBounds()

# To play with bigger sd

#automaticParam <- parasiteLoad::getParamBounds(model = "negbin",
#                                               data = Genus.samples,
#                                               response = "")


#automaticParam["mysdUB"] <- 50

#parasiteLoad::analyse(data = alphaDiv, response = "Chao1", model = "normal", group = "SexAl",
#                      myparamBounds = automaticParam)
