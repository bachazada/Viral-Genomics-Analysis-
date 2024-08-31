rm(list=ls())
setwd("C:/Users/Bacha Zada/Desktop/tttt2/")
#### exercise 1 ####
print("my first script")

myFirstF <- function(n)
  {
    n*n
}
x <- runif(100)
y <- myFirstF(x)
xy <- list(x,y)
save(xy,file = "xy.rdat")
load("xy.rdat")
saveRDS(object = xy,file = "xy.rds")
new_object <- readRDS("xy.rds")

#### exercise 2 ####
#install.packages("ape")
#install.packages("rentrez")

library("ape")
library("rentrez")
seq0 <- read.GenBank('AY274119.3', species.names = TRUE, as.character = TRUE)
seq0$AY274119.3
seq <- paste(seq0$AY274119.3, collapse = "")
seq <- toupper(seq)

library(Biostrings)
dna1 <- DNAString("TTGATATGGCCCTTATAA")
class(dna1)
translate(dna1)
nchar(seq)
GENETIC_CODE[["TTG"]]
GENETIC_CODE[["ATA"]]
GENETIC_CODE[["TGG"]]
GENETIC_CODE[["CCC"]]
GENETIC_CODE[["TTA"]]
GENETIC_CODE[["TAA"]]

seqBioStr <- DNAString(seq)
prot_seq <- translate(seqBioStr)
saveRDS(object = prot_seq, file = "SARS_CoV_2_AA.rds")
#BiocManager::install("ORFik")
library("ORFik")

s <- findORFs(
  seq,
  #startCodon = startDefinition(1),
  #stopCodon = stopDefinition(1),
  #longestORF = TRUE,
  minimumLength = 100
)
s
seq0 <- seq0$AY274119.3
str(s)
s@unlistData@start[8]
s@unlistData@width[8]
codingregion <- seq0[s@unlistData@start[8]:(s@unlistData@start[8]+s@unlistData@width[8])]
codingregion <- paste(codingregion, collapse = "")
codingregion <- toupper(codingregion)

codingregion <- DNAString(codingregion)
codingregion_translated <- translate(codingregion)
str(codingregion_translated)
codingregion_translated <- as.character(codingregion_translated)
nchar(codingregion_translated)

saveRDS(object = codingregion_translated, "own_sequence.rds")



#### retrieve virus data on your own ####
dat = c('Bovine CoV1',           'AAL57308.1',
  'Bovine CoV2'   ,        'AAK83356.1',
  'Human CoV OC43'   ,     'NP_937950.1',
  'Porcine HEV3'      ,    'AAL80031',
  'Murine HV1'         ,   'YP_209233.1',
  'Murine HV2'          ,  'NP_045300.1',
  'IBV3'                 , 'NP_040831.1',
  'Porcine PEDV'          ,'NP_598310.1',
  'Canine CoV1'     ,      'S41453',
  'Feline CoV4'      ,     'BAA06805',
  'Human Sars CoV'    ,    'AAP41037.1',# % genbank nucleotide id 'AY274119.3';
  'Palm civet'         ,   'AAV49723')
other_viruses <- matrix(dat, ncol=2, nrow=12, byrow=TRUE)

library("ape")#, 
library("rentrez")
prot_seqs <- list()
DNA_seqs <- list()
for(i in 1:nrow(other_viruses)){
  temp <- rentrez::entrez_fetch(id = other_viruses[i,2],
                                          db = "protein", 
                                          rettype = "fasta")
  temp2 <- strsplit(temp,"\n")
  temp3 <- paste(temp2[[1]][2:length(temp2[[1]])], collapse = "")
  prot_seqs[[i]] <- temp3
}
names(prot_seqs) <- other_viruses[,1]
saveRDS(object = prot_seqs, file = "viruses.rds")



#### exercise 3 ####
prot_seqs <- readRDS("viruses.rds")
prot_SARS <- prot_seqs$`Human Sars CoV`
nchar(prot_SARS)

codingregion_translated <- gsub("*","",codingregion_translated, fixed = TRUE)
nchar(prot_SARS)
nchar(codingregion_translated)
prot_SARS==codingregion_translated #Yippie!


library("msa")
seqs_for_msa <- c(prot_seqs[[1]],prot_seqs[[2]],prot_seqs[[3]],prot_seqs[[4]],
                  prot_seqs[[5]],prot_seqs[[6]],prot_seqs[[7]],prot_seqs[[8]],
                  prot_seqs[[9]],prot_seqs[[10]],prot_seqs[[11]],prot_seqs[[12]])
names(seqs_for_msa) <- (other_viruses[,1])
seq_alignment <- msa(seqs_for_msa, type="protein")
msaPrettyPrint(seq_alignment, file="test.pdf", output="pdf",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="rasmol",
               verbose=FALSE, askForOverwrite=FALSE)

library(phangorn)

ph.dat <- as.phyDat(seq_alignment,type = "AA")
dists <- dist.ml(ph.dat, model="JC69")
tree1 <- NJ(dists)
plot(tree1)


#### exercise 4 ####
rm(list=ls())
#install.packages("seqinr")
library(seqinr)
dat <- read.fasta("sars_spike.fasta")

seqs <- NULL
nam <- NULL
for(i in 1:length(dat)){
  seqs <- c(seqs, toupper(paste(dat[[i]], collapse="")))
  nam <- c(nam, attributes(dat[[i]])$Annot)
}
names(seqs) <- nam
v <- c(1,2,3,4)
names(v) <- c("A","B","C","D")
v["A"]
library("msa")

seq_alignment <- msa(seqs, type="dna")
msaPrettyPrint(seq_alignment, file="test_2.pdf", output="pdf",
               showNames="left", showNumbering="none", showLogo="top",
               showConsensus="bottom", logoColors="rasmol",
               verbose=FALSE, askForOverwrite=FALSE)


library(phangorn)

ph.dat <- as.phyDat(seq_alignment,type = "dna")
dists <- dist.ml(ph.dat, model="JC69")
#dists <- dist.ml(ph.dat, model="F81")
#dists <- dist.logDet(ph.dat)

tree1 <- NJ(x = dists)
pdf("tree.pdf")
plot(tree1,"phylogram", main="NJ")
dev.off()

fit <- upgma(dists)
plot(fit)

fit <- wpgma(dists)
pdf("tree_wpgma.pdf")
plot(fit,"phylogram", main="NJ")
dev.off()


