############################################
#Purpose: Code to do compile gene lists
#Author: Pichai Raman
#Date: 5/10/2018
############################################

#Call libraries
library("tidyverse")

myCancerGeneList <- list();

############################################
#COSMIC Data - 
#From Here : https://cancer.sanger.ac.uk/census#cl_search
cosmicData <- read.csv("../data/raw/CancerGenes/CancerGeneCensus/Census_allThu May 10 17_20_12 2018.csv")
cosmicGenes <- unique(cosmicData[,1]); # 716 Unique Genes
myCancerGeneList[["COSMIC"]] <- as.character(cosmicGenes)

############################################
#FoundationOne Panel, MyCancerGenome, oncoKB, MskImpact, Cancer Commons
#From Here : http://www.dgidb.org/downloads (http://www.dgidb.org/data/genes.tsv)
############################################
dgiDBGenes <- read.delim("../data/raw/CancerGenes/DGIDB/genes.tsv")

#OncoKB
oncoKBGenes <- dgiDBGenes %>% filter(gene_claim_source=="OncoKB") %>% select(gene_name)
oncoKBGenes <- unique(oncoKBGenes[,1]); # 50 Genes
myCancerGeneList[["OncoKB"]] <- as.character(oncoKBGenes)

#Foundation Medicine
fmGenes <- dgiDBGenes %>% filter(gene_claim_source=="FoundationOneGenes") %>% select(gene_name)
fmGenes <- unique(fmGenes[,1]); # 243 Genes
myCancerGeneList[["FM"]] <- as.character(fmGenes)

#MSK Genes
mskGenes <- dgiDBGenes %>% filter(gene_claim_source=="MskImpact") %>% select(gene_name)
mskGenes <- unique(mskGenes[,1]); # 350 Genes
myCancerGeneList[["MskImpact"]] <- as.character(fmGenes)

#MyCancerGenome Genes
mcgGenes <- dgiDBGenes %>% filter(gene_claim_source=="MyCancerGenome") %>% select(gene_name)
mcgGenes <- unique(mcgGenes[,1]); # 181 Genes
myCancerGeneList[["MyCancerGenome"]] <- as.character(mcgGenes)

#MyCancerGenome CLinical Trials Genes
mcgCTGenes <- dgiDBGenes %>% filter(gene_claim_source=="MyCancerGenomeClinicalTrial") %>% select(gene_name)
mcgCTGenes <- unique(mcgCTGenes[,1]); # 79 Genes
myCancerGeneList[["MyCancerGenomeClinTrial"]] <- as.character(mcgCTGenes)

#Cancer Commons Genes
ccGenes <- dgiDBGenes %>% filter(gene_claim_source=="CancerCommons") %>% select(gene_name)
ccGenes <- unique(ccGenes[,1]); # 48 Genes
myCancerGeneList[["CancerCommons"]] <- as.character(ccGenes)

#Cancer Genome Interpreter
cgiGenes <- dgiDBGenes %>% filter(gene_claim_source=="CGI") %>% select(gene_name)
cgiGenes <- unique(cgiGenes[,1]); # 118 Genes
myCancerGeneList[["CancerGenomeInt"]] <- as.character(cgiGenes)

#Caris Molecular Intelligence
cmiGenes <- dgiDBGenes %>% filter(gene_claim_source=="CarisMolecularIntelligence") %>% select(gene_name)
cmiGenes <- unique(cmiGenes[,1]); # 60 Genes
myCancerGeneList[["CarisMolIntel"]] <- as.character(cmiGenes)

#Jackson Labs
ckbGenes <- dgiDBGenes %>% filter(gene_claim_source=="CKB") %>% select(gene_name)
ckbGenes <- unique(ckbGenes[,1]); # 86 Genes
myCancerGeneList[["JacksonLabs"]] <- as.character(ckbGenes)

#NCI 
nciGenes <- dgiDBGenes %>% filter(gene_claim_source=="NCI") %>% select(gene_name)
nciGenes <- unique(nciGenes[,1]); # 1037 Genes
myCancerGeneList[["NCI"]] <- as.character(nciGenes)

#CIVIC 
cvcGenes <- dgiDBGenes %>% filter(gene_claim_source=="CIViC") %>% select(gene_name)
cvcGenes <- unique(cvcGenes[,1]); # 192 Genes
myCancerGeneList[["CIVIC"]] <- as.character(cvcGenes)

############################################

############################################
#cBioPortal
#From Here : https://pedcbioportal.org/cancer_gene_list.jsp
############################################

cbioPortalGenes <- unique(read.delim("../data/raw/CancerGenes/cbioPortal/genes.txt", header=F)[,1])
myCancerGeneList[["cBioPortal"]] <- as.character(cbioPortalGenes)

############################################
#Pediatric Landscape Papers
############################################

#TARGET St. Jude 
#From here https://www.nature.com/articles/nature25795#supplementary-information
#Supplementary table 2 

targetSTJude <- unique(read.delim("../data/raw/CancerGenes/pediatricPublications/TARGET_STJUDE.txt", header=F)[,1])
myCancerGeneList[["PEDStJude"]] <- as.character(targetSTJude)

#Pfister Paper
#From here https://www.nature.com/articles/nature25480#supplementary-information
#Supplementary table 20

pfister <- unique(read.delim("../data/raw/CancerGenes/pediatricPublications/PFISTER_PEDLANDSCAPE.txt", header=F)[,1])
myCancerGeneList[["PEDPfister"]] <- as.character(pfister)

############################################
#TCGA Landscape Paper
############################################

#Comprehensive Characterization of Cancer Driver Genes and Mutations
#From Here : https://www.cell.com/cell/fulltext/S0092-8674(18)30237-X
#Supplementary table 1

tcgaBailey <- unique(read.delim("../data/raw/CancerGenes/TCGA/TCGA_BAILEY.txt")[,1])
myCancerGeneList[["TCGABailey"]] <- as.character(tcgaBailey)


#Given a set of lists this function will create a matrix

getOverlap <- function(fullList, y)
{
	out <- as.numeric(fullList%in%y)
	return(out);	
}

createMatrix <- function(x)
{
	allGenes <- unique(as.character(unlist(x)));
	out <- do.call(rbind, lapply(x, FUN=getOverlap, fullList=allGenes))
	out <- data.frame(t(out));
	rownames(out) <- allGenes;
	return(out);

}

CancerGeneMatrix <- createMatrix(myCancerGeneList)
CancerGeneMatrix[,"Count"] <- rowSums(CancerGeneMatrix)
write.table(CancerGeneMatrix, "../data/processed/CancerGeneList.txt", sep="\t", row.names=T)

barplot(table(CancerGeneMatrix[,"Count"]))
colSums(CancerGeneMatrix)












