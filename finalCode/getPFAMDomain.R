#####################################################
#Purpose: Obtain PFAM ID's associated with output from Fusion Prioritization Pipeline       
#Author: Pichai Raman
#Date: June 25, 2018
#####################################################



#This function will get the 
#Get PFAMgenomic coordinates from UCSC
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/pfamDesc.txt.gz
# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ucscGenePfam.txt.gz
# PFAM export (map gene and pfam) from http://www.ensembl.org/biomart/martview/2a436835bd0ca1f076c993bc1e062a5e


#Will create a table for other functions to use
pfamDesc <- read.delim("/home/nick/Desktop/Fusion_prioritization/ignore/data/pfamDesc.txt", header=F);
colnames(pfamDesc) <- c("PFAM_ID", "PFAM_SD", "PFAM_DESC");
pfamLoc <- read.delim("/home/nick/Desktop/Fusion_prioritization/ignore/data/ucscGenePfam.txt", header=F)
colnames(pfamLoc) <- c("BIN", "CHROM", "CHROM_START", "CHROM_END", "NAME", "SCORE",
						"STRAND", "THICK_START", "THICK_END", "RESERVED", "BLOCK_COUNT",
						 "BLOCK_SIZES", "CHROM_STARTS");
pfamDescLoc <- merge(pfamDesc, pfamLoc, by.x="PFAM_SD", by.y="NAME", all.y=T)
pfamDescLoc[,"PFAM_ID_DESC"] <- paste(pfamDescLoc[,"PFAM_ID"], pfamDescLoc[,"PFAM_DESC"], sep=": ")
pfamDescLocFusion <- unique(pfamDescLoc[,c("PFAM_ID", "PFAM_ID_DESC", "CHROM", "CHROM_START", "CHROM_END")]);
pfamDescLocFusion[,"CHROM"] <- gsub("chr", "", pfamDescLocFusion[,"CHROM"]);

pfamDataBioMart <- read.delim("/home/nick/Desktop/Fusion_prioritization/ignore/data/Biomart_PFAM_Export.txt", header=T);
colnames(pfamDataBioMart) <- c("start_position", "end_position", "pfam" ,"pfam_start", "pfam_end", "hgnc_symbol");

#Support function to collapse a tall skinny list
#so that data is short and comma separated (serialized)
collapseData <- function(geneDomain=NULL)
{
	out <- geneDomain %>% group_by(hgnc_symbol) %>% summarise(PFAM=paste(PFAM_ID_DESC, collapse="; "))
	out <- data.frame(out);
	return(out);
}


#Input's 
getPFAMDomain <- function(starFusionCombo=NULL)
{
	#Call library
	require("tidyverse")

	#Get data for H Gene
	hGenes <- as.character(starFusionCombo[,"H_Gene"]);
	pfamData <- pfamDataBioMart[pfamDataBioMart[,"hgnc_symbol"]%in%hGenes,]
	pfamData <- na.omit(pfamData);
	
	#Add all PFAM Domains
	pfamTmp <- unique(pfamData[,c("hgnc_symbol", "pfam")]);
	pfamTmp <- unique(merge(pfamTmp, unique(pfamDescLoc[,c("PFAM_ID", "PFAM_ID_DESC")]), by.x="pfam", by.y="PFAM_ID"));
	hGeneAllDomains <- collapseData(pfamTmp);
	colnames(hGeneAllDomains)[2] <- "H_Gene_PFAM_All"
	starFusionCombo <- merge(starFusionCombo, hGeneAllDomains, by.x="H_Gene", by.y="hgnc_symbol", all.x=T);

	#####################################
	#Add Only overlapping domains
	#####################################

	#Pull out STAR Fusion 
	starFusionComboTmp <- starFusionCombo[,c("H_Gene", "Left_Breakpoint_Chr", "Left_Breakpoint_Pos")]
	pfamTmp <- unique(pfamData[,c("hgnc_symbol", "start_position")]);
	pfamTmp <- pfamTmp %>% group_by(hgnc_symbol) %>% summarise(start_position=min(start_position))
	pfamTmp <- data.frame(pfamTmp);
	starFusionComboTmp <- na.omit(merge(starFusionComboTmp, pfamTmp, by.x="H_Gene", by.y="hgnc_symbol", all.x=T))
	starFusionComboTmp <- merge(starFusionComboTmp, unique(pfamData[,c("hgnc_symbol", "pfam")]), by.x="H_Gene", by.y="hgnc_symbol")
	starFusionComboTmp <- merge(starFusionComboTmp, pfamDescLocFusion, 
						  by.x=c("pfam", "Left_Breakpoint_Chr"),
						  by.y=c("PFAM_ID", "CHROM"));

	starFusionComboTmp[,"inside"] <- ifelse((starFusionComboTmp[,"CHROM_START"]>starFusionComboTmp[,"start_position"])&(starFusionComboTmp[,"CHROM_END"]<starFusionComboTmp[,"Left_Breakpoint_Pos"]), 1, 0);
	starFusionComboTmp <- starFusionComboTmp[starFusionComboTmp[,"inside"]==1,]
	starFusionComboTmp2 <- unique(starFusionComboTmp[,c("H_Gene", "PFAM_ID_DESC", "Left_Breakpoint_Pos")]);
	starFusionComboTmp2 <- starFusionComboTmp2 %>% group_by(H_Gene, Left_Breakpoint_Pos) %>% summarise(PFAM=paste(PFAM_ID_DESC, collapse="; "))
	starFusionComboTmp2 <- data.frame(starFusionComboTmp2);
	colnames(starFusionComboTmp2)[3] <- "H_Gene_PFAM_IN_FUSION"
	starFusionCombo <- merge(starFusionCombo, starFusionComboTmp2, by=c("H_Gene", "Left_Breakpoint_Pos"), all.x=T)

	#Get data for T Gene
	tGenes <- as.character(starFusionCombo[,"T_Gene"]);
	pfamData <- pfamDataBioMart[pfamDataBioMart[,"hgnc_symbol"]%in%tGenes,]
	pfamData <- na.omit(pfamData);
	
	#Add all PFAM Domains
	pfamTmp <- unique(pfamData[,c("hgnc_symbol", "pfam")]);
	pfamTmp <- unique(merge(pfamTmp, unique(pfamDescLoc[,c("PFAM_ID", "PFAM_ID_DESC")]), by.x="pfam", by.y="PFAM_ID"));
	tGeneAllDomains <- collapseData(pfamTmp);
	colnames(tGeneAllDomains)[2] <- "T_Gene_PFAM_All"
	starFusionCombo <- merge(starFusionCombo, tGeneAllDomains, by.x="T_Gene", by.y="hgnc_symbol", all.x=T);

	#####################################
	#Add Only overlapping domains
	#####################################

	#Pull out STAR Fusion 
	starFusionComboTmp <- starFusionCombo[,c("T_Gene", "Right_Breakpoint_Chr", "Right_Breakpoint_Pos")]
	pfamTmp <- unique(pfamData[,c("hgnc_symbol", "end_position")]);
	pfamTmp <- pfamTmp %>% group_by(hgnc_symbol) %>% summarise(end_position=max(end_position))
	pfamTmp <- data.frame(pfamTmp);
	starFusionComboTmp <- na.omit(merge(starFusionComboTmp, pfamTmp, by.x="T_Gene", by.y="hgnc_symbol", all.x=T))
	starFusionComboTmp <- merge(starFusionComboTmp, unique(pfamData[,c("hgnc_symbol", "pfam")]), by.x="T_Gene", by.y="hgnc_symbol")
	starFusionComboTmp <- merge(starFusionComboTmp, pfamDescLocFusion, 
						  by.x=c("pfam", "Right_Breakpoint_Chr"),
						  by.y=c("PFAM_ID", "CHROM"));

	starFusionComboTmp[,"inside"] <- ifelse((starFusionComboTmp[,"CHROM_START"]>starFusionComboTmp[,"Right_Breakpoint_Pos"])&(starFusionComboTmp[,"CHROM_END"]<starFusionComboTmp[,"end_position"]), 1, 0);
	starFusionComboTmp <- starFusionComboTmp[starFusionComboTmp[,"inside"]==1,]
	starFusionComboTmp2 <- unique(starFusionComboTmp[,c("T_Gene", "PFAM_ID_DESC", "Right_Breakpoint_Pos")]);
	starFusionComboTmp2 <- starFusionComboTmp2 %>% group_by(T_Gene, Right_Breakpoint_Pos) %>% summarise(PFAM=paste(PFAM_ID_DESC, collapse="; "))
	starFusionComboTmp2 <- data.frame(starFusionComboTmp2);
	colnames(starFusionComboTmp2)[3] <- "T_Gene_PFAM_IN_FUSION"
	starFusionCombo <- merge(starFusionCombo, starFusionComboTmp2, by=c("T_Gene", "Right_Breakpoint_Pos"), all.x=T)

	return(starFusionCombo);

}



#Example
#starFusionCombo <- read.delim("../data/processed/star_fusion_combo.tsv");
#starFusionCombo2 <- getPFAMDomain(starFusionCombo);

