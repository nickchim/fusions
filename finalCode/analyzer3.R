library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)
#library(biomaRt)


source("/home/nick/Desktop/Fusion_prioritization/finalCode/getPFAMDomain.R")

#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
input <- read.delim("/home/nick/Desktop/Fusion_prioritization/processed/star_fusion_combo.tsv")

#Cancer Gene and Fusion Prep
geneList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/processed/CancerGeneList.tsv")
fuseList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/processed/FusionList2.txt")

input$Cancerous_Gene <- F
input$Cancerous_Gene_Symbol <- NA
input$Cancerous_Fusion <- F
input$Cancerous_Fusion_Symbol <- NA
input$Fusion_Cancer_Type <- NA

geneRep <- 1
fuseRep <- 1
repAdjustedGene <- geneList[geneList$Count >= geneRep,]
repAdjustedFuse <- fuseList[fuseList$Count >= fuseRep,]

#Target Prep
input$Targetable <- FALSE
input$Targetable_gene <- NA
input$Drug_name <- NA
input$Drug_chembl_id <- NA

Target_genes <- read_tsv("/home/nick/Desktop/Fusion_prioritization/processed/Target_List.tsv")

#PFAM Domain
input <- getPFAMDomain(input)
#input <- getPFAMDomain(input) #Need to fix this, have to run twice in order for it to work
# input <- as.tibble(input)
# input <- subset(input, select=-c(H_Gene_PFAM_All.x,T_Gene_PFAM_All.x))
# colnames(input)[colnames(input)=="H_Gene_PFAM_All.y"] <- "H_Gene_PFAM_All"
# colnames(input)[colnames(input)=="T_Gene_PFAM_All.y"] <- "T_Gene_PFAM_All"

#Domain Prep
input$Cancer_Domains <- NA
domainList <-read_tsv("/home/nick/Desktop/Fusion_prioritization/processed/Domain_List.tsv")

#Big Apply Function
bigApply <- function(x) {
  #Cancer G+F
  gene1 <- x["H_Gene"]
  gene2 <- x["T_Gene"]
  fuse <- paste(gene1,gene2,sep="--")
  
  gene1In <- gene1 %in% repAdjustedGene$Gene
  gene2In <- gene2 %in% repAdjustedGene$Gene
  fuseIn <- fuse %in% repAdjustedFuse$Fusions
  
  x["Cancerous_Gene"] <- gene1In || gene2In
  x["Cancerous_Fusion"] <- fuseIn
  
  if (gene1In && gene2In) {
    x["Cancerous_Gene_Symbol"] <- paste(gene1,gene2,sep=",")
  }else if (gene1In) {
    x["Cancerous_Gene_Symbol"] <- gene1
  }else if (gene2In) {
    x["Cancerous_Gene_Symbol"] <- gene2
  }
  
  if (fuseIn) {
    x["Cancerous_Fusion_Symbol"] <- paste(gene1,gene2,sep="--")
    
    cancers <- unlist(str_split(fuseList[fuseList$Fusions == paste(gene1,gene2,sep="--"),"Cancer_Types"],","))
    indices <- which(fuseList[fuseList$Fusions == paste(gene1,gene2,sep="--"),-c(1,2,ncol(fuseList))] != 0) + 2
    for (i in 1:length(cancers)) {
      if (!is.na(cancers[i])) {
        if (i == 1) {
          x["Fusion_Cancer_Type"] <- paste(cancers[i],indices[i],sep=":")
        }else {
          x["Fusion_Cancer_Type"] <- paste(x["Fusion_Cancer_Type"],paste(cancers[i],indices[i],sep=":"),sep=",")
        }
      }
    }
  }
  
  #Target
  if (x["H_Gene"] %in% Target_genes$Trgt_Genes) {
    x["Targetable"] <- TRUE
    
    index <- Target_genes$Trgt_Genes == x["H_Gene"]
    
    x["Targetable_gene"] <- x["H_Gene"]
    x["Drug_name"] <- Target_genes$Drug_name[index]
    x["Drug_chembl_id"] <- Target_genes$Drug_chembl_id[index]
  }
  
  if (x["T_Gene"] %in% Target_genes$Trgt_Genes) {
    x["Targetable"] <- TRUE
    
    index <- Target_genes$Trgt_Genes == x["T_Gene"]
    
    if (is.na(x["Targetable_gene"])) {
      x["Targetable_gene"] <- x["T_Gene"]
      x["Drug_name"] <- Target_genes$Drug_name[index]
      x["Drug_chembl_id"] <- Target_genes$Drug_chembl_id[index]
    }else {
      x["Targetable_gene"] <- paste(x["Targetable_gene"],x["T_Gene"],sep=",")
      x["Drug_name"] <- paste(x["Drug_name"],Target_genes$Drug_name[index],sep=",")
      x["Drug_chembl_id"] <- paste(x["Drug_chembl_id"],Target_genes$Drug_chembl_id[index],sep=",")
    }
    
  }
  
  #Domain Checker
  doms <- c()
  domList <- str_split(unlist(str_split(x["H_Gene_PFAM_IN_FUSION"],";")),":")
  for (i in 1:length(domList)) {
    doms <- c(doms,domList[[i]][1])
  }
  
  for (dom in doms) {
    
    if (is.na(dom)) {
      next
    }
    
    if (dom %in% domainList$PF) {
      # print(dom)
      if (is.na(x["Cancer_Domains"] )) {
        x["Cancer_Domains"] <- dom
      }else {
        x["Cancer_Domains"] <- paste(x["Cancer_Domains"],dom,sep=",")  
      }
    }
  }
  
  doms <- c()
  domList <- str_split(unlist(str_split(x["T_Gene_PFAM_IN_FUSION"],";")),":")
  for (i in 1:length(domList)) {
    doms <- c(doms,domList[[i]][1])
  }
  
  for (dom in doms) {
    if (is.na(dom)) {
      next
    }
    
    if (dom %in% domainList$PF) {
      #print(dom)
      if (is.na(x["Cancer_Domains"] )) {
        x["Cancer_Domains"] <- dom
      }else {
        x["Cancer_Domains"] <- paste(x["Cancer_Domains"],dom,sep=",")  
      }
    }
  }
  
  if (!is.na(x["Cancer_Domains"])) {
    list <- str_split(x["Cancer_Domains"],",")[[1]]
    x["Cancer_Domains"] <- paste(levels(factor(list)),collapse=",")
  }
  return(x)
}

input <- apply(input, MARGIN = 1, FUN = bigApply)
input <- t(input)
input <- as.tibble(input)
input <- input[,c(5,3,1,6,4,7,8,2,9,as.numeric(10:32))]

input <- add_column(input,H_Gene_PFAM_NOT_IN_FUSION=NA,.before = "T_Gene_PFAM_All")
input <- add_column(input,T_Gene_PFAM_NOT_IN_FUSION=NA,.before = "Cancer_Domains")

domainFunc <- function(x) {
  allFuse <- unlist(str_split(unlist(str_split(x["H_Gene_PFAM_All"],":")),"; "))
  allFuse <- allFuse[seq(1,length(allFuse),by=2)]
  inFuse <- unlist(str_split(unlist(str_split(x["H_Gene_PFAM_IN_FUSION"],":")),"; "))
  inFuse <- inFuse[seq(1,length(inFuse),by=2)]

  notInFuse <- allFuse[which(!allFuse %in% inFuse)]
  if (length(notInFuse) > 0) {x["H_Gene_PFAM_NOT_IN_FUSION"] <- paste(notInFuse,collapse=",")}



  allFuse <- unlist(str_split(unlist(str_split(x["T_Gene_PFAM_All"],":")),"; "))
  allFuse <- allFuse[seq(1,length(allFuse),by=2)]
  inFuse <- unlist(str_split(unlist(str_split(x["T_Gene_PFAM_IN_FUSION"],":")),"; "))
  inFuse <- inFuse[seq(1,length(inFuse),by=2)]


  notInFuse <- allFuse[which(!allFuse %in% inFuse)]
  if (length(notInFuse) > 0) {x["T_Gene_PFAM_NOT_IN_FUSION"] <- paste(notInFuse,collapse=",")}

  return(x)
}
# 
# 
input <- apply(input,FUN=domainFunc,MARGIN = 1)
input <- as.data.frame(t(input))


input$Tier <- 4

tierize <- function(x,input) {
  cancerousGene <- as.logical(x["Cancerous_Gene"])
  cancerousFusion <- as.logical(x["Cancerous_Fusion"])
  spanning <- FALSE
  if (!is.na(x["FC_Spanning_Pairs"]) && !is.na(x["FC_Spanning_Unique_Reads"])) {
    spanning <- x["FC_Spanning_Pairs"] > 3 || x["FC_Spanning_Unique_Reads"] > 3
  }else if (!is.na(x["Star_Spanning_Fragcount"]) && !is.na(x["Star_Junction_Readcount"])) {
    spanning <- x["Star_Spanning_Fragcount"] > 3 || x["Star_Junction_Readcount"] > 3
  }

  star <- x["Catcher"] == "Star"
  fusion <- x["Catcher"] == "Fusion_Catcher"
  starAndFusion <- FALSE
  if (sum(x["Fusion_Name"] == input[,"Fusion_Name"]) > 1) {
    if (any(input[x["Fusion_Name"] == input[,"Fusion_Name"],"Catcher"] != x["Catcher"])) {
      starAndFusion <- TRUE
    }
  }
  targetable <- as.logical(x["Targetable"])


  if (cancerousGene & spanning & starAndFusion & targetable) {
    tier <- "1A"
  }else if (cancerousFusion & spanning & starAndFusion & targetable) {
    tier <- "1B"
  }else if (cancerousGene & spanning & targetable) {
    tier <- "2A"
  }else if (cancerousFusion & spanning & targetable) {
    tier <- "2B"
  }else if (cancerousGene & targetable) {
    tier <- "3A"
  }else if (cancerousFusion & targetable) {
    tier <- "3B"
  }else {
    tier <- "4"
  }

  x["Tier"] <- tier
  return(x)
}


input <- t(apply(input, FUN = tierize, MARGIN = 1,input))

View(input)
write.table(input, "/home/nick/Desktop/Fusion_prioritization/processed/star_fusion_analyzed.tsv", sep="\t", row.names = FALSE)
