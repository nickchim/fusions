library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)

#From DGIDB
cats <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/categories.tsv")
ints <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/interactions.tsv")

cats_genes <- levels(factor(cats[cats$category=="DRUGGABLE GENOME","entrez_gene_symbol"]$entrez_gene_symbol))
ints_genes <- levels(factor(ints$gene_name))

combo_genes <- levels(factor(cats_genes,ints_genes))

target_genes <- tibble(Trgt_Genes=combo_genes,Drug_name=NA,Drug_chembl_id=NA)

for (i in 1:nrow(target_genes)) {
  val <- target_genes$Trgt_Genes[i]
  if (val %in% ints_genes) {
    list <- ints$gene_name==val
    list[is.na(list)] <- FALSE
    
    list1 <- ints[list,"drug_name"]$drug_name
    list1 <- list1[!is.na(list1)]
    
    list2 <- ints[list,"drug_chembl_id"]$drug_chembl_id
    list2 <- list2[!is.na(list2)]
    
    if (length(list1) > 0) {
      target_genes$Drug_name[i] <- paste(levels(factor(list1)),collapse=",")
    }
    
    if (length(list2) > 0) {
      target_genes$Drug_chembl_id[i] <- paste(levels(factor(list2)),collapse=",")
    }
    
  }
  
  if (i %% 500 == 0) {
    print(paste(i,"/",nrow(target_genes)))
  }
}

write.table(target_genes, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/Target_List.tsv", sep="\t", row.names = FALSE)
