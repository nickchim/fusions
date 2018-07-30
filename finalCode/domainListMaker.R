library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)

#colOI stands for column of interest
addDomains <- function(domList,input,colOI) {
  output <- domList
  PF_dom <- input[str_detect(input[,colOI], "PF") & !is.na(input[,colOI]),colOI]
  IPR_dom <- input[str_detect(input[,colOI], "IPR") & !is.na(input[,colOI]),colOI]
  if (nrow(PF_dom) > nrow(IPR_dom)) {
    for (i in 1:nrow(PF_dom)) {
      dom <- as.character(PF_dom[i,1])
    
      if (!(dom %in% output$PF)) {

        output[nrow(output)+1,"PF"] <- dom
        output[nrow(output),"IPR"] <- NA
      }
    }
  }else {
    for (i in 1:nrow(IPR_dom)) {
      dom <- as.character(IPR_dom[i,1])

      
      if (!(dom %in% output$IPR)) {

        output[nrow(output)+1,"IPR"] <- dom
        output[nrow(output),"PF"] <- NA
      }
    }
  }
  print(output)
  return(output)
}

input1 <- read_excel("/home/nick/Desktop/Fusion_prioritization/Data/pfam_domains.XLS")
input2 <- read_excel("/home/nick/Desktop/Fusion_prioritization/Data/interpro_domains.XLSX")

Domain_List <- tibble(PF = NA,IPR = NA)

Domain_List <- addDomains(Domain_List,input1,"ACCESSION")
Domain_List <- Domain_List[-1,]

Domain_List <- addDomains(Domain_List,input2,"INTERPRO_ID")

write.table(Domain_List, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/Domain_List.tsv", sep="\t", row.names = FALSE)