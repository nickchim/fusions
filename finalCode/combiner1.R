library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tibble)

#Input from Star Catcher must be named "star_data.tsv" and that from Fusion Catcher "fusion_catcher_data.tsv"
star_input <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/star_data.tsv")
catcher_input <- read_tsv("/home/nick/Desktop/Fusion_prioritization/Data/fusion_catcher_data.tsv")


#Cleaning star BreakPoints
r_raw <- unlist(str_split(star_input$RightBreakpoint,"chr"))
l_raw <- unlist(str_split(star_input$LeftBreakpoint,"chr"))
rowCounter<-1
for (i in 1:length(r_raw)) {
    if (i %% 2 == 0) {
      star_input[rowCounter,"RightBreakpoint"] = r_raw[i]
      rowCounter<-rowCounter+1
    }
}
rowCounter<-1
for (i in 1:length(l_raw)) {
  if (i %% 2 == 0) {
    star_input[rowCounter,"LeftBreakpoint"] = l_raw[i]
    rowCounter<-rowCounter+1
  }
}

#Cleaning Catcher Fusion_Type
rowCounter<-1
for (shift in catcher_input$Predicted_effect) {
  if (shift == "in-frame") {
    catcher_input[rowCounter,"Predicted_effect"] <- "INFRAME"
  }else if(shift == "out-of-frame") {
    catcher_input[rowCounter,"Predicted_effect"] <- "FRAMESHIFT"
  }else {
    catcher_input[rowCounter,"Predicted_effect"] <- toupper(shift)
  }
  rowCounter<-rowCounter+1
}

star_fusions_factor <- factor(star_input$`#FusionName`)
catcher_fusions_factor <- factor(paste(catcher_input$`Gene_1_symbol(5end_fusion_partner)`,catcher_input$`Gene_2_symbol(3end_fusion_partner)`,sep="--"))


#output <- tibble(Fusion_Name=levels(factor(c(levels(star_fusions_factor),levels(catcher_fusions_factor)))))

#sigColNames_star <- c("LeftBreakpoint","RightBreakpoint","JunctionReadCount","SpanningFragCount","FUSION_CDS","FUSION_TRANSL","PROT_FUSION_TYPE")

output <- tibble(Fusion_Name=NA)
output$H_Gene <- NA
output$T_Gene <-NA
output$Left_Breakpoint <-NA
output$Left_Breakpoint_Chr <-NA
output$Left_Breakpoint_Pos <-NA
output$Left_Breakpoint_Str <-NA
output$Right_Breakpoint<- NA
output$Right_Breakpoint_Chr<- NA
output$Right_Breakpoint_Pos<- NA
output$Right_Breakpoint_Str<- NA
output$Star_Junction_Readcount <- NA
output$Star_Spanning_Fragcount <- NA
output$FC_Common_Mapping_Reads <- NA
output$FC_Spanning_Pairs <- NA
output$FC_Spanning_Unique_Reads <- NA
output$Fusion_Transcript_Sequence <- NA
output$Fusion_Protein_Sequence <- NA
output$Fusion_Type <- NA

#The Big For Loop: Loops through each fusion and checks for repeats
for (fuse in levels(factor(c(levels(star_fusions_factor),levels(catcher_fusions_factor))))) {
  
  if (fuse %in% star_fusions_factor && fuse %in% catcher_fusions_factor) {
    
    index_star <- which(star_fusions_factor %in% fuse)
    index_catcher <- which(catcher_fusions_factor %in% fuse)
    if (length(index)>1) {
      #Checking for duplicates in Star
      sequences_star <- star_input[index_star[1:length(index_star)],"FUSION_TRANSL"]
      index_star <- index_star[!duplicated(sequences_star)]
      
      #Checking for duplicates in Catcher
      sequences_catcher <- catcher_input[index_catcher[1:length(index_catcher)],"Fusion_sequence"]
      index_catcher <- index_catcher[!duplicated(sequences_catcher)]
      for (i in index_star) {
        newRowNum <-nrow(output)+1
        output[newRowNum,"Fusion_Name"] <- fuse
        output[newRowNum,"Left_Breakpoint"] <- star_input[i,"LeftBreakpoint"]
        output[newRowNum,"Right_Breakpoint"] <- star_input[i,"RightBreakpoint"]
        output[newRowNum,"Star_Junction_Readcount"] <- star_input[i,"JunctionReadCount"]
        output[newRowNum,"Star_Spanning_Fragcount"] <- star_input[i,"SpanningFragCount"]
        output[newRowNum,"Fusion_Transcript_Sequence"] <- star_input[i,"FUSION_CDS"]
        output[newRowNum,"Fusion_Protein_Sequence"] <- star_input[i,"FUSION_TRANSL"]
        output[newRowNum,"Fusion_Type"] <- star_input[i,"PROT_FUSION_TYPE"]
        output[newRowNum,"Catcher"] <- "Star"
      }
      
      for (i in index_catcher) {
        newRowNum <-nrow(output)+1
        output[newRowNum,"Fusion_Name"] <- fuse
        output[newRowNum,"Left_Breakpoint"] <- catcher_input[i,"Fusion_point_for_gene_1(5end_fusion_partner)"]
        output[newRowNum,"Right_Breakpoint"] <- catcher_input[i,"Fusion_point_for_gene_2(3end_fusion_partner)"]
        output[newRowNum,"FC_Common_Mapping_Reads"] <- catcher_input[i,"Counts_of_common_mapping_reads"]
        output[newRowNum,"FC_Spanning_Pairs"] <- catcher_input[i,"Spanning_pairs"]
        output[newRowNum,"FC_Spanning_Unique_Reads"] <- catcher_input[i,"Spanning_unique_reads"]
        output[newRowNum,"Fusion_Transcript_Sequence"] <- catcher_input[i,"Fusion_sequence"]
        output[newRowNum,"Fusion_Protein_Sequence"] <- catcher_input[i,"Predicted_fused_transcripts"]
        output[newRowNum,"Fusion_Type"] <- catcher_input[i,"Predicted_effect"]
        output[newRowNum,"Catcher"] <- "Fusion_Catcher"
      }
    }else {
      newRowNum <-nrow(output)+1
      output[newRowNum,"Fusion_Name"] <- fuse
      output[newRowNum,"Left_Breakpoint"] <- star_input[index[1],"LeftBreakpoint"]
      output[newRowNum,"Right_Breakpoint"] <- star_input[index[1],"RightBreakpoint"]
      output[newRowNum,"Star_Junction_Readcount"] <- star_input[index[1],"JunctionReadCount"]
      output[newRowNum,"Star_Spanning_Fragcount"] <- star_input[index[1],"SpanningFragCount"]
      output[newRowNum,"Fusion_Transcript_Sequence"] <- star_input[index[1],"FUSION_CDS"]
      output[newRowNum,"Fusion_Protein_Sequence"] <- star_input[index[1],"FUSION_TRANSL"]
      output[newRowNum,"Fusion_Type"] <- star_input[index[1],"PROT_FUSION_TYPE"]
      output[newRowNum,"Catcher"] <- "Star, Fusion_Catcher"
    }
    
  }else if (fuse %in% star_fusions_factor) {
    
    index <- which(star_fusions_factor %in% fuse)
    
    if (length(index)>1) {
      sequences <- star_input[index[1:length(index)],"FUSION_TRANSL"]
      index <- index[!duplicated(sequences)]
      for (i in index) {
        newRowNum <- nrow(output)+1
        output[newRowNum,"Fusion_Name"] <- fuse
        output[newRowNum,"Left_Breakpoint"] <- star_input[i,"LeftBreakpoint"]
        output[newRowNum,"Right_Breakpoint"] <- star_input[i,"RightBreakpoint"]
        output[newRowNum,"Star_Junction_Readcount"] <- star_input[i,"JunctionReadCount"]
        output[newRowNum,"Star_Spanning_Fragcount"] <- star_input[i,"SpanningFragCount"]
        output[newRowNum,"Fusion_Transcript_Sequence"] <- star_input[i,"FUSION_CDS"]
        output[newRowNum,"Fusion_Protein_Sequence"] <- star_input[i,"FUSION_TRANSL"]
        output[newRowNum,"Fusion_Type"] <- star_input[i,"PROT_FUSION_TYPE"]
        output[newRowNum,"Catcher"] <- "Star"
      }
    }else {
      newRowNum <-nrow(output)+1
      output[newRowNum,"Fusion_Name"] <- fuse
      output[newRowNum,"Left_Breakpoint"] <- star_input[index[1],"LeftBreakpoint"]
      output[newRowNum,"Right_Breakpoint"] <- star_input[index[1],"RightBreakpoint"]
      output[newRowNum,"Star_Junction_Readcount"] <- star_input[index[1],"JunctionReadCount"]
      output[newRowNum,"Star_Spanning_Fragcount"] <- star_input[index[1],"SpanningFragCount"]
      output[newRowNum,"Fusion_Transcript_Sequence"] <- star_input[index[1],"FUSION_CDS"]
      output[newRowNum,"Fusion_Protein_Sequence"] <- star_input[index[1],"FUSION_TRANSL"]
      output[newRowNum,"Fusion_Type"] <- star_input[index[1],"PROT_FUSION_TYPE"]
      output[newRowNum,"Catcher"] <- "Star"
    }
  }else {
    
    index <- which(catcher_fusions_factor %in% fuse)
    if (length(index)>1) {
      sequences <- catcher_input[index[1:length(index)],"Fusion_sequence"]
      index <- index[!duplicated(sequences)]
      for (i in index) {
        newRowNum <- nrow(output)+1
        output[newRowNum,"Fusion_Name"] <- fuse
        output[newRowNum,"Left_Breakpoint"] <- catcher_input[i,"Fusion_point_for_gene_1(5end_fusion_partner)"]
        output[newRowNum,"Right_Breakpoint"] <- catcher_input[i,"Fusion_point_for_gene_2(3end_fusion_partner)"]
        output[newRowNum,"FC_Common_Mapping_Reads"] <- catcher_input[i,"Counts_of_common_mapping_reads"]
        output[newRowNum,"FC_Spanning_Pairs"] <- catcher_input[i,"Spanning_pairs"]
        output[newRowNum,"FC_Spanning_Unique_Reads"] <- catcher_input[i,"Spanning_unique_reads"]
        output[newRowNum,"Fusion_Transcript_Sequence"] <- catcher_input[i,"Fusion_sequence"]
        output[newRowNum,"Fusion_Protein_Sequence"] <- catcher_input[i,"Predicted_fused_transcripts"]
        output[newRowNum,"Fusion_Type"] <- catcher_input[i,"Predicted_effect"]
        output[newRowNum,"Catcher"] <- "Fusion_Catcher"
      }
    }else {
      newRowNum <-nrow(output)+1
      output[newRowNum,"Fusion_Name"] <- fuse
      output[newRowNum,"Left_Breakpoint"] <- catcher_input[index[1],"Fusion_point_for_gene_1(5end_fusion_partner)"]
      output[newRowNum,"Right_Breakpoint"] <- catcher_input[index[1],"Fusion_point_for_gene_2(3end_fusion_partner)"]
      output[newRowNum,"FC_Common_Mapping_Reads"] <- catcher_input[index[1],"Counts_of_common_mapping_reads"]
      output[newRowNum,"FC_Spanning_Pairs"] <- catcher_input[index[1],"Spanning_pairs"]
      output[newRowNum,"FC_Spanning_Unique_Reads"] <- catcher_input[index[1],"Spanning_unique_reads"]
      output[newRowNum,"Fusion_Transcript_Sequence"] <- catcher_input[index[1],"Fusion_sequence"]
      output[newRowNum,"Fusion_Protein_Sequence"] <- catcher_input[index[1],"Predicted_fused_transcripts"]
      output[newRowNum,"Fusion_Type"] <- catcher_input[index[1],"Predicted_effect"]
      output[newRowNum,"Catcher"] <- "Fusion_Catcher"
    }
  }
}

#Creating Left and Right Gene Vars
fusions_split <- unlist(str_split(output$Fusion_Name,"-"))
rowCounter <- 1
for (i in 1:length(fusions_split)) {
  if (i%%2!=0) {
    output[rowCounter,"H_Gene"] <- fusions_split[i]
  }else {
    output[rowCounter,"T_Gene"] <- fusions_split[i]
    rowCounter<-rowCounter+1
  }
}
output <- output[-1,]

#Cleaning up "."s and changing them to NA
dotClean <- function(x) {
  loci <- str_locate(x,"\\.")[,1] == 1
  loci[is.na(loci)] <- FALSE
  x[!is.na(str_locate(x,"\\."))[,1] & loci] <- NA
  return (x)
}
output <- apply(output,FUN=dotClean,MARGIN=2)

output<-as.tibble(output)

#Adding Breakpoint Chromosone Position and Strand Columns
LBP_List <-str_split(output$Left_Breakpoint,":")
rowCounter<-1
for (item in LBP_List) {
  output[rowCounter,"Left_Breakpoint_Chr"] <- item[1]
  output[rowCounter,"Left_Breakpoint_Pos"] <- item[2]
  output[rowCounter,"Left_Breakpoint_Str"] <- item[3]
  rowCounter<-rowCounter+1
}

RBP_List <-str_split(output$Right_Breakpoint,":")
rowCounter<-1
for (item in RBP_List) {
  output[rowCounter,"Right_Breakpoint_Chr"] <- item[1]
  output[rowCounter,"Right_Breakpoint_Pos"] <- item[2]
  output[rowCounter,"Right_Breakpoint_Str"] <- item[3]
  rowCounter<-rowCounter+1
}

output <- subset(output, select=-c(Right_Breakpoint,Left_Breakpoint))

write.table(output, "/home/nick/Desktop/Fusion_prioritization/Data/Processed/star_fusion_combo.tsv", sep="\t", row.names = FALSE)