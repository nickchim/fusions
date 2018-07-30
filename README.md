# Fusion_prioritization


Pipeline to prioritize gene fusions coming from the STAR-Fusion workflow.


Processed Data:
CancerGeneList.tsv - Compilation of cancerous genes and a record of the studies which deemed them cancerous.
Data from:
Cosmic - https://cancer.sanger.ac.uk/census#cl_search
DgiDB -  http://www.dgidb.org/downloads (http://www.dgidb.org/data/genes.tsv)
cBio - https://pedcbioportal.org/cancer_gene_list.jsp
TARGET - https://www.nature.com/articles/nature25795#supplementary-information
Pfister - https://www.nature.com/articles/nature25480#supplementary-information
Comprehensive Characterization of Cancer Driver Genes and Mutations - https://www.cell.com/cell/fulltext/S0092-8674(18)30237-X


Domain_List.tsv - Compilation of cancerous domains.
Taken from:
http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002797#pcbi.1002797.s009
http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0004805

FusionList2.txt - Compilation of canerous fusions.
Sources:
CHIMERDB30 - https://academic.oup.com/nar/article/45/D1/D784/2605708
Jackson Lab - https://www.jax.org/clinical-genomics/clinical-offerings/fusionseq 
Cosmic - https://cancer.sanger.ac.uk/cosmic/fusion
FusionCancer - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4517624/ 
ChiTaRS - http://chitars.bioinfo.cnio.es/
TICdb - http://www.unav.es/genetica/TICdb/
ConjoinG - https://metasystems.riken.jp/conjoing/

Target_List.tsv - Compilation of targetable genes and the applicable medication based on the DgiDB database
Source: http://www.dgidb.org/downloads

star_fusion_combo.tsv - The output of the combination of data sets from Fusion and Star Catchers. Outputted from the R script "star_catcher_combiner_small.R"/"star_catcher_combiner_big.R"
star_fusion_analyzed - the output of the anaylisis of star_fusion_combo.tsv. Outputted from the R script star_catcher_analyzer3.R


Code:
CancerGeneCompilation.R - Composes CancerGeneList.txt

combiner1.R - Takes in two files, one from Star, one from Fusion Catcher, and combines them into one coherent file with standarized data. 

analyzer3.R - Takes in output file of star_catcher_combiner, analyzes and outputs a new file with several more calculated variables.

domainListMaker.R - Composes Domain_List.tsv

fusionCollecter2.R - Composes FusionList2.txt

Target_Gene_Maker.R - Composes Target_List.tsv

getPFAMDomain.R - Used in anaylyzer3.R for domain anaylisis
