################################################################################
# File: Pax6Epi_PCO_Build_Master_Tables.R					                             #
# Author: Adam Faranda			                                     				       #
# Created: June 27, 2019						                                           #
# Purpose: Combine differential expression data for pairwise contrasts between #
#          wildtype and Pax6 Het Lens epithelium with DE data from             #
#          Lens Epithelium at 6, 24 and 48 hours PCS.  Data from a total of    #
#          five pariwise contrasts is concatenated into a single "tidy" table  #
#                                                                              #
################################################################################

# Setup Workspace
library('openxlsx')
library('dplyr')
wd<-getwd()
data_dir<-paste(wd, "data_files", sep="/")

# Function to select one row from a set of duplicates.  For any
# duplicate symbol, the record with the greatest absolute logFC is
# retained and all others are discarded
uniqueMaxLfc<-function(df, idc="MGI.symbol", lfc="logFC", fdr="FDR"){
  names(df)[grep(idc, names(df))]<-"idc"
  names(df)[grep(lfc, names(df))]<-"lfc"
  names(df)[grep(fdr, names(df))]<-"fdr"
  
  df<-df %>% 
    group_by(idc) %>%
    top_n(1, abs(lfc))  %>% 
    group_by(idc) %>%
    filter(row_number(lfc) == 1)
  
  names(df)[grep("idc", names(df))]<-idc
  names(df)[grep("lfc", names(df))]<-lfc
  names(df)[grep("fdr", names(df))]<-fdr
  return (data.frame(df, stringsAsFactors=F))
}

# Load in data files
pax6DEG<-read.xlsx(paste(data_dir,'Cuffdiff1_max.xlsx', sep="/"))
dbi<-list()
for(f in list.files(data_dir, pattern="_expressedTags-all.txt")){
  c<-gsub("_expressedTags-all.txt","", f)
  file<-paste(data_dir,f, sep="/")
  dbi[[c]]<-read.table(
      file, sep="\t",
      header= T, quote="", stringsAsFactors = F
  )
  names(dbi[[c]])[grep('GeneID$', names(dbi[[c]]))]<-'MGI.symbol'
  dbi[[c]]$Lab<-"DBI"
  dbi[[c]]$experiment<-"PCO"
  dbi[[c]]<-uniqueMaxLfc(dbi[[c]])
}

dna<-list()
for(f in list.files(data_dir, pattern="WT0vs")){
  c<-gsub(".xlsx","", f)
  file<-paste(data_dir,f, sep="/")
  dna[[c]]<-read.xlsx(file)
  names(dna[[c]])[grep('^gene$', names(dna[[c]]))]<-'MGI.symbol'
  names(dna[[c]])[grep('q_value', names(dna[[c]]))]<-'FDR'
  names(dna[[c]])[grep('sample_1', names(dna[[c]]))]<-'Group_1'
  names(dna[[c]])[grep('sample_2', names(dna[[c]]))]<-'Group_2'
  names(dna[[c]])[grep('value_1', names(dna[[c]]))]<-'Avg1'
  names(dna[[c]])[grep('value_2', names(dna[[c]]))]<-'Avg2'
  dna[[c]]$Lab<-"DNA"
  dna[[c]]$description<-NA
  dna[[c]]$ensembl_gene_id<-NA
  dna[[c]]$experiment<-"PCO"
  dna[[c]]$logFC<-log2(dna[[c]]$Avg2/dna[[c]]$Avg1)
  dna[[c]]<-dna[[c]] %>% filter(status=="OK")
}

# Fix pax6DEG Column Headers
names(pax6DEG)[grep('^gene$', names(pax6DEG))]<-'MGI.symbol'
names(pax6DEG)[grep('q_value', names(pax6DEG))]<-'FDR'
names(pax6DEG)[grep('sample_1', names(pax6DEG))]<-'Group_1'
names(pax6DEG)[grep('sample_2', names(pax6DEG))]<-'Group_2'
names(pax6DEG)[grep('value_1', names(pax6DEG))]<-'Avg1'
names(pax6DEG)[grep('value_2', names(pax6DEG))]<-'Avg2'
pax6DEG$Lab<-"DNA"
pax6DEG$description<-NA
pax6DEG$ensembl_gene_id<-NA
pax6DEG$experiment<-'Pax6Epi'
pax6DEG$logFC<-log2(pax6DEG$Avg2/pax6DEG$Avg1)
pax6DEG<-pax6DEG %>% filter(status=="OK")

# Validate MGI.Symbol as primary key in each data set
x<-list(pax6DEG, dbi[[1]], dbi[[2]], dna[[1]], dna[[2]])
for (i in 1:5){
  print(
    paste("Rows in ",i, ":", nrow(x[[i]]), 
          "Unique MGI.Symbol: ", length(unique(x[[i]]$MGI.symbol)))
  )
}

#
ag_master<-data.frame(stringsAsFactors = F)
ss_master<-data.frame(stringsAsFactors = F)

for (d in names(dbi)){
  contrast<-d
  
  # Load in DEGs
  dg<-dbi[[d]]
  dg$Lab<-"DBI"
  samples<-names(dg)[grep("_GeneCount.cpm", names(dg))]
  
  print(paste("############ FILTER DEG: ",contrast, "##########"))	
  
  # Setup for RPKM Filtering Criteria
  samples<-gsub("_GeneCount.cpm", "_RPKM", samples)
  print(paste("RPKM Filtering Cols:", paste(samples, collapse=" "), sep=" "))
  dg$Avg1<-apply(dg[,samples[1:3]], 1, mean)
  dg$Avg2<-apply(dg[,samples[4:6]], 1, mean)
  
  dg$Group_1<-unlist(strsplit(contrast,"_vs_"))[1]
  dg$Group_2<-unlist(strsplit(contrast,"_vs_"))[2]
  
  # Get Statistically Significant Genes
  ds<-dg %>% filter(abs(logFC) > 1 & FDR < 0.05) 

  # Reorder Columns For readability, drop sample specific columns
  cols<-c("MGI.symbol",
          "ensembl_gene_id", 
          "description", "Lab", "experiment",
          "logFC", "FDR", "Group_1", 
          "Group_2", "Avg1", "Avg2"
  )
  ag_master<-bind_rows(
    ag_master,
    dg<-dg[,cols]
  )
  
  ss_master<-bind_rows(
    ss_master,
    ds<-ds[,cols]
  )
}


for (d in names(dna)){
  contrast<-d
  
  # Load in DEGs
  dg<-dna[[d]]
  
  print(paste("############ FILTER DEG: ",contrast, "##########"))	
  
  # Fix Condition labels to match DBI labels
  dg$Group_1<-gsub("W_h0", "WT_0_Hour", dg$Group_1)
  dg$Group_2<-gsub("W_h6", "WT_6_Hour",  
                   gsub("W_h24", "WT_24_Hour", dg$Group_2)
              )
  
  # Get Statistically Significant Genes
  ds<-dg %>% filter(abs(logFC) > 1 & FDR < 0.05) 
  
  # Reorder Columns For readability, drop sample specific columns
  cols<-c("MGI.symbol",
          "ensembl_gene_id", 
          "description", "Lab", "experiment",
          "logFC", "FDR", "Group_1", 
          "Group_2", "Avg1", "Avg2"
  )
  ag_master<-bind_rows(
    ag_master,
    dg<-dg[,cols]
  )
  
  ss_master<-bind_rows(
    ss_master,
    ds<-ds[,cols]
  )
  
}




ds<-pax6DEG %>% filter(abs(logFC) > 1 & FDR < 0.05) 

# Reorder Columns For readability, drop sample specific columns
cols<-c("MGI.symbol",
        "ensembl_gene_id", 
        "description", "Lab", # "experiment",
        "logFC", "FDR", "Group_1", 
        "Group_2", "Avg1", "Avg2"
)
ag_master<-bind_rows(
  ag_master,
  pax6DEG<-pax6DEG[,cols]
)

ss_master<-bind_rows(
  ss_master,
  ds<-ds[,cols]
)

print(
  paste(
    "Row Count Verification (expect 0):",
    nrow(ag_master) - nrow(pax6DEG) - nrow(dbi[[1]]) -nrow(dbi[[2]]) -nrow(dna[[1]]) -nrow(dna[[2]])
  )
)
save(ag_master, ss_master, dbi, dna, pax6DEG, file = "master_tables.Rdata")