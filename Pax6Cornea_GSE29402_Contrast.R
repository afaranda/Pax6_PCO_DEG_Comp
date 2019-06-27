################################################################################
# File: Pax6Epi_PCO_BuildMaster_Tables.R					                             #
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



joinDegLists<-function(df1, df2){
  allResults<-as.data.frame(
    inner_join(
      gseDEG %>%
        filter(MGI.symbol %in% intersect(gseDEG$MGI.symbol, pax6DEG$MGI.symbol)) %>%
        dplyr::select(
          Human_Gene=Gene.Symbol,
          Human_Ensembl_ID=Gene.stable.ID,
          Human_Description=Gene.description,
          GSE29402_Cornea_Avg_Intensity=Cornea_Mean,
          GSE29402_Conjuncitva_Avg_Intensity=Conjunctiva_Mean,
          GSE29402_logFC=logFC,
          GSE29402_FDR=FDR,
          Mouse_Gene=MGI.symbol,
          Mouse_Ensembl_ID=Gene.stable.ID.1,
          Mouse_Description=MGI.description
        ),
      pax6DEG %>%
        filter(MGI.symbol %in% intersect(gseDEG$MGI.symbol, pax6DEG$MGI.symbol)) %>%
        dplyr::select(
          Mouse_Gene=MGI.symbol,
          Pax6_Heterozygote_Avg_FPKM=value_1,
          Pax6_WT_Littermate_Avg_FPKM=value_2,
          Pax6Cornea_logFC=logFC, 
          Pax6Cornea_FDR=FDR
        ),
      by='Mouse_Gene'
      
    ),
    stringsAsFactors=F
  )
}






allResults<-as.data.frame(
  inner_join(
    gseDEG %>%
      filter(MGI.symbol %in% intersect(gseDEG$MGI.symbol, pax6DEG$MGI.symbol)) %>%
      dplyr::select(
        Human_Gene=Gene.Symbol,
        Human_Ensembl_ID=Gene.stable.ID,
        Human_Description=Gene.description,
        GSE29402_Cornea_Avg_Intensity=Cornea_Mean,
        GSE29402_Conjuncitva_Avg_Intensity=Conjunctiva_Mean,
        GSE29402_logFC=logFC,
        GSE29402_FDR=FDR,
        Mouse_Gene=MGI.symbol,
        Mouse_Ensembl_ID=Gene.stable.ID.1,
        Mouse_Description=MGI.description
      ),
    pax6DEG %>%
      filter(MGI.symbol %in% intersect(gseDEG$MGI.symbol, pax6DEG$MGI.symbol)) %>%
      dplyr::select(
        Mouse_Gene=MGI.symbol,
        Pax6_Heterozygote_Avg_FPKM=value_1,
        Pax6_WT_Littermate_Avg_FPKM=value_2,
        Pax6Cornea_logFC=logFC, 
        Pax6Cornea_FDR=FDR
      ),
    by='Mouse_Gene'
    
  ),
  stringsAsFactors=F
)

bioResults<-allResults %>%
  filter(Pax6_Heterozygote_Avg_FPKM > 2 | Pax6_WT_Littermate_Avg_FPKM > 2) %>%
  filter(abs(Pax6_Heterozygote_Avg_FPKM - Pax6_WT_Littermate_Avg_FPKM) > 2) %>%
  filter(
    GSE29402_Cornea_Avg_Intensity > minExp  |
    GSE29402_Conjuncitva_Avg_Intensity > minExp
  )

  
# Build Excel Workbook
wb<-loadWorkbook(
        list.files( pattern='Contrast_Description.xlsx')
  )



tables<-list(
  `Stat Sig Intersection`=allResults,
  `Bio Sig Intersection`=bioResults,
  `SS Down Cornea Up P6`=allResults %>% filter(GSE29402_logFC < 0, Pax6Cornea_logFC > 0 ),
  `SS Down Cornea Down P6`=allResults %>% filter(GSE29402_logFC < 0, Pax6Cornea_logFC < 0 ),
  `SS Up Cornea Up P6`=allResults %>% filter(GSE29402_logFC > 0, Pax6Cornea_logFC > 0 ),
  `SS Up Cornea Down P6`=allResults %>% filter(GSE29402_logFC > 0, Pax6Cornea_logFC < 0 ),
  `BS Down Cornea Up P6`=bioResults %>% filter(GSE29402_logFC < 0, Pax6Cornea_logFC > 0 ),
  `BS Down Cornea Down P6`=bioResults %>% filter(GSE29402_logFC < 0, Pax6Cornea_logFC < 0 ),
  `BS Up Cornea Up P6`=bioResults %>% filter(GSE29402_logFC > 0, Pax6Cornea_logFC > 0 ),
  `BS Up Cornea Down P6`=bioResults %>% filter(GSE29402_logFC > 0, Pax6Cornea_logFC < 0 )
)


# Tabulate Statistically Significant, and Biologically Significant Genes from
# GSE29402
degTable.GSE29402<-data.frame(
  criteria = c("Statistically Significant", "Biologically Significant"),
  Total = c(
      nrow(gseDEG), nrow(gseDEG %>% filter(
        Cornea_Mean > minExp  |
        Conjunctiva_Mean > minExp
      )
    )
  ),
  Up =  c(
      nrow(gseDEG %>% filter(logFC > 0)), 
      nrow(gseDEG %>% filter(
        Cornea_Mean > minExp  |
        Conjunctiva_Mean > minExp
      ) %>% filter(logFC > 0)
    )
  ),
  Down =  c(
    nrow(gseDEG %>% filter(logFC < 0)), 
    nrow(gseDEG %>% filter(
        Conjunctiva_Mean> minExp  |
        Cornea_Mean > minExp
      ) %>% filter(logFC < 0)
    )
  )
)
writeData(wb, 1, degTable.GSE29402, startCol=2, startRow=32,colNames=F) # Add DEG counts to main page

# Tabulate Statistically Significant, and Biologically Significant Genes from
# the pax6 Experiment
degTable.pax6<-data.frame(
  criteria = c("Statistically Significant", "Biologically Significant"),
  Total = c(
    nrow(pax6DEG), 
    nrow(pax6DEG %>% 
           filter(value_1 > 2  | value_2 > 2) %>%
           filter(abs(value_1 - value_2) > 2)
    )
  ),
  Up =  c(
    nrow(pax6DEG %>% filter(logFC > 0)), 
    nrow(pax6DEG %>% 
           filter(value_1 > 2  | value_2 > 2) %>%
           filter(abs(value_1 - value_2) > 2) %>%
           filter(logFC > 0)
    )
  ),
  Down =  c(
    nrow(pax6DEG %>% filter(logFC < 0)), 
    nrow(pax6DEG %>% 
           filter(value_1 > 2  | value_2 > 2) %>%
           filter(abs(value_1 - value_2) > 2) %>%
           filter(logFC < 0)
    )
  )
)
writeData(wb, 1, degTable.pax6, startCol=2, startRow=37,colNames=F) # Add DEG counts to main page


# Tabulate Statistically Significant Intersection between the two data sets
stat.Intersect<-data.frame(
  U<-c(nrow(tables[[5]]), nrow(tables[[3]])),
  D<-c(nrow(tables[[6]]), nrow(tables[[4]]))
)
writeData(wb, 1, stat.Intersect, startCol=3, startRow=42,colNames=F)

# Tabulate Biologically Significant Intersection between the two data sets
biol.Intersect<-data.frame(
  U<-c(nrow(tables[[9]]), nrow(tables[[7]])),
  D<-c(nrow(tables[[10]]), nrow(tables[[8]]))
)
writeData(wb, 1, biol.Intersect, startCol=3, startRow=46,colNames=F)

tx<-createStyle(numFmt="@")
for(i in names(tables)){
  addWorksheet(wb,i)
  print(nrow(tables[[i]]))
  cells<-expand.grid(row=nrow(tables[[i]]), col=grep("_Gene", names(tables[[i]])))
  addStyle(wb, i, rows=cells$row, cols=cells$col, style=tx)
  writeData(wb, i, tables[[i]])
}
saveWorkbook(wb, paste("Pax6Cornea_GSE29402_Contrast_Results.xlsx", sep="_"), overwrite=T)

print(sessionInfo())


### NewBorn GEO Contrast 


