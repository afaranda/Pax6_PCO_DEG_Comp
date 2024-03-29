################################################################################
# File: Pax6Epi_PCO_Write_Excel.R                                              #
# Author: Adam Faranda                                                         #
# Created: June 28, 2019                                                       #
# Purpose: Write overlapping differential expression results to spreadsheets   #
#                                                                              #
################################################################################

# Setup Workspace
options(echo=T)
library('openxlsx')
library('dplyr')
wd<-getwd()
data_dir<-paste(wd, "data_files", sep="/")
out_dir<-paste(wd, "analysis_results", sep="/")
source('Pax6Epi_PCO_Compare_Overlap.R')

master<-list.files(pattern="master_tables")
if (length(master) > 0) { 
  load(master[1])
} else {
  source("Pax6Epi_PCO_Build_Master_Tables.R")
}

print(data_dir)
# Fix group variable from Pax6 Experiment
ss_master$Group_1<-gsub("PAX6LEplusplus", "plusplus", ss_master$Group_1)
ag_master$Group_1<-gsub("PAX6LEplusplus", "plusplus", ag_master$Group_1)
ss_master$Group_2<-gsub("PAX6LEplusminus", "plusminus", ss_master$Group_2)
ag_master$Group_2<-gsub("PAX6LEplusminus", "plusminus", ag_master$Group_2)

# Build annotation table
an<-ss_master %>% 
  group_by(MGI.symbol) %>% 
  filter(row_number(MGI.symbol) == 1) %>%
  select(MGI.symbol, ensembl_gene_id, description)
an<-data.frame(an, stringsAsFactors = F)

# list of filter condtions to select comparison subsets of ss_master
evaluationList<-list(
  c("DNA", "WT_6_Hour"),
  c("DNA", "WT_24_Hour"),
  c("DBI", "WT_24_Hour"),
  c("DBI", "WT_48_Hour")
)


# First Contrast is Wildtype vs Pax6 Heterozygous lens epithelium
dg1<-ss_master %>% filter(Lab == "DNA", Group_2 == "plusminus")
e<-evaluationList[1]


# Iterate through a list of second contrasts from the PCO data set
for(e in evaluationList){
  
  dg2<-ss_master %>% filter(Lab == e[1], Group_2 == e[2])
  
  # Query joins two pairwise contrasts
  allResults<-query(dg1, dg2)
  bioResults<-bioSig(allResults)
  
  # Get individual deg summary tables
  degTable.dg1<-degSummary(dg1)
  degTable.dg2<-degSummary(dg2)
  
  # Get intersection tables and directional subsets
  stat.tables<-subsetTables(
    allResults, descname = T, annot = an, unlog = T,
    Contrast_1 = "WTvP6", 
    Contrast_2= paste(
      "WT0v", gsub("_","",gsub("Hour","",e[2])), sep=""
    )
  )
  bio.tables<-subsetTables(
    bioResults, descname = T, annot = an, stat=F, unlog=T,
    Contrast_1 = "Pax6_LE", Contrast_2= paste(
      "WT0v", gsub("_","",gsub("_","",e[2])), sep=""
    )                  
  )
  # Get contingency tables for intersection
  stat.inx<-tabulateOverlap(stat.tables, rename = T)
  bio.inx <-tabulateOverlap(bio.tables, rename = T)
  

  print(degTable.dg1)
  print(degTable.dg2)

  print(stat.inx)
  print(bio.inx)
  print(
    paste(
      e[1], e[2], "All Results:", nrow(allResults),
      "Bio Results:", nrow(bioResults) 
    )
  )
  
  wb<-loadWorkbook(
    paste(data_dir,
          list.files(data_dir,pattern='Comparisons.xlsx'),
          sep="/"
    )
  )
  
  # Delete "Contrasts tab from directional subsets
  stat.tables["Contrasts"]<-NULL
  bio.tables["Contrasts"]<-NULL
  
  # Tabulate Statistically Significant, and Biologically Significant Genes from each contrast
  writeData(wb, 1, degTable.dg1, startCol=2, startRow=32,colNames=F) # Add DEG counts to main page
  writeData(wb, 1, degTable.dg2, startCol=2, startRow=37,colNames=F) # Add DEG counts to main page
  
  # Tabulate Statistically Significant Intersection between the two data sets
  writeData(wb, 1, stat.inx, startCol=3, startRow=42,colNames=F, rowNames = F)
  
  # Tabulate Biologically Significant Intersection between the two data sets
  writeData(wb, 1, bio.inx, startCol=3, startRow=46,colNames=F, rowNames = F)
  
  tx<-createStyle(numFmt="TEXT")
  
  # Add statistically significant directional subsets to workbook object
  for(i in names(stat.tables)){
    addWorksheet(wb,i)
    print(nrow(stat.tables[[i]]))
    cells<-expand.grid(
      row=nrow(stat.tables[[i]]), 
      col=grep("MGI.symbol", 
               names(stat.tables[[i]])
      )
    )
    addStyle(wb, i, rows=cells$row, cols=cells$col, style=tx)
    writeData(wb, i, stat.tables[[i]])
  }
  
  # Add biologically significant directional subsets to workbook object
  for(i in names(bio.tables)){
    addWorksheet(wb,i)
    print(nrow(bio.tables[[i]]))
    cells<-expand.grid(
      row=nrow(bio.tables[[i]]),
      col=grep("MGI.symbol",
               names(bio.tables[[i]])
      )
    )
    addStyle(wb, i, rows=cells$row, cols=cells$col, style=tx)
    writeData(wb, i, bio.tables[[i]])
  }

  # Save workbooks to file
  fn<-paste("WTvP6", paste(e[1],"_",
    "WT0v", gsub("_","",gsub("Hour","",e[2])), sep=""
  ), sep="_")
  fn<-paste(fn, ".xlsx", sep="")
  #Only uncomment the following line to re-generate spreadsheets!!!!
  #saveWorkbook(wb, paste(outdir_dir,fn, sep="/"), overwrite=T)

  # Print Venn Diagarams for
  for(i in c("All", "Up", "Dn")){
      fp<-gsub(".xlsx", paste("_",i,"_Venn.png", sep=""), fn)
      if(i == "Up"){
        sets<-list(
          g1 = (dg1 %>% filter(logFC > 0))$MGI.symbol,
          g2 = (dg2 %>% filter(logFC > 0))$MGI.symbol
        )
      } else if (i == "Dn"){
        sets<-list(
          g1 = (dg1 %>% filter(logFC < 0))$MGI.symbol,
          g2 = (dg2 %>% filter(logFC < 0))$MGI.symbol
        )
      } else {
        sets<-list(
          g1 = dg1$MGI.symbol,
          g2 = dg2$MGI.symbol
        )
      }
      names(sets)<-c(
        paste("WTvP6", i,sep="_"),
        paste(
          e[1],"_","WT0v",
          gsub("_","",gsub("Hour","",e[2])),
          "_",i,sep=""
        )
      )
      png(filename = paste(out_dir, fp, sep="/"))
      venn(sets, cexsn=1, cexil=1.2, zcolor="style")
      dev.off()
  }

  print(fn)
}
print(sessionInfo())



