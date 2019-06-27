################################################################################
# File: Pax6Cornea_GSE29402_Contrast.R					                               #
# Author: Adam Faranda			                                     				       #
# Created: May 17, 2019							                                           #
# Purpose: Compare Corneal Genes from GSE29402 to Genes that are differentially#
#          regulated in a pax6 heterozygous mutant mouse corneas               #
################################################################################

# Setup Workspace
library('openxlsx')
wd<-getwd()
source('Process_GSE29402.R')
gseDEG<-deg
pax6DEG<-read.xlsx('Pax6_Cornea_DEG.xlsx')

# filter pax6DEG for significant genes
names(pax6DEG)[grep('^gene$', names(pax6DEG))]<-'MGI.symbol'
names(pax6DEG)[grep('q_value', names(pax6DEG))]<-'FDR'
pax6DEG$experiment<-'Pax6Cornea'
pax6DEG$logFC<-log2(pax6DEG$value_1/pax6DEG$value_2)
pax6DEG<-pax6DEG %>% filter(abs(logFC) > 1, FDR < 0.05)

# filter DEG from GSE29402 -- this is a very crude filter.
minExp<-log(0.01*(2^max(gseDEG$AveExpr)),2)
names(gseDEG)[grep('adj.P.Val', names(gseDEG))]<-'FDR'
# names(gseDEG)[grep('Conjunctiva_Mean', names(gseDEG))]<-'value_1'
# names(gseDEG)[grep('Cornea_Mean', names(gseDEG))]<-'value_2'
gseDEG<- gseDEG %>%
	 filter(abs(logFC)> 1, FDR < 0.05)
	 	
gseDEG<-gseDEG %>%
	       group_by(Gene.Symbol) %>%
	       filter(abs(logFC) == max(abs(logFC)))
gseDEG$experiment<-'GSE29402'


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




