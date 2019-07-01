################################################################################
# File: Pax6Epi_PCO_Compare_Overlap_Tables.R			                             #
# Author: Adam Faranda			                                     				       #
# Created: June 27, 2019						                                           #
# Purpose: Compare differential expression results from pairwise contrasts     #
#          between wildtype and Pax6 Het Lens epithelium with DE results from  #
#          Lens Epithelium at 6, 24 and 48 hours PCS.                          #
#                                                                              #
################################################################################

# Setup Workspace
library('openxlsx')
library('dplyr')
wd<-getwd()
master<-list.files(pattern="master_tables")
if (length(master) > 0) { 
  load(master[1])
} else {
  source("Pax6Epi_PCO_Build_Master_Tables.R")
}

# Function recieves a data frame and a list of name pairs, returns a data frame
# with names updated based on the list
dfRename<-function(df, name_pairs){
  if (!is.null(names(name_pairs))){
    for(np in names(name_pairs)){
      names(df)[grep(name_pairs[np], names(df))]<-np
    }
  }
  df
}

# dfSubname<-function(df, name_pairs){
#   if (!is.null(names(name_pairs))){
#     for(np in names(name_pairs)){
#       names(df)<-gsub(name_pairs[np], np, names(df))
#     }
#   }
#   df
# }

dfSubname<-function(df, name_pairs){
  if (!is.null(names(name_pairs))){
    for(np in names(name_pairs)){
      r<-paste(name_pairs[np], "$",sep="")
      names(df)<-gsub(r, np, names(df))
    }
  }
  df
}


# Validate MGI.symbol as primary key for each contrast (NOT FOR GLOBAL TABLE)

print("########## Confirm Unique Gene IDs in each Contrast (All)#########")
for( i in unique(ag_master$Group_1)){
  for(j in unique(ag_master$Group_2)){
    for(l in unique(ag_master$Lab)){
      val<-ag_master %>% filter(Group_1 ==i, Group_2==j, Lab==l)
      if(nrow(val) > 0){
        print(
          paste(i,j,
                "Rows:", nrow(val),
                "Unique_Genes:", length(unique(val$MGI.symbol))      
          )
        )
      }
    }
  }
}

# Query for Pairwise overlap between two Deglists -- assumes both tables
# have the same column headers, and that "MGI.symbol" is a unique id in both
# lists
query<-function(
  dg1=ss_master, dg2=ss_master, id_col="MGI.symbol"
){
 cols<-c("logFC", "FDR", "Group_1", "Group_2", "Avg1", "Avg2")
 dg1<-dg1[,c(id_col, cols)]
 dg2<-dg2[,c(id_col, cols)]
 
 for(c in cols){
   names(dg1)[grep(c, names(dg1))]<-paste("dg1.",c, sep="")
   names(dg2)[grep(c, names(dg1))]<-paste("dg2.",c, sep="")
 }
 out<-inner_join(dg1, dg2, by=id_col)
 out
}


# Function takes a set of DEGs tablutes Total, Up and Down Genes
# For statistical and biological significance
degSummary<-function(
  df, lfc_min=1, fdr_max=0.05, 
  Avg1="Avg1", Avg2="Avg2",
  lfc="logFC", fdr="FDR"
){
  names(df)[grep(lfc, names(df))]<-"lfc"
  names(df)[grep(fdr, names(df))]<-"fdr"
  names(df)[grep(Avg1, names(df))]<-"Avg1"
  names(df)[grep(Avg2, names(df))]<-"Avg2"
  
  dg<-data.frame(
    criteria = c("Statistically Significant", "Biologically Significant"),
    Total = c(
      nrow(df %>% filter(abs(lfc) > lfc_min, fdr < fdr_max)), 
      nrow(
        df %>%
          filter(Avg1 > 2 | Avg2 > 2) %>%
          filter(abs(Avg1 - Avg2) > 2) %>% 
          filter(abs(lfc) > lfc_min, fdr < fdr_max)
      )
    ),
    Up = c(
      nrow(df %>% filter(lfc > lfc_min, fdr < fdr_max)), 
      nrow(
        df %>%
          filter(Avg1 > 2 | Avg2 > 2) %>%
          filter(abs(Avg1 - Avg2) > 2) %>%
          filter(lfc > lfc_min, fdr < fdr_max)
      )
    ),
    Down = c(
      nrow(df %>% filter(lfc < -lfc_min, fdr < fdr_max)), 
      nrow(
        df %>%
          filter(Avg1 > 2 | Avg2 > 2) %>%
          filter(abs(Avg1 - Avg2) > 2) %>%
          filter(lfc < -lfc_min, fdr < fdr_max)
      )
    )
  )
  dg
}

# Extract directional subsets of statistically significant genes
subsetTables<-function(
  df,                       # Data frame with a joined pair of results
  id_col="MGI.symbol",      # Unique Identifier for this gene
  Contrast_1 = "LE",        # Name of the first contrast in df
  Contrast_2 = "PCO",       # Name of the second contrast in df
  dg1="dg1",                # Prefix for contrast 1
  dg2="dg2",                # Prefix for contrast 2
  lfc="logFC",              # Column with log 2 fold change values
  fdr="FDR",                # Column with FDR values
  g1 = "Group_1",           # Column with Group_1 label
  g2 = "Group_2",           # Column with Group_2 label
  a1 = "Avg1",              # Column with average values for Group_1
  a2 = "Avg2",              # Column with average values for Group_2
  stat = T,                 # Whether to use 'Stat' or 'Bio' naming scheme
  unlog = T,                # Whether to report absolute or log2 fold changes
  descname = F,             # Use original, or descriptive attribute names
  annot = NULL,             # Optionally provide table (keyed on ID)
  dropGroup = T             # Drop or keep group label columns
){
  # Standardize column headers
  cols<-c(lfc=lfc, fdr=fdr, g1=g1, g2=g2, a1=a1, a2=a2)
  df<-dfSubname(df, cols)
  
  prf<-c(dg1=dg1, dg2=dg2)
  df<-dfSubname(df, prf)
  
  if(unlog){
    df$dg1.lfc<-ifelse(df$dg1.lfc >= 0, 2^df$dg1.lfc, -1/2^df$dg1.lfc)
    df$dg1.lfc<-ifelse(df$dg1.lfc >= 0, 2^df$dg1.lfc, -1/2^df$dg1.lfc)
  }
  if(!is.null(annot)){
    if(any(grepl(id_col, names(annot)))){
       if(length(unique(annot[,id_col])) == nrow(annot)){
        acols<-setdiff(names(annot), id_col)
        dcols<-setdiff(names(df), id_col)
        df<-left_join(df, annot, by=id_col)
        df<-data.frame(df, stringsAsFactors = F)
        df<-df[, c(id_col, acols, dcols)]
      }
    }
  }
  # Subset results
  if (stat){
    tables<-list(
      `Stat Sig Intersection`=df,
      `SS Dn C1 Up C2` = df %>% filter(dg1.lfc < 0, dg2.lfc > 0 ),
      `SS Dn C1 Dn C2` = df %>% filter(dg1.lfc < 0, dg2.lfc < 0 ),
      `SS Up C1 Up C2` = df %>% filter(dg1.lfc > 0, dg2.lfc > 0 ),
      `SS Up C1 Dn C2` = df %>% filter(dg1.lfc > 0, dg2.lfc < 0 )
    )
    
  # Use if the data tables submitted via df are biologically significant
  } else {
    tables<-list(
      `Bio Sig Intersection`=df,
      `BS Dn C1 Up C2`= df %>% filter(dg1.lfc < 0, dg2.lfc > 0 ),
      `BS Dn C1 Dn C2`= df %>% filter(dg1.lfc < 0, dg2.lfc < 0 ),
      `BS Up C1 Up C2`= df %>% filter(dg1.lfc > 0, dg2.lfc > 0 ),
      `BS Up C1 Dn C2`= df %>% filter(dg1.lfc > 0, dg2.lfc < 0 )
    )
  }
  
  # Rename Contrasts
  names(tables)<-gsub("C1", Contrast_1, names(tables))
  names(tables)<-gsub("C2", Contrast_2, names(tables))
  if(descname){
    cols<-c(
      dg1.g1 = paste(Contrast_1, df$dg1.g1[1], sep="_"),
      dg1.g2 = paste(Contrast_1, df$dg1.g2[1], sep="_"),
      dg2.g1 = paste(Contrast_2, df$dg2.g1[1], sep="_"),
      dg2.g2 = paste(Contrast_2, df$dg2.g2[1], sep="_"),
      dg1.a1 = paste(Contrast_1, df$dg1.g1[1], "Avg", sep="_"),
      dg1.a2 = paste(Contrast_1, df$dg1.g2[1], "Avg", sep="_"),
      dg2.a1 = paste(Contrast_2, df$dg2.g1[1], "Avg", sep="_"),
      dg2.a2 = paste(Contrast_2, df$dg2.g2[1], "Avg", sep="_"),
      dg1.lfc = ifelse(
        unlog, paste(Contrast_1, "Fold_Change",sep="_"),
        paste(Contrast_1, lfc,sep="_")
      ),
      dg2.lfc=ifelse(
        unlog, paste(Contrast_2, "Fold_Change",sep="_"),
        paste(Contrast_2, lfc,sep="_")
      ),
      dg1.fdr = paste(Contrast_1, fdr,sep="_"),
      dg2.fdr = paste(Contrast_2, fdr,sep="_")
    )
    
    for (t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        names(sloc)<-cols
        keep<-setdiff(
          names(tables[[t]]), 
          sloc[c(grep("\\.g1$", sloc), grep("\\.g2$", sloc))]
        )
        sloc<-sloc[-c(grep("\\.g1$", sloc), grep("\\.g2$", sloc)) ]
        print(keep)
        tables[[t]]<-tables[[t]][, keep]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        
      } else {
        sloc<-names(cols)
        names(sloc)<-cols
        tables[[t]]<-dfSubname(tables[[t]], sloc)
      }
    }
    
  } else {
    for(t in names(tables)){
      if(dropGroup){
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        
        drop<-c(
          grep("\\.g1$", names(tables[[t]])), 
          grep("\\.g2$", names(tables[[t]]))
        )
        keep<-names(tables[[t]])[-drop]
        print(keep)
        sloc<-sloc[-c(grep("g1$", sloc), grep("g2$", sloc))]
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
        
      } else {
        sloc<-names(cols)
        frp<-names(prf)
        names(sloc)<-cols
        names(frp)<-prf
        tables[[t]]<-dfSubname(tables[[t]], sloc)
        tables[[t]]<-dfSubname(tables[[t]], frp)
      }
    }
  }
  tables<-append(
    tables, list(Contrasts=c(Contrast_1=Contrast_1, Contrast_2=Contrast_2))
  )
  tables
}

# Tabulate Directional Intersections between the two data sets
# Recieves set of tables generated by "subsetTables()" returns a
# data frame with row counts for each directional intersect
tabulateOverlap<-function(tables, rename=F){
  stat.Intersect<-data.frame(
    C1_UP=c(nrow(tables[[4]]), nrow(tables[[5]])),
    C1_Down=c(nrow(tables[[2]]), nrow(tables[[3]]))
  )
  row.names(stat.Intersect)<-c("C2_Up", "C2_Down")
  if(rename){
    names(stat.Intersect)<-gsub(
      "C1", tables[['Contrasts']][1], names(stat.Intersect)
    )
    row.names(stat.Intersect)<-gsub(
      "C2", tables[['Contrasts']][2], row.names(stat.Intersect)
    )
  }
  stat.Intersect
}

# Biosig Filter -- filter a pair of deg-lists for biologically significant
# genes in both contrasts

bioSig<-function(ar){
  ar %>%
    filter(dg1.Avg1 > 2 | dg1.Avg2 > 2) %>%
    filter(abs(dg1.Avg1 - dg1.Avg2) > 2) %>%
    filter(dg2.Avg1 > 2 | dg2.Avg2 > 2) %>%
    filter(abs(dg2.Avg1 - dg2.Avg2) > 2)
}


