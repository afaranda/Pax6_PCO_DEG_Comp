# Pax6_PCO_DEG_Comp
Integrated analysis comparing genes differentially expressed in five pairwise
contrasts.

## Contrasts Evaluated
  1. Lens Epithelium: Pax6 heterozygotes vs wildtype littermate control
   + Data file "Cuffdiff1_max.xlsx" sent by DNA Link on May 28, 2019
   + Data imported from sheet "gene_exp.diff", other sheets deleted
   + Deleted Multi-line header section -- first line is now column headers
   + After filtering on "status == OK" -- attribute "gene" is primary key
    
  2. Lens Epithelium: 0 hours PCS vs 6 hours PCS (DNA_Link)
   + Data file "WT0vsWT6.xlsx" sent by DNA Link on Sept 26, 2018
   + Data imported from sheet "gene_exp.diff", other sheets deleted
   + Deleted Multi-line header section -- first line is now column headers
   + After filtering on "status == OK" -- attribute "gene" is primary key
    
  3. Lens Epithelium: 0 hours PCS vs 24 hours PCS (DNA_Link)
   + Data file "WT0vsWT6.xlsx" sent by DNA Link on Sept 26, 2018
   + Data imported from sheet "gene_exp.diff", other sheets deleted
   + Deleted Multi-line header section -- first line is now column headers
   + After filtering on "status == OK" -- attribute "gene" is primary key
    
  4. Lens Epithelium: 0 hours PCS vs 24 hours PCS (DBI)
   + Data file: WT_0_Hour_vs_WT_48_Hour_expressedTags-all.txt
   + File is output from "Fibronectin_Analysis" run on June 27, 2019
   + attribute "GeneID" is primary key after duplicate filtering
    
  5. Lens Epithelium: 0 hours PCS vs 48 hours PCS (DBI) 
   + Data file: WT_0_Hour_vs_WT_48_Hour_expressedTags-all.txt
   + File is output from "Fibronectin_Analysis" run on June 27, 2019
   + attribute "GeneID" is primary key after duplicate filtering

## Notes
  + Data files from DBI (Contrasts 4 & 5) were filtered for duplicates. For
    any present gene with valid duplicate measurements, the record with the
    maximum absolute fold change was selected as representative. In both 
    contrasts, the gene Epb4 was the only one with duplicates.  None of the
    duplicate entries indicated fold change conflicts, and the entry with
    the maximum absolute fold change was statistically significant
