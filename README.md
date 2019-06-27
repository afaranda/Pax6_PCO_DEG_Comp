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
    + Data file: WT_0_Hour_vvs_WT_48_Hour_expressedTags-all.txt
    + File is output from "Fibronectin_Analysis" run on June 27, 2019
    
  5. Lens Epithelium: 0 hours PCS vs 48 hours PCS (DBI) 
    + Data file: WT_0_Hour_vvs_WT_48_Hour_expressedTags-all.txt
    + File is output from "Fibronectin_Analysis" run on June 27, 2019
    
## R Scripts

