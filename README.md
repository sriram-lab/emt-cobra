# Modeling the metabolic changes of the epithelial-to-mesenchymal transition 
This repo contains code that analyzes the metabolic changes occuring in the epithelial-to-mesenchymal transition (EMT). There are three main operations performed in this analysis:

  1. Meta-analysis of multiple -omics modalities (Bulk RNA-Seq, single-cell RNA-Seq, Proteomics)
  2. Differential gene expression analysis across several datasets
  3. Constraint-based metabolic modeling for hypothesis generation  
  
**Performed so far:**
  * Differentially expressed genes and proteins for 2 bulk RNASeq datasets, 2 proteomics datasets, and 1 single-cell RNASeq dataset.
  * Parameter sensitivity analysis for these different models under 5 different parameters
  * Differential reaction sensitivity analysis (DRSA) for these datasets

  
**To do:**
  - [ ] Create a comprehensive EMT database
  - [ ] Create a pipeline that can analyzes and perform differential expression analysis / enrichment for metabolic genes
  - [ ] Tie with either literature curation and/or pathway level enrichment
  
