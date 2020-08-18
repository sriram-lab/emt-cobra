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
  
## Installation
MATLAB, R, and Python scripts are available in this repository. The R programming language was used mainly to perform bioinformatics analysis. MATLAB was used for constraint-based metabolic modeling. Python was used to build the interactive explorer. 
  
### Dependencies
There are several software dependencies required to run different analyses for each language. A shortlist of each dependency is located in each subdirectory as the file name `requirements.txt`. 
  
## Usage
Notebooks are available that walk through each step of the analyses in more detail, and the README documents the sequence to run each file, as well as a description of what each livescript does. An object-oriented implementation written in Python is also available, if this analysis needs to be integrated into a single pipeline.
  
## Contributing
Contributions are always welcome! If there are any bugs/questions, please raise them in the `Issues` tab. For adding new documentation, feature enhancements, and other new additions, feel free to submit a pull request, and our team will get to it within 1-2 business days.
  
## License
This work is submitted under the MIT license
