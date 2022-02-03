# Constraint-based modeling identifies cell-state specific metabolic vulnerabilities during the epithelial to mesenchymal transition  

## Summary

This repository contains the code from the paper Constraint-based modeling identifies cell-state specific metabolic vulnerabilities during the epithelial to mesenchymal transition by Campit, S.E., Keshamouni, V.G., and Chandrasekaran, S. 

**Key analyses contained in notebooks:**

  1. Data preprocessing for transcriptomics, proteomics, single-cell transcriptomics, CERES Score data, and other =omics datasets.
  2. Constraint-based metabolic reconstruction and analysis code for simulating metabolic fluxes and growth resulting from gene and reaction knockout. 
  3. Statistical analyses for assessing differences between groups. 
  
## Programming languages used in this analysis

  * MATLAB version R2020b Update 4
  * R version 4.03
  * Python version 3.8.6

## Usage
Three programming languages (Python / R / MATLAB) were used, based on availability of scientific libraries and strengths in specific tasks. Thus, we would recommend the following workflow to perform the entire analysis end-to-end. We will point to specific directories and scripts that are numbered by usage.

  1. Exploratory data analysis and general understanding of data distributions: `notebooks/r/01_EDA/*.Rmd`
  2. Preprocessing bulk -omics data for COBRA: `notebooks/r/02_DifferentialExpression/*.Rmd`
  3. Preprocessing single-cell omics data for COBRA: `notebooks/r/03_Preprocess/*.Rmd`
  4. Performing MAGIC data imputation for single-cell COBRA analysis: `notebooks/python/magic.ipynb`
  5. Constraint-based reconstruction and analysis for bulk -omics data: `notebooks/matlab/01_bulk_analysis/RECON1/*.mlx`
  6. Constraint-based reconstruction and analysis for single-cell -omics data: `notebooks/matlab/02_single_cell_analysis/recon1_scCOBRA.mlx`
  7. Generating FBA-UMAP profiles: `notebooks/r/05_Embeddings/*.Rmd`
  8. Statistical analyses: Google Colab notebooks can be found [here](https://drive.google.com/drive/folders/1kCNsrULvzgaTEH3387mAx7KbB_dJSO4p).

Note that there are additional QA/QC scripts and notebooks available as well.
  
## Contributing
Contributions to make this analysis better, more robust, and easier to follow are greatly appreciated. Here are the steps we ask of you:
  1. Fork the project
  2. Create a new branch
  3. Make your changes
  4. Commit your changes
  5. Push to the branch
  6. Open a pull request

## License
Released via GPL GNU License . See `LICENSE` for more information.

&copy; 2022 The Regents of the University of Michigan
## Contact
Chandrasekaran Research Group - https://systemsbiologylab.org/

Contact: csriram@umich.edu