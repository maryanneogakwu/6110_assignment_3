# 6110_assignment_3
Shotgun metagenomics Analysis of Human Gut Microbiomes in severe non-reponsive celiac disease patients (NRCD) vs asymptomatic celiac pateients (control)


## Dependencies and Packages
| Tool | Version | Source | Purpose |
|------|---------|--------|---------|
| Aspera-cli | 4.20.0 | Bioconda | High speed Raw data download |
| Fastqc | 0.12.1 | Bioconda | Raw reads QC analyser |
| Multiqc | 1.33 |Bioconda/noarch| Summarizing Fastqc reports in one document |
| DESeq2 | 1.48.2 | BioC 3.21 | Differential Gene Expression Analysis | 
| Kraken2 | 2.17.1 | Bioconda| Taxonomic classification |
| Braken | 3.1 | Bioconda | Computation of species abundance |
| Kraken-biom | Smdabdoub | Conversion of kraken files to BIOM | 1.0.1 |
| Phyloseq | BioC 3.21 | 1.52.0 | Handling high-throughput microbiome census data | 
| Ggplot2 | 4.0.2 | CRAN | Data Visualisation |
| Vegan | 2.7 - 3 | CRAN | Community Ecoclogy Package |
| Biomformat | 1.36.0 | BioC 3.21 | Interface package for BIOM file format |


## References
Bracken: Bayesian Reestimation of Abundance with Kraken (n.d.). John Hopkins University: Center for Computational Biology https://ccb.jhu.edu/software/bracken/

Marcos-Zambrano, L. J., Lacruz-Pleguezuelos, B., Aguilar-Aguilar, E., Marcos-Pasero, H., Valdés, A., Loria-Kohen, V., Cifuentes, A., De Molina, A. R., Diaz-Ruiz, A., Pancaldi, V., & De Santa Pau, E. C. (2025). Microbiome gut community structure and functionality are associated with symptom severity in non-responsive celiac disease patients undergoing a gluten-free diet. mSystems, 10(7), e0014325. https://doi.org/10.1128/msystems.00143-25

Francavilla A, Ferrero G, Pardini B, Tarallo S, Zanatto L, Caviglia GP, Sieri S,Grioni S, Francescato G, Stalla F, Guiotto C, Crocella L, Astegiano M,Bruno M, Calvo PL, Vineis P, Ribaldone DG, Naccarati A. 2023. Gluten-freediet affects fecal small non-coding RNA profiles and microbiomecomposition in celiac disease supporting a host-gut microbiotacrosstalk. Gut Microbes 15:2172955. https://doi.org/10.1080/19490976.2023.2172955

Dabdoub, SM (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software]. Available at https://github.com/smdabdoub/kraken-biom.
