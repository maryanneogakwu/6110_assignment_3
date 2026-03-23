# 6110_assignment_3
Shotgun metagenomics Analysis of Human Gut Microbiomes in severe non-reponsive celiac disease patients (NRCD) vs asymptomatic celiac pateients (control)

----

## Introduction
Celiac disease is an autoimmune disease, often characterized by the production of autoantibodies linked to celiac disease (CD) under inflammatory conditions, triggered by gluten ingestion in genetically predisposed individuals [1]. The disease manifests in varying degrees of inflammation in the small intestine, leading to  gastrointestinal and axra-intestinal symptoms [4]. Gluten is the protein found in Wheat, Barley, Rye, Spelt and Kamut, it serves the function of creating the structure and texture of baked goods, as it forms a network of strands that trap air bubbles, allowing dough to rise [2]. It has been theorized that he production of new grain variants due to technological rather than nutritional reasons may have influenced the observed increase in the number of CD diagnoses in recent years [3]. Non-responsive Celiac Disease (NRCD) occurs when celiac patients that follow a strict gluten free diet continue to experience and exhibit symptoms linked to CD [5]. NRCD can be categorised into groups;  
Primary: there is an initial lack of response to a gluten free diet and inflammation continues,   Secondary: symptoms are consistent with CD reoccur following an initial period of normalization while maintaining a gluten free diet [6].  
Shotgun metagenomics was carried out in the study because it provides species‑level and functional resolution that 16S rRNA sequencing simply cannot achieve, which is essential for understanding the microbiome disruptions seen in non‑responsive celiac disease. By sequencing all microbial DNA in the samples, the researchers could precisely identify which bacterial species were depleted or enriched, quantify shifts in metabolic pathways. Using this method allowed the study to investigate broad taxonomic patterns and uncover the specific microbial and functional signatures that distinguish high‑NRCD from treated celiac controls, strengthening the evidence that NRCD is linked to a dysbiotic, metabolically impaired gut ecosystem.

---
## Methods
### 1. Data Acquistion
Data from the Celiac study by Zambrano et al. (2025) [4], under the identifier PRJEB65879, contained 39 NRCD patients, of which 14 were classified as "high-NRCD" within cluster 2, the control group used in the study was obtained from publicly available data of CD patients in a study by  Francavilla et al. (2023), which contained 48 tCD control samples, under the identifier PRJNA904924.[7].Three samples from each group were randomly selected and high-NRCD and tCD samples were downloaded from the ENA browser.  
```
#Downloading the seleted high-NRCD samples from ENA
for run in ERR12025082 ERR12025102 ERR12025074; do
    echo "Getting URLs for $run..."
    urls=$(wget -qO- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp&format=tsv" \
           | tail -1 | cut -f2 | tr ';' '\n')
    
    # Download each file
    for url in $urls; do
        wget -c -P ~assignment_03/metagenomics/high_nrcd "ftp://${url}"
    done
    echo "Done $run"
done
```
```
#Downloading the tCD controls from the ENA database browser
for run in SRR22402259 SRR22402265 SRR22402328; do
    echo "$(date): Downloading $run..."
    urls=$(wget -qO- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp&format=tsv" \
           | tail -1 | cut -f2 | tr ';' '\n')
    
    for url in $urls; do
        wget -c "ftp://${url}"
    done
    echo "$(date): Done $run"
done
```
#### Table ofthe selected high-NRCD and tCD samples

| Accension | Sample Alias | Cluster |  
|-----------|--------------|---------|
| ERR12025082 |	1807 |	high-NRCD |
| ERR12025102	| 1837 | high-NRCD |
| ERR12025074 |	1797 | high-NRCD |   
| SRR22402259	| Cii071 | tCD_control |
| SRR22402265 |	Cii045 | tCD_control |
| SRR22402328	| CMa003 |tCD_control |


### 2. Quality Control
The analysis involved pairwise comparison of the 2 groups. Poor QC can lead to biological conclusions derived from technical artifacts. Trimming the low quality reads can lead to inflated species counts and false diversity measures, and adapter contamination can cause false positives as taxa will be incorrectly placed between groups. 
```
fastqc ~/binf_6110/assignment_03/metagenomics/high_nrcd/*.fastq.gz \
       ~/binf_6110/assignment_03/metagenomics/tcd/*.fastq.gz \
       --outdir ~/binf_6110/assignment_03/metagenomics/qc/raw \
       --threads 4
```
Based on the fastqc report, there is no need to trim the sequences. The samples overall have good quality scores, no adapter contamination and no N content issues. The samples all failed on the per base sequence content , likely because the first few bases of Illumina reads tend to show biased composition due to random hexamer priming during library preparation, this won’t have any effect on the Kraken2 classification and diversity analyses. Sample ERR12025102 is flagged for having a duplication percentage of 66% compared to others with 75-80%. This indicates that more reads are unique, showing higher diversity. Overall the tCD control samples have more reads compared to the high-NRCD samples, likely because the reads are from different studies and different Illumina runs. Deeper sequencing detects rarer species, which can inflate alpha diversity in the tCD group independently. This will be c=accounted for by subsampling to equal depth usuing rarefaction before running the diversity analysis.

### 

---
## Results


----
## Discussion

---

## Conclusion

---
## Dependencies and Packages
| Tool | Version | Source | Purpose |
|------|---------|--------|---------|
| Aspera-cli | 4.20.0 | Bioconda | High speed Raw data download |
| Fastqc | 0.12.1 | Bioconda | Raw reads QC analyser |
| Multiqc | 1.33 |Bioconda/noarch| Summarizing Fastqc reports in one document |
| DESeq2 | 1.48.2 | BioC 3.21 | Differential Gene Expression Analysis | 
| Kraken2 | 2.17.1 | Bioconda| Taxonomic classification |
| Braken | 3.1 | Bioconda | Computation of species abundance |
| Kraken-biom | 1.0.1| Smdabdoub | Conversion of kraken files to BIOM |
| Phyloseq | 1.52.0 | BioC 3.21 | Handling high-throughput microbiome census data | 
| Ggplot2 | 4.0.2 | CRAN | Data Visualisation |
| Vegan | 2.7 - 3 | CRAN | Community Ecoclogy Package |
| Biomformat | 1.36.0 | BioC 3.21 | Interface package for BIOM file format |
| R | 4.5.1 | CRAN | Statistical analysis |
| SRA Toolkit | 3.3.0 | NCBI | Raw data download |


## References
- [1] Caio, G., Volta, U., Sapone, A., Leffler, D. A., De Giorgio, R., Catassi, C., & Fasano, A. (2019). Celiac disease: a comprehensive current review. BMC medicine, 17(1), 142. https://doi.org/10.1186/s12916-019-1380-z

- [2] Martinez, J. (2025, October 6). Unraveling the role of gluten in baked Goods: Structure, function, and Alternatives - TheKitchenPursuits. TheKitchenPursuits. https://thekitchenpursuits.com/does-gluten-give-baked-goods-structure/

 - [3] de Lorgeril M, Salen P. Gluten and wheat intolerance today: are modern wheat strains involved? Int J Food Sci Nutr. 2014;65:577–581. doi: 10.3109/09637486.2014.886185.

- [4] Marcos-Zambrano, L. J., Lacruz-Pleguezuelos, B., Aguilar-Aguilar, E., Marcos-Pasero, H., Valdés, A., Loria-Kohen, V., Cifuentes, A., De Molina, A. R., Diaz-Ruiz, A., Pancaldi, V., & De Santa Pau, E. C. (2025). Microbiome gut community structure and functionality are associated with symptom severity in non-responsive celiac disease patients undergoing a gluten-free diet. mSystems, 10(7), e0014325. https://doi.org/10.1128/msystems.00143-25

- [5] Leonard MM, Cureton P, Fasano A. (2017). Indications and use of the gluten contamination elimination diet for patients with non-responsive celiac disease. Nutrients 9:1129. https://doi.org/10.3390/nu9101129

- [6] Leffler DA, Dennis M, Hyett B, Kelly E, Schuppan D, Kelly CP. (2007). Etiologies and predictors of diagnosis in nonresponsive celiac disease. Clin Gastroenterol Hepatol 5:445–450. https://doi.org/10.1016/j.cgh.2006.12.006

- [7] Francavilla A, Ferrero G, Pardini B, Tarallo S, Zanatto L, Caviglia GP, Sieri S,Grioni S, Francescato G, Stalla F, Guiotto C, Crocella L, Astegiano M,Bruno M, Calvo PL, Vineis P, Ribaldone DG, Naccarati A. (2023). Gluten-freediet affects fecal small non-coding RNA profiles and microbiomecomposition in celiac disease supporting a host-gut microbiotacrosstalk. Gut Microbes 15:2172955. https://doi.org/10.1080/19490976.2023.2172955

- [8] Laurichi. (n.d.). GitHub - laurichi13/NRCD_analysis. GitHub. https://github.com/laurichi13/NRCD_analysis/tree/main


- Bracken: Bayesian Reestimation of Abundance with Kraken (n.d.). John Hopkins University: Center for Computational Biology https://ccb.jhu.edu/software/bracken/




- Dabdoub, SM (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software]. Available at https://github.com/smdabdoub/kraken-biom.
