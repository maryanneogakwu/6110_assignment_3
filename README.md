# 6110_assignment_3

**Shotgun metagenomics Analysis of Human Gut Microbiomes in severe non-reponsive celiac disease patients (NRCD) vs asymptomatic celiac pateients (control)**

Author: Maryanne Ogakwu

Dataset:   PRJEB65879 - Zambrano et al. (2025) and PRJNA904924 -Francavilla et al. (2023)

----

## Introduction
Celiac disease is an autoimmune disease, often characterized by the production of autoantibodies linked to celiac disease (CD) under inflammatory conditions, triggered by gluten ingestion in genetically predisposed individuals [1]. The disease manifests in varying degrees of inflammation in the small intestine, leading to  gastrointestinal and axra-intestinal symptoms [4]. Gluten is the protein found in Wheat, Barley, Rye, Spelt and Kamut, it serves the function of creating the structure and texture of baked goods, as it forms a network of strands that trap air bubbles, allowing dough to rise [2]. It has been theorized that he production of new grain variants due to technological rather than nutritional reasons may have influenced the observed increase in the number of CD diagnoses in recent years [3]. Non-responsive Celiac Disease (NRCD) occurs when celiac patients that follow a strict gluten free diet continue to experience and exhibit symptoms linked to CD [5]. NRCD can be categorised into groups;  
Primary: there is an initial lack of response to a gluten free diet and inflammation continues,   Secondary: symptoms are consistent with CD reoccur following an initial period of normalization while maintaining a gluten free diet [6].  
Shotgun metagenomics was carried out in the study because it provides species‑level and functional resolution that 16S rRNA sequencing simply cannot achieve, which is essential for understanding the microbiome disruptions seen in non‑responsive celiac disease. By sequencing all microbial DNA in the samples, the researchers could precisely identify which bacterial species were depleted or enriched, quantify shifts in metabolic pathways. Using this method allowed the study to investigate broad taxonomic patterns and uncover the specific microbial and functional signatures that distinguish high‑NRCD from treated celiac controls, strengthening the evidence that NRCD is linked to a dysbiotic, metabolically impaired gut ecosystem.  
This project aims to compare species‑level taxonomic composition and functional potential between high‑NRCD and treated celiac (tCD) controls, and to identify taxa and pathways associated with persistent symptoms in NRCD.

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
The analysis involved pairwise comparison of the 2 groups. Quality Control was carried out usingg Fastqc (v 0.12.1). Poor QC can lead to biological conclusions derived from technical artifacts. Trimming the low quality reads can lead to inflated species counts and false diversity measures, and adapter contamination can cause false positives as taxa will be incorrectly placed between groups. 
```
fastqc ~/binf_6110/assignment_03/metagenomics/high_nrcd/*.fastq.gz \
       ~/binf_6110/assignment_03/metagenomics/tcd/*.fastq.gz \
       --outdir ~/binf_6110/assignment_03/metagenomics/qc/raw \
       --threads 4
```
Based on the fastqc report, there is no need to trim the sequences. The samples overall have good quality scores, no adapter contamination and no N content issues. The samples all failed on the per base sequence content , likely because the first few bases of Illumina reads tend to show biased composition due to random hexamer priming during library preparation, this won’t have any effect on the Kraken2 classification and diversity analyses. Sample ERR12025102 is flagged for having a duplication percentage of 66% compared to others with 75-80%. This indicates that more reads are unique, showing higher diversity. Overall the tCD control samples have more reads compared to the high-NRCD samples, likely because the reads are from different studies and different Illumina runs. Deeper sequencing detects rarer species, which can inflate alpha diversity in the tCD group independently. This will be c=accounted for by subsampling to equal depth usuing rarefaction before running the diversity analysis.

### 3. Classification using Kraken2
Classification was carried out using Kraken2 (v 2.17.1) as opposed to MetaPhlan that was used in the original study because Kraken2 provides higher sensitivity and broader taxonomic coverage, making it better fro detecting the full microbial community in shotgun metagenomics. Kraken2 and DRAC Nibi Cluster was used to run computationally intensive Kraken2 and Bracken analyses with 16 CPU cores and 32GB RAM, significantly reducing processing time compared to local computation.
```
for run in ERR12025082 ERR12025102 ERR12025074 \
           SRR22402259 SRR22402265 SRR22402328; do
    echo "$(date): Running Kraken2 on $run..."
    if [[ $run == ERR* ]]; then
        OUTDIR=$OUT/high_nrcd
        INDIR=$READS/high_nrcd
    else
        OUTDIR=$OUT/tcd
        INDIR=$READS/tcd
    fi
    kraken2 --db $DB \
            --paired \
            --gzip-compressed \
            --threads 16 \
            --output ${OUTDIR}/${run}.kraken \
            --report ${OUTDIR}/${run}.report \
            ${INDIR}/${run}_1.fastq.gz \
            ${INDIR}/${run}_2.fastq.gz
    echo "$(date): Done $run"
done
```
| Flag              | Purpose                                                                 |
|-------------------|-------------------------------------------------------------------------|
| `--db`            | Specifies the path to your Kraken2 database                             |
| `--paired`        | Indicates that the input consists of paired-end reads                   |
| `--gzip-compressed` | Allows Kraken2 to read `.fastq.gz` files directly without unzipping    |
| `--threads 4`     | Uses 4 CPU threads to speed up classification                           |
| `--output`        | Writes the per-read classification results to a file                    |
| `--report`        | Generates a summary report used for downstream diversity analyses       |

#### Abundance reestimation using Bracken 
 Bracken (v 3.1) was used to reestimate species-level abundances from Kraken2 output by statistically redistributing reads assigned to higher taxonomic levels down to species level.
```
for run in ERR12025082 ERR12025102 ERR12025074; do
    bracken -d $DB \
            -i $OUT/high_nrcd/${run}.report \
            -o $OUT/high_nrcd/${run}.bracken \
            -r 150 -l S -t 10
done

for run in SRR22402259 SRR22402265 SRR22402328; do
    bracken -d $DB \
            -i $OUT/tcd/${run}.report \
            -o $OUT/tcd/${run}.bracken \
            -r 150 -l S -t 10
done
```

##### Data conversion using Kraken-biom
Report files, generated with Kraken were combined and converted using Kraken-biom (v 1.0.1). 
`-fmt json` — biom files come in two formats, JSON and HDF5. The biomformat R package reads JSON more reliably, so specifying this avoids compatibility issues when you import into R.
`--min S` - tells kraken-biom to assign reads at below min rank will be recorded as being assigned to the min rank level.
```
kraken-biom \
  high_nrcd/ERR12025082.report \
  high_nrcd/ERR12025102.report \
  high_nrcd/ERR12025074.report \
  tcd/SRR22402259.report \
  tcd/SRR22402265.report \
  tcd/SRR22402328.report \
  --min S \
  --fmt json \
  -o combined4.biom
```

### 4. Diversity Abundance
Alpha diversity, beta diversity, and differential abundance each reveal a different layer of biological insight in the high-NRCD vs tCD dataset. Phyloseq (1.52.0) was used to store and manage the combined taxonomic abundance data, sample metadata and taxonomy table in a single unified object for downstream diversity analyses. Biomformat (v 1.36.0) was used to import the combined BIOM file into R and convert it into a format compatible with the phyloseq package.

#### Alpha Diversity
Alpha Diversity was carried out to measure within‑sample richness and evenness, to determine whether high‑NRCD patients have a depleted or imbalanced microbiome compared to the tCD controls.
```
plot_richness(physeq_rare,
              x = "Cluster",
              color = "Cluster",
              measures = c("Observed", "Shannon", "Simpson")) +
  geom_boxplot(alpha = 0.5) +
  labs(title = "Alpha Diversity: high-NRCD vs tCD",
       x = "Group",
       y = "Diversity") +
  theme_bw()
```

#### Beta Diversity
Beta diversity was carried out to evaluate between‑sample differences, showing whether the two groups form distinct microbial community structures which is essential for demonstrating that high-NRCD is associated with a unique dysbiotic signature. The following three tests help give the complete picture on abundance differences, community structure differences, and compositional turnover.

##### Bray-Curtis PCoA
This was carried out to capture abundance‑based differences, to observe shifts in dominant taxa separate high‑NRCD from tCD controls
```
ord.pcoa.bray <- ordinate(physeq_rare, method = "PCoA", distance = "bray")
plot_ordination(physeq_rare, ord.pcoa.bray,
                color = "Cluster",
                title = "Bray-Curtis PCoA: high-NRCD vs tCD") +
  geom_point(size = 4) +
  theme_bw()
```
##### NMDS with Bray-Curtis
NMDS with Bray–Curtis was carried out to provide a non‑linear, rank‑based view of the same distances, which is useful because it can reveal group separation even when the data don’t fit linear assumptions.
```
ord.nmds.bray <- ordinate(physeq_rare, method = "NMDS", distance = "bray")
plot_ordination(physeq_rare, ord.nmds.bray,
                color = "Cluster",
                title = "Bray-Curtis NMDS: high-NRCD vs tCD") +
  geom_point(size = 4) +
  theme_bw()
```
##### Jaccard PCoA
Jaccard PCoA was carried out to investigate the presence or absence of taxa within the groups and whether the groups differ in which taxa are present at all, regardless of abundance.
```
ord.pcoa.jaccard <- ordinate(physeq_rare, method = "PCoA", distance = "jaccard")
plot_ordination(physeq_rare, ord.pcoa.jaccard,
                color = "Cluster",
                title = "Jaccard PCoA: high-NRCD vs tCD") +
  geom_point(size = 4) +
  theme_bw()
```
#### Differential Abundance

Differential abundance was carried using DESeq2 (1.48.2) to identify the specific taxa driving the differences, identifying which species or phyla are enriched or depleted in high‑NRCD. DESeq2 was used over ANCOMBC because it is better for small datasets because as it has negative‑binomial model gives stronger power to detect real differences. It handles low‑abundance taxa and uneven sequencing depth more robustly than ANCOMBC.
ANCOM‑BC is more conservative, so DESeq2 will typically identify more biologically meaningful shifts in your NRCD vs tCD samples [11]. ANCOMBC was attempted for differential abundance analysis but excluded due to a package dependency conflict with CVXR.
```
ggplot(as.data.frame(res), 
       aes(x = log2FoldChange, 
           y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), size = 2) +
  geom_vline(xintercept = 0, color = "red") +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "blue"),
                     name = "Significant (padj<0.05)") +
  labs(title = "Differential Abundance: tCD vs high-NRCD",
       x = "Log2 Fold Change",
       y = "-log10(adjusted p-value)") +
  theme_bw()
```
---
## Results
Extra plots and discussions of the plots can be found in the Results directory and in the corresponding plot folders. These were not included in the README.md file for the sake of brevity.

### Alpha Diversity
<img width="1000" height="500" alt="Alpha Diversity Plot" src="https://github.com/user-attachments/assets/a7ea9d5d-2281-49af-b990-0d4d11e960bd" />

**Figure 1.** This plot compares within-sample miv=crobial diversity between the high-NRCD patients and the tCD controls using 3 metrics: Observed richness, Shannon diversity and Simpson diversity. Across all measures, the tCD control group shows consistently higher microbial diversity counts, with tighter distributions and higher medians. The high-NRCD group displays reduced richness and eveness, indicating that individuals with Non-responsive Celiac Disease (NRCD) have a microbiome that is significantly less diverse and potentially more dysbiotic. These findings support the paper's observed results, in which high- NRCD patients exhibit depleted microbial ecosystem characterized by lower alpha diversity. These can be interpreted as a sign of imparied microbial resilience and altered gut ecology.

#### Relative Abundance Test: Top 10 Species
<img width="1000" height="700" alt="Relative Abundance species level" src="https://github.com/user-attachments/assets/99eace8e-c9b0-4450-b8b7-7882d5833464" />

**Figure 2.** This plot revels a clear compositional shift between the 2 groups, as the high-NRCD samples are shown to have a dominant subset of species,This aligns with the paper's findings that link NRCS to a dysbiotic microbiome and overrepresentation of pro-inflammatory Bacteroides species. In the high-NRCD group, Bacteroides cellulosilyticus (salmon/coral) is the dominant species across all three samples, particularly in ERR12025082 where it accounts for over 50% of the classified reads. Faecalibacterium prausnitzii (green) and Bacteroides sunii (purple) are also consistently present across the high-NRCD samples. The tCD controls show a more balanced distribution of species with no single taxon dominating to the same extent. The samples in tCD can also be observed to have higher ralative abundances of Bacteroidota  (Bacteroides dorei (olive green) and Bacteroides vulgatus (pink)) compared to the high-NRCD group. Bacteroides dorei and Bacteroides vulgatus are common gut Bacteroides species that have been associated with inflammatory signatures in several studies; enrichment of these taxa in NRCD may reflect a pro‑inflammatory community configuration [14]. Bacteroides cellulosilyticus is a polysaccharide‑degrading Bacteroides species that ferments complex dietary fibers into short‑chain fatty acids (SCFAs) such as acetate and propionate, and shifts in its abundance have been linked to altered gut metabolic states and inflammation [12]. Overall, the high-NRCD group appears to have a more uneven community composition dominated by fewer species, while the tCD control group shows a comparatively more balanced distribution among the top taxa, suggesting potential differences in microbial community structure between the two groups.

#### Beta Diversity: Bray-Curtis PCoA
<img width="937" height="392" alt="PCoA Bray-Curtis" src="https://github.com/user-attachments/assets/72997b44-88b9-4e43-852f-76530a8c6706" />

**Figure 3** This plot represents how similar or different the microbiomes are between the 2 groups. Each dot frepresents a sample the distance between dots reflects how different their microbial communities are. The red high‑NRCD samples cluster together on one side, while the cyan tCD controls cluster on the other, meaning the two groups have noticeably different overall microbiome compositions. This indicates that NRCD microbiomes resemble each other more closely than the the tCD controls, which supports the paper's conclusions of NRCD microbiomes being distinct dysbiotic community structure dominated by specific Bacteroides species and reduced SCFA‑producers, while tCD controls have a more balanced microbiome.

#### Beta Diversity 2: Jaccard PCoA
<img width="1000" height="500" alt="PCoA Jaccard" src="https://github.com/user-attachments/assets/233be297-cd31-447e-ac31-99e84c62a542" />

**Figure 4** This plot represents how different the microbiomes are between the two groups based on the presence or absence of the species, disregarding how much abundance they have. The dots that sit closer together share more of the same species. The red high‑NRCD samples cluster separately from the cyan tCD controls, meaning the two groups differ not just in abundance patterns but also in the overall taxa composition they contain. This shows that NRCD and tCD samples don’t just differ in how much of each species they, have but also in which species show up in the first place, reinforcing that the two groups have distinct microbial communities. This aligns with the paper’s observation that NRCD is associated with the loss of beneficial taxa like Faecalibacterium prausnitzii and the appearance or dominance of inflammation‑linked species.

----
## Discussion

Across all analyses, the microbial profiles of high‑NRCD patients consistently diverged from those of tCD controls, reinforcing the dysbiosis pattern described in the original study. The Bray–Curtis PCoA and NMDS plots both showed clear separation between groups, indicating that NRCD is associated with broad shifts in microbial abundance rather than random variation. The Jaccard PCoA added an important complementary perspective by demonstrating that the two groups also differ in species presence and absence, suggesting that NRCD involves the loss of beneficial taxa and the emergence of taxa not typically found in treated celiac individuals. These community level differences were further supported by the differential‑abundance volcano plot, which highlighted specific taxa enriched in each condition, mirroring the paper’s findings of reduced SCFA‑producing Firmicutes and increased Bacteroides associated species in NRCD. Together, these results suggest that NRCD is characterized by a stable but altered microbial community structure that may contribute to persistent symptoms despite adherence to a gluten‑free diet.

Taken together, these findings highlight the importance of examining NRCD not only as a clinical condition but also as a disorder with a strong microbial component. The consistency of the high‑NRCD samples across multiple beta‑diversity metrics suggests that once dysbiosis is established, it may form a stable and resilient community structure that persists despite dietary treatment. This stability could help explain why some patients continue to experience symptoms even with strict gluten avoidance. In contrast, the slightly more variable profiles in tCD controls may reflect a healthier, more adaptable microbial ecosystem. These observations underscore the need to consider the microbiome as a potential contributor to disease persistence and as a target for therapeutic intervention.  

This study is limited by a small, cross‑study sample set and differences in sequencing depth and library preparation between the source datasets, which can introduce batch effects and reduce power to detect low‑abundance taxa. Although there were maesures taken to control for depth by rarefaction and used sensitive classifiers (Kraken2/Bracken) to maximize taxonomic resolution, residual biases may remain and could influence differential abundance results. The use of Kraken2 increases sensitivity but may also change classification specificity at low abundance; therefore, key findings should be validated with independent methods like targeted qPCR, strain‑resolved metagenomics, or metabolomics, and in larger, longitudinal cohorts to establish temporal directionality and causal links.  

Future studies should move beyond taxonomic profiling to explore the functional consequences of this dysbiosis, such as shifts in metabolic pathways, SCFA production, or immune‑modulating microbial metabolites.  Faecalibacterium prausnitzii is a major butyrate‑producing commensal with well‑documented anti‑inflammatory effects on the gut mucosa; its depletion is commonly associated with impaired epithelial health and chronic intestinal inflammation [13]. Longitudinal sampling would help determine whether these microbial signatures precede NRCD onset or arise as a consequence of chronic inflammation. Additionally, integrating metagenomics with host transcriptomics or metabolomics could clarify how microbial changes interact with mucosal immune responses. Finally, intervention studies, such as targeted probiotics, microbiota directed foods, or fecal microbiota transplantation—may help determine whether modifying the gut microbiome can restore microbial balance and improve clinical outcomes in NRCD.

---

## Conclusion

Overall, the combined beta‑diversity and differential‑abundance analyses demonstrate that high‑NRCD patients harbor a distinctly altered gut microbiome compared to tCD controls, consistent with the dysbiosis signature described in the original study. While the two groups appear similar at broad taxonomic levels, finer‑scale analyses reveal clear shifts in community structure, species presence, and key functional taxa. These findings highlight the potential role of the microbiome in the persistence of symptoms in NRCD and underscore the need for deeper functional and longitudinal investigations to clarify how microbial alterations contribute to disease mechanisms and therapeutic outcomes.

---
## Dependencies and Packages
| Tool | Version | Source | Purpose |
|------|---------|--------|---------|
| SRA Toolkit | 3.3.0 | NCBI | Raw data download |
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

- [9] Bracken: Bayesian Reestimation of Abundance with Kraken (n.d.). John Hopkins University: Center for Computational Biology https://ccb.jhu.edu/software/bracken/

- [10] Dabdoub, SM (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software]. Available at https://github.com/smdabdoub/kraken-biom.

- [11] DESeq2 vs. ANCOM results. (2020). QIIME 2 Forum. https://forum.qiime2.org/t/deseq2-vs-ancom-results/1738

- [12] Dinić, M., Đokić, J., Jakovljević, S. et al. Insight into immunoregulatory and neuromodulatory capability of Bacteroides cellulosilyticus and Bacteroides xylanisolvens human gut microbiota isolates. Sci Rep 15, 38058 (2025). https://doi.org/10.1038/s41598-025-21839-0
- Dabdoub, SM (2016). kraken-biom: Enabling interoperative format conversion for Kraken results (Version 1.2) [Software]. Available at https://github.com/smdabdoub/kraken-biom.

- [13] Lenoir, M., Martín, R., Torres-Maravilla, E., Chadi, S., González-Dávila, P., Sokol, H., … Bermúdez-Humarán, L. G. (2020). Butyrate mediates anti-inflammatory effects of Faecalibacterium prausnitzii in intestinal epithelial cells through Dact3. Gut Microbes, 12(1). https://doi.org/10.1080/19490976.2020.1826748

- [14] Singh V, West G, Fiocchi C, Cominelli F, Good CE, Jacobs MR, Rodriguez-Palacios A. 2024. Genomes of Bacteroides ovatus, B. cellulosilyticus, B. uniformis, Phocaeicola vulgatus, and P. dorei isolated from gut cavernous fistulous tract micropathologies in Crohn's disease. Microbiol Resour Announc 13:e01152-23.
https://doi.org/10.1128/mra.01152-23 
