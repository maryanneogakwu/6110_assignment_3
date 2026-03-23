# *****************************************************************
# 6110 Genomic Methods for Bioinformatics 
#Shotgun Metagenomics Analysis
#Comparing high-NRCD vs tCD patients
#Maryanne Ogakwu
# *****************************************************************

#------------------- Install packages --------------
#if (!require("BiocManager")) install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("ANCOMBC")
#install.packages("ggplot2")
#install.packages("vegan")
#install.packages("biomformat")

#Load libraries
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(biomformat)

# ***********************************************************************
# -------------- Part 1 - Import Bracken/Kraken2 data -------------------
# ***********************************************************************

#Braken species reports for all 6 samples (high-nrcd and tCD) were combined into one BIOM file using Kraken-Biom on the DRAC Nibi Cluster
#Loading the combined BIOM file
setwd("C:/Users/marya/OneDrive/Desktop/6110_assign_3/data")
biom_data <- read_biom("combined4.biom")
biom_data

#Checking properties of the biom-data object 
colnames(biom_data(biom_data))
rownames(biom_data(biom_data))
dim(biom_data(biom_data))
physeq <- import_biom(biom_data)

# ******************************************************************
# ---------------- Part 2 - Create metadata table ------------------
# ******************************************************************

metadata <- read.table("metadata_table_R.txt",
                       header = TRUE,
                       sep = "\t",
                       skip = 1,
                       row.names = 1,    
                       strip.white = TRUE)
head(metadata)
colnames(metadata)
rownames(metadata)
# ********************************************************
# ------------- Part 3 - Build phyloseq object -----------
# ********************************************************
#Making sure the sample names in the phyloseq object and the metadata file match
sample_names(physeq)
rownames(metadata)

#Sample names dont match completely, so they will be renamed in the physeq object to match the metadata file
sample_names(physeq) <- gsub("_bracken_species", "", sample_names(physeq))
sample_names(physeq)

sample_data(physeq) <- sample_data(metadata)
physeq

#Checking how many taxa and samples are in the object
ntaxa(physeq)
nsamples(physeq)
sample_data(physeq)
# ************************************************************
# --------------- Part 4 - Rarefy to even depth --------------
# ************************************************************

# The tCD samples have more reads than high-NRCD samples (about 39-49M vs 17-20M reads), so to correct this I will rarefy to the depth of your lowest sample

# Check sample read counts
sample_sums(physeq)

# Rarefy to minimum depth
set.seed(123)  
physeq_rare <- rarefy_even_depth(physeq,
                             sample.size = min(sample_sums(physeq)),
                             rngseed = 123,
                             replace = FALSE)
physeq_rare

# ***************************************************************
# ------------------ Part 5 - Alpha diversity -------------------
# ***************************************************************

#Plot alpha diversity
plot_richness(physeq_rare,
              x = "Cluster",
              color = "Cluster",
              measures = c("Observed", "Shannon", "Simpson")) +
  geom_boxplot(alpha = 0.5) +
  labs(title = "Alpha Diversity: high-NRCD vs tCD",
       x = "Group",
       y = "Diversity") +
  theme_bw()

#Statistical test for alpha diversity
#Extracting diversity measures
alpha_div <- estimate_richness(physeq_rare,
                               measures = c("Observed", "Shannon", "Simpson"))
alpha_div$Group <- metadata[rownames(alpha_div), "Group"]
head(alpha_div)
colnames(alpha_div)

#Adding the Cluster column from metadata
alpha_div$Cluster <- metadata[rownames(alpha_div), "Cluster"]

#Wilcoxon test (non-parametric, matches paper's approach)
wilcox.test(Shannon ~ Cluster, data = alpha_div)
wilcox.test(Simpson ~ Cluster, data = alpha_div)
wilcox.test(Observed ~ Cluster, data = alpha_div)

# ************************************************************
# -------------- Part 6 - Relative abundance plots -----------
# ************************************************************
head(tax_table(physeq))
tax_table(physeq)

#The taxanomic table contains columns named "Rank 1 - Rank 7", so it needs to be adjusted to contain the standard classification names "Kindom, Phylum, etc"
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", 
                                 "Order", "Family", "Genus", "Species")
head(tax_table(physeq))

#Removing the prefixes that show up in the taxonomy table
tax_table(physeq) <- gsub(".__", "", tax_table(physeq))
head(tax_table(physeq))

#Updating my rarefied object
colnames(tax_table(physeq_rare)) <- c("Kingdom", "Phylum", "Class",
                                      "Order", "Family", "Genus", "Species")
tax_table(physeq_rare) <- gsub(".__", "", tax_table(physeq_rare))
head(tax_table(physeq_rare))

#Convert to relative abundance
physeq_rel <- transform_sample_counts(physeq_rare, function(x) x / sum(x))
#Phylum level plot
physeq_phy <- tax_glom(physeq_rel, taxrank = "Phylum")
# Melt for plotting
df_phy <- psmelt(physeq_phy)

#````````
#Relative abundance plot 1
ggplot(df_phy, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Cluster, scales = "free_x") +
  labs(title = "Phylum-level Relative Abundance",
       y = "Relative Abundance",
       x = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#````````
#Relative abundance plot 2
# Species level plot - top 10 only 
ps_rel_species <- tax_glom(physeq_rel, taxrank = "Species")
top10a <- names(sort(taxa_sums(ps_rel_species), decreasing = TRUE))[1:10]
ps_top10a <- prune_taxa(top10a, ps_rel_species)
df_top10a <- psmelt(ps_top10a)

ggplot(df_top10a, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Cluster, scales = "free_x") +
  labs(title = "Top 10 Species Relative Abundance",
       y = "Relative Abundance",
       x = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Plot top 10 most abundant taxa
top10 <- names(sort(taxa_sums(physeq_rel), decreasing = TRUE))[1:10]
ps_top10 <- prune_taxa(top10, physeq_rel)
df_top10 <- psmelt(ps_top10)

ggplot(df_top10, aes(x = Sample, y = Abundance, fill = OTU)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Cluster, scales = "free_x") +
  labs(title = "Top 10 Species Relative Abundance",
       y = "Relative Abundance",
       x = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# *************************************************************
# -------------- Part 7 - Beta diversity ---------------------
# *************************************************************

# Bray-Curtis PCoA
ord.pcoa.bray <- ordinate(physeq_rare, method = "PCoA", distance = "bray")
plot_ordination(physeq_rare, ord.pcoa.bray,
                color = "Cluster",
                title = "Bray-Curtis PCoA: high-NRCD vs tCD") +
  geom_point(size = 4) +
  theme_bw()

# NMDS with Bray-Curtis
ord.nmds.bray <- ordinate(physeq_rare, method = "NMDS", distance = "bray")
plot_ordination(physeq_rare, ord.nmds.bray,
                color = "Cluster",
                title = "Bray-Curtis NMDS: high-NRCD vs tCD") +
  geom_point(size = 4) +
  theme_bw()

# Jaccard PCoA (presence/absence only)
ord.pcoa.jaccard <- ordinate(physeq_rare, method = "PCoA", distance = "jaccard")
plot_ordination(physeq_rare, ord.pcoa.jaccard,
                color = "Cluster",
                title = "Jaccard PCoA: high-NRCD vs tCD") +
  geom_point(size = 4) +
  theme_bw()

# **************************************************************
# ------------------- Part 8 - PERMANOVA -----------------------
# **************************************************************
# Testing whether community composition differs significantly between high-NRCD and tCD groups

metadata_df <- as(sample_data(physeq_rare), "data.frame")

# Bray-Curtis PERMANOVA
adonis2(phyloseq::distance(physeq_rare, method = "bray") ~ Cluster,
        data = metadata_df)

# Jaccard PERMANOVA
adonis2(phyloseq::distance(physeq_rare, method = "jaccard") ~ Cluster,
        data = metadata_df)

# *************************************************************
# --------- Part 9 - Differential abundance with DESeq2 -----
# *************************************************************
#Identifying which species are significantly different between high-NRCD and tCD groups

#Convert phyloseq to DESeq2 object
dds <- phyloseq_to_deseq2(physeq_rare, ~ Cluster)

#Run DESeq2
dds <- DESeq(dds, test = "Wald", fitType = "parametric")

#Get results
res <- results(dds, contrast = c("Cluster", "high-NRCD", "tCD_control"))
res <- res[order(res$padj), ]

#View significant results
sig <- subset(res, padj < 0.05)
sig

#Plot log fold changes
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

