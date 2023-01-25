# Note importing BioC pkgs after dplyr requires explicitly using dplyr::select()

library(dplyr)
library(DESeq2)

# The below code is an example that was provided by the author - can use the DESeq2 portion to apply to your own RNAseq expression dataset

# Which data do you want to use? Let's use the sailfish counts.
# browseURL("http://dx.doi.org/10.6084/m9.figshare.1601975")
# countDataURL = "http://files.figshare.com/2439061/GSE37704_featurecounts.csv"
countDataURL = "http://files.figshare.com/2600373/GSE37704_sailfish_genecounts.csv"

# Import countdata
countData = read.csv(countDataURL, row.names=1) %>% 
  dplyr::select(-length) %>% 
  as.matrix()

# Filter data where you only have 0 or 1 read count across all samples.
countData = countData[rowSums(countData)>1, ]
head(countData)

# Import metadata
colData = read.csv("http://files.figshare.com/2439060/GSE37704_metadata.csv", row.names=1)
colData

# Set up the DESeqDataSet Object and run the DESeq pipeline
dds = DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design=~condition)
dds = DESeq(dds)
dds

# Next, get results for the HoxA1 knockdown versus control siRNA, and reorder them by p-value. 
# Call summary on the results object to get a sense of how many genes are up or down-regulated at FDR 0.1.
res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
res = res[order(res$pvalue),]
summary(res)

# Since we mapped and counted against the Ensembl annotation, our results only have information about Ensembl gene IDs. 
# But, our pathway analysis downstream will use KEGG pathways, and genes in KEGG pathways are annotated with Entrez gene IDs. 
# The author wrote an R package for doing this offline the dplyr way (https://github.com/stephenturner/annotables), 
# but the canonical Bioconductor way to do it is with the AnnotationDbi and organism annotation packages. Here we’re using the 
# organism package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi database package (“db”) using Entrez Gene IDs (“eg”) as primary keys. 
# To see what all the keys are, use the columns function.

library("AnnotationDbi")
library("org.Hs.eg.db")
columns(org.Hs.eg.db)

res$symbol = mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez = mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
res$name =   mapIds(org.Hs.eg.db,
                     keys=row.names(res), 
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
head(res, 10)

