# Load in required packages
library(pathview)
library(gage)
library(gageData)

# Need to use kegg.sets.hs, kegg.gs, and kegg.gs.dise for the analysis as these all include differening pathways
# This will require running the same analysis multiple times

# Load kegg databases into R environment
data(kegg.sets.hs)
data(kegg.gs)
data(kegg.gs.dise)

# Need to convert ENSEMBLIDs to ENTREZIDs for the KEGG pathway analysis - KEGG stores information in Entrez ID format
differential_expression_res$entrez = mapIds(org.Hs.eg.db,
                                            keys=row.names(res), 
                                            column="ENTREZID",
                                            keytype="ENSEMBL",
                                            multiVals="first")

# The gage() functon requires a named vector of fold changes, where the names of the values are the Entrez gene IDs
foldchanges = differential_expression_res$log2FoldChange
names(foldchanges) = differential_expression_res$entrez
head(foldchanges)

# Run the pathway analysis - use 'same.dir=TRUE' to get seperate lists for upregulated and downregulated

# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE) # gsets = any kegg database you are currently using
# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

# Process the results to pull out the top 5 upregulated pathways, then further process that just to get the IDs
# We’ll use these KEGG pathway IDs downstream for plotting

# Get the pathways
keggrespathways_up = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways_up

# Do the same for getting the top 5 downregulated pathways
keggrespathways_down = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways_down

# The pathview() function in the pathview package makes the plots
# Let’s write a function so we can loop through and draw plots for the top 5 pathways we created above

# Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# Plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))

# We can also do a similar process with gene ontology (GO) as above
# Can create pathway object and plot as well
             
# Start with this and then repeat the steps where the up/down regulated pathway
# Objects were created and then plot the reults
data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]
gobpres = gage(foldchanges, gsets=gobpsets, same.dir=TRUE)
lapply(gobpres, head)
