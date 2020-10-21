
setwd("/Users/laura/Google Drive/KEGG-SIF/no-indirect/comunidades-2level-withNames-isolatedNetworks/")

source("../bin/modularity_zscore.R")

mn.nodes <- "MN/MN.ncbi.pathways.nodes"
mn.edges <- "MN/MN.ncbi.pathways.edges"
mn.pdf <- "MN/Modularity-MN.pdf"

mn.q <- modularity_analysis(mn.nodes, mn.edges, mn.pdf)

rn.nodes <- "RN/RN.ncbi.pathways.nodes"
rn.edges <- "RN/RN.ncbi.pathways.edges"
rn.pdf <- "RN/Modularity-RN.pdf"

rn.q <- modularity_analysis(rn.nodes, rn.edges, rn.pdf)

ppi.nodes <- "PPI/PPI.ncbi.pathways.nodes"
ppi.edges <- "PPI/PPI.ncbi.pathways.edges"
ppi.pdf <- "PPI/Modularity-PPI.pdf"

ppi.q <- modularity_analysis(ppi.nodes, ppi.edges, ppi.pdf)
