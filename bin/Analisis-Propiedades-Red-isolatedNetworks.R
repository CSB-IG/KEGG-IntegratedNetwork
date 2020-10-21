library(plotly)
library(ggplot2)

source("bin/undirected_plots.R")
source("bin/directed_plots.R")

############## PPI
############## UNDIRECTED
############## 

ppi.nodes.all.file <- c("sif-isolatedNetworks/PPI.ncbi.pathways.nodes.csv")
ppi.nodes.gc.nodes.file <- c("sif-isolatedNetworks/PPI.ncbi.pathways.BigComponent.nodes.csv")
ppi.nodes.gc.edges.file <- c("sif-isolatedNetworks/PPI.ncbi.pathways.BigComponent.edges.csv")

ppi.pdf.all.file <- c("Analisis-Propiedades-isolatedNetworks/PPI-All.pdf")
ppi.pdf.gc.file <- c("Analisis-Propiedades-isolatedNetworks/PPI-BigComponent.pdf")

plot_all_network(ppi.nodes.all.file, ppi.pdf.all.file)
plot_gcc_network(ppi.nodes.gc.nodes.file, ppi.nodes.gc.edges.file, ppi.pdf.gc.file) 


############## RN
############## 
############## 

rn.nodes.all.file <- c("sif-isolatedNetworks/RN.ncbi.pathways.nodes.csv")
rn.nodes.gc.nodes.file <- c("sif-isolatedNetworks/RN.ncbi.pathways.BigComponent.nodes.csv")
rn.nodes.gc.edges.file <- c("sif-isolatedNetworks/RN.ncbi.pathways.BigComponent.edges.csv")

rn.pdf.all.file <- c("Analisis-Propiedades-isolatedNetworks/RN-All.pdf")
rn.pdf.gc.file <- c("Analisis-Propiedades-isolatedNetworks/RN-BigComponent.pdf")

plot_all_network_directed(rn.nodes.all.file, rn.pdf.all.file)
plot_gcc_network_directed(rn.nodes.gc.nodes.file, rn.nodes.gc.edges.file, rn.pdf.gc.file) 


############## MN
############## 
############## 

mn.nodes.all.file <- c("sif-isolatedNetworks/MN.ncbi.pathways.nodes.csv")
mn.nodes.gc.nodes.file <- c("sif-isolatedNetworks/MN.ncbi.pathways.BigComponent.nodes.csv")
mn.nodes.gc.edges.file <- c("sif-isolatedNetworks/MN.ncbi.pathways.BigComponent.edges.csv")

mn.pdf.all.file <- c("Analisis-Propiedades-isolatedNetworks/MN-All.pdf")
mn.pdf.gc.file <- c("Analisis-Propiedades-isolatedNetworks/MN-BigComponent.pdf")

plot_all_network_directed(mn.nodes.all.file, mn.pdf.all.file)
plot_gcc_network_directed(mn.nodes.gc.nodes.file, mn.nodes.gc.edges.file, mn.pdf.gc.file) 

