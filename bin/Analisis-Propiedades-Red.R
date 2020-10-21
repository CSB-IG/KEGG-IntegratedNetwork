library(plotly)
library(ggplot2)

setwd("/Users/laura/Google Drive/KEGG-SIF/no-indirect/Analisis-Propiedades/")

######### TOPOLOGY ALL NETWORK -MINUS ISOLATED NODES
######### 
######### 


pdf("Analisis-Propiedades/Topological-Metrics-All.pdf", onefile = T)

nodes.all <- read.table("sif-allKEGG/AllKEGG.ncbi.filter.header.txt.Comunidades.node.csv", sep=",", header=T, stringsAsFactors = F)

###### LOG-LOG

table.k.log <- table(log10(nodes.all$Degree))
max.table.k.log <- max(as.numeric(names(table.k.log)))


y <- as.vector(log10(table.k.log))
x <- as.numeric(names(table.k.log))

data <- data.frame(x=x, y=y)
ggplot(data = data, aes(x = x, y = y)) + 
  geom_point(color='blue') +
  geom_smooth(method = "lm", se = FALSE, col="black") + 
  labs( x = "log10(k)", y = "log10(Frequency)") +
  ylim(0, max(data$y)) + 
  theme_bw(base_size = 25)

dev.off()




############# 
############# 
############# 
############# 




pdf("Analisis-Propiedades/Topological-Metrics.pdf", onefile = T)
############# 
############# EDGES
############# 
############# 

edges <- read.table("sif-allKEGG/AllKEGG.ncbi.filter.header.BigComponent.edges.csv", sep=",", header=T, stringsAsFactors = F)


dfP <- as.data.frame(edges$NPATHWAYS)

pdf("Analisis-Propiedades/Cumulative-per-edge.pdf")
ggplot(dfP, aes(edges$NPATHWAYS)) + stat_ecdf(geom = "point", size = 1, color='blue') +
  xlab( "Number of pathways per edge") +
  ylab("Cumulative frequency") +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=0.3) +
  theme_bw(base_size = 25)
dev.off()
  
dfE <- as.data.frame(log10(edges$EdgeBetweenness))

pdf("Analisis-Propiedades/Cumulative-edge-betwenness.pdf")
ggplot(dfE, aes(log10(edges$EdgeBetweenness))) + stat_ecdf(geom = "point", size = 0.3, color='blue') +
  xlab("log10(Edge Betweenness)") +
  ylab("Cumulative frequency") +
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=0.3) + 
  theme_bw(base_size = 25)
dev.off()


dfEBV <- data.frame(logEB = log10(edges$EdgeBetweenness), vias = edges$NPATHWAYS)

pdf("Analisis-Propiedades/Pathways-betwenness.pdf")
ggplot(dfEBV, aes(x=logEB, y=vias)) + 
  geom_point(size = 0.5, color='blue')+
  labs( x = "log10(Edge Betweenness)", y = "Number of pathways") + 
  theme_bw(base_size = 25)
dev.off()

plot(edges$EdgeBetweenness, edges$NPATHWAYS, xlab ="log10(Edge Betweenness)", ylab = "Number of pathways", cex = 0.5)

cor.test(log10(edges$EdgeBetweenness), edges$NPATHWAYS, alternative = "two.sided", method="spearman")


#### EDGES THAT BELONG TO ONE PATHWAY
sum(edges$NPATHWAYS == 1)/nrow(edges)

#### EDGES THAT BELONG TO MORE THAN 10 PATHWAYS
sum(edges$NPATHWAYS > 10)/nrow(edges)

#### EDGES THAT BELONG TO MORE THAN 20 PATHWAYS
sum(edges$NPATHWAYS > 20)/nrow(edges)



############# 
############# NODES
############# 
############# 

nodes <- read.table("sif-allKEGG/AllKEGG.ncbi.filter.header.BigComponent.nodes.csv", sep=",", header=T, stringsAsFactors = F)

table.k <- table(nodes$Degree)
max.table.k <- max(as.numeric(names(table.k)))

###### LOG-LOG

table.k.log <- table(log10(nodes$Degree))
max.table.k.log <- max(as.numeric(names(table.k.log)))

y <- as.vector(log10(table.k.log))
x <- as.numeric(names(table.k.log))

data <- data.frame(x=x, y=y)
pdf("Analisis-Propiedades/Degree-Distribution.pdf")
ggplot(data = data, aes(x = x, y = y)) + 
  geom_point(color='blue') +
  geom_smooth(method = "lm", se = FALSE, col="black") + 
  labs( x = "log10(k)", y = "log10(Frequency)") +
  ylim(0, max(data$y)) + 
  theme_bw(base_size = 25)
dev.off()


pdf("Analisis-Propiedades/Average-Shortest-Path.pdf")
ggplot(as.data.frame(nodes$AverageShortestPathLength), aes(x=nodes$AverageShortestPathLength)) + 
  geom_histogram(color="darkblue", fill="white") + 
  labs( x = "Average Shortest Path Length", y = "Frequency") +
  theme_bw(base_size = 25)
dev.off()


pdf("Analisis-Propiedades/core-periphery-GCC.pdf")
ggplot(as.data.frame(nodes$ClosenessCentrality), aes(x=nodes$ClosenessCentrality)) + 
  geom_histogram(color="darkblue", fill="white", bins = 50) + 
  labs( x = "Closeness Centrality", y = "Frequency") + 
  theme_bw(base_size = 25)
dev.off()

