library(igraph)
library(ggplot2)


nodes <- read.table(file = "comunidades-2level-withNames/Analisis-Comunidades/AllKEGG.ncbi.filter.nodes", header=F, stringsAsFactors = F)
names(nodes) <- c("name", "community")

edges <- read.table(file = "comunidades-2level-withNames/Analisis-Comunidades/AllKEGG.ncbi.filter.edges", header=F, stringsAsFactors = F)
names(edges) <- c("from", "to", "npath", "path")

net <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)

modularity.net <- modularity(net, nodes$community)
### 0.6824

modularity.permutate <- sapply(seq(1,1000), function(x, net, nodes){
  net.permute <- permute(net, permutation= sample(seq(1,length(V(net))), length(V(net))))
  modularity(net.permute, nodes$community)
}, net=net, nodes=nodes, simplify = T)


z.test = function(x,mu,popvar){
  one.tail.p <- NULL
  z.score <- (mean(x)-mu)/(popvar/sqrt(length(x)))
  one.tail.p <- pnorm(abs(z.score),lower.tail = FALSE)
  cat(" z =",z.score,"\n","one-tailed probability =", one.tail.p,"\n",
      "two-tailed probability =", 2*one.tail.p )
}

z.test(modularity.net, mean(modularity.permutate), sd(modularity.permutate))
#z = 1425.736
#one-tailed probability = 0 
#two-tailed probability = 0

####### HISTOGRAM RANDOM

pdf ("comunidades-2level-withNames/Analisis-Comunidades/Q-value.pdf")
ggplot(as.data.frame(modularity.permutate), aes(x=modularity.permutate)) + 
  geom_histogram(color="darkblue", fill="lightblue", binwidth = 0.0005) + 
  geom_vline(aes(xintercept=modularity.net),
             color="black",  size=0.5, linetype="dashed")  + 
  labs( x = "Q value", y = "Frequency") + 
  theme_bw(base_size = 22) 
dev.off()

pdf("comunidades-2level-withNames/Analisis-Comunidades/Q-value-Dist.pdf")
ggplot(as.data.frame(modularity.permutate), aes(x=modularity.permutate)) + 
  geom_histogram(color="darkblue", fill="lightblue") + 
  labs( x = "", y = "") + 
  theme_bw(base_size = 32)
dev.off()

############# SCATTERPLOT - TAMA??O COMUNIDAD - NUMERO DE VIAS ENRIQUECIDAS

summary <- read.table("../2-level-summaryPerVia.txt", header=T, sep="\t", stringsAsFactors = F)

comunidades <- colnames(summary)[2:322]

comunidades.size <- sapply(comunidades, function(x, summary){
  sum(summary[,x]) 
}, summary = summary, simplify = T)

enriq <- read.table("comunidades-2level-withNames/2-level-summaryPerVia_ALLKEGG_enrichment_NoRedundantes.txt", header=T, sep="\t", stringsAsFactors = F)

comunidades.vias <- sapply(comunidades, function(x, enriq){
  nrow(subset(enriq, Comunidad == x))
}, enriq = enriq)


#### CORRELACION
df <- data.frame(size = comunidades.size, vias = comunidades.vias)

pdf("comunidades-2level-withNames/Analisis-Comunidades/CommunitySize.pdf")
ggplot(df, aes(x=size, y=vias)) + 
  geom_point()+
  geom_smooth(method=lm) +
  labs( x = "Community size", y = "Number of enriched pathways") + 
  theme_bw(base_size = 22)
dev.off()

sum(comunidades.vias == 0)
#74

sum(comunidades.vias == 1)
#163

cor(comunidades.size, comunidades.vias, method="spearman")
cor.test(comunidades.size, comunidades.vias, alternative = "two.sided", method="spearman")

#### CORRELACION -95% TRIMMED DISTRIBUTION
size.quantile <- quantile(df$size, probs=seq(0,1,0.025))
df.trimmed <- subset(df, size >= size.quantile[2] & size < size.quantile[40]) 

pdf("comunidades-2level-withNames/Analisis-Comunidades/CommunitySize-Trimmed.pdf")
ggplot(df.trimmed, aes(x=size, y=vias)) + 
  geom_point()+
  geom_smooth(method=lm) +
  labs( x = "Community size", y = "Number of enriched pathways") + 
  theme_bw()
dev.off()

cor(comunidades.size, comunidades.vias, method="spearman")
cor.test(comunidades.size, comunidades.vias, alternative = "two.sided", method="spearman")
