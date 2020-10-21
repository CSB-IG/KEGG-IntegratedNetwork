######### TOPOLOGY ALL NETWORK -MINUS ISOLATED NODES
############# 
############# 

plot_all_network_directed <- function(nodes.file, pdf.file){
  pdf(pdf.file, onefile = T)
  
  nodes.all <- read.table(nodes.file, sep=",", header=T, stringsAsFactors = F)

  ######### LOG-LOG DEGREE DISTRIBUTION  
  
  table.k.log <- table(log10(nodes.all$Indegree + nodes.all$Outdegree ))
  max.table.k.log <- max(as.numeric(names(table.k.log)))
  min.table.k.log <- min(as.numeric(names(table.k.log)))
  
  partitions <- seq(floor(min.table.k.log), ceiling(max.table.k.log), length.out=4)
  y <- as.vector(log10(table.k.log))
  x <- as.numeric(names(table.k.log))
  #plot(log10(table.k.log), type="p", xaxt="n", xlab = "log10(k)", ylab="log10(Frequency)", col="blue")
  #axis(side = 1, at = partitions, labels = partitions)
  #abline(lm(y ~ x  ), lwd = 4)
  
  data <- data.frame(x=x, y=y)
  
  print(ggplot(data = data, aes(x = x, y = y)) + 
          geom_point(color='blue') +
          geom_smooth(method = "lm", se = FALSE, col="black") + 
          labs( x = "log10(k)", y = "log10(Frequency)") +
          ylim(0, max(data$y)) + 
          theme_bw(base_size = 25))
  
  dev.off()
}


######### TOPOLOGY GIANT CONNECTED COMPONENT
############# 
############# 

plot_gcc_network_directed <- function(nodes.file, edges.file, pdf.file){
  pdf(pdf.file, onefile = T)
  
  edges <- read.table(edges.file, sep=",", header=T, stringsAsFactors = F)
  nodes <- read.table(nodes.file, sep=",", header=T, stringsAsFactors = F)
  
  dfP <- as.data.frame(edges$NPATHWAYS)
  
  print(ggplot(dfP, aes(edges$NPATHWAYS)) + stat_ecdf(geom = "point", size = 1, color='blue') +
          xlab( "Number of pathways per edge") +
          ylab("Cumulative frequency") +
          geom_hline(yintercept = 1, linetype="dashed", 
                     color = "black", size=0.3) +
          theme_bw(base_size = 25))
  
  dfE <- as.data.frame(log10(edges$EdgeBetweenness))
  
  print(ggplot(dfE, aes(log10(edges$EdgeBetweenness))) + stat_ecdf(geom = "point", size = 0.3, color='blue') +
          xlab("log10(Edge Betweenness)") +
          ylab("Cumulative frequency") +
          geom_hline(yintercept = 1, linetype="dashed", color = "black", size=0.3) + 
          theme_bw(base_size = 25))
  
  
  dfEBV <- data.frame(logEB = log10(edges$EdgeBetweenness), vias = edges$NPATHWAYS)
  
  print(ggplot(dfEBV, aes(x=logEB, y=vias)) + 
          geom_point(size = 0.5, color='blue')+
          labs( x = "log10(Edge Betweenness)", y = "Number of pathways") + 
          theme_bw(base_size = 25))
  
  plot(edges$EdgeBetweenness, edges$NPATHWAYS, xlab ="log10(Edge Betweenness)", ylab = "Number of pathways", cex = 0.5)
  
  cor.test(log10(edges$EdgeBetweenness), edges$NPATHWAYS, alternative = "two.sided", method="spearman")
  
  
  ###### LOG-LOG DEGREE DISTRIBUTION
  
  table.k.log <- table(log10(nodes$Indegree + nodes$Outdegree))
  max.table.k.log <- max(as.numeric(names(table.k.log)))
  min.table.k.log <- min(as.numeric(names(table.k.log)))
  
  
  partitions <- seq(floor(min.table.k.log), ceiling(max.table.k.log), length.out=4)
  
  y <- as.vector(log10(table.k.log))
  x <- as.numeric(names(table.k.log))
  #plot(log10(table.k.log), type="p", xaxt="n", xlab = "log10(k)", ylab="log10(Frequency)", col="blue")
  #axis(side = 1, at = partitions, labels = partitions)
  #abline(lm(y ~ x  ), lwd = 4)
  
  data <- data.frame(x=x, y=y)
  
  print(ggplot(data = data, aes(x = x, y = y)) + 
    geom_point(color='blue') +
    geom_smooth(method = "lm", se = FALSE, col="black") + 
    labs( x = "log10(k)", y = "log10(Frequency)") +
    ylim(0, max(data$y)) + 
    theme_bw(base_size = 25))
  
  ####### SHORTEST-PATH-LENGTH
  
  print(ggplot(as.data.frame(nodes$AverageShortestPathLength), aes(x=nodes$AverageShortestPathLength)) + 
          geom_histogram(color="darkblue", fill="white") + 
          labs( x = "Average Shortest Path Length", y = "Frequency") +
          theme_bw(base_size = 25))
  
  ####### CORE-PERIPHERY
  
  print(ggplot(as.data.frame(nodes$ClosenessCentrality), aes(x=nodes$ClosenessCentrality)) + 
          geom_histogram(color="darkblue", fill="white", bins = 50) + 
          labs( x = "Closeness Centrality", y = "Frequency") + 
          theme_bw(base_size = 25))
  
  dev.off()
}
