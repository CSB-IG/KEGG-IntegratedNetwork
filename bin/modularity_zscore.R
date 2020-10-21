

library(igraph)
library(ggplot2)

modularity_analysis <- function(nodes.file, edges.file, pdf.file){
  nodes <- read.table(file = nodes.file, header=F, stringsAsFactors = F)
  names(nodes) <- c("name", "community")
  
  edges <- read.table(file = edges.file, header=F, stringsAsFactors = F)
  names(edges) <- c("from", "to", "npath", "path")
  
  net <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
  
  modularity.net <- modularity(net, nodes$community)
  
  modularity.net
  
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
  ####### HISTOGRAM RANDOM
  
  pdf (pdf.file)
  print(ggplot(as.data.frame(modularity.permutate), aes(x=modularity.permutate)) + 
          geom_histogram(color="darkblue", fill="lightblue", binwidth = 0.0005) + 
          geom_vline(aes(xintercept=modularity.net),
                     color="black",  size=0.5, linetype="dashed")  + 
          labs( x = "Q value", y = "Frequency") + 
          theme_bw())
  
  print(ggplot(as.data.frame(modularity.permutate), aes(x=modularity.permutate)) + 
          geom_histogram(color="darkblue", fill="lightblue") + 
          labs( x = "", y = "") + 
          theme_bw(base_size = 22))
  
  dev.off() 
  
  return(modularity.net)
}

