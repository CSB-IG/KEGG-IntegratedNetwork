
library(igraph)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

data.dir <- args[1]
pdf.out <- args[2]
type.element <- args[3]
limit <- as.numeric(args[4])
limit.by <- as.numeric(args[5])


pdf(pdf.out, onefile = T) 


result.all  <- as.data.frame(t(sapply(seq(0, limit, limit.by), function(x){
  file = paste(paste(paste(data.dir, type.element, sep="/"), x, sep=""), "final.txt", sep=".")
  net.matrix <- as.matrix(read.table(file, header=F, stringsAsFactors = F))
  graph <- graph_from_edgelist(net.matrix, directed = FALSE)
  vias.number <- length(V(graph))
  interacciones.number <- gsize(graph)
  degree.median <- median(degree(graph, v = V(graph), loops = FALSE, normalized = TRUE))
  path.mean <- mean_distance(graph, directed = FALSE, unconnected = TRUE)
  componentes.no <- components(graph)$no
  result <- c(vias.number, interacciones.number, degree.median, path.mean, componentes.no)
  return(result)
}, simplify = T)))

names(result.all) <- c("vias", "interacciones", "degree", "path.mean", "components")
vias.total <- result.all$vias[1]
result.all$isolated <- c(vias.total)-result.all$vias
result.all$removed <- seq(0, limit, limit.by)

#plot(result.all$removed, result.all$components)
#plot(result.all$removed, result.all$vias)
#plot(result.all$removed, result.all$interacciones)
#plot(result.all$removed, result.all$degree)
#plot(result.all$removed, result.all$path.mean)


combined.data <- data.frame(removed=rep(result.all$removed, 5), 
                            variable=c(rep("components", nrow(result.all)), rep("pathways", nrow(result.all)), rep("edges", nrow(result.all)), 
                                       rep("mean.degree", nrow(result.all)), rep("shortest.path", nrow(result.all))), 
                            values=c(result.all$components, result.all$vias, result.all$interacciones,
                                     result.all$degree, result.all$path.mean))

ggplot(combined.data, aes(x = removed, y = values)) +
  geom_line(aes(color = variable)) +
  facet_grid(variable ~ ., scales = "free_y") +
  theme(strip.text.y = element_blank(),  
        axis.text.x = element_text(size = 17), 
        axis.text.y = element_text(size = 17),  
        legend.position = "none", 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        panel.spacing = unit(1.1, "lines")) 

dev.off()
