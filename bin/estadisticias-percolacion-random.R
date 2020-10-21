
library(igraph)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

data.dir <- args[1]
pdf.out <- args[2]
type.element <- args[3]
limit <- as.numeric(args[4])
limit.by <- as.numeric(args[5])
rep.number <- as.numeric(args[6])

#data.dir <- "Analisis-Percolacion-nodesRandom/percolacion-results"
#pdf.out <- "Analisis-Percolacion-nodesRandom/Estadisticas-percolacion-random.pdf"
#type.element <- "nodes_rep_"
#limit <- as.numeric("6750")
#limit.by <- as.numeric("10")
#rep.number <- as.numeric("50")

pdf(pdf.out, onefile = T) 


result.all  <- as.data.frame(t(sapply(seq(0, limit, by = limit.by), function(x, rep.number){
  sapply(seq(1,rep.number), function(rep, x){
    file = paste(paste(paste(paste(data.dir, type.element, sep="/"), rep, sep=""), x, sep="_"), "txt", sep=".")
    if (file.size(file) == 0){
      vias.number <- 0
      interacciones.number <- 0
      degree.median <- 0
      path.mean <- 0
      componentes.no <- 0
    }else{
      net.matrix <- as.matrix(read.table(file, header=F, stringsAsFactors = F))
      graph <- graph_from_edgelist(net.matrix, directed = FALSE)
      vias.number <- length(V(graph))
      interacciones.number <- gsize(graph)
      degree.median <- median(degree(graph, v = V(graph), loops = FALSE, normalized = TRUE))
      path.mean <- mean_distance(graph, directed = FALSE, unconnected = TRUE)
      componentes.no <- components(graph)$no
    }
    result <- c(vias.number, interacciones.number, degree.median, path.mean, componentes.no)
    return(result)
  }, x=x, simplify = T)
}, rep.number = rep.number, simplify = T)))



names(result.all) <- c("vias", "interacciones", "degree", "path.mean", "components")

vias <- apply(as.matrix(result.all[,seq(1, by=5, length.out = rep.number)]), 1, median)
interacciones <- apply(as.matrix(result.all[,seq(2, by=5, length.out = rep.number)]), 1, median)
degree.median <- apply(as.matrix(result.all[,seq(3, by=5, length.out = rep.number)]), 1, median)
path.mean <- apply(as.matrix(result.all[,seq(4, by=5, length.out = rep.number)]), 1, median)
componentes <- apply(as.matrix(result.all[,seq(5, by=5, length.out = rep.number)]), 1, median)

result.median <- data.frame(vias, interacciones, degree.median, path.mean, componentes)

vias.total <- result.median$vias[1]
result.median$isolated <- c(vias.total)-result.median$vias
result.median$removed <- seq(0, limit, limit.by)

#plot(result.median$removed, result.median$componentes.no)
#plot(result.median$removed, result.median$vias)
#plot(result.median$removed, result.median$interacciones)
#plot(result.median$removed, result.median$degree)
#plot(result.median$removed, result.median$path.mean)


combined.data <- data.frame(removed=rep(result.median$removed, 5), 
                            variable=c(rep("components", nrow(result.median)), rep("pathways", nrow(result.median)), rep("edges", nrow(result.median)), 
                                       rep("mean.degree", nrow(result.median)), rep("shortest.path", nrow(result.median))), 
                            values=c(result.median$componentes, result.median$vias, result.median$interacciones,
                                     result.median$degree, result.median$path.mean))

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
