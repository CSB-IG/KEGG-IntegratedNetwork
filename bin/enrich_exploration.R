setwd("/Users/laura/Google Drive/KEGG-SIF/no-indirect/")

data.file <- "comunidades-2level-withNames/2-level-summaryPerVia.txt"
data.enri <- "comunidades-2level-withNames/2-level-summaryPerVia-enrichment.txt"


dataComu <- read.delim(file = data.file)
dataComu$ID -> rownames(dataComu)

dataComu <- dataComu[,-1] #Eliminar la columna de nombres repetida

dataEnri <- read.delim(file = data.enri , header = FALSE)

#### Funcion
getgenes <- function(i, df, comparacion){
  col1 = as.character(comparacion[i,2])
  col2 = as.character(comparacion[i,3])
  logical <- df[,col1] & df[,col2]
  genes <- rownames(df)[logical]
  return(genes)
}
####

genesPerComp <- lapply(1:nrow(dataEnri), 
                       getgenes, 
                       df = dataComu, 
                       comparacion = dataEnri)

names.comp <- as.vector(dataEnri[,2]) #Nombres de las comparaciones, columna 2 de los resultados del enriquecimiento de comunidades. 


names(genesPerComp)[1:length(genesPerComp)] <- names.comp #Rango segun el numero de comparaciones

#################################################################################
########## SEPARACION DE LAS COMUNIDADES ENRIQUECIDAS EN MAS DE UNA VIA
#################################################################################

names.comp <- as.vector(dataEnri[,2])

comuEn.M2 <- as.list(unique(names.comp[duplicated(names.comp)])) #Nombres de comunidades enriquecidas mas de una vez.
#length(comuEn.M2)
# [1] 95 KEGGSIF
# [1] 34 Apoptosis

#### Función
getlistsPerComuRep <- function( i, lists, listnames){
  logici <- names(lists)  == listnames[i] #Vector logico por nombre
  list_perComuRep <- lists[logici]
}
####

genesLComEnri <- lapply(1:length(genesPerComp), 
                        FUN = getlistsPerComuRep, 
                        lists =genesPerComp, 
                        listnames = comuEn.M2)


genesLComEnri <- genesLComEnri[1:length(comuEn.M2)] #Rango del total de comunidades enriquecidas mas de una vez

###### Asignar nombres completos 
###### Son los de la columna 1 del archivo dataEnri

asign.names <- function(i, list, dataEnri){
  Lnames <- unique(names(list[[i]]))
  namesdf <- split(dataEnri$V1, dataEnri$V2)
  names(list[[i]]) <- unlist(namesdf[Lnames])
  return(list[[i]])
  
}

genesLComEnri.NC <- lapply(1:length(genesLComEnri), asign.names, list = genesLComEnri, dataEnri = dataEnri)

################# FUNCIONES
########################
########################
########################

pareado <- function(activo, prueba, lista, id.thr){
  longactivo <- length(lista[[activo]])
  longprueba <- length(lista[[prueba]])
  longintersecto <- length(intersect(lista[[activo]], lista[[prueba]]))
  if (longintersecto > id.thr*longprueba){
    res <- c(activo)
  } else {
    res <- c(activo, prueba)
  }
  return(res)
}

contra1 <- function(activos, comparacion, lista, id.thr){
  res.pareados <- lapply(activos, 
                         pareado, 
                         prueba = comparacion, 
                         lista = lista, 
                         id.thr = id.thr)
  nuevos.activos <- unlist(res.pareados)
  if (sum(nuevos.activos == comparacion) == length(activos)){
    return(c(activos, comparacion))
  }else{
    return(activos)
  }
  
}

comp.gandalla <- function(activos, posiblesAgandallados, lista, id.thr){
  for (x in 1:length(posiblesAgandallados) ){
    activos <- contra1(activos, 
                       comparacion=posiblesAgandallados[x], 
                       lista=lista, 
                       id.thr = id.thr)
  }
  return(activos)
}

obtener_ganadores <- function(j, id.thr){
  j.order <- j[order(sapply(j, length, simplify = T), decreasing=T)]
  ganadores.posicion <- comp.gandalla(activos = c(1) , 
                                      posiblesAgandallados = c(2:length(j.order)), 
                                      lista = j.order, 
                                      id.thr = id.thr)
  ganadores.Enriquecimiento <- j.order[ganadores.posicion]
}

########################
########################
########################

Vias.EnriquecidasPerCom.Ganadoras.0.5 <- lapply(genesLComEnri.NC, obtener_ganadores, 0.5)
Vias.EnriquecidasPerCom.Ganadoras.0.6 <- lapply(genesLComEnri.NC, obtener_ganadores, 0.6)
Vias.EnriquecidasPerCom.Ganadoras.0.7 <- lapply(genesLComEnri.NC, obtener_ganadores, 0.7)
Vias.EnriquecidasPerCom.Ganadoras.0.8 <- lapply(genesLComEnri.NC, obtener_ganadores, 0.8)
Vias.EnriquecidasPerCom.Ganadoras.0.9 <- lapply(genesLComEnri.NC, obtener_ganadores, 0.9)
Vias.EnriquecidasPerCom.Ganadoras.1.0 <- lapply(genesLComEnri.NC, obtener_ganadores, 1.0)


hist(sapply(Vias.EnriquecidasPerCom.Ganadoras.0.5, length), breaks = seq(0.5, 25.5), 
     main = "", xlab = "Number of enriched pathways per module\nIntersection size=0.5" )
hist(sapply(Vias.EnriquecidasPerCom.Ganadoras.0.6, length), breaks = seq(0.5, 25.5), 
     main = "", xlab = "Number of enriched pathways per module\nIntersection size=0.6" )
hist(sapply(Vias.EnriquecidasPerCom.Ganadoras.0.7, length), breaks = seq(0.5, 25.5), 
     main = "", xlab = "Number of enriched pathways per module\nIntersection size=0.7" )
hist(sapply(Vias.EnriquecidasPerCom.Ganadoras.0.8, length), breaks = seq(0.5, 25.5), 
     main = "", xlab = "Number of enriched pathways per module\nIntersection size=0.8" )
hist(sapply(Vias.EnriquecidasPerCom.Ganadoras.0.9, length), breaks = seq(0.5, 25.5), 
     main = "", xlab = "Number of enriched pathways per module\nIntersection size=0.9" )
hist(sapply(Vias.EnriquecidasPerCom.Ganadoras.0.9, length), breaks = seq(0.5, 25.5), 
     main = "", xlab = "Number of enriched pathways per module\nIntersection size=1.0" )

nodes.shared <- unlist(sapply(genesLComEnri.NC, names))
number.nodes <- unlist(sapply(genesLComEnri.NC, function(x){
  sapply(x, length)
}))


hist(number.nodes)
