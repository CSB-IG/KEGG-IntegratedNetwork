######################################################
######################################################
########## ELIMINAR LAS VIAS ENRIQUECIDAS REDUNDANTES
################## DATOS REALES
######################################################


args = commandArgs(trailingOnly=TRUE)

data.file <- args[1]
data.enri <- args[2]
# "comunidades-2level-withNames/total.PCom.R"
#total.Pcom <- args[3]
# "comunidades-2level-withNames/total.PVia.R"
#total.Pvia <- args[4]
# "comunidades-2level-withNames/tablascontinguencuaEnriKEGG.R"
#tablascontingencia <- args[5]
file.out = args[3]

#"~/Documentos/posgradoBIOQUIMICA/APOPTOSIS_KEGG-SIF/con-Ecrel/comunidades-2level-withNames/2-level-summaryPerVia.Apoptosis.txt"

#Archivos directorio: comunidades-2level-withNames
dataComu <- read.delim(file = data.file)
#dataComu <- read.delim(file = "~/Documentos/posgradoBIOQUIMICA/APOPTOSIS_KEGG-SIF/con-Ecrel/comunidades-2level-withNames/2-level-summaryPerVia.Apoptosis.txt")

dataComu$ID -> rownames(dataComu)

dataComu <- dataComu[,-1] #Eliminar la columna de nombres repetida
#"~/Documentos/posgradoBIOQUIMICA/APOPTOSIS_KEGG-SIF/con-Ecrel/comunidades-2level-withNames/2-level-summaryPerVia.APOP-enrichment.txt"
dataEnri <- read.delim(file = data.enri , header = FALSE)
#dataEnri <- read.delim(file = "~/Documentos/posgradoBIOQUIMICA/APOPTOSIS_KEGG-SIF/con-Ecrel/comunidades-2level-withNames/2-level-summaryPerVia.APOP-enrichment.txt" , header = FALSE)

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

pareado <- function(activo, prueba, lista){
  longactivo <- length(lista[[activo]])
  longprueba <- length(lista[[prueba]])
  longintersecto <- length(intersect(lista[[activo]], lista[[prueba]]))
  if (longintersecto > 0.7*longprueba){
    res <- c(activo)
  } else {
    res <- c(activo, prueba)
  }
  return(res)
}

contra1 <- function(activos, comparacion, lista){
  res.pareados <- lapply(activos, 
                         pareado, 
                         prueba = comparacion, 
                         lista = lista)
  nuevos.activos <- unlist(res.pareados)
  if (sum(nuevos.activos == comparacion) == length(activos)){
    return(c(activos, comparacion))
  }else{
    return(activos)
  }
  
}

comp.gandalla <- function(activos, posiblesAgandallados, lista){
  for (x in 1:length(posiblesAgandallados) ){
    activos <- contra1(activos, 
                       comparacion=posiblesAgandallados[x], 
                       lista=lista)
  }
  return(activos)
}

obtener_ganadores <- function(j){
  j.order <- j[order(sapply(j, length, simplify = T), decreasing=T)]
  ganadores.posicion <- comp.gandalla(activos = c(1) , 
                                      posiblesAgandallados = c(2:length(j.order)), 
                                      lista = j.order)
  ganadores.Enriquecimiento <- j.order[ganadores.posicion]
}

########################
########################
########################

Vias.EnriquecidasPerCom.Ganadoras <- lapply(genesLComEnri.NC, obtener_ganadores)

######### Generar la tabla resumen de las comunidades y vias enriquecidas ganadoras
#######################
#######################
#######################

#Objetos a utilizar creados en enriquecimiento_comunidades-via.R y genes_ComPerVia.R. 
#1 genesPerComp, lista con los genes de todos los enriquecimiento de cada comunidad por vias; AQUI. 
#2 total.Pcom, total.Pvia, total de elementos en cada comunidad y en cada via
load(file = "comunidades-2level-withNames/total.PCom.R")
load(file = "comunidades-2level-withNames/total.PVia.R")

#3 mat.cont.obs, en este objeto tenemos las tablas de contingencia
load(file = "comunidades-2level-withNames/tablascontinguencuaEnriKEGG.R")
#Ahora necesito obtener los nombres de los enriquecimientos de las ganadoras. 
# Nombres de todos los enriquecimientos, Apoptosis-219, ALLKEGG 511
names.comp.via <- as.vector(dataEnri[,1]) #Nombres de las comparaciones, columna 1 de los resultados del enriquecimiento de comunidades. 
names(genesPerComp)[1:length(genesPerComp)] <- names.comp.via #Rango segun el numero de comparaciones

names.all.comvia <- as.vector(names(genesPerComp))
# Nombres de todos los repetidos, Apoptosis-143, ALLKEGG-372
names.allrep.comvia <- as.vector(unlist(lapply(genesLComEnri.NC, names)))
# Nombres no repetidos y unicos, Apoptosis-76, ALLKEGG-139
names.Norep.comvia <- setdiff(names.all.comvia, names.allrep.comvia)
##### Lista con todos los enriquecimientos ganadores 
namescom.vias.ganadoras <- unlist(lapply(Vias.EnriquecidasPerCom.Ganadoras, names))

# Nombres con enriquecimientos no redundantes y los únicos, unir names.Norep.comvia y namescom.vias.ganadoras
############################################ GANADORAS
enriquecimientos.NoRedundantes <- c(namescom.vias.ganadoras, names.Norep.comvia)
############################################
######## Extraer las GANADORAS de la tabla resultado del enriquecimiento.
library(dplyr)
tabla.summ.Enri.NoRed <- data.frame(filter(dataEnri, 
                                           dataEnri$V1 %in% enriquecimientos.NoRedundantes))

######  Crear los columnas extras
########## NUMERO DE ELEMENTOS EN CADA COMUNIDAD
# nombres de comunidades
zz <- unique(as.character(tabla.summ.Enri.NoRed$V2))
# extraer numero de elementos de cada comunidad
nCom <- total.Pcom[zz]
#Convertir a dataframe para despues hacer la union
t1 <- data.frame(cname = names(nCom), ncomu = unname(nCom))
# Unir tabla inicial con el numero de elementos de cada comunidad
new.Summ.Com <- dplyr::inner_join(tabla.summ.Enri.NoRed, t1, by = c("V2" = "cname"))

########## NUMERO DE ELEMENTOS EN CADA VIA
zv <- unique(as.character(tabla.summ.Enri.NoRed$V3))
nVia <- total.Pvia[zv]
t2 <- data.frame(vname = names(nVia), nvia = unname(nVia))
new.Summ.ComVia <- dplyr::inner_join(new.Summ.Com, t2, by = c("V3" = "vname"))

#NUMERO DE ELEMENTOS DE LA COMUNIDAD QUE ESTAN EN LA VIA (INTERSECCION)
x <- mat.cont.obs[enriquecimientos.NoRedundantes]
totalCV <- unlist(lapply(x, FUN = function(x) { 
  x[1,1]
  }))
t3 <- data.frame(cvname = names(totalCV), nComVia = unname(totalCV))
new.Summ.ComVia.CV <- dplyr::inner_join(new.Summ.ComVia, t3, by = c("V1" = "cvname")) 

#ELEMENTOS DE LA INTERSECCION
y <- genesPerComp[enriquecimientos.NoRedundantes]
t4 <- data.frame(CVnames = names(y), elementsCV = NA)
t4$elementsCV <- sapply(y, paste, collapse = ",")
new.Summ.ComVia.CV.E <- dplyr::inner_join(new.Summ.ComVia.CV, t4, by = c("V1" = "CVnames")) 

colnames(new.Summ.ComVia.CV.E) <- c("ID", "Comunidad", "Via", "pValue", "ncomu", "nvia", "nComVia", "ElementsComVia")

#write.table(new.Summ.ComVia.CV.E, file = "~/Documentos/posgradoBIOQUIMICA/APOPTOSIS_KEGG-SIF/con-Ecrel/comunidades-2level-withNames/2-level-summaryPerVia_Apop_enrichment_NoRedundantes.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(new.Summ.ComVia.CV.E, file = file.out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
