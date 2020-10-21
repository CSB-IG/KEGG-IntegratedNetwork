##### Fisher
args = commandArgs(trailingOnly=TRUE)

data.file <- args[1]
ncomunidades <- as.numeric(args[2])
pvalue <- as.numeric(args[3])
file.out = args[4]

#ncomunidades


data <- read.table(file = data.file,
                   header = TRUE,
                   sep = "\t"
)

#ncomunidades <- 218
### Total de moleculas en cada via y en cada comunidad

total.Pvia <- apply(data[, (ncomunidades+2):(ncol(data)-1)], 2, sum)#Marginales Total de cada via
total.Pcom <- apply(data[,2:(ncomunidades+1)], 2, sum)#Marginales Total de cada comunidad
save(total.Pcom, file = "comunidades-2level-withNames/total.PCom.R")
save(total.Pvia, file = "comunidades-2level-withNames/total.PVia.R")
total.pvia.m <- matrix(data = rep(total.Pvia, ncomunidades), #Num de repeticiones = num de comunidades 
                       nrow = length(total.Pcom), 
                       ncol = length(total.Pvia),
                       byrow = TRUE)
#dim(total.pvia.m)
# [1] 218 185

total.Pcom.m <- matrix(data = rep(total.Pcom, length(total.Pvia)), #num de repeticiones = num de vias
                       nrow = length(total.Pcom), 
                       ncol = length(total.Pvia),
                       byrow = FALSE)
# dim(total.Pcom.m)
# [1] 218 185

####Cuantos estan en cada comunidad y en cada via

via.comu <- list()
for (i in 2:(ncomunidades+1)){ #rango de comunidades 
  x <- as.vector(apply(data[,i]&data[, (ncomunidades+2):(ncol(data)-1) ],2,sum)) #rango de vias
  via.comu[[i-1]] <- x
}

#via.comu 
via.comu.m <- matrix(data = unlist(via.comu), ncol = length(total.Pvia), nrow = ncomunidades, byrow = TRUE)
#ncols segun el numero de vias y nrow segun el numero de comunidades. 

### Cuantos de TOTAL de cada via(185) que no estan en cada comunidad(218)

via.nocomun.m <- total.pvia.m - via.comu.m

#dim(via.nocomun.m)
# [1] 218 185

### TOTAL
total <- nrow(data)
#2515

### Cuantos de TOTAL de cada comunidad(4) que no estan en cada via

comunidad.novia.m <- total.Pcom.m - via.comu.m
# dim(comunidad.novia.m)
# [1] 218 185

### Cuantos del total no estan en la via y no estan en la comunidad

novia.nocomunidad.m <- total - via.comu.m - via.nocomun.m - comunidad.novia.m
#dim(novia.nocomunidad.m)
# [1] 218 185

####Tabla de continguencia para cada comunidad para cada via, conteos observados.

n.comparaciones <- ncomunidades * length(total.Pvia)
mat.cont.obs <- vector("list", n.comparaciones)

for (i in 1:n.comparaciones) {
  mat.cont.obs[[i]] <- matrix(NA, nrow=2, ncol=2)
}
for (k in 1:n.comparaciones) {
  mat.cont.obs[[k]][1,1] = t(via.comu.m)[k]
  mat.cont.obs[[k]][1,2] = t(comunidad.novia.m)[k]
  mat.cont.obs[[k]][2,1] = t(via.nocomun.m)[k]
  mat.cont.obs[[k]][2,2] = t(novia.nocomunidad.m)[k]
}

###Asignar los nombres correspondientes a casa tabla de continguencia
namescom <- colnames(data[2:(ncomunidades+1)])
namesvias <- colnames(data[(ncomunidades+2):(ncol(data)-1)])

namescom.long <- rep(namescom, each=length(total.Pvia))
namesvias.long <- rep(namesvias, times= ncomunidades)
names.m.test <- paste(namescom.long, namesvias.long,  sep=".")

#names.m.test <- as.vector(apply(expand.grid(namescom, namesvias), 1, paste, collapse="."))
# [1] "X1.hsa00010"   "X1.hsa00020"   "X1.hsa00030"   "X1.hsa00040"   "X1.hsa00051"
names(mat.cont.obs) <- names.m.test

save(mat.cont.obs, file = "comunidades-2level-withNames/tablascontinguencuaEnriKEGG.R")

#Prueba fisher.test, calculara todos los estadisticos...

tests <- sapply(mat.cont.obs, function(x){
  fisher.test(x, alternative = "greater")$p.value
  }, simplify = T)
# length(tests)
# [1] 96726
# tests$X1.hsa04933$observed #Verificacion de la prueba echa con anterioridad. 
pvalue.seg <- pvalue / n.comparaciones
test.sig <- tests[tests < pvalue.seg ]
#### Quedarnos solo son las comunidades que estan enriquecidas en ciertas vias de acuerdo 
#### a un p.value menos a 1e-7

res <- data.frame(comunidad = namescom.long[tests < pvalue.seg],
                  via = namesvias.long[tests < pvalue.seg], 
                  pvalor = test.sig)

write.table(x = res, file=file.out, sep="\t", quote=F, col.names = F)
