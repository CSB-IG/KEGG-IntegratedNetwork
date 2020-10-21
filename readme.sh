####### KGML DOWNLOAD 
## get_kgml.sh

####### CONF DOWNLOAD
## sh bin/get_all_conf.sh


####### GET INTERACTIONS FROM KGML
mkdir sif-perKGML
sh bin/run_fromKGML.sh

####### CONCATENATE ALL KGMLS
mkdir sif-allKEGG
cat sif-perKGML/*sif | uniq > sif-allKEGG/AllKEGG.txt

####### DICTIONARY WITH MOLECULE TYPE 
python3.7 bin/genes_compoundFromSif.py -s sif-allKEGG/AllKEGG.txt -o diccionario/diccionario_molecule_type.txt


####### DICTIONARY WITH GENE NAME AND ENSEMBL ID 
python bin/FromDICT_toENSEMBL.py --sif sif-allKEGG/AllKEGG.txt --mol diccionario/diccionario_molecule_type.txt --ensembl anotacion/biomart-ensembl-genes96-GRCh38.txt --ncbi anotacion/gene_result_NCBI_GRCh38.p12.txt --out diccionario/diccionario_ncbi_ensembl.txt


####### REMOVE REDUNDANT ECREL INTERACTION
####### REMOVE PCREL INDIRECT-EFFECT INTERACTIONS
python3.7 bin/Filter_existent_ecrel.py --sif sif-allKEGG/AllKEGG.txt --out sif-allKEGG/AllKEGG.direct.txt


####### TRANSFORM SIF ALL KEG TO INCLUDE GENE NAMES 
python3.7 bin/Transform_Syf_ID_NCBI.py --sif sif-allKEGG/AllKEGG.direct.txt  --dict diccionario/diccionario_ncbi_ensembl.txt --out sif-allKEGG/AllKEGG.ncbi.txt


####### FILTER REPEATED INTERACTIONS 
####### COUNT NUMBER OF PATHWAYS PER INTERACTION
python3.7 bin/SimplifySyf_UniquePathway.py -s sif-allKEGG/AllKEGG.ncbi.txt -o sif-allKEGG/AllKEGG.ncbi.filter.txt


####### INCLUDE A COLUMN WITH ALL PATHWAYS NAMES FOR VISUALIZATION
python3.7 bin/SimplifySyf_UniquePathway_WithHSA.py -s sif-allKEGG/AllKEGG.ncbi.txt -o sif-allKEGG/AllKEGG.ncbi.filter.pathway.txt


####### CREATE EDGE FILE FOR VISUALIZATION
####### ADD HEADER FOR VISUALIZATION
awk -F "\t" '{print $1 " (interaction) " $3 "\t" $6}' sif-allKEGG/AllKEGG.ncbi.filter.pathway.txt > sif-allKEGG/AllKEGG.ncbi.filter.pathway.interactions.txt
cp ../focus-interaction-type/sif-allKEGG/header.ncbi.filter sif-allKEGG/
cat sif-allKEGG/header.ncbi.filter sif-allKEGG/AllKEGG.ncbi.filter.txt > sif-allKEGG/AllKEGG.ncbi.filter.header.txt


####### MODULE ANALYSIS
####### GENERATE NET FILE REQUIREF BY INFOMAP
mkdir comunidades-2level-withNames
python bin/FromSIF_To_NET.py -s sif-allKEGG/AllKEGG.ncbi.filter.txt -o comunidades-2level-withNames/AllKEGG.ncbi.filter.net


####### RUN INFOMAP
mkdir comunidades-2level-withNames/2-level/
/Users/laura/Desktop/Infomap/Infomap comunidades-2level-withNames/AllKEGG.ncbi.filter.net comunidades-2level-withNames/2-level/ --two-level --undirected  --include-self-links  --num-trials 1000 --seed 1 --clu --map --tree --bftree

###### GET GENES PER MODULE 
mkdir comunidades-2level-withNames/comunidades-genes/
python3.7 bin/Extraer_Nodos_Comunidad.py --map comunidades-2level-withNames/2-level/AllKEGG.ncbi.filter.map --out comunidades-2level-withNames/comunidades-genes/

###### GENERATE FILE FOR VISUALIZATION OF MODULES
mkdir comunidades-2level-withNames/Analisis-Comunidades/
python bin/Escribir-archivos-nodos-edges.py --map comunidades-2level-withNames/2-level/AllKEGG.ncbi.filter.map --kegg sif-allKEGG/AllKEGG.ncbi.filter.txt --out-nodes comunidades-2level-withNames/Analisis-Comunidades/AllKEGG.ncbi.filter.nodes --out-edges comunidades-2level-withNames/Analisis-Comunidades/AllKEGG.ncbi.filter.edges

###### IMPORT ON CYTOSCAPE: 
## NOMBRE DE COLUMA = ALL-PATHS
## NOMBRE DE TABLA: sif-allKEGG/AllKEGG.ncbi.filter.pathway.interactions.txt

## FILTER BIG COMPONENT
## ANALYZE NETWORK
## EXPORT NODES AND EDGES FILES
#  sif-allKEGG/AllKEGG.ncbi.filter.header.BigComponent.edges.csv
#  sif-allKEGG/AllKEGG.ncbi.filter.header.BigComponent.nodes.csv


####### FUNCTIONAL ENRICHMENT
####### ENRICHMENT STATISTICS
python3.7 bin/Tabla_Comunidades_ncbi.py -m comunidades-2level-withNames/2-level/AllKEGG.ncbi.filter.map -k sif-allKEGG/AllKEGG.ncbi.txt -o comunidades-2level-withNames/2-level-summaryPerVia.txt
Rscript bin/enriquecimiento_comunidades-vias.R comunidades-2level-withNames/2-level-summaryPerVia.txt 321 0.01 comunidades-2level-withNames/2-level-summaryPerVia-enrichment.txt
Rscript bin/getComuViasUniq.R comunidades-2level-withNames/2-level-summaryPerVia.txt comunidades-2level-withNames/2-level-summaryPerVia-enrichment.txt comunidades-2level-withNames/2-level-summaryPerVia_ALLKEGG_enrichment_NoRedundantes.txt



## ADD COMMUNITY INFO TO CYTOSCAPE
## NOMBRE DE COLUMA = COMUNIDADES
## NOMBRE DE TABLA: comunidades-2level-withNames/Analisis-Comunidades/AllKEGG.ncbi.filter.nodes
## ANALYZE NETWORK
## EXPORT NODES
# sif-allKEGG/AllKEGG.ncbi.filter.header.txt.Comunidades.node.csv


####### STROGEN SIGNALING PATHWAY
#######

mkdir hsa04915
grep "hsa04915" sif-allKEGG/AllKEGG.ncbi.filter.pathway.txt > hsa04915/AllKEGG.ncbi.filter.pathway.hsa04915.txt
grep "hsa04915" sif-isolatedNetworks/PPI.ncbi.pathways.txt > hsa04915/PPI.ncbi.pathways.hsa04915.txt
grep "hsa04915" sif-isolatedNetworks/MN.ncbi.pathways.txt > hsa04915/MN.ncbi.pathways.hsa04915.txt
grep "hsa04915" sif-isolatedNetworks/RN.ncbi.pathways.txt  > hsa04915/RN.ncbi.pathways.hsa04915.txt


cat hsa04915/header-all-paths.txt hsa04915/AllKEGG.ncbi.filter.pathway.hsa04915.txt > hsa04915/AllKEGG.ncbi.filter.pathway.hsa04915.header.txt
cat hsa04915/header-all-paths.txt hsa04915/MN.ncbi.pathways.hsa04915.txt > hsa04915/MN.ncbi.pathways.hsa04915.header.txt
cat hsa04915/header-all-paths.txt hsa04915/RN.ncbi.pathways.hsa04915.txt > hsa04915/RN.ncbi.pathways.hsa04915.header.txt
cat hsa04915/header-all-paths.txt hsa04915/PPI.ncbi.pathways.hsa04915.txt > hsa04915/PPI.ncbi.pathways.hsa04915.header.txt



####### MODULES PLOTS
Rscript bin/analisis-comunidad.R

####### PLOTS TOPOLOGICAL PROPERTIES METABOLISM NETWORK AND GCC 
Rscript bin/Analisis-Propiedades-Red.R

####### POWER LAW FIT METABOLISM NETWORK AND GCC
python3.7 bin/power_law_fit.py --network sif-allKEGG/AllKEGG.ncbi.filter.header.txt.Comunidades.node.csv
python3.7 bin/power_law_fit.py --network sif-allKEGG/AllKEGG.ncbi.filter.header.BigComponent.nodes.csv


####### PERCOLATION ANALYSIS WAS DONE ON A SERVER
####### BINARIES ARE FOUND ON percolacion/bin 
####### mkfiles ARE FOUND ON percolacion/mk

####### PERCOLATION PLOTS
####### PERCOLATION EDGES
Rscript bin/estadisticias-percolacion.R percolacion/results-edges percolacion/percolacion-edges.pdf edges_rep_ 48000 100

####### PERCOLATION NODES DEGREE
Rscript bin/estadisticias-percolacion.R percolacion/results-nodesDegree percolacion/percolacion-nodesDegree.pdf nodes_rep_ 6060 20

####### PERCOLATION NODES RANDOM
Rscript bin/estadisticias-percolacion-random.R percolacion/results-nodesRandom percolacion/percolacion-nodesRandom.pdf nodes_rep_ 6750 50 30





####### GENERATE BIOLOGICAL RELEVANT NETWORKS
####### PROTEIN-PROTEIN INTERACTIO NETWORK: GR + PP
####### REGULATORY NETWORK: PD
####### METABOLIC NETWORK: RT, ONLY KEEP COMPOUNDS, ADD REVERSIBLE REACTIONS

mkdir sif-isolatedNetworks

python3.7 bin/Isolate_Network_PPI.py --sif sif-allKEGG/AllKEGG.ncbi.txt --out sif-isolatedNetworks/PPI.ncbi.pathways.txt
##ADD HEADER
#NODEA	INTERACTION	NODEB	NPATHWAYS	PATHWAY	ALL-PATH-LIST
## ANALYZE NETWORK - UNDIRECTED
## EXPORT NODES AND EDGES FILE

python3.7 bin/Isolate_Network_RN.py --sif sif-allKEGG/AllKEGG.ncbi.txt --out sif-isolatedNetworks/RN.ncbi.pathways.txt
##ADD HEADER
#NODEA	INTERACTION	NODEB	NPATHWAYS	PATHWAY	ALL-PATH-LIST
## ANALYZE NETWORK - DIRECTED
## EXPORT NODES AND EDGES FILE


python3.7 bin/Isolate_Network_MetabolicNet.py --dict diccionario/diccionario_ncbi_ensembl.txt --sif sif-allKEGG/AllKEGG.ncbi.txt --out sif-isolatedNetworks/MN.ncbi.pathways.txt
##ADD HEADER
#NODEA	INTERACTION	NODEB	NPATHWAYS	PATHWAY	ALL-PATH-LIST
## ANALYZE NETWORK
## EXPORT NODES AND EDGES FILE


####### POWER LAW FOR NON-ISOLATED NODES
python3.7 bin/power_law_fit.py --network sif-isolatedNetworks/PPI.ncbi.pathways.nodes.csv
python3.7 bin/power_law_fit.py --network sif-isolatedNetworks/PPI.ncbi.pathways.BigComponent.nodes.csv

python3.7 bin/power_law_fit_directed.py --network sif-isolatedNetworks/MN.ncbi.pathways.nodes.csv
python3.7 bin/power_law_fit_directed.py --network sif-isolatedNetworks/MN.ncbi.pathways.BigComponent.nodes.csv

python3.7 bin/power_law_fit_directed.py --network sif-isolatedNetworks/RN.ncbi.pathways.nodes.csv
python3.7 bin/power_law_fit_directed.py --network sif-isolatedNetworks/RN.ncbi.pathways.BigComponent.nodes.csv


####### MODULES

####### GENERATE NET FILES
mkdir comunidades-2level-withNames-isolatedNetworks

python bin/FromSIF_To_NET.py -s sif-isolatedNetworks/PPI.ncbi.pathways.txt -o comunidades-2level-withNames-isolatedNetworks/PPI.ncbi.pathways.net
python bin/FromSIF_To_NET.py -s sif-isolatedNetworks/MN.ncbi.pathways.txt -o comunidades-2level-withNames-isolatedNetworks/MN.ncbi.pathways.net
python bin/FromSIF_To_NET.py -s sif-isolatedNetworks/RN.ncbi.pathways.txt -o comunidades-2level-withNames-isolatedNetworks/RN.ncbi.pathways.net


####### RUN INFOMAP
mkdir comunidades-2level-withNames-isolatedNetworks/PPI/ 
/Users/laura/Desktop/Infomap/Infomap comunidades-2level-withNames-isolatedNetworks/PPI.ncbi.pathways.net comunidades-2level-withNames-isolatedNetworks/PPI/ --two-level --undirected  --include-self-links  --num-trials 1000 --seed 1 --clu --map --tree --bftree
python bin/Escribir-archivos-nodos-edges.py --map comunidades-2level-withNames-isolatedNetworks/PPI/PPI.ncbi.pathways.map --kegg sif-isolatedNetworks/PPI.ncbi.pathways.txt --out-nodes comunidades-2level-withNames-isolatedNetworks/PPI/PPI.ncbi.pathways.nodes --out-edges comunidades-2level-withNames-isolatedNetworks/PPI/PPI.ncbi.pathways.edges
######## IMPORT  NODE FILE TO CYTOSCAPE  
######## SELECT MODULE NUMBER AS COLOR LABEL 

mkdir comunidades-2level-withNames-isolatedNetworks/MN/
/Users/laura/Desktop/Infomap/Infomap comunidades-2level-withNames-isolatedNetworks/MN.ncbi.pathways.net comunidades-2level-withNames-isolatedNetworks/MN/ --two-level --directed  --include-self-links  --num-trials 1000 --seed 1 --clu --map --tree --bftree
python bin/Escribir-archivos-nodos-edges.py --map comunidades-2level-withNames-isolatedNetworks/MN/MN.ncbi.pathways.map --kegg sif-isolatedNetworks/MN.ncbi.pathways.txt --out-nodes comunidades-2level-withNames-isolatedNetworks/MN/MN.ncbi.pathways.nodes --out-edges comunidades-2level-withNames-isolatedNetworks/MN/MN.ncbi.pathways.edges
######## IMPORT  NODE FILE TO CYTOSCAPE  
######## SELECT MODULE NUMBER AS COLOR LABEL 

mkdir comunidades-2level-withNames-isolatedNetworks/RN/
/Users/laura/Desktop/Infomap/Infomap comunidades-2level-withNames-isolatedNetworks/RN.ncbi.pathways.net comunidades-2level-withNames-isolatedNetworks/RN/ --two-level --directed  --include-self-links  --num-trials 1000 --seed 1 --clu --map --tree --bftree
python bin/Escribir-archivos-nodos-edges.py --map comunidades-2level-withNames-isolatedNetworks/RN/RN.ncbi.pathways.map --kegg sif-isolatedNetworks/RN.ncbi.pathways.txt --out-nodes comunidades-2level-withNames-isolatedNetworks/RN/RN.ncbi.pathways.nodes --out-edges comunidades-2level-withNames-isolatedNetworks/RN/RN.ncbi.pathways.edges
######## IMPORT  NODE FILE TO CYTOSCAPE  
######## SELECT MODULE NUMBER AS COLOR LABEL 


######## TOPOLOGICAL PROPERTIES OF ISOLATED NETWORKS
Rscript bin/Analisis-Propiedades-Red-isolatedNetworks.R



######## COUNTNG NUMBER OF INERACTION TYPES
########

## PROTEIN-PROTEIN INTERACTIONS
grep "\tgr\t" sif-allKEGG/AllKEGG.ncbi.txt  | cut -f1-3 | sort -u | wc -l
grep "\tpp\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-g\t" | cut -f1-3 | sort -u | wc -l
grep "\tpp\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tgr-g\t" | cut -f1-3 | sort -u | wc -l
grep "\tpp\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-gr\t" | cut -f1-3 | sort -u | wc -l
grep "\tpp\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tgr-gr\t" | cut -f1-3 | sort -u | wc -l

## REGULATORY INTERACTIONS
grep "\tpd\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-g\t" | cut -f1-3 | sort -u | wc -l
grep "\tpd\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tgr-g\t" | cut -f1-3 | sort -u | wc -l
grep "\tpd\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-gr\t" | cut -f1-3 | sort -u | wc -l
grep "\tpd\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tgr-gr\t" | cut -f1-3 | sort -u | wc -l

## ENZYMATIC REACTIONS
grep "\trt\t" sif-allKEGG/AllKEGG.ncbi.txt | grep "\treversible\t" | cut -f1-3 | sort -u | wc -l
grep "\trt\t" sif-allKEGG/AllKEGG.ncbi.txt | grep "\tirreversible\t" | cut -f1-3 | sort -u | wc -l
grep "\trt\t" sif-allKEGG/AllKEGG.ncbi.txt | cut -f1-3 | sort -u | wc -l

## METABOLIC INTERACTIONS
grep "\tpc\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-c\t" | cut -f1-3 | sort -u | wc -l
grep "\tpc\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tc-g\t" | cut -f1-3 | sort -u | wc -l
grep "\tpc\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tc-gr\t" | cut -f1-3 | sort -u | wc -l
grep "\tpc\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tgr-c\t" | cut -f1-3 | sort -u | wc -l
grep "\tpc\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tc-c\t" | cut -f1-3 | sort -u | wc -l
grep "\tpc\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-g\t" | cut -f1-3 | sort -u | wc -l
grep "\tpc\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-gr\t" | cut -f1-3 | sort -u | wc -l
grep "\tpc\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tgr-g\t" | cut -f1-3 | sort -u | wc -l

## SUCCESIVE REACTIONS
grep "\tec\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-g\t" | cut -f1-3 | sort -u | wc -l
grep "\tec\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tgr-g\t" | cut -f1-3 | sort -u | wc -l
grep "\tec\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tg-gr\t" | cut -f1-3 | sort -u | wc -l
grep "\tec\t" sif-allKEGG/AllKEGG.ncbi.txt  | grep "\tgr-gr\t" | cut -f1-3 | sort -u | wc -l

