
# KEGG INTEGRATED NETWORK REPOSITORY

This repository hosts all code and resources needed to generate the results in the research paper "The large scale structure of human metabolism reveals resilience via extensive signaling crosstalk"

Some important directories and files are:
 - readme.sh. Commands required to generate the final integrated network of human metabolism from KEGG data and all associated results
 - AllKEGG.ncbi.filter.pathway.txt. This is the integrated network of human metabolism build from KEGG data
 - anotation. Contains resources required to do the name standarization and conversion
 - bin. Contains all source codes
 - sif-allKEGG. Contains different representations of the integrated network:
  - AllKEGG.ncbi.txt: one row per interaction, multiple rows will contain the same interactionif the interactions is present in several pathways 
  - AllKEGG.ncbi.filter.pathway.txt: one row per unique interactions, the pathways in which that interactions is present are summarized in column 6.


