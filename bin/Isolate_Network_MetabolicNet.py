import os
import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sif", help="sif file", dest='sif',required=True, metavar="sif.txt")
parser.add_argument("-d", "--dict", help="Dictionary file", dest='dict',required=True, metavar="diccionario_ncbi_ensembl.txt")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="sif_simple.txt")
#parser.add_argument("-r", "--resultsPATH", help="Path for strelka results.\n Output dirs mut be named $sample_$RunType.\n Output files must have the columns NORMAL and TUMOR", dest='resultsPATH',required=True, metavar="resultsStrelkaPATH")

args = parser.parse_args()

sif_file = args.sif
dict_file = args.dict
out_file = args.out

out = open(out_file, 'w')


###### READ DICTIONARY FILE
###### TYPE PER MOLECULE
######
molecules = {}
with open(dict_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        type = info[1]
        name = info[3]
        molecules[name] = type



nodes_pathway = {}

####### DIRECTED NETWORK
####### CREATE ID PER INTERACTION
####### ASSOCIATE PATHWAYS TO EACH ID

with open(sif_file) as handle:
    lines = handle.readlines()
    for line in lines:
        if "NODEA" in line:
            continue

        line = line.rstrip()
        info = line.split("\t")
        nodea = info[0]
        interaction  = info[1]
        nodeb = info[2]
        type_reaction = info[3]
        pathway = info[4]

        if interaction == "rt":
            
            id_interaction = nodea + ":" + nodeb

            if id_interaction not in nodes_pathway:
                nodes_pathway[id_interaction] = []

            if pathway not in nodes_pathway[id_interaction]:
                nodes_pathway[id_interaction].append(pathway)

            ###### ADD ANOTHER INTERACTIONS FOR REVERSIBLE REACTIONS

            if type_reaction == "reversible":
                id_interaction_r = nodeb + ":" + nodea

                if id_interaction_r not in nodes_pathway:
                    nodes_pathway[id_interaction_r] = []

                if pathway not in nodes_pathway[id_interaction_r]:
                    nodes_pathway[id_interaction_r].append(pathway)


####### IMPRIMIR PATHWAYS
####### CALCULAR NUMERO DE PATHWAYS

for interaction in nodes_pathway:
    nodes = interaction.split(":")
    pathway_number = len(nodes_pathway[interaction])

    all_paths = ""
    for path in nodes_pathway[interaction]:
        all_paths = all_paths + path + ","

    out.write(nodes[0] + "\tinteraction\t" + nodes[1] + "\t" + str(pathway_number) )

    if pathway_number == 1:
        out.write("\t" + str(nodes_pathway[interaction][0]) + "\t" + str(nodes_pathway[interaction][0]) + "\n")
    else:
        out.write("\tNA\t" + all_paths + "\n")
