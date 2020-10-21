import os
import argparse
import re
import numpy as np
import random

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--edges", help="edges file", dest='edges',required=True, metavar="edges.csv")
parser.add_argument("-n", "--nodes", help="nodes file", dest='nodes',required=True, metavar="nodes.csv")
parser.add_argument("-k", "--k", help="Number of interactions",dest='k',required=True, metavar="k")
parser.add_argument("-s", "--seed", help="Seed for random generator",dest='seed',required=True, metavar="seed")
parser.add_argument("-o", "--out", help="Output dir",dest='out',required=True, metavar="Output_dir/")
#parser.add_argument("-r", "--resultsPATH", help="Path for strelka results.\n Output dirs mut be named $sample_$RunType.\n Output files must have the columns NORMAL and TUMOR", dest='resultsPATH',required=True, metavar="resultsStrelkaPATH")

args = parser.parse_args()

edges_file = args.edges
nodes_file = args.nodes
k = int(args.k)
seed = int(args.seed)
dir_out = args.out
#out = open(out_file, 'w')


############################
############################ ALL VS ALL
############################

def all_all(pathway, vias):
    #print(pathway)
    all_paths = pathway.split(",")
    for i in range(0, len(all_paths)):
        for j in range(i+1, len(all_paths)):
            id_vias = all_paths[i] + "-" + all_paths[j]
            id_vias2 = all_paths[j] + "-" + all_paths[i]
            if id_vias not in vias and id_vias2 not in vias and id_vias != id_vias2:
                vias[id_vias] = 1
                #print (id_vias)
    return(vias)


############################
############################ MOLECULES PER INTERACTION
############################

def molecules(interaction):
    interaction_simple = interaction.replace(" (interaction) ", "/")
    mol = interaction_simple.split("/")
    return(mol)

############################
############################ FUNCION PARA CALCULAR LA RED DE VIAS, A PARTIR DE UN VECTOR DE INTERACCIONES Y VIAS
############################

def get_red_vias(interactions, path, removed, iteration, out_dir):
    out_file =  out_dir + "/nodes_rep_" + str(iteration) + "_" + str(removed) + ".txt"
    out = open(out_file, 'w')
    vias_connected = {}
    for xa in range(0, len(interactions)):
        for xb in range(xa+1, len(interactions)):
            interaction_a = interactions[xa]
            interaction_b = interactions[xb]

            pathway_a = path[xa]
            pathway_b = path[xb]

            ######### SI HAY MAS DE UN PATHWAY POR INTERACCION CONECTAR LOS PATHWAY

            if "," in pathway_a:
                vias_connected = all_all(pathway_a, vias_connected)

            if "," in pathway_b:
                vias_connected = all_all(pathway_b, vias_connected)


            ########## CONECTAR TODOS LOS PATHWAY SI LAS INTERACCIONES COMPARTEN MOLECULAS

            mol_a_vector = molecules(interaction_a)
            mol_b_vector = molecules(interaction_b)

            for mola in mol_a_vector:
                if mola in mol_b_vector:
                    all_paths_a = pathway_a.split(',')
                    all_paths_b = pathway_b.split(',')
                    for i in range(0, len(all_paths_a)):
                        for j in range(0, len(all_paths_b)):
                            id_vias = all_paths_a[i] + "-" + all_paths_b[j]
                            id_vias2 = all_paths_b[j] + "-" + all_paths_a[i]
                            if id_vias not in vias_connected and id_vias2 not in vias_connected and id_vias != id_vias2:
                                vias_connected[id_vias] = 1

    for id_vias in vias_connected:
        vias = id_vias.split("-")
        out.write(vias[0] + "\t" + vias[1] + "\n")

    out.close()

############################
############################
############################

####### ALMACENAR LOS VALORES DE DEGREE POR CADA NODO

nodes_degree = []
nodes_names = []

with open(nodes_file) as handle:
    lines = handle.readlines()
    for line in lines:
        if "AverageShortestPathLength" in line:
            continue

        line = line.rstrip()
        info = line.split('","')
        degree = float(info[5].replace('"', ""))
        node = info[8].replace('"', "")
        nodes_names.append(node)


random.seed(seed)
random.shuffle(nodes_names)

#### SE ELIMINAN LOS PRIMEROS K NODOS
new_nodes = nodes_names[k:len(nodes_names)]

####### ALMACENAR LAS INTERACCIONES EN LAS QUE PARITCIPAN EL RESTO DE LOS NODOS EN LA RED
####### ALMACENAR LA INFORMACION DE LAS VIAS POR EDGE

edges_path = []
edges_interaction = []
with open(edges_file) as handle:
    lines = handle.readlines()
    for line in lines:
        if "EdgeBetweenness" in line:
            continue

        line = line.rstrip()
        info = line.split('","')
        path = info[1].replace('"', "")
        interaction = info[9].replace('"', "")

        if path.endswith(","):
            path = path[:-1]

        genes_interaction_array = interaction.split(" ")
        gene1 = genes_interaction_array[0]
        gene2 = genes_interaction_array[2]

        if gene1 in new_nodes and gene2 in new_nodes:
            edges_path.append(path)
            edges_interaction.append(interaction)
        

#print (edges_path)
#print (edges_interaction)
######## OBTENER LA RED DE VIAS PARA LAS INTERACCIONES SELECCIONADAS

get_red_vias(edges_interaction, edges_path, k, seed, dir_out)



    
