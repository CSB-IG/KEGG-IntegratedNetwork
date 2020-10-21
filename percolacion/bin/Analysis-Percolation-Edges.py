import os
import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-e", "--edges", help="edges file", dest='edges',required=True, metavar="edges.csv")
parser.add_argument("-k", "--k", help="Number of interactions",dest='k',required=True, metavar="k")
parser.add_argument("-o", "--out", help="Output dir",dest='out',required=True, metavar="Output_dir/")
#parser.add_argument("-r", "--resultsPATH", help="Path for strelka results.\n Output dirs mut be named $sample_$RunType.\n Output files must have the columns NORMAL and TUMOR", dest='resultsPATH',required=True, metavar="resultsStrelkaPATH")

args = parser.parse_args()

edges_file = args.edges
k = int(args.k)
dir_out = args.out
#out = open(out_file, 'w')

edges_between = []
edges_path = []
edges_interaction = []



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

def get_red_vias(interactions, path, rep, out_dir):
    out_file =  out_dir + "/edges_rep_" + str(rep) + ".txt"
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
				#print(interaction_a + "\t" + interaction_b + "\t" + id_vias)
				#print(pathway_a + "\t" + pathway_b  + "\t" + id_vias)

    for id_vias in vias_connected:
        vias = id_vias.split("-")
        out.write(vias[0] + "\t" + vias[1] + "\n")


    out.close()

############################
############################
############################


####### ALMACENAR LOS VALORES DE EDGE BETWEENNESS
####### ALMACENAR LA INFORMACION DE LAS VIAS POR EDGE

with open(edges_file) as handle:
    lines = handle.readlines()
    for line in lines:
        if "EdgeBetweenness" in line:
            continue

        line = line.rstrip()
        info = line.split('","')
        between = float(info[2].replace('"', ""))
        path = info[1].replace('"', "")
        interaction = info[9].replace('"', "")

        if path.endswith(","):
            path = path[:-1]


        edges_between.append(between)
        edges_path.append(path)
        edges_interaction.append(interaction)

print("READ EDGES")        

edges_between_sorted = sorted(edges_between, reverse=True)

thr_edge = edges_between_sorted[k]
new_interaction = []
new_path = []
for x in range(0, len(edges_between)):
	if edges_between[x] <= thr_edge:
		new_interaction.append(edges_interaction[x])
		new_path.append(edges_path[x])
#print(edges_between_sorted)
#print(edges_path_sorted)
#print(edges_interaction_sorted)

print("SORT EDGES")

###### ELIMINAR THE TOP X INTERACCIONES
#print (k)
#print (edges_between_sorted[1:k])
#print (edges_path_sorted[1:k])
#print (edges_interaction_sorted[1:k])

#new_interacion = edges_interaction_sorted[k:len(edges_interaction_sorted)]
#new_path = edges_path_sorted[k:len(edges_path_sorted)]

print(len(edges_path))
print(len(edges_interaction))

print(len(new_interaction))
print(len(new_path))


#get_red_vias(edges_interaction, edges_path, k, dir_out)

get_red_vias(new_interaction, new_path, k, dir_out)



    
