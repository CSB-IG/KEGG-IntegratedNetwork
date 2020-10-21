import os
import argparse
import re
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sif", help="sif file", dest='sif',required=True, metavar="sif.txt")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="sif_simple.txt")
#parser.add_argument("-r", "--resultsPATH", help="Path for strelka results.\n Output dirs mut be named $sample_$RunType.\n Output files must have the columns NORMAL and TUMOR", dest='resultsPATH',required=True, metavar="resultsStrelkaPATH")

args = parser.parse_args()

sif_file = args.sif
out_file = args.out

out = open(out_file, 'w')

nodes=[]
nodes_matrix = {}
nodes_isolated = []
nodes_isolated_path = {}
nodes_connected = []



####### CREAR LISTA DE NODOS AISLADOS
####### CREAR UNA MATRIZ CON LAS INTERACCIONES Y UN VECTOR DE PATHAWAYS POR INTERACCION, MATRIX CUADRADA
####### cCONTAR EL  UMERO DE PATHWAYS EN LOS QUE EXISTE UN NODO AISLADO

with open(sif_file) as handle:
    lines = handle.readlines()
    for line in lines:
        if "NODEA" in line:
            continue

        line = line.rstrip()
        info = line.split("\t")
        nodea = info[0]
        nodeb = info[2]
        pathway = info[4]

        if nodeb == "":
            if nodea not in nodes_isolated:
                nodes_isolated.append(nodea)
                
            if nodea not in nodes_isolated_path:
                nodes_isolated_path[nodea] = []
            nodes_isolated_path[nodea].append(pathway)
            #print (nodea + " " + pathway)

        else:
            if nodea in nodes_matrix:
                if nodeb in  nodes_matrix[nodea]:
                    nodes_matrix[nodea][nodeb].append(pathway)
                else:
                    nodes_matrix[nodea][nodeb] = []
                    nodes_matrix[nodea][nodeb].append(pathway)
            else:
                nodes_matrix[nodea] = {}
                nodes_matrix[nodea][nodeb] = []
                nodes_matrix[nodea][nodeb].append(pathway)

            if nodea not in nodes:
                nodes.append(nodea)

            if nodeb not in nodes:
                nodes.append(nodeb)

            nodes_connected.append(nodea)
            nodes_connected.append(nodeb)


####### SUMAR LOS PATHWAYS PARA INTERACCIONES DEL TIPO HOLA-ADIOS Y ADIOS-HOLA

for x in range(0, len(nodes)):
    for y in range(x, len(nodes)):
        pathway_number = 0
        pathways = []
        pathways_espejo = []

        if nodes[x] in nodes_matrix:
            if nodes[y] in nodes_matrix[nodes[x]]:
                pathways = nodes_matrix[nodes[x]][nodes[y]]

        if nodes[y] in nodes_matrix:
            if nodes[x] in nodes_matrix[nodes[y]]:
                pathways_espejo = nodes_matrix[nodes[y]][nodes[x]]

        pathways_total = np.concatenate([pathways,pathways_espejo])
        pathway_number = len(set(pathways_total))

            
        if pathway_number > 0:
            all_paths = ""
            for path in set(pathways_total):
                all_paths = all_paths + path + ","

            out.write(nodes[x] + "\tinteraction\t" + nodes[y] + "\t" + str(pathway_number) )

            if pathway_number == 1:
                out.write("\t" + str(pathways_total[0]) + "\t" + str(pathways_total[0]) + "\n")
            else:
                out.write("\tNA\t" + all_paths + "\n")


###### AGREGAR LOS NODOS AISLADOS

for node in nodes_isolated:
    if node not in nodes_connected:
        pathway_number = int(len(set(nodes_isolated_path[node])))

        all_paths = ""
        for path in set(nodes_isolated_path[node]):
            all_paths = all_paths + path + ","

        out.write(node + "\t\t\t" + str(pathway_number))
        if pathway_number == 1:
            out.write("\t" + str(nodes_isolated_path[node][0]) + "\t" + str(nodes_isolated_path[node][0]) + "\n")
        else:
            out.write("\tNA\t" + all_paths + "\n")

