import os
import argparse
import re

##### python3.7 bin/Tabla_Comunidades.py -m COMUNIDADES/output3/AllKEGG.map -k sif-allKEGG/AllKEGG.txt -o prueba-out3

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--map", help="Map community file",dest='map',required=True, metavar="Community.map")
parser.add_argument("-k", "--kegg", help="KEGG all sif file",dest='kegg',required=True, metavar="AllKEGG.ncbi.filter.txt")
parser.add_argument("-o", "--out-nodes", help="Output nodes file",dest='out_nodes',required=True, metavar="Output.nodes.txt")
parser.add_argument("-p", "--out-edges", help="Output edges file",dest='out_edges',required=True, metavar="Output.edges.txt")


args = parser.parse_args()

map_file = args.map
kegg_file = args.kegg
out_nodes_file = args.out_nodes
out_edges_file = args.out_edges

out_nodes = open(out_nodes_file, "w")
out_edges = open(out_edges_file, "w")

######### GENERAR UN DICCIONARIO CON LOS NOMBRES DE LAS COMUNIDADES
######### GENERAR UN DICCIONARIO CON TODAS LAS MOLECULAS QUE PERTENECEN A CADA COMUNIDAD
start = 0
start_comunidades = 0

comunidades_nodos = {}

with open(map_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        if "*Modules " in line:
            start_comunidades = 1
            continue
        if "*Nodes" in line: 
            start = 1
            start_comunidades = 0
            continue
        if "*Links" in line:
            break


        ######## CUANDO SE ESTAN RECORRIENDO LOS NODOS
        if start == 1 and start_comunidades == 0:
            info = line.split(" ")
            comunidad_id = info[0].split(":")[0]
            node_id = info[1].replace('"', '')

            comunidades_nodos[node_id] = comunidad_id


nodes = []
with open(kegg_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        nodea = info[0]
        nodeb = info[2]
        npath = info[3]
        path = info[4]

        if nodea in comunidades_nodos and nodeb in comunidades_nodos:
            out_edges.write(nodea + "\t" + nodeb + "\t" + npath + "\t" + path + "\n")

            if nodea not in nodes:
                out_nodes.write(nodea + "\t" + comunidades_nodos[nodea] + "\n")
                nodes.append(nodea)

            if nodeb not in nodes:
                out_nodes.write(nodeb + "\t" + comunidades_nodos[nodeb] + "\n")
                nodes.append(nodeb)

