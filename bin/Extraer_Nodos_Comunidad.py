import os
import argparse
import re

##### python3.7 bin/Tabla_Comunidades.py -m COMUNIDADES/output3/AllKEGG.map -k sif-allKEGG/AllKEGG.txt -o prueba-out3

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--map", help="Map community file",dest='map',required=True, metavar="Community.map")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="Output.txt")


args = parser.parse_args()

map_file = args.map
out_dir = args.out


######### GENERAR UN DICCIONARIO CON LOS NOMBRES DE LAS COMUNIDADES
######### GENERAR UN DICCIONARIO CON TODAS LAS MOLECULAS QUE PERTENECEN A CADA COMUNIDAD
start = 0
start_comunidades = 0

comunidades_names = {}
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

        ######## CUANDO SE ESTAN RECORRIENDO LAS COMUNIDADES
        if start == 0 and start_comunidades == 1:
            info = line.split(" ")
            comunidad_id = info[0]
            gene_name = info[1].split(",")[0].replace('"', '')
            comunidades_names[comunidad_id] = gene_name
            comunidades_nodos[comunidad_id] = []

        ######## CUANDO SE ESTAN RECORRIENDO LOS NODOS
        if start == 1 and start_comunidades == 0:
            info = line.split(" ")
            comunidad_id = info[0].split(":")[0]
            node_id = info[1].replace('"', '')

            comunidades_nodos[comunidad_id].append(node_id)

       

######## IMPIRMIR ARCHIVO POR COMUNIDAD

for comunidad_id in comunidades_names:
    out_file = out_dir + comunidades_names[comunidad_id] + ".txt"
    out = open(out_file, "w")

    for gene in comunidades_nodos[comunidad_id]:
        out.write(gene + "\n")
