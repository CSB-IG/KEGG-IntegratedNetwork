import os
import argparse
import re

##### python3.7 bin/Tabla_Comunidades.py -m COMUNIDADES/output3/AllKEGG.map -k sif-allKEGG/AllKEGG.txt -o prueba-out3

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--map", help="Map community file",dest='map',required=True, metavar="Community.map")
parser.add_argument("-k", "--kegg", help="All Kegg file",dest='kegg',required=True, metavar="AllKEGG.txt")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="Output.txt")


args = parser.parse_args()

map_file = args.map
kegg_file = args.kegg
out_file = args.out

out = open(out_file, 'w')

vias = {}
comunidades = {}

vias_total = {}
comunidades_total = {}

######## GENERAR UN VECTOR POR MOLECULAS CON TODAS LAS VIAS A LAS QUE PERTENECE
with open(kegg_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        name1 = info[0]
        name2 = info[2]
        via = info[4]

        if name2 == "":
            continue
            
        if name1 in vias:
            vias[name1].append(via)
        else:
            vias[name1] = []
            vias[name1].append(via)

        if name2 in vias:
            vias[name2].append(via)
        else:
            vias[name2] = []
            vias[name2].append(via)

        vias_total[via] = 1

######### GENERAR UN DICCIONARIO CON LA COMUNIDAD A LA QUE PERTENECE CADA MOLECULA
start = 0
start_comunidades = 0
comunidades_names = {}
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
        if start == 1 and start_comunidades == 0:
            info = line.split(" ")
            comunidad_id = info[0].split(":")[0]
            node_id = info[1].replace('"', '')
            comunidades[node_id] = comunidad_id

            comunidades_total[comunidad_id] = 1

        if start == 0 and start_comunidades == 1:
            info = line.split(" ")
            comunidad_id = info[0]
            gene_name = info[1].split(",")[0].replace('"', '')
            comunidades_names[comunidad_id] = gene_name

######## GENERAR TABLA RESUMEN
######## IMPRIMIR TODAS LAS COMUNIDADES, SEGUIDAS DE TODAS LAS VIAS (COLUMNAS)
######## RENGLONES = MOLECULAS
######## 1 = ESA MOLECULA PERTENECE A ESA COMUNIDAD O VIA
######## 0 = ESA MOLECULA NO PERTENECE A ESA COMUNIDAD O VIA

out.write("ID\t")
for comunidad in comunidades_total:
    out.write(comunidades_names[comunidad] + "\t")
 
for via_ask in vias_total:
    out.write(via_ask + "\t")
out.write("ID\n")


for molecule in comunidades:

    out.write(molecule  + "\t")

    for comunidad in comunidades_total:
        if comunidades[molecule] == comunidad:
            out.write("1\t")
        else:
            out.write("0\t")

    for via_ask in vias_total:
        if via_ask in vias[molecule]:
            out.write("1\t")
        else:
            out.write("0\t")

    out.write(molecule + "\n")


