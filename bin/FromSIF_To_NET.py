import os
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sif", help="Network sif file", dest='sif',required=True, metavar="all_kegg.unique.sif")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="out.net")
#parser.add_argument("-r", "--resultsPATH", help="Path for strelka results.\n Output dirs mut be named $sample_$RunType.\n Output files must have the columns NORMAL and TUMOR", dest='resultsPATH',required=True, metavar="resultsStrelkaPATH")

args = parser.parse_args()

sif_file = args.sif
out_file = args.out

nodes = {}

edges = []
cont = 1
cont_edges = 0

out = open(out_file, 'w')

with open(sif_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        nodea = info[0]
        nodeb = info[2]

        if "NODEA" in line:
            continue

        if nodeb == "":
            continue

	if nodea == nodeb:
		continue

        if nodea not in nodes.values():
            nodes[cont] = nodea
            cont = cont + 1

        if nodeb not in nodes.values():
            nodes[cont] = nodeb
            cont = cont + 1

        cont_nodea = nodes.keys()[nodes.values().index(nodea)]
        cont_nodeb = nodes.keys()[nodes.values().index(nodeb)]

        edges.append(str(cont_nodea) + "-" + str(cont_nodeb))
        cont_edges = cont_edges + 1

print (cont)
out.write("*Vertices " + str(cont-1) + "\n" )

for i in range(1, cont):
	out.write(str(i) + " \"" + nodes[i] + "\"\n")


out.write("*Edges " + str(cont_edges) + "\n" )

for edge in edges:
	elements = edge.split("-")
	out.write(elements[0] + " " + elements[1] + " 1.0\n")


