import os
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sif", help="Dictionary of sif file", dest='sif',required=True, metavar=".txt")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="out.sif")
#parser.add_argument("-r", "--resultsPATH", help="Path for strelka results.\n Output dirs mut be named $sample_$RunType.\n Output files must have the columns NORMAL and TUMOR", dest='resultsPATH',required=True, metavar="resultsStrelkaPATH")

args = parser.parse_args()

sif_file = args.sif
out_file = args.out

out = open(out_file, 'w')

genes = []
compuestos = []

with open(sif_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        nodea = info[0]
        nodeb = info[2]

#### ALMACENAR INFORMACION SOBRE GENES Y COMPUESTOS
        match_compound = re.search('^[CG]\d{5}$', nodea)
        #match_compound = re.search('hola', nodea)

        if match_compound:
            compuestos.append(nodea)
        else:
            genes.append(nodea)

        if nodeb == "":
            continue

        match_compound = re.search('^[CG]\d{5}$', nodeb)
        #match_compound = re.search('hola', nodeb)

        if match_compound:
            compuestos.append(nodeb)
        else:
            genes.append(nodeb)


#### ACCEDER A LOS ELEMENTOS UNICOS DE LOS ARREGLOS GENES Y COMPUESTOS

for gene in set(genes):
    out.write(gene + "\tgene\n")

for compuesto in set(compuestos):
    out.write(compuesto + "\tcompuesto\n")
