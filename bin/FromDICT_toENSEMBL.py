import os
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sif", help="sif file", dest='sif',required=True, metavar="sif.txt")
parser.add_argument("-m", "--mol", help="Molecule type dictionary file", dest='mol',required=True, metavar="diccionario_molecule_type.txt")
parser.add_argument("-e", "--ensembl", help="ENSEMBL gene name file", dest='ensembl',required=True, metavar="biomart_ensembl.txt")
parser.add_argument("-n", "--ncbi", help="NCBI gene name file", dest='ncbi',required=True, metavar="gene_ncbi.txt")
parser.add_argument("-o", "--out", help="Dictionary per name",dest='out',required=True, metavar="dict_name.txt")

args = parser.parse_args()

sif_file = args.sif
molecule_file = args.mol
ensembl_file = args.ensembl
ncbi_file = args.ncbi

out_file = args.out

out = open(out_file, 'w')

mol_type = {}
ncbi_names = {}
ensembl = {}

###### LEER LOS NOMBRES DE LAS MOLECULAS Y TIPO DE MOLECULA DE ACUERDO A KEGG
######
with open(molecule_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        mol_type[info[0]] = info[1]


###### LEER LOS GENE SYMBOLS Y ALIASES DE ACUERDO A NCBI
######
with open(ncbi_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        if info[0] != "9606":
        	continue
        symbol = info[5].replace(".", "").replace(" ", "_")
        names = info[6].replace(".", "").replace(" ", "_")
        ncbi_names[symbol] = names

##### LEER LOS GENE SYMBOLS Y ENSEMBL ID
#####

with open(ensembl_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        ensembl_id = info[0]
        name = info[1]
        ensembl[name] = ensembl_id



for molecule in mol_type:
	if mol_type[molecule] == "compuesto":
		out.write(molecule + "\t" + mol_type[molecule] + "\t" + molecule + "\t" + molecule + "\t" + molecule + "\n")
		#print(molecule + "\t" + mol_type[molecule] + "\t" + mol_id[molecule] + "\t" + molecule + "\t" + molecule + "\n")

	else:
		gene_blank = molecule.replace(" ", "")
		gene_array = gene_blank.split(",")

		###### SI LA PRIMERA POSICION CORRESPONDE AL NOMBRE DEL GENE
		if gene_array[0] in ncbi_names:
			if gene_array[0] in ensembl:
				ensembl_id = ensembl[gene_array[0]]
			else:
				ensembl_id = "NA"

			match = 0
			if len(gene_array) > 1:
				for x in range(1, len(gene_array)):
					gene_alias = gene_array[x]
					if gene_alias in ncbi_names[gene_array[0]]:
						match = match + 1
				match_ratio = float(match)/float(len(gene_array)-1)
				if match_ratio > 0.5:
					out.write(molecule + "\t" + mol_type[molecule]  + "\t" + gene_array[0] + "\t" + gene_array[0].replace("-", "") + "\t" + ensembl_id + "\n")
					#print(molecule + "\t" + mol_type[molecule] + "\t" + mol_id[molecule] + "\t" + gene_array[0] + "\t" + ensembl_id + "\n")
				#else:
				#	print(molecule + "\t" + mol_type[molecule] + "\t" + mol_id[molecule] + "\t" + gene_array[0] + "\t" + ensembl_id + "\n")
			
			###### SI NO EXISTEN SINONIMOS
			else:
				out.write(molecule + "\t" + mol_type[molecule]  + "\t" + gene_array[0] + "\t" + gene_array[0].replace("-", "") + "\t" + ensembl_id + "\n")
				#print(molecule + "\t" + mol_type[molecule] + "\t" + mol_id[molecule] + "\t" + gene_array[0] + "\t" + ensembl_id + "\n")
	

		##### BUSCAR SI ALGUNA POSICION CORRESPONDE AL NCBI GENE SYMBOL
		else:
			found = 0
			for x in range(1, len(gene_array)):
				if gene_array[x] in ncbi_names:
					exclude = x
					if gene_array[x] in ensembl:
						ensembl_id = ensembl[gene_array[x]]
					else:
						ensembl_id = "NA"

					for y in range(0, len(gene_array)):
						if y != exclude:
							gene_alias = gene_array[y]
							if gene_alias in ncbi_names[gene_array[exclude]]:
								match = match + 1
					match_ratio = float(match)/float(len(gene_array)-1)
					if match_ratio > 0.5:
						found = 1
						out.write(molecule + "\t" + mol_type[molecule]  + "\t" + gene_array[exclude] + "\t" + gene_array[exclude].replace("-", "") + "\t" + ensembl_id + "\n")
						#print(molecule + "\t" + mol_type[molecule] + "\t" + mol_id[molecule] + "\t" + gene_array[exclude] + "\t" + ensembl_id + "\n")

			###### NOT ALIASES OR MAIN GENE NAME FOUND
			if found == 0:
				loci = "LOC" + molecule
				if loci in ncbi_names:
					if loci in ensembl:
						ensembl_id = ensembl[loci]
					else:
						ensembl_id = "NA"

					out.write(molecule + "\t" + mol_type[molecule]  + "\t" + loci + "\t" + loci.replace("-", "") + "\t" + ensembl_id + "\n")
					#print(molecule + "\t" + mol_type[molecule] + "\t" + mol_id[molecule] + "\t" + loci + "\t" + ensembl_id + "\n")					
				#else:
					#print ("NOT-FOUND")
					#print (molecule)





