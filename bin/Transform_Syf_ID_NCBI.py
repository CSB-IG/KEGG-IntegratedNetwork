import os
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sif", help="sif file", dest='sif',required=True, metavar="sif.txt")
parser.add_argument("-d", "--dict", help="Dictionary file",dest='dict',required=True, metavar="dictionary.txt")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="sif_simple.txt")

args = parser.parse_args()

sif_file = args.sif
dict_file = args.dict
out_file = args.out

out = open(out_file, 'w')

nodes = {}


with open(dict_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        nodes[info[0]] = info[3]


with open(sif_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()

        #print (line)
        if "NODEA" in line:
            out.write(line + "\n")
            continue

        #print (line)
        info = line.split("\t")
        nodea = info[0]
        int_type = info[1]
        nodeb = info[2]
        number = info[3]
        path = info[4]


        if nodeb == "":
            if nodea in nodes:
                out.write(nodes[nodea] + "\t" + int_type + "\t" + nodeb + "\t" + number + "\t" + path + "\n")
        else:
            if nodea in nodes and nodeb in nodes:
                out.write(nodes[nodea] + "\t" + int_type + "\t" + nodes[nodeb] + "\t" + number + "\t" + path +  "\n")



