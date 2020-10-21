import os
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--sif", help="sif file", dest='sif',required=True, metavar="sif.txt")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="sif_simple.txt")
#parser.add_argument("-r", "--resultsPATH", help="Path for strelka results.\n Output dirs mut be named $sample_$RunType.\n Output files must have the columns NORMAL and TUMOR", dest='resultsPATH',required=True, metavar="resultsStrelkaPATH")

args = parser.parse_args()

sif_file = args.sif
out_file = args.out

out = open(out_file, 'w')

###### CREATE DICTIONARY WITH ENZYMES AND PRODUCTS
###### USE PRODUCTIO DICTIONARY TO CONNECT WITH SUCCESIVE ENZYMES
reaction_product = {}
reaction_sustrate = {}

with open(sif_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")

        if (len(info) < 6):
            continue

        nodea = info[0]
        interaction  = info[1]
        nodeb = info[2]
        subtype = info[5]
        
        if interaction == "rt":
            if subtype == "ep":
                if nodea not in reaction_product:
                    reaction_product[nodea] = []
                reaction_product[nodea].append(nodeb)

            #### CHECK IF THE COMPOUND IS CONNECTED TO OTHER ENZYME
            if subtype == "se":
                if nodeb not in reaction_sustrate:
                    reaction_sustrate[nodeb] = []
                reaction_sustrate[nodeb].append(nodea)

#print(reaction_sustrate["UGT1A7"])
#print(reaction_product["CYP2E1"])

#### OBTAIN ALL SUCCESIVE ENZYMES FROM REACTIONS
enzymes_suc = {}
for e1 in reaction_product:
    for e1_product in reaction_product[e1]:
        for e2 in reaction_sustrate:
            for e2_sustrate in reaction_sustrate[e2]:
                if e1_product == e2_sustrate:
                    if e1 not in enzymes_suc:
                        enzymes_suc[e1] = []
                    if e2 not in enzymes_suc:
                        enzymes_suc[e2] = []
                    enzymes_suc[e1].append(e2)
                    enzymes_suc[e2].append(e1)


#### OBTAIN ALL ENZYMES WITH SHARED PRODUCTS FROM REACTIONS
for e1 in reaction_sustrate:
    for e1_sustrate in reaction_sustrate[e1]:
        for e2 in reaction_sustrate:
            for e2_sustrate in reaction_sustrate[e2]:
                if e1_sustrate == e2_sustrate and e1 != e2:
                    if e1 not in enzymes_suc:
                        enzymes_suc[e1] = []
                    if e2 not in enzymes_suc:
                        enzymes_suc[e2] = []
                    enzymes_suc[e1].append(e2)
                    enzymes_suc[e2].append(e1)


#### OBTAIN ALL ENZYMES WITH SHARED SUSTRATES FROM REACTIONS
for e1 in reaction_product:
    for e1_product in reaction_product[e1]:
        for e2 in reaction_product:
            for e2_product in reaction_product[e2]:
                if e1_product == e2_product and e1 != e2:
                    if e1 not in enzymes_suc:
                        enzymes_suc[e1] = []
                    if e2 not in enzymes_suc:
                        enzymes_suc[e2] = []
                    enzymes_suc[e1].append(e2)
                    enzymes_suc[e2].append(e1)


#print(enzymes_suc["CYP2E1"])
#print(enzymes_suc["UGT1A7"])



##### FILTER PREVIOUSLY CONSIDERED SUCCESIVE ENZYMES
with open(sif_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")

        if (len(info) < 6):
            out.write(line + "\n")
            continue

        nodea = info[0]
        interaction  = info[1]
        nodeb = info[2]
        subtype = info[5]

        if interaction == "pc":
            if subtype == "indirect effect" or subtype == "NA":
                continue
            out.write(line + "\n")
        elif interaction == "ec":
            if nodea in enzymes_suc:
                if nodeb in enzymes_suc[nodea]:
                    continue
            if nodeb in enzymes_suc:
                if nodea in enzymes_suc[nodeb]:
                    continue
                ### ONLY WILL BE WRITTEN IF NEITHER A-B OR B-A ARE SUCCESIVE ENZYMES
                out.write(line + "\n")
            else:
                out.write(line + "\n")
        else:
            out.write(line + "\n")


