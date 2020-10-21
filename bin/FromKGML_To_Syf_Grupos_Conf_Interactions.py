import os
import argparse
import re


############### FUNCTIONS 
#########################
#relation_type + "\t" + id_dict[id2] + "\t" + subtype + "\t" + hsa
def write_interactions(synonims_function, id1, id2, relation_type, subtype, hsa, out, effect):

    ##### SOLO ESCRIBE LA RELACION CUANDO NINGUNO DE LOS DOS ES DRUGA D######
    ##### ESRCIBIR LA RELACION DE ID1 CON ID2
    match_id1 = re.search('^D\d{5}$', id1)
    match_id2 = re.search('^D\d{5}$', id2)
    if not (match_id1 or match_id2):
        out.write(id1 + "\t" + relation_type + "\t" + id2+ "\t" + subtype + "\t" + hsa + "\t" + effect + "\n")

    ##### EL ALIAS A UTILIZAR SERA EL QUE ESTE EN SINONIMOS, PROVENIENTE DEL ARCHIVO CONF
    split_id1 = id1.split(", ")
    split_id2 = id2.split(", ")

    id_id1 = None
    id_id2 = None

    for x in range(0, len(split_id1)):
        if split_id1[x] in synonims_function:
            id_id1 = split_id1[x]

    for x in range(0, len(split_id2)):
        if split_id2[x] in synonims_function:
            id_id2 = split_id2[x]
    
    ##### ESRCIBIR LAS RELACIONES DE KEY-ID1 CON CADA SINONIMO DE ID2
    if id_id2 is not None:
        for sin_id2 in synonims_function[id_id2]:
            match_id1 = re.search('^D\d{5}$', id1)
            match_id2 = re.search('^D\d{5}$', sin_id2)
            if not (match_id1 or match_id2):
                out.write(id1 + "\t" + relation_type + "\t" + sin_id2 + "\t" + subtype + "\t" + hsa + "\t" + effect+ "\n")

    ###### ESCRIBIR LAS RELACIONES DE SIN-ID1 CON KEY-ID2
    if id_id1 is not None:
        for sin_id1 in synonims_function[id_id1]:
            match_id1 = re.search('^D\d{5}$', sin_id1)
            match_id2 = re.search('^D\d{5}$', id2)
            if not (match_id1 or match_id2):
                out.write(sin_id1 + "\t" + relation_type + "\t" + id2 + "\t" + subtype + "\t" + hsa + "\t" + effect+ "\n")

    ###### ESCRIBIR LAS RELACIONES DE SIN-ID1 CON SIN-ID2
    if id_id1 is not None and id_id2 is not None:
        for sin_id1 in synonims_function[id_id1]:
            for sin_id2 in synonims_function[id_id2]:
                match_id1 = re.search('^D\d{5}$', sin_id1)
                match_id2 = re.search('^D\d{5}$', sin_id2)
                if not (match_id1 or match_id2):
                    out.write(sin_id1 + "\t" + relation_type + "\t" + sin_id2 + "\t" + subtype + "\t" + hsa + "\t" + effect+ "\n")



    



parser = argparse.ArgumentParser()
parser.add_argument("-n", "--hsa", help="Network kgml file", dest='kgml',required=True, metavar="network.kgml")
parser.add_argument("-c", "--conf", help="Configuration kgml file", dest='conf',required=True, metavar="network.conf")
parser.add_argument("-o", "--out", help="Output file",dest='out',required=True, metavar="out.sif")
#parser.add_argument("-r", "--resultsPATH", help="Path for strelka results.\n Output dirs mut be named $sample_$RunType.\n Output files must have the columns NORMAL and TUMOR", dest='resultsPATH',required=True, metavar="resultsStrelkaPATH")

args = parser.parse_args()

kgml_file = args.kgml
conf_file = args.conf
out_file = args.out

match_hsa = re.search('(hsa\d+)\.', kgml_file)
hsa = match_hsa.group(1)

print (hsa)

out = open(out_file, 'w')

id_dict={}
type_dict = {}

####### almacena los ids de los elementos usados
elements_used = {}

####### alamcena los nombres de los elementos usados , el mismo componente puede tener varios ids
names_elements_used = {}
brite= {}

found_gene = 0
found_compound = 0
found_group = 0
found_relation = 0

found_pprel = 0
found_pcrel = 0
found_ecrel = 0

found_reaction = 0

synonims = {}

with open(conf_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()
        info = line.split("\t")
        sinonym_column = info[2]
        #print(line)
        ######## SOLO LAS LINEAS DE GENES SE ANALIZAN PARA OBTENER SINONIMOS
        #### PARA CADA GENE SE GENERA UN VECTOR DE SINONIMOS

        if ("rect" in line or "line" in line) and "/dbget-bin/www_bget?hsa:" in line:
            sinonym_array = sinonym_column.split(",")
            for x in range(0, len(sinonym_array)):
                ### cada sinomimo es la llave
                match = re.search('\((.+)\)', sinonym_array[x])
                if match:
                    key = match.group(1)
                else:
                    key = sinonym_array[x].replace(" ", "")

                if key in synonims:
                    continue

                synonims[key] = []
                #print (key)
                for y in range(0, len(sinonym_array)):
                    #### el resto de los elementos se agregan como un vector asociado a la llave
                    if x != y:
                    ###### SI ESA ENTRADA NO TIENE UN NOMBRE DE MOLECULA, SALTARLA
                        match = re.search('\((.+)\)', sinonym_array[y])
                        if match:
                            key_specific = match.group(1)
                        else:
                            key_specific = sinonym_array[y].replace(" ", "")
                        synonims[key].append(key_specific)

        if "circ" in line:
            sinonym_array = sinonym_column.split(", ")
            #print (sinonym_array[0])
            match = re.search('^(\w+)', sinonym_array[0])
            key = match.group(1)
            synonims[key] = []
            for i in range(1, len(sinonym_array)):
                #print (sinonym_array[i])
                match = re.search('(\w+)', sinonym_array[i])
                key_specific = match.group(1)
                synonims[key].append(key_specific)


#print(synonims)


with open(kgml_file) as handle:
    lines = handle.readlines()
    for line in lines:
        line = line.rstrip()

        #### ALMACENAR INFORMACION SOBRE GENES
        if "entry" in line and "type=\"gene\"" in line :
        	m = re.search('id="(\d+)"', line)
        	id = m.group(1)
        	found_gene = 1
        	continue
        if "name" in line and found_gene == 1:
        	m = re.search('name="(.+)" fgcolor', line)
        	name_genes =  m.group(1)
        	name_genes_final = name_genes.replace(".", "")
        	id_dict[id] = name_genes_final
        	type_dict[id] = "gene"
        	found_gene = 0
        	continue

        #### ALMACENAR INFROMACION SOBRE COMPUESTOS
        if "entry" in line and "type=\"compound\"" in line :
        	m = re.search('id="(\d+)"', line)
        	id = m.group(1)
        	found_compound = 1
        	continue
        if "name" in line and found_compound == 1:
        	m = re.search('name="(.+)" fgcolor', line)
        	name_cpd =  m.group(1)
        	name_cpd_final = name_cpd.replace(".", "")
        	id_dict[id] = name_cpd_final
        	type_dict[id] = "cpd"
        	found_compound = 0
        	continue

        #### ALMACENAR INFROMACION SOBRE GRUPOS
        if "entry" in line and "type=\"group\"" in line :
        	m = re.search('id="(\d+)"', line)
        	id = m.group(1)
        	found_group = 1
        	continue
        if "component" in line and found_group == 1:
        	m = re.search('id="(\d+)"', line)
        	name_comp =  m.group(1)
        	if id in id_dict:
        		id_dict[id].append(name_comp)
        	else:
        		id_dict[id] = []
        		type_dict[id] = "group"
        		id_dict[id].append(name_comp)
        		continue
        if "</entry>" in line and found_group == 1:
            found_group = 0
            for x in range(0, len(id_dict[id])):
                for y in range(x+1, len(id_dict[id])):
                    write_interactions(synonims, id_dict[id_dict[id][x]], id_dict[id_dict[id][y]], "gr", "group", hsa, out, "group")
                    #out.write(id_dict[id_dict[id][x]] + "\tgr\t" + id_dict[id_dict[id][y]] + "\tgroup\t" + hsa + "\n")
                    elements_used[id_dict[id][x]] = 1
                    elements_used[id_dict[id][y]] = 1

                    #### nombres de los elementos usados
                    names_elements_used[id_dict[id_dict[id][x]]] = 1
                    names_elements_used[id_dict[id_dict[id][y]]] = 1


 		#### ALMACENAR INFROMACION SOBRE BRITES
        if "entry" in line and "type=\"brite\"" in line :
        	m = re.search('id="(\d+)"', line)
        	id = m.group(1)
        	brite[id] = 1
        	continue

        ######## OBTENER LAS RELACIONES
        ### Obtener relaciones PPrel y GErel

        ### Para las relaciones con grupos, establecer una relaciones entre todos los elementos involucrados
        if "relation" in line and ("PPrel" in line or "GErel" in line):
        	m1 = re.search('entry1="(\d+)"', line)
        	id1 = m1.group(1)
        	m2 = re.search('entry2="(\d+)"', line)
        	id2 = m2.group(1)
        	found_pprel = 1
        	subtype ="NA"
        	if "PPrel" in line:
        		relation_type = "pp"
        	else:
        		relation_type = "pd"

        	if id1 in brite or id2 in brite:
        		found_pprel = 0

        	continue
        if "subtype" in line and found_pprel == 1 and subtype == "NA":
        	m = re.search('name="([\w|/| ]+)"', line)
        	#print (line)
        	subtype = m.group(1)
        if "</relation>" in line and found_pprel == 1:
            if type_dict[id1] == "gene" and type_dict[id2] == "gene":
                write_interactions(synonims, id_dict[id1], id_dict[id2], relation_type, "g-g", hsa, out, subtype)
        		#out.write(id_dict[id1] + "\t" + relation_type + "\t" + id_dict[id2] + "\t" + subtype + "\t" + hsa + "\n")
                elements_used[id1] = 1
                elements_used[id2] = 1
                #### nombres de los elementos usados
                names_elements_used[id_dict[id1]] = 1
                names_elements_used[id_dict[id2]] = 1

            if type_dict[id1] == "gene" and type_dict[id2] == "group":
                for id_group in id_dict[id2]:
                    write_interactions(synonims, id_dict[id1], id_dict[id_group], relation_type, "g-gr", hsa, out, subtype)
        			#out.write(id_dict[id1] + "\t" + relation_type + "\t" + id_dict[id_group] + "\t" + subtype + "\t" + hsa + "\n")
                    elements_used[id1] = 1
                    elements_used[id_group] = 1
                    names_elements_used[id_dict[id1]] = 1
                    names_elements_used[id_dict[id_group]] = 1

            if type_dict[id1] == "group" and type_dict[id2] == "gene":
                for id_group in id_dict[id1]:
                    write_interactions(synonims, id_dict[id_group], id_dict[id2], relation_type, "gr-g", hsa, out, subtype)
        			#out.write(id_dict[id2] + "\t" + relation_type + "\t" + id_dict[id_group] + "\t" + subtype + "\t" + hsa + "\n")
                    elements_used[id2] = 1
                    elements_used[id_group] = 1
                    names_elements_used[id_dict[id2]] = 1
                    names_elements_used[id_dict[id_group]] = 1

            if type_dict[id1] == "group" and type_dict[id2] == "group":
                for x in range(0, len(id_dict[id1])):
                    for y in range(0, len(id_dict[id2])):
                        write_interactions(synonims, id_dict[id_dict[id1][x]], id_dict[id_dict[id2][y]], relation_type, "gr-gr", hsa, out, subtype)
        				#out.write(id_dict[id_dict[id1][x]] + "\t" + relation_type + "\t" + id_dict[id_dict[id2][y]]+ "\t" + subtype + "\t" + hsa + "\n")
                        elements_used[id_dict[id1][x]] = 1
                        elements_used[id_dict[id2][y]] = 1
                        names_elements_used[id_dict[id_dict[id1][x]]] = 1
                        names_elements_used[id_dict[id_dict[id2][y]]] = 1

            found_pprel = 0
            continue

        ### Obtener relaciones PCrel

        if "relation" in line and "PCrel" in line:
        	#print (line)
        	m1 = re.search('entry1="(\d+)"', line)
        	id1 = m1.group(1)
        	m2 = re.search('entry2="(\d+)"', line)
        	id2 = m2.group(1)
        	found_pcrel = 1
        	subtype ="NA"

        	if id1 in brite or id2 in brite:
        		found_pcrel = 0

        	continue
        if "subtype" in line and found_pcrel == 1 and subtype == "NA":
        	#print (line)
        	m = re.search('name="([\w|/| ]+)"', line)
        	subtype = m.group(1)
        if "</relation>" in line and found_pcrel == 1:
            if type_dict[id1] == "gene" and type_dict[id2] == "cpd":
                write_interactions(synonims, id_dict[id1], id_dict[id2], "pc", "g-c", hsa, out, subtype)
                #out.write(id_dict[id1] + "\tpm\t" + id_dict[id2] + "\t" + subtype + "\t" +hsa + "\n")
            if  type_dict[id1] == "cpd" and type_dict[id2] == "gene":
                write_interactions(synonims, id_dict[id1], id_dict[id2], "pc", "c-g", hsa, out, subtype)
                #out.write(id_dict[id1] + "\tmp\t" + id_dict[id2] + "\t" + subtype + "\t" +hsa + "\n")

             
            ##### ESTE TIPO DE INTERACCION PCREL ES GENE-GEGE, EN EL SIF FINAL SE CONFUNDIRA CON PPREL 
            if type_dict[id1] == "gene" and type_dict[id2] == "gene":
                write_interactions(synonims, id_dict[id1], id_dict[id2], "pc", "g-g", hsa, out, subtype)
                #out.write(id_dict[id1] + "\tpp\t" + id_dict[id2] + "\t" + subtype + "\t" +hsa + "\n")

            if type_dict[id1] == "cpd" and type_dict[id2] == "cpd":
                write_interactions(synonims, id_dict[id1], id_dict[id2], "pc", "c-c", hsa, out, subtype)
                #out.write(id_dict[id1] + "\tmm\t" + id_dict[id2] + "\t" + subtype + "\t" +hsa + "\n")

            if type_dict[id1] == "cpd" and type_dict[id2] == "group":
                for id_group in id_dict[id2]:
                    write_interactions(synonims, id_dict[id1], id_dict[id_group], "pc", "c-gr", hsa, out, subtype)
                    #out.write(id_dict[id1] + "\tmp\t" + id_dict[id_group] + "\t" + subtype + "\t" + hsa + "\n")
                    elements_used[id_group] = 1
                    names_elements_used[id_dict[id_group]] = 1

            if type_dict[id1] == "group" and type_dict[id2] == "cpd":
                for id_group in id_dict[id1]:
                    write_interactions(synonims, id_dict[id2], id_dict[id_group], "pc", "gr-c", hsa, out, subtype)
                    #out.write(id_dict[id1] + "\tmp\t" + id_dict[id_group] + "\t" + subtype + "\t" + hsa + "\n")
                    elements_used[id_group] = 1
                    names_elements_used[id_dict[id_group]] = 1

            ##### ESTE TIPO DE INTERACCION PCREL ES GENE-GEGE, EN EL SIF FINAL SE CONFUNDIRA CON PPREL 
            if type_dict[id1] == "gene" and type_dict[id2] == "group":
                for id_group in id_dict[id2]:
                    write_interactions(synonims, id_dict[id1], id_dict[id_group], "pc", "g-gr", hsa, out, subtype)
                    #out.write(id_dict[id1] + "\tpp\t" + id_dict[id_group] + "\t" + subtype + "\t" + hsa + "\n")
                    elements_used[id_group] = 1
                    names_elements_used[id_dict[id_group]] = 1

            if type_dict[id1] == "group" and type_dict[id2] == "gene":
                for id_group in id_dict[id1]:
                    write_interactions(synonims, id_dict[id2], id_dict[id_group], "pc", "gr-g", hsa, out, subtype)
                    #out.write(id_dict[id1] + "\tpp\t" + id_dict[id_group] + "\t" + subtype + "\t" + hsa + "\n")
                    elements_used[id_group] = 1
                    names_elements_used[id_dict[id_group]] = 1



            #print(id1)
            #print(id2)
            #print(id_dict[id1])
            #print(type_dict[id2])
            #print(id_dict[id2])

            elements_used[id1] = 1
            elements_used[id2] = 1

            if type_dict[id1] is not "group":
                names_elements_used[id_dict[id1]] = 1

            if type_dict[id2] is not "group":
                names_elements_used[id_dict[id2]] = 1

            #elements_used[id_dict[id1]] = 1
           # elements_used[id_dict[id2]] = 1

            found_pcrel = 0
            continue


        ### Obtener relaciones ECrel
        if "relation" in line and "ECrel" in line:
            #print (line)
            m1 = re.search('entry1="(\d+)"', line)
            id1 = m1.group(1)
            m2 = re.search('entry2="(\d+)"', line)
            id2 = m2.group(1)
            found_ecrel = 1
            subtype = "NA"
            relation_type = "ec"

            if id1 in brite or id2 in brite:
                found_ecrel = 0

            continue

        if "subtype" in line and found_ecrel == 1 and subtype == "NA":
            #print (line)
            m = re.search('name="([\w|/| ]+)"', line)
            subtype = m.group(1)
        if "</relation>" in line and found_ecrel == 1:


            if type_dict[id1] == "gene" and type_dict[id2] == "gene":
                write_interactions(synonims, id_dict[id1], id_dict[id2], relation_type, "g-g", hsa, out, subtype)
                #out.write(id_dict[id1] + "\t" + relation_type + "\t" + id_dict[id2] + "\t" + subtype + "\t" + hsa + "\n")
                elements_used[id1] = 1
                elements_used[id2] = 1
                names_elements_used[id_dict[id1]] = 1
                names_elements_used[id_dict[id2]] = 1

            if type_dict[id1] == "gene" and type_dict[id2] == "group":
                for id_group in id_dict[id2]:
                    write_interactions(synonims, id_dict[id1], id_dict[id_group], relation_type, "g-gr", hsa, out, subtype)
                    #out.write(id_dict[id1] + "\t" + relation_type + "\t" + id_dict[id_group] + "\t" + subtype + "\t" + hsa + "\n")
                    elements_used[id1] = 1
                    elements_used[id_group] = 1
                    names_elements_used[id_dict[id1]] = 1
                    names_elements_used[id_dict[id_group]] = 1


            if type_dict[id1] == "group" and type_dict[id2] == "gene":
                for id_group in id_dict[id1]:
                    write_interactions(synonims, id_dict[id2], id_dict[id_group], relation_type, "gr-g", hsa, out, subtype)
                    #out.write(id_dict[id2] + "\t" + relation_type + "\t" + id_dict[id_group] + "\t" + subtype + "\t" + hsa + "\n")
                    elements_used[id2] = 1
                    elements_used[id_group] = 1
                    names_elements_used[id_dict[id2]] = 1
                    names_elements_used[id_dict[id_group]] = 1


            if type_dict[id1] == "group" and type_dict[id2] == "group":
                for x in range(0, len(id_dict[id1])):
                    for y in range(0, len(id_dict[id2])):
                        write_interactions(synonims, id_dict[id_dict[id1][x]], id_dict[id_dict[id2][y]], relation_type, "gr-gr", hsa, out, subtype)
                        #out.write(id_dict[id_dict[id1][x]] + "\t" + relation_type + "\t" + id_dict[id_dict[id2][y]]+ "\t" + subtype + "\t" + hsa + "\n")
                        elements_used[id_dict[id1][x]] = 1
                        elements_used[id_dict[id2][y]] = 1
                        names_elements_used[id_dict[id_dict[id1][x]]] = 1
                        names_elements_used[id_dict[id_dict[id2][y]]] = 1


            found_ecrel = 0
            continue

        #### OBTENER REACCIONES

        if "<reaction " in line and "</reaction>" not in line:
            #print(line)
            m = re.search('id="(\d+)"', line)
            enzyme = m.group(1)
            m = re.search('type="(\w+)"', line)
            subtype = m.group(1)
            substrates = []
            products = []
            found_reaction = 1
            continue
        if found_reaction == 1 and "</reaction>" not in line:
            m = re.search('substrate id="(\d+)"', line)
            if m:
                substrates.append(m.group(1))
            m = re.search('product id="(\d+)"', line)
            if m:
                products.append(m.group(1))
            continue
        if "</reaction>" in line and found_reaction == 1:
            for participant in substrates:
                write_interactions(synonims, id_dict[participant], id_dict[enzyme], "rt", subtype, hsa, out, "se")
        		#out.write(id_dict[enzyme] + "\trt\t" + id_dict[participant] + "\t" + subtype + "\t" +hsa + "\n")
                if subtype == "reversible":
                    write_interactions(synonims, id_dict[enzyme], id_dict[participant], "rt", subtype, hsa, out, "ep")

                elements_used[enzyme] = 1
                elements_used[participant] = 1
                names_elements_used[id_dict[enzyme]] = 1
                names_elements_used[id_dict[participant]] = 1

            for participant in products:
                write_interactions(synonims, id_dict[enzyme], id_dict[participant], "rt", subtype, hsa, out, "ep")
                #out.write(id_dict[enzyme] + "\trt\t" + id_dict[participant] + "\t" + subtype + "\t" +hsa + "\n")
                if subtype == "reversible":
                    write_interactions(synonims, id_dict[participant], id_dict[enzyme], "rt", subtype, hsa, out, "se")
                
                elements_used[enzyme] = 1
                elements_used[participant] = 1
                names_elements_used[id_dict[enzyme]] = 1
                names_elements_used[id_dict[participant]] = 1

            found_reaction = 0
            continue


#print(elements_used)

for id_specific in id_dict:
    if type_dict[id_specific] != "group" and id_specific not in elements_used and id_dict[id_specific] not in names_elements_used:
        #print(id_specific)
        #print(id_dict[id_specific])
        match = re.search('^D\d{5}$', id_dict[id_specific])
        if not match:
            out.write(id_dict[id_specific] + "\t\t\t\t" + hsa + "\t\n")

        id_sin = id_dict[id_specific].split(",")[0]
        if id_sin in synonims:
            for sin in synonims[id_sin]:
                match = re.search('^D\d{5}$', sin)
                if not match:
                    out.write(sin + "\t\t\t\t" + hsa + "\t\n")

