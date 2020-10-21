awk '{print "wget http://rest.kegg.jp/get/" $0 "/kgml -O " $0 ".kegg.kgml" }' pathways.names.hsa.txt > get_all_kegg.sh
sh get_all_kegg.sh