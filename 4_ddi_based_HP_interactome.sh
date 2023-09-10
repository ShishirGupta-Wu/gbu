#!/bin/bash

main(){
	 set_variables
#	 getting_latest_pfam_interactions
#     annotate_domains
#	 create_ddi_pathogen_interactome
#	 create_ddi_host_interactome
	 create_ddi_HP_interactome
	 }

set_variables(){
	pathogen=afum
	host=Hsap
	Source_DIR=/home/shg29ny/hpidb_23_08_20
	PFAM_DIR=/home/shg29ny/pfam_assignment/PfamScan
}

getting_latest_pfam_interactions(){
	cd $Source_DIR
	mkdir ddi
	cd ddi
	printf '\n%s' "Getting latest interacting domains from i-Pfam database"	
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.1/database_files/pfamA_interactions.txt.gz
	gunzip pfamA_interactions.txt.gz
	head pfamA_interactions.txt
	# removing repeated interactions
	awk '!(SEEN[$1,$2]++) && !(($2,$1) in SEEN)' pfamA_interactions.txt > pfamA_interactions_1
	echo removing self-domian interactions
	awk '($1 != $2)' pfamA_interactions_1 > interactors_final.txt
}

annotate_domains(){
    cd $Source_DIR
	cp inparanoid_pathogen/"$pathogen".fa inparanoid_pathogen/"$pathogen".fasta
	cp inparanoid_pathogen/"$pathogen".fasta $PFAM_DIR
	cp inparanoid_host/"$host".fa inparanoid_host/"$host".fasta
	cp inparanoid_host/"$host".fasta $PFAM_DIR
	printf '\n%s' "Annotating domains in Pathogen Sequences"
#	perl pfam_scan.pl -fasta "$pathogen".fasta -dir ~/pfam_assignment/ > "$pathogen"_domain_result
	printf '\n%s' "Annotating domains in Host Sequences"
#	perl pfam_scan.pl -fasta "$host".fasta -dir ~/pfam_assignment/ > "$host"_domain_result	
}	

create_ddi_pathogen_interactome(){
	cd $Source_DIR/ddi
	mkdir pathogen_domain_based_ppi
	cd pathogen_domain_based_ppi
	cp $PFAM_DIR/"$pathogen"_domain_result .
	cp "$pathogen"_domain_result "$pathogen"_domain_result1
	# deleting first 28 lines
	sed -i '1,28d' "$pathogen"_domain_result1
	awk {'print $1,$6'} "$pathogen"_domain_result1 > "$pathogen"_domains.txt
	cut -f1 -d"." "$pathogen"_domains.txt > "$pathogen"_domains_final.txt
	cut -f2 -d"|" "$pathogen"_domains_final.txt > "$pathogen"_domains_final2.txt
	awk {'print $2,$1'} "$pathogen"_domains_final2.txt > "$pathogen"_domains_final3.txt
	cp "$pathogen"_domains_final3.txt domains_final3.txt
	ln -s $Source_DIR/ddi/interactors_final.txt .
	ln -s ../../Interactome_scripts/pfam_filter.pl .
	perl pfam_filter.pl
	cp pfam_assigned_to_species_pairs pfam_assigned_to_"$pathogen"_pairs
	awk '!(SEEN[$2,$3]++) && !(($3,$2) in SEEN)'  pfam_assigned_to_"$pathogen"_pairs > "$pathogen"_ddi_interactions 
	awk {'print $2,"+",$3'} "$pathogen"_ddi_interactions  > preFinal_ddi_"$pathogen"_PPI.sif
	sed 's/"$pathogen"|//g' preFinal_ddi_"$pathogen"_PPI.sif > pre2Final_interactome.sifddi
	sed 's/_/\./g' pre2Final_interactome.sifddi > Final_ddi_"$pathogen"_PPI.sif
	awk {'print $1'} Final_ddi_"$pathogen"_PPI.sif > t1
	awk {'print $3'} Final_ddi_"$pathogen"_PPI.sif > t2
	cat t1 t2 > t3
	sort t3 | uniq > Final_ddi_"$pathogen"_PPI_IDs.txt
	rm -rf t1 t2 t3
	printf '\n%s' "Looking head of pathogen interactome"
	head Final_ddi_"$pathogen"_PPI.sif
}

create_ddi_host_interactome(){
	cd $Source_DIR/ddi
	mkdir host_domain_based_ppi
	cd host_domain_based_ppi
	cp $PFAM_DIR/"$host"_domain_result .
	cp "$host"_domain_result "$host"_domain_result1
	# deleting first 28 lines
	sed -i '1,28d' "$host"_domain_result1
	awk {'print $1,$6'} "$host"_domain_result1 > "$host"_domains.txt
	cut -f1 -d"." "$host"_domains.txt > "$host"_domains_final.txt
	cut -f2 -d"|" "$host"_domains_final.txt > "$host"_domains_final2.txt
	awk {'print $2,$1'} "$host"_domains_final2.txt > "$host"_domains_final3.txt
	cp "$host"_domains_final3.txt domains_final3.txt
	ln -s $Source_DIR/ddi/interactors_final.txt .
	ln -s ../../Interactome_scripts/pfam_filter.pl .
	perl pfam_filter.pl
	cp pfam_assigned_to_species_pairs pfam_assigned_to_"$host"_pairs
	awk '!(SEEN[$2,$3]++) && !(($3,$2) in SEEN)'  pfam_assigned_to_"$host"_pairs > "$host"_ddi_interactions 
	awk {'print $2,"+",$3'} "$host"_ddi_interactions  > preFinal_ddi_"$host"_PPI.sif
	sed 's/"$host"|//g' preFinal_ddi_"$host"_PPI.sif > pre2Final_interactome.sifddi
	sed 's/_/\./g' pre2Final_interactome.sifddi > Final_ddi_"$host"_PPI.sif
	awk {'print $1'} Final_ddi_"$host"_PPI.sif > t1
	awk {'print $3'} Final_ddi_"$host"_PPI.sif > t2
	cat t1 t2 > t3
	sort t3 | uniq > Final_ddi_"$host"_PPI_IDs.txt
	rm -rf t1 t2 t3
	printf '\n%s' "Looking head of host interactome"
	head Final_ddi_"$host"_PPI.sif
}

create_ddi_HP_interactome(){
	cd $Source_DIR/ddi
	mkdir Host-Pathogen_domain_based_ppi
	cd Host-Pathogen_domain_based_ppi
	cp ../pathogen_domain_based_ppi/"$pathogen"_domains_final3.txt .
	cp ../host_domain_based_ppi/"$host"_domains_final3.txt .
	sed -i 's/ / '$pathogen'|/g' "$pathogen"_domains_final3.txt
	sed -i 's/ / '$host'|/g' "$host"_domains_final3.txt
	cat "$pathogen"_domains_final3.txt "$host"_domains_final3.txt > domains_final3.txt
	ln -s $Source_DIR/ddi/interactors_final.txt .
	ln -s ../../Interactome_scripts/pfam_filter.pl .
	perl pfam_filter.pl
	cp pfam_assigned_to_species_pairs pfam_assigned_to_"$host"_"$pathogen"_pairs
	awk '/'$host'/&&/'$pathogen'/' pfam_assigned_to_"$host"_"$pathogen"_pairs > "$host"_"$pathogen"_DDI_interactome.txt
	awk {'print $2, $3'} "$host"_"$pathogen"_DDI_interactome.txt > "$host"_"$pathogen"_DDI_interactions_only_redundant
	awk '!(SEEN[$1,$2]++) && !(($2,$1) in SEEN)'  "$host"_"$pathogen"_DDI_interactions_only_redundant > FINAL_"$host"_"$pathogen"_DDI_interactome_nonredundant
	sort  FINAL_"$host"_"$pathogen"_DDI_interactome_nonredundant > FINAL_"$host"_"$pathogen"_DDI_interactome_nonredundant_sort
	grep "^"$pathogen"" FINAL_"$host"_"$pathogen"_DDI_interactome_nonredundant_sort > FINAL_"$host"_"$pathogen"_DDI_interactome_nonredundant_sort_"$pathogen"
	grep "^"$host"" FINAL_"$host"_"$pathogen"_DDI_interactome_nonredundant_sort  > FINAL_"$host"_"$pathogen"_DDI_interactome_nonredundant_sort_"$host"
	awk '!(SEEN[$1,$2]++) && !(($2,$1) in SEEN)'  FINAL_"$host"_"$pathogen"_DDI_interactome_nonredundant_sort_"$host" > Final_HP_DDI_INTERACTOME_nonredundant
	echo "Duplicate Data Count: $(sort Final_HP_DDI_INTERACTOME_nonredundant | uniq -dc | awk '{count+=$1} END {print count}')"
	echo "Total Data Count: $(wc -l < Final_HP_DDI_INTERACTOME_nonredundant)"
	cut -f1 -d" " Final_HP_DDI_INTERACTOME_nonredundant > "$pathogen"_list
	sort "$pathogen"_list | uniq > "$pathogen"_list_unique
	cut -f2 -d" " Final_HP_DDI_INTERACTOME_nonredundant > "$host"_list
	sort "$host"_list | uniq > "$host"_list_unique
}

main
