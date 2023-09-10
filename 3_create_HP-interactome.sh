#!/bin/bash


main(){
	 create_interactome
	 }

create_interactome(){

	##ADJUST 
	pathogen_species=afum
	host_species=Hsap

	printf '\n%s' "Parsing seed for Species $pathogen_species and $host_species"
	cp inparanoid_pathogen/seed_sql_"$pathogen_species"_pathogenhpidb/FINAL_seed_"$pathogen_species"_pathogenhpidb_inparanoid_sql.txt .
	cp inparanoid_host/seed_sql_"$host_species"_hosthpidb/FINAL_seed_"$host_species"_hosthpidb_inparanoid_sql.txt .
	cat FINAL_seed_"$pathogen_species"_pathogenhpidb_inparanoid_sql.txt FINAL_seed_"$host_species"_hosthpidb_inparanoid_sql.txt > all_vaues_to_map_0
	awk {'print $2, $1'} all_vaues_to_map_0 > all_vaues_to_map

	ln -s Interactome_scripts/mapping_ortholog_to_interaction_file.pl .
	perl mapping_ortholog_to_interaction_file.pl
	awk {'print $2, $3'} primary_interactome_HP > primary_interactome_HP_v1
	awk '!(SEEN[$1,$2]++) && !(($2,$1) in SEEN)' primary_interactome_HP_v1 > Human_"$pathogen_species"_primary_interactome_HP_v2

	mkdir interactome_pathogen-human
	mv all_vaues_to_map_0 all_vaues_to_map primary_interactome_HP primary_interactome_HP_v1 interactome_pathogen-human/
    rm -rf Final*
	mv Human_"$pathogen_species"_primary_interactome_HP_v2 interactome_pathogen-human/
}

main