#!/bin/bash
if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` [Please Enter the four letter abbrivation of species: run as script.sh cdip]"
  exit 0
fi

# for awk passing variables
#var1="$1"
#var2="$2" pathogenhpidb


printf "\n"
echo Parsing_seed for Species "$1"

# fasta extension is must 

cp sqltable."$1".fa-pathogenhpidb.fasta sqltable."$1".fasta-pathogenhpidb.fasta


# scr_seed_new_sql.sh "$1" pathogenhpidb

cp seed_sql_"$1"_pathogenhpidb/FINAL_seed_"$1"_pathogenhpidb_inparanoid_sql.txt .

cp /home/shg29ny/PROJECTS_New/HPIDB_PHISTO_NEW_DATA/hpidb_new/Human_host_interolog/seed_sql_Hsap_hosthpidb/FINAL_seed_Hsap_hosthpidb_inparanoid_sql.txt .


cat FINAL_seed_"$1"_pathogenhpidb_inparanoid_sql.txt FINAL_seed_Hsap_hosthpidb_inparanoid_sql.txt > all_vaues_to_map_0

awk {'print $2, $1'} all_vaues_to_map_0 > all_vaues_to_map

cp ../HPIDB_interaction_file_nrd.txt_notab .

cp ../mapping_ortholog_to_interaction_file.pl .



perl mapping_ortholog_to_interaction_file.pl


awk {'print $2, $3'} primary_interactome_HP > primary_interactome_HP_v1
awk '!(SEEN[$1,$2]++) && !(($2,$1) in SEEN)' primary_interactome_HP_v1 > Human_"$1"_primary_interactome_HP_v2


mkdir pathogen-human
mv all_vaues_to_map_0 pathogen-human/
mv all_vaues_to_map pathogen-human/

mv primary_interactome_HP pathogen-human/
mv primary_interactome_HP_v1 pathogen-human/
mv Human_"$1"_primary_interactome_HP_v2 pathogen-human/
