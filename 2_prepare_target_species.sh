#!/bin/bash

#### Currently using the filtered ENSMBL human sequences
#### For human host no filtering is required
#### Create folder Input_files and keep host.fasta and pathogen.fasta in the folder. The folder slould be in the work directory.

main(){
#	 set_variables
#	 get_scripts
#         edit_host_file
#	 edit_pathogen_file
#	 run_host_inparanoid
#         run_pathogen_inparanoid
	 clean_up
	 }


set_variables(){
#     eval "$(conda shell.bash hook)"
#     conda activate /opt/anaconda/envs/shared_env
#     countlines="perl -lne 'END { print $. }'"
#     countseq="grep -c '>'"
      dos2unix Interactome_scripts/*.sh
      chmod +x Interactome_scripts/*.sh
      chmod +x Interactome_scripts/*.pl
}


get_scripts(){
    mkdir bin
    cp Interactome_scripts/orthomclAdjustFasta bin
    wait
    ln -s bin/* .
}

edit_host_file(){
 	printf '%s\n' "############# Processing host #############"
	printf '%s\n' "Not required for human host"
	printf '%s\n' "Skipping.."
	cp Input_files/Hsap.fa .
	### ADJUST: Human is formatted do directly using it here
	cp Input_files/Hsap.fa .
#	perl orthomclAdjustFasta Hsap host.fasta 1
	sed -i 's/\./_/g' Hsap.fa
	sed -i '/^[^>]/s/[U]/X/g' Hsap.fa
	head Hsap.fa
	mv Hsap.fa inparanoid_host/
}	 
	 
edit_pathogen_file(){
 	printf '%s\n' "############# Processing pathogen #############"
	### ADJUST: Use species name as First letter of Genus and Three letters of Species
	cp Input_files/pathogen.fasta .
	perl orthomclAdjustFasta afum pathogen.fasta 1
	sed -i 's/\./_/g' afum.fasta
	sed -i '/^[^>]/s/[U]/X/g' afum.fasta
	head afum.fasta
	mv afum.fasta afum.fa
	mv afum.fa inparanoid_pathogen/
}

run_host_inparanoid(){
 	printf '%s\n' "############# Running host Inparanoid #############"
	cd inparanoid_host
	ln -s ../Interactome_scripts/iparanoid_metazoa/* .
	ln -s ../Interactome_scripts/scr_seed_new_sql.sh .
	perl inparanoid.pl Hsap.fa hosthpidb.fasta
	bash scr_seed_new_sql.sh Hsap hosthpidb

}

run_pathogen_inparanoid(){
 	printf '%s\n' "############# Running pathogen Inparanoid #############"
	cd inparanoid_pathogen
	ln -s ../Interactome_scripts/inparanoid_4_1/* .
	ln -s ../Interactome_scripts/scr_seed_new_sql.sh .
	perl inparanoid.pl afum.fa pathogenhpidb.fasta
	cp sqltable.afum.fa-pathogenhpidb.fasta sqltable.afum.fasta-pathogenhpidb.fasta
	bash scr_seed_new_sql.sh afum pathogenhpidb
}



clean_up(){
	rm -rf orthomclAdjustFasta
	rm -rf bin/
	chmod -R 700 Interactome_scripts
}



main
