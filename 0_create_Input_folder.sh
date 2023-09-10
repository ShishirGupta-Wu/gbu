#!/bin/bash

main(){
	set_variables
	create_folder
	create_goodProteins_fasta
	create_goodProteins_CDS_fasta
	get_orthogroup_file
}

set_variables(){
	##Set orthofinder input fasta path
	SOURCE_FOLDER=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/
	ALL_PROTEIN_CONTAINING_FOLDER=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/1_sequence_editor/Output_edited_sequences_for_orthology
	ALL_CDS_CONTAINING_FOLDER=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/1_sequence_editor_CDS/Output_edited_sequences_for_orthology
	##Input_folder will be created
	PAML_SCRIPT_INPUT_FOLDER=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/Input_folder
	ORTRHOFINDER_ORTHOGROUP_FILE=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/2_busco_and_orthofinder/good_quality_proteomes/OrthoFinder/Results_Mar14/Orthogroups/Orthogroups.txt
}

create_folder(){
	cd $SOURCE_FOLDER
	mkdir -p Input_folder
}

create_goodProteins_fasta(){
	printf '%s\n' && printf '\e[1;34m%-6s\e[m' "############# Creating goodProteins.fasta #############" && printf "\n" 
	cd $ALL_PROTEIN_CONTAINING_FOLDER
	 
	if [[ ! -e goodProteins.fasta ]]; then
		cat *.fa > goodProteins.fasta
	fi

	PROTEIN_COUNTS="$(grep -c "^>" goodProteins.fasta)"
	printf "\n"
	echo "$PROTEIN_COUNTS"

	echo goodProteins.fasta contains total "$PROTEIN_COUNTS" PROTEIN sequences; printf "\n"
	cp goodProteins.fasta $PAML_SCRIPT_INPUT_FOLDER
}

create_goodProteins_CDS_fasta(){
	printf '%s\n' && printf '\e[1;34m%-6s\e[m' "############# Creating goodProteins_CDS.fasta #############" && printf "\n" 
	cd $ALL_CDS_CONTAINING_FOLDER
	 
	if [[ ! -e goodProteins_CDS.fasta ]]; then
		cat *.fa > goodProteins_CDS.fasta
	fi

	CDS_COUNTS="$(grep -c "^>" goodProteins_CDS.fasta)"
	printf "\n"
	echo "$CDS_COUNTS"

	echo goodProteins_CDS.fasta contains total "$CDS_COUNTS" CDS sequences; printf "\n"
	cp goodProteins_CDS.fasta $PAML_SCRIPT_INPUT_FOLDER
}

get_orthogroup_file(){
	printf '%s\n' && printf '\e[1;34m%-6s\e[m' "############# Getting OrthoFinder Orthogroup file, part of file looks like #############" && printf "\n" 
	cp $ORTRHOFINDER_ORTHOGROUP_FILE $PAML_SCRIPT_INPUT_FOLDER
	printf "\n"
	tail $PAML_SCRIPT_INPUT_FOLDER/Orthogroups.txt
	
	printf '\n%s\n' && printf '\e[1;34m%-6s\e[m' "############# Here is the contents of input folder #############" && printf "\n" 
	cd $PAML_SCRIPT_INPUT_FOLDER
	ls -ltrh
}	
main
