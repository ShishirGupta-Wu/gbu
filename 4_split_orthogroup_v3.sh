#!/bin/bash
## Scripts folder and Input_folder is needed

main(){
	set_variables
#	get_orthogroup_file
#	process_orthofinder_orthogroups
#	create_list_of_species
#	split_orthogroups
#	split_singleton_and_inparalogs
#	run_count_matrix
#	 clean_up
}

set_variables(){
	Project_DIR=alicia_p
	SOURCE_FOLDER=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/4_split_orthofinder_clusters
	## Verify orthofinder results folder name
	ORTHOGROUP_FILE=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/3_Busco_and_Orthofinder_analysis/good_quality_proteomes/OrthoFinder/Results_Nov14/Orthogroups/Orthogroups.txt

}

get_orthogroup_file(){
	cd $SOURCE_FOLDER
	mkdir -p Input_folder
	cd Input_folder
	cp $ORTHOGROUP_FILE .
}

process_orthofinder_orthogroups(){
	cd $SOURCE_FOLDER
	chmod +x $SOURCE_FOLDER/Scripts/*
	Scripts/dos2unix Input_folder/Orthogroups.txt
	sed 's/\./_/g' Input_folder/Orthogroups.txt > orthogroups_edited.txt
	printf "\nOrthogroups processed! Here is how the orthogroup looks\n"
	tail -n 5 orthogroups_edited.txt
}

create_list_of_species(){

	printf "\nCreating protein count of each species in orthofider orthogroup file:\n"
	cut -f2 -d":" orthogroups_edited.txt > temp1
	perl -pe 's/ /\n/g' temp1 > temp2
	cut -f1 -d"|" temp2 | sort > temp3
	sed -i '/^$/d' temp3
	sort temp3 | uniq -c > initial_OG_all_species_counts.txt
	cat initial_OG_all_species_counts.txt
	printf "\n"

	printf "\nCreating list of species:\n"
	sort temp3 | uniq > species.txt
	cat species.txt
	printf "\n"; echo Total number of species: $(wc -l species.txt); printf "\n"
	
	rm -rf temp*

}

split_orthogroups(){

	printf "\nSplitting orthofider orthogroup file:\n"
	OUTFILE=(singleton_and_inparalogs.txt real_orthogroups.txt)
	: >${OUTFILE[0]} 2>${OUTFILE[1]}

	while read line
	  do
			echo "$line" >> ${OUTFILE[$(( $(grep -oE " [^|]*|" <<<$line | sort | uniq | wc -l) > 1))]}
	  done < orthogroups_edited.txt

	printf "\nOrthogroups processed! Here is how the singleton_and_inparalogs looks\n"	  
	head singleton_and_inparalogs.txt
	printf "\nHere is how the real_orthogroups looks\n"
	tail real_orthogroups.txt

	printf "\nSplitting orthofider orthogroup file using perl to cross-check :\n"	
	# It will not correctly if your species name has "-" symbol (solution- replace it to "_")
#	perl Scripts/split_orthogroup.pl
	
	
	printf "\nCounting proteins in real orthogroup file:\n"
	IFS=$'\n' read -d '' -r -a speciesarray < species.txt
	j=0
	for i in "${speciesarray[@]}"
	do
			  echo 'processing... '"$i"
			  grep "${speciesarray[$j]}" real_orthogroups.txt > "${speciesarray[$j]}"_species_group
			  printf "\n"
	let j++
	done 

	wc -l *_species_group | cat > species_groups_in_real_OG.txt
	cut -f2 -d ":" real_orthogroups.txt | sed 's/\s\+/\n/g' | sed '/^$/d' | cut -f1 -d "|" | sort | uniq -c > real_OG_protein_counts.txt

}

split_singleton_and_inparalogs(){

	printf "\nSplitting Singletons and Inparalogs in two files:\n"
	[ -e singleton.txt ] && rm singleton.txt
	[ -e inparalogs.txt ] && rm inparalogs.txt

	while read line
	do
			wdCount=$(echo $line | wc -w)
			if [[ $wdCount == 2 ]]
			then
					printf '%s\n' "$line" >> singleton.txt
			else
					printf '%s\n' "$line" >> inparalogs.txt
			fi
	done < singleton_and_inparalogs.txt
	
	cut -f1 -d "|" singleton.txt | cut -f2 -d " " | sort | uniq -c > singleton_protein_counts.txt
	cut -f1 -d "|" inparalogs.txt | cut -f2 -d " " | sort | uniq -c > inparalogs_group_counts.txt
	cut -f2 -d ":" inparalogs.txt | sed 's/\s\+/\n/g' | sed '/^$/d' | cut -f1 -d"|" | sort | uniq -c  > inparalogs_protein_counts.txt
	
	printf "\nWriting all generated counts all_counts_log file:\n"
	more *_counts.txt | cat > all_counts_log.txt
	
}

run_count_matrix() {

	printf "\nGenerating counts for real_orthogroups file:\n"	
	Scripts/countmatrix.sh real_orthogroups
	
	printf "\nGenerating counts for singleton_and_inparalogs file:\n"	
	Scripts/countmatrix.sh singleton_and_inparalogs
	
	printf "\nGenerating counts for inparalogs file:\n"	
	Scripts/countmatrix.sh inparalogs
	
	printf "\nGenerating counts for singleton file:\n"	
	Scripts/countmatrix.sh singleton
	
}

clean_up(){
		mkdir -p results_groups results_countmatrix results_others
		
		mv *_group *_log.txt *_protein_counts* inparalogs_group_counts.txt species_groups_in_real_OG.txt species.txt results_others
		mv *_count_matrix* initial_OG_all_species_counts.txt results_countmatrix
		mv orthogroups_edited.txt real_orthogroups.txt inparalogs.txt singleton_and_inparalogs.txt singleton.txt results_groups
		
		cd results_countmatrix
		
		printf "\nGenerating counts for orthofinder orthogroup file:\n"	
		grep -v Orthogroup singleton_and_inparalogs_count_matrix.tsv > temp_OG
		cat real_orthogroups_count_matrix.tsv temp_OG > Initial_ofdr_OG_count_matrix.tsv
		rm -rf temp_OG
		
}

main