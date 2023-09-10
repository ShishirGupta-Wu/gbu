#!/bin/bash

## This version of script does need "Input_folder". 
## Input_folder should contain goodProteins_CDS.fasta, goodProteins.fasta and orthoMCL's output orthogroups.txt. Will be better to keep the same name for input files as written here (optional). Set FULL Path in set_variables()
## It also want list of protein of interest from the orthogrops. Here vir_gene_list.txt is needed which should be in Input_folder.
## copy_number_variation() can be change if working on few species. This function need to be updated using automated construction of species array.

main(){
	set_variables
	color_output
	make_execulables
	process_nucleotide_sequences
	process_amino_acid_sequences
##	process_orthomcl_orthogroups
    process_orthofinder_orthogroups
	copy_number_variation 
	gene_counts
	cleanup
}

set_variables(){
##  SOURCE_FOLDER will be the folder where jobs will be executed
	SOURCE_FOLDER=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis
	ALL_CDS_CONTAINING_FILE=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/Input_folder/goodProteins_CDS.fasta
	ALL_PROTEINS_CONTAINING_FILE=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/Input_folder/goodProteins.fasta
	ORTHOGROUP_FILE=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/Input_folder/Orthogroups.txt
	REQUIRED_CLUSTERS=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/Input_folder/all_effector.txt
}

make_execulables(){
	cd github_uploaded_scripts/
	chmod +x *
	chmod 700 *
	cd ..
	
	if [[ ! -e dos2unix ]]; then
			ln -s github_uploaded_scripts/dos2unix .
	fi
	
    printf '\n'"%s\n" "${blu}Processing sequences:${end}" 
}

process_nucleotide_sequences(){
	cd $SOURCE_FOLDER/Input_folder
	sed -i 's/\./_/g' $ALL_CDS_CONTAINING_FILE
	printf "\nNucleotide sequence processed! Here is how the sequence looks\n"
	head -n 5 $ALL_CDS_CONTAINING_FILE
	cd $SOURCE_FOLDER
}

process_amino_acid_sequences(){
	cd $SOURCE_FOLDER/Input_folder
	sed -i 's/\./_/g' $ALL_PROTEINS_CONTAINING_FILE
	printf "\nAmino acid sequence processed! Here is how the sequence looks\n"
	head -n 5 $ALL_PROTEINS_CONTAINING_FILE
	cd $SOURCE_FOLDER
}

process_orthomcl_orthogroups(){
	cd $SOURCE_FOLDER/Input_folder
	sed 's/\./_/g' $ORTHOGROUP_FILE > orthogroups_edited.txt
	printf "\nOrthogroups processed! Here is how the orthogroup looks\n"
	tail -n 5 orthogroups_edited.txt
	cd $SOURCE_FOLDER
	cp $SOURCE_FOLDER/Input_folder/orthogroups_edited.txt .	
}

process_orthofinder_orthogroups(){
	cd $SOURCE_FOLDER/Input_folder
	sed 's/\./_/g' $ORTHOGROUP_FILE > orthogroups_edited.txt
	printf "\nOrthogroups processed! Here is how the orthogroup looks\n"
	tail -n 5 orthogroups_edited.txt
	cd $SOURCE_FOLDER
	cp $SOURCE_FOLDER/Input_folder/orthogroups_edited.txt .
}

copy_number_variation(){
	dos2unix $REQUIRED_CLUSTERS
	sed -i 's/\./_/g' $REQUIRED_CLUSTERS
	printf "\nExtracting clusters what you want! Please wait... .. ... \n"	
	grep -w -f $REQUIRED_CLUSTERS orthogroups_edited.txt > parsed_clusters
	printf "\nParsed_clusters!\nFollowing is proteins not found or grepped (if avalible) in orthogroups\n"
	grep -vf $REQUIRED_CLUSTERS parsed_clusters

	ls -ltrh
    ###ADJUST SPECIES NAME ABBREVIATIONS ###  
	awk '{c1=gsub(var1,x);c2=gsub(var2,x);c3=gsub(var3,x);c4=gsub(var4,x);c5=gsub(var5,x);c6=gsub(var6,x);c7=gsub(var7,x);c8=gsub(var8,x);print $1,var1"="c1,var2"="c2,var3"="c3,var4"="c4,var5"="c5,var6"="c6,var7"="c7,var8"="c8}' var1="8L" var2="BR-32"  var3="EBR"  var4="LMK"  var5="Mfg5"  var6="Mfrc123"  var7="Mfrg269"  var8="Mlax316" parsed_clusters > counts_parse_clusters.txt
	
	printf "\nCopy numbers of proteins in parsed clusters:\n"
	head -n 5 counts_parse_clusters.txt
	wc -l counts_parse_clusters.txt
}

gene_counts(){
	sed 's/[^|]//g' parsed_clusters | awk '{ print length }' > temp_gene_counts
	awk '{print $1}' parsed_clusters > temp_cluster_ids
	paste temp_cluster_ids temp_gene_counts > temp_total_gene_in_parsed_clusters
	sed 's/\s\+/ /g' temp_total_gene_in_parsed_clusters > total_gene_in_parsed_clusters.txt
	printf "\nRemoving Clusters with less than 2 species\n"
	awk '($2 > 2) {print}' total_gene_in_parsed_clusters.txt > parsed_cluster_with_2+species
	comm -23 total_gene_in_parsed_clusters.txt parsed_cluster_with_2+species > rejected_clusters.txt
	printf "\nNumber of deleted clusters\n"
	perl -lne 'END { print $. }' rejected_clusters.txt
	awk '{print $1}' parsed_cluster_with_2+species > temp_header
	sed 's/[:]//g' temp_header > parsed_clusters_id_Final
	rm -rf temp*
	grep -w -f parsed_clusters_id_Final parsed_clusters > parsed_clusters_Final
	tail -n 3 parsed_clusters_Final > parsed_clusters_tail_temp
	
    # printf "%'d\n Number of extracted clusters" $(wc -l < parsed_clusters_Final )
	Extracted_Cluster_COUNTS="$(wc -l < parsed_clusters_Final)"
	printf "\n"
	echo Number of extracted clusters: "$Extracted_Cluster_COUNTS"

	
	### access some cluster with duplications ###ADJUST### species here
#	grep 'afum=3' counts_parse_clusters.txt | tail  -n 2 | cut -f1 -d":" > temp_dup_tail_id
#	grep -w -f temp_dup_tail_id parsed_clusters > temp_parsed_clusters_2dups
#	cat parsed_clusters_tail_temp temp_parsed_clusters_2dups > parsed_clusters_tail
	### use parsed_clusters_tail to test next script
	
    ###ADJUST - REMOVE SPECIES NAME ABBREVIATIONS ### 	
	sed -e 's/://g' -e 's/8L=//g' -e 's/BR-32=//g' -e 's/EBR=//g' -e 's/LMK=//g' -e 's/Mfg5=//g' -e 's/Mfrc123=//g' -e 's/Mfrg269=//g' -e 's/Mlax316=//g' counts_parse_clusters.txt > parsed_cluster_counts_for_heatmap.txt
    
	###ADJUST - ADD FIRST LINE for HEATMAP INPUT FILE### 		
	(echo "Orthogroup 8L BR-32 EBR LMK Mfg5 Mfrc123 Mfrg269 Mlax316" && cat parsed_cluster_counts_for_heatmap.txt) > temp && mv temp parsed_cluster_counts_for_heatmap.txt
	
	printf "\nCopy numbers of proteins in parsed clusters:\n"
	head -n 5 parsed_cluster_counts_for_heatmap.txt
	
	awk '{for(i=3;i<=NF;i++)if($i!=$(i-1))next}1' parsed_cluster_counts_for_heatmap.txt > extracted_cluster_with_same_counts.txt
	printf "\nPrinting extracted cluster with same counts:\n"
	head extracted_cluster_with_same_counts.txt
	
	cut -f1 -d' ' extracted_cluster_with_same_counts.txt > temp
	grep -v -f temp parsed_cluster_counts_for_heatmap.txt > gene_loss_gain_in_extracted_clusters.txt
	printf "\nPrinting extracted cluster with possible gene gain or loss:\n"
	head gene_loss_gain_in_extracted_clusters.txt
	
	## Getting single copy ortholog cluster in extracted clusters
	awk '{for(i=3;i<=NF;i++)if($i!="1")next}1' extracted_cluster_with_same_counts.txt > SCO.txt
	printf "\nPrinting head of single copy ortholog cluster in extracted clusters:\n"
	head SCO.txt
	cut -f1 -d' ' SCO.txt > SCO_cluster_ids.txt
	wc -l SCO.txt
	
	grep -w -f SCO_cluster_ids.txt parsed_clusters_Final > SCO_cluster_Final

	
	printf "\nFile prepared for next steps. For testing the scripts use parsed_clusters_tail. If sussessful, then use parsed_clusters_Final or SCO_cluster_Final\n"
}

cleanup(){
    dir_name=preprocessing
    if [ -d "$dir_name" ]; then
        echo "Removing $dir_name"
        rm -rf "$dir_name"
    elif [ -f "$dir_name" ]; then
        echo "File with this name already exists, not a directory."
        exit
    fi
    if mkdir "$dir_name"; then
        echo "Clean directory created: $dir_name"
			mv parsed_* "$dir_name"
			mv temp* "$dir_name"
			mv *.txt "$dir_name"
			mv SCO_* "$dir_name"
			printf "\nResults moved to preprocessing folder \n"
        return 0
    else
        echo "Creating directory failed: $dir_name"
        return 1
    fi 

	printf '\n'"%s\n" "If all goes fine run script ${red}2_set_up_project_folders.sh${end}"
	printf '\n'
}

color_output(){
	red=$'\e[1;31m'
	grn=$'\e[1;32m'
	yel=$'\e[1;33m'
	blu=$'\e[1;34m'
	mag=$'\e[1;35m'
	cyn=$'\e[1;36m'
	end=$'\e[0m'
	# printf "%s\n" "Text in ${red}your text${end}, white and ${blu}your text${end}."
}


main

