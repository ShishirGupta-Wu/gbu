#!/bin/bash

## conda activate ncbi_datasets
## Following python modules are needed: pandas numpy sys matplotlib.pyplot
## please edit file pyscript_iqr_filter.py mentioned in run_iqr_filter()
## please edit thresold in file OG_filter.py in create_assessment_plot()
#conda install -c anaconda scipy
#conda install -c anaconda pandas
#pip install matplotlib

main(){
	 set_variables
#	 correct_files
#	 get_good_and_bad_clusters
#	 create_assessment_plot
#	 run_iqr_filter
	 run_count_proteins
	 }

set_variables(){
		Project_DIR=alicia_p
		Source_Folder=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/5_filter_unreliable_orthogroups
		Gene_Count_File=$Source_Folder/Input_folder/real_orthogroups_count_matrix.tsv
		ORTHOGROUP_FILE=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/3_Busco_and_Orthofinder_analysis/good_quality_proteomes/OrthoFinder/Results_Nov14/Orthogroups/Orthogroups.txt
		## set decided threshold as per cluster_loss_chart.pdf and Final_removed_cluster_counts_for_plot.txt
		threshold=50
		
}

correct_files(){
	## Coverting windows opened file into Unix format to avoid ^M error
	chmod +x Scripts/dos2unix
	Scripts/dos2unix $Gene_Count_File
	grep -v Orthogroup $Gene_Count_File > tmp
	awk -v OFS="\t" '$1=$1' tmp > tmp2
	sort tmp2 > Gene_counts_edited.txt
	printf "\n"	
	tail Gene_counts_edited.txt
	rm -rf tmp tmp2
	printf "\n"
}

get_good_and_bad_clusters(){
	cd $Source_Folder
	
		awk -v FS="\t"  '{ for(k = 2; k <= NF; k++) if($k > '$threshold') { print > "bc_'$threshold'.txt" }}' $Source_Folder/Gene_counts_edited.txt
		sort bc_$threshold.txt | uniq > bad_clusters_$threshold.txt
		rm -rf bc_$threshold.txt
		comm -23 Gene_counts_edited.txt bad_clusters_$threshold.txt > good_clusters_$threshold.txt
		
		wc -l good_clusters_$threshold.txt
		head good_clusters_$threshold.txt
		printf "\n"		
		wc -l bad_clusters_$threshold.txt
		head bad_clusters_$threshold.txt

	    printf "\nUse bad_clusters_$threshold.txt for max-median analysis to re-add few clusters!!\n"

}

create_assessment_plot(){
	printf "\n Mapping plot with python!!\n"
	cp Input_folder/real_orthogroups_count_matrix.tsv .
	mv real_orthogroups_count_matrix.tsv Gene_counts_all.txt
	
	##Set threshold in OG_filter.py script
	python3 Scripts/OG_filter.py
	
	mkdir -p Temp
	mv filtered_og_* Temp
}

run_iqr_filter(){
	##Edit orthogroup lines in the script pyscript_iqr_filter.py  before running this function
	##Adjust threshold file name for good and bad clusters accoding to set threshold
	
	printf "\n Adding few clusters from bad file to good file by applying interquartile range (IQR) filtering, if essential\n"	
	python3 Scripts/pyscript_iqr_filter.py
	
	printf "\n Result saved in Final_filtered_gene_count.tsv!!\n"
	
	mv Gene_counts_edited.txt Intermediate_files/
	
	grep -v Orthogroup Final_filtered_gene_count.tsv | cut -f1  | sort > Final_filtered_good_OG_ids.txt
	awk '{print $1}' Gene_counts_all.txt | grep -v Orthogroup > all_real_OG_Ids.txt
	comm -3 all_real_OG_Ids.txt Final_filtered_good_OG_ids.txt > Final_filtered_bad_OG_ids.txt
	
	##Adjust
	awk '{print $1}' bad_clusters_50.txt | sort > temp
    comm -3 temp Final_filtered_bad_OG_ids.txt > readded_into_good_OG_by_IQR.txt
	rm -rf temp
	
	grep -f Final_filtered_good_OG_ids.txt $ORTHOGROUP_FILE > Final_filtered_good_OGs.txt
	cut -f2 -d":" Final_filtered_good_OGs.txt | sed 's/ //' > Final_filtered_good_OGs_ClusterVenn.txt
	wc -l *.txt

}

run_count_proteins() {

	printf "\nGenerating protein counts for good orthogroups file:\n"
	##Adjust
	cp /scratch-rz/users/shg29ny/PROJECT_22/alicia_p/4_split_orthofinder_clusters/results_others/species.txt .

	IFS=$'\n' read -d '' -r -a speciesarray < species.txt
	j=0
	for i in "${speciesarray[@]}"
	do
			  echo 'processing... '"$i"
			  grep "${speciesarray[$j]}" Final_filtered_good_OGs.txt > "${speciesarray[$j]}"_species_group
			  printf "\n"
	let j++
	done 

	wc -l *_species_group | cat > species_groups_in_Final_filtered_OG.txt
	cut -f2 -d ":" Final_filtered_good_OGs.txt | sed 's/\s\+/\n/g' | sed '/^$/d' | cut -f1 -d "|" | sort | uniq -c > Final_filtered_OG_protein_counts.txt
	
	head -n4 species_groups_in_Final_filtered_OG.txt Final_filtered_OG_protein_counts.txt

}

main
