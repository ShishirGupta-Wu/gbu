#!/bin/bash
## Scripts folder is needed


main(){
	 set_variables
#	 get_file
#	 correct_files
#	 filter_clusters
	 create_R_plot
	 }

set_variables(){
		Project_DIR=alicia_p
		Source_Folder=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/5_filter_unreliable_orthogroups
		Gene_Count_File=/scratch-rz/users/shg29ny/PROJECT_22/alicia_p/4_split_orthofinder_clusters/results_countmatrix/real_orthogroups_count_matrix.tsv
		## declare an array variable (cluster threshold here)
		## If threshold is 5 it means the orthogroup which has any species having more more than 6 copies of gene will be treated as bad cluster
		declare -ga thresholdarray=('5' '10' '20' '30' '40' '50' '100' '150' '200' '250' '300' '350' '400' '450' '500' '550' '600' '650' '700');
}


get_file(){

	mkdir -p $Source_Folder/Input_folder
	cd $Source_Folder/Input_folder
	cp $Gene_Count_File .
	printf "\nIf file copied correctly you can see here how the real_orthogroups_count_matrix looks\n"
	tail -n 5 real_orthogroups_count_matrix.tsv
}		
		


correct_files(){
	cd $Source_Folder
	## Coverting windows opened file into Unix format to avoid ^M error
	chmod +x Scripts/dos2unix
	Scripts/dos2unix $Source_Folder/Input_folder/real_orthogroups_count_matrix.tsv
	grep -v Orthogroup $Source_Folder/Input_folder/real_orthogroups_count_matrix.tsv > tmp
	awk -v OFS="\t" '$1=$1' tmp > tmp2
	sort tmp2 > Gene_counts_edited.txt
	rm -rf tmp tmp2
	printf "\n"
}

filter_clusters(){
	cd $Source_Folder
	mkdir Intermediate_files
	file="removed_cluster_counts"
	if [ -f "$file" ] ; then
		rm "$file"
	fi
	
	j=0
	for i in "${thresholdarray[@]}"
	do
		    echo 'processing...data with threshold '"$i"
			awk -v FS="\t"  '{ for(k = 2; k <= NF; k++) if($k > '${thresholdarray[$j]}') { print > "bc_'${thresholdarray[$j]}'.txt" }}' $Source_Folder/Gene_counts_edited.txt
			sort bc_${thresholdarray[$j]}.txt | uniq > bad_clusters_${thresholdarray[$j]}.txt
			rm -rf bc_${thresholdarray[$j]}.txt
			
			comm -23 Gene_counts_edited.txt bad_clusters_${thresholdarray[$j]}.txt > good_clusters_${thresholdarray[$j]}.txt

			perl -lne 'END { print $. }' bad_clusters_${thresholdarray[$j]}.txt > bad_cluster_count_${thresholdarray[$j]}.txt
			awk '{print $0,'${thresholdarray[$j]}'}' bad_cluster_count_${thresholdarray[$j]}.txt >> removed_cluster_counts

			sort removed_cluster_counts | uniq > temp
			awk '{print $2,"\t",$1}' temp | sort -n > Final_removed_cluster_counts_for_plot.txt
			sed -i '1i\Threshold	Deleted_clusters' Final_removed_cluster_counts_for_plot.txt
			rm -rf temp


		  printf "\n"
	let j++
	done 
		  printf "\nSee the output file for plotting\n"
	      cat Final_removed_cluster_counts_for_plot.txt
		  echo Analysis done!
		  rm -rf removed_cluster_counts
		  mv bad* Intermediate_files
		  mv good* Intermediate_files
}

create_R_plot(){
		Rscript Scripts/threshold_set.R
}


main
