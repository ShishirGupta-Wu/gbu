#!/bin/bash

#### script needs one pre-existing folder: Input_files with formatted canonical fasta files
#### folder database_files will be created
#### set_variables will be ON all the times
#### for plot problems login with ssh -X wbbiIP : Not in remote PC   (like ssh -X wbbi201)
#### for plot problems if using putty please enable X11 option https://superuser.com/questions/119792/how-to-use-x11-forwarding-with-putty : Rscript may give error for busco plot
# conda activate kinfin_env2 (for new orthofinder and busco)

main(){
#	download_database
	set_variables
	run_busco
	create_busco_plot
	remove_bad_quality_proteome
	run_orthofinder
	 }
	 
download_database(){
	## Go to https://busco-data.ezlab.org/v5/data/lineages/
	## Select database according to your species and ADJUST WGET below
	mkdir database_files && cd "$_"
	wget https://busco-data.ezlab.org/v5/data/lineages/liliopsida_odb10.2020-09-10.tar.gz
	tar -xvzf liliopsida_odb10.2020-09-10.tar.gz
	rm -rf liliopsida_odb10.2020-09-10.tar.gz && ls
}

set_variables(){
	Project_DIR=monoc_p
	SOURCE_FOLDER=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/3_Busco_and_Orthofinder_analysis
	INPUT_FOLDER=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/2_Sequence_editor/Output_edited_sequences_for_orthology
	MY_BUSCO_DATABASE=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/3_Busco_and_Orthofinder_analysis/database_files/liliopsida_odb10
#	CONFIG_FILE=/home/shg29ny/anaconda3/envs/kinfin_env2/share/busco/config.ini
	DATABASE_NAME=liliopsida_odb10
		
}



run_busco(){
    cd $SOURCE_FOLDER
	for file in $INPUT_FOLDER/*.fa
	do
	   echo running $file 
	   tempvar=$file
	   echo $tempvar | sed 's/\.fa//g' | sed 's#/# #g' | awk '{print $(NF)}' > temp
	   species=`cat temp`
	   echo "$species"
	   
	   ##ADJUST Database
	 busco -i $file -l $MY_BUSCO_DATABASE -o assessment_$species -m prot -c 50 -f 
	done  
	rm -rf temp
   
}	 

create_busco_plot(){
    mkdir busco_result
	cp assessment_*/short_summary* busco_result/
	echo Generating plotting Rscript....
	## busco_plot has a bug and it is not generating plot so just produce script here
 generate_plot.py -wd busco_result --no_r 
    cp busco_result/busco_figure.R busco_figure_v1.R
	mkdir busco_plots
	sed -i 's/busco_result/busco_plots/g' busco_figure_v1.R
	sed -i 's/\.png/\.pdf/g' busco_figure_v1.R
	
#	conda deactivate in 201
	Rscript busco_figure_v1.R
}

remove_bad_quality_proteome(){
    mkdir -p good_quality_proteomes bad_quality_proteomes
	cd busco_result
	#ADJUST: percentage threshold of missing gene above which you think genome is bad 
    threshold=20
	for file in *.txt
	  do
		echo running $file 
		tempvar=$file
		cat $tempvar | awk 'NR==8' | sed 's/%,n:/@/g' | cut -f1 -d"@" | sed 's/%,M/@M/g' | cut -f2 -d"@" | sed 's/M://g' > temp
		Missing_gene_percentage=`cat temp`
		echo Missing gene percentage: "$Missing_gene_percentage"
		echo Percentage threshold of missing gene above which genome is bad : "$threshold"
		
			if (( $(echo $Missing_gene_percentage $threshold | awk '{if ($1 > $2) print 1;}') )); then
				echo In $file missing gene SCO are more than 20% : bad quality genome; printf "\n"
				
				##ADJUST FILE NAME FROM busco_result FOLDER
				ls $file | sed 's/short_summary\.specific\.'$DATABASE_NAME'\.assessment_//g' | sed 's/\.txt/\.fa/g' > bad_genome
				proteome_file=`cat bad_genome`
				cp $INPUT_FOLDER/$proteome_file ../bad_quality_proteomes
				printf "Copied $proteome_file to bad_quality_proteomes \n \n"
				rm -rf bad_genome
			else 
				echo $file missing genes SCO are less than 20% : good quality genome;
				ls $file | sed 's/short_summary\.specific\.'$DATABASE_NAME'\.assessment_//g' | sed 's/\.txt/\.fa/g' > good_genome
				proteome_file=`cat good_genome`
				cp $INPUT_FOLDER/$proteome_file ../good_quality_proteomes
				printf "Copied $proteome_file to good_quality_proteomes \n \n"
				rm -rf good_genome
			fi
				rm -rf temp
	   done  
	
	printf "Following proteomes have good quality: \n"
    ls ../good_quality_proteomes
	printf "\nFollowing proteomes have bad quality: \n"
	ls ../bad_quality_proteomes
}

run_orthofinder(){
#	printf "\nRunning Ortofinder on all proteomes ... \n"
#	orthofinder -f $INPUT_FOLDER -S diamond
#	printf "Orthofinder results are save in a new folder under $INPUT_FOLDER/Orthofinder \n"

	printf "\nRunning Ortofinder on all GOOD quality proteomes ... \n"
	orthofinder -f $SOURCE_FOLDER/good_quality_proteomes -S diamond
	printf "Ortofinder results are save in a new folder under good_quality_proteomes/Ortofinder \n"
}

main

