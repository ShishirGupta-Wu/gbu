#!/bin/bash

## run it in Cluster 201 to use 80 cores for t-coffee
## github_uploaded_scripts and preprocessing folder (i.e., output of collect_data_script_v1.sh) is needed in the work directory
## ADJUST set_variables, get_files
## ADJUST replace afum to desires species 8L in prepare_MLtreedir()
## ADJUST replace afum| to desires species 8L_ ONLY if working on branch site model in prepare_MLtreedir()
## On keeping it like 8L| sed will not work and branch will not mark with # which is required for site model

main(){
	 set_variables
#	 get_files
#	 get_scripts
#	 prepare_clusters
#	 create_sequence_files
#	 run_tcoffee
#	 process_CDS_files
#	 run_pal2nal
#	 run_gblocks
#	 fasta_alignment_to_phylip
#	 run_phyml
#	 cleanup
#	 prepare_MLtreedir
#	 prepare_phylipdir
#	 move_final_together
}

set_variables(){
	SOURCE_FOLDER=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis
	ALL_CDS_CONTAINING_FILE=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/Input_folder/goodProteins_CDS.fasta
	ALL_PROTEINS_CONTAINING_FILE=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/Input_folder/goodProteins.fasta
	## Use parsed_clusters_tail, parsed_clusters_Final or SCO_cluster_Final
	CLUSTER_TO_PROCESS=/scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/preprocessing/parsed_clusters_Final
	#### full path is needed
}

get_files(){
   sed -i 's/|/_/g' $ALL_CDS_CONTAINING_FILE
   sed -i 's/|/_/g' $ALL_PROTEINS_CONTAINING_FILE
   sed -i 's/|/_/g' $CLUSTER_TO_PROCESS
   cp $ALL_CDS_CONTAINING_FILE datafiles/
   cp $ALL_PROTEINS_CONTAINING_FILE datafiles/
   cp $CLUSTER_TO_PROCESS datafiles/
   wait
   ln -s datafiles/* .

}

get_scripts(){
    cp github_uploaded_scripts/run_tcoffee.pl bin/
    cp github_uploaded_scripts/run_pal2nal.pl bin/
    cp github_uploaded_scripts/run_gblocks.pl bin/
	cp github_uploaded_scripts/pal2nal.pl bin/
	cp github_uploaded_scripts/run_PhyML.pl bin/
	cp github_uploaded_scripts/Gblocks bin/
	cp github_uploaded_scripts/faSomeRecords bin/
	wait
	ln -s bin/* .
	
#	cp github_uploaded_scripts/ElConcatenero.py output_dir/  #failed to copy
#	cp github_uploaded_scripts/ElParsito.py output_dir/      #failed to copy
}

prepare_clusters(){
	## Check if orthogroups has group tag? Replace the Group tags of clusters ##ADJUST##
	sed -e s/'Group[^ ]*'/''/g $CLUSTER_TO_PROCESS > no_others_no_duplicates_rmv_spc_no_prefix.txt
#	sed -e s/'Group[^ ]*'/''/g parsed_clusters_tail > no_others_no_duplicates_rmv_spc_no_prefix.txt
	perl -pne 'BEGIN {$\="\n"}s/ /\n/g' no_others_no_duplicates_rmv_spc_no_prefix.txt > no_others_no_duplicates_rmv_spc_no_prefix_in_colomn.txt # Convert the rows into columns and separate the groups by blank line
	cp no_others_no_duplicates_rmv_spc_no_prefix_in_colomn.txt ids_in_all_cluster_blank.txt
	cat -s ids_in_all_cluster_blank.txt > ids_in_all_cluster.txt # replace mutliple empty lines with a single empty line in bash
	# Split the file into clusters when blank line is used as a record separator #IMPORTANT# Do it from cluster otherwise it will not generate all files
	awk '{print $0 > "clu_no"NR".out"}' FS="\n+" RS= ids_in_all_cluster.txt
	
	printf "\nSorting the required cluster sequence file name to see last file .. \n"
	ls *.out | sort -k2 -to -n  | tail 
	
	# SEE the number of cluster here for furter loops and change the loop below accordingly
	# last cluster number should be the exactly number which is the result of following command
	perl -lne 'END { print $. }' no_others_no_duplicates_rmv_spc_no_prefix.txt &>> loopend.txt  
}


create_sequence_files(){
	local loopend=$(head -n 1 loopend.txt)
	echo "Value of loopend is : $loopend"
	for ((i=1;i<="$loopend";i++)); do faSomeRecords $ALL_PROTEINS_CONTAINING_FILE clu_no$i.out clus_grp_seq$i.fa; done
	for ((i=1;i<="$loopend";i++)); do sed -i '/^$/d' clus_grp_seq$i.fa; done
	mv *.fa protein_cluster_dir/
}

run_tcoffee(){
    cd protein_alignment_dir
	sed -i -e '/[^>]/s/[U]/X/g'
	cd ..
	# Running t-coffee for creating protein alignments
	perl run_tcoffee.pl -i protein_cluster_dir
	## Move created protein alignments in previous step to protein_alignment_directory
	mv *.aln protein_alignment_dir/
	mv clus_grp_seq* protein_alignment_dir/
	
	## check in log, there should not be FATAL, if yes then re-run
}

process_CDS_files(){
	## Here CDS fasta  goodProteins_CDS.fasta is needed. Use get_datafiles function
	local loopend=$(head -n 1 loopend.txt)
	for ((i=1;i<="$loopend";i++)); do faSomeRecords $ALL_CDS_CONTAINING_FILE clu_no$i.out clus_grp_seq$i.fa; done
	for ((i=1;i<="$loopend";i++)); do sed -i '/^$/d' clus_grp_seq$i.fa; done
    mv *.fa gene_cluster_dir/
}

run_pal2nal(){
	perl run_pal2nal.pl -p protein_alignment_dir -d gene_cluster_dir -c output_dir 2>&1 | tee pal2nal.txt
}

#If you get a gblocks error saying Permission denied, just do chmod -R 755 path/to/Gblocks(folder)

run_gblocks(){
	cd output_dir
	for file in *; do mv "${file}" "${file/.codon.aln/}"; done #change the extension to .fa	
	for file in *seq*; do mv $file "$file".fa; done
	cd ..
	perl run_gblocks.pl -i output_dir
    cd output_dir/
	mv *.fa.fa ../gblocks_results/
	mv *.htm ../gblocks_results/
    cd ../gblocks_results/
	rename "s/fa.fa/fa/" *.fa.fa
	cd ..
}

fasta_alignment_to_phylip(){

	## In case of argparse error use 108: try few times
    cd output_dir
    ln -s $SOURCE_FOLDER/github_uploaded_scripts/ElConcatenero.py .
	ln -s $SOURCE_FOLDER/github_uploaded_scripts/ElParsito.py .
	ls *.py
	#rename "s/\.codon.fa/_bp.codon.fa/" *.codon.fa 
	###testing - rename "s/_bp.codon.fa/\.codon.fa/" *.codon.fa
	for i in *.fa; do python ElConcatenero.py  -if fasta -of phylip -in $i -o "$i"_phy; done
	rm -rf *.File
	rm -rf *.pyc
	for i in *.phy; do sed -i 's/_bp /_bp  /g' $i;done  ##Add one or two extra spaces to the sequence name: to overcome the PAML error: Make sure to separate the sequence from its name by 2 or more spaces. 
	echo fasta_alignment_to_phylip Analysis done! Results are saved in folder results_Phylip_Files
###Add one or two extra spaces to the sequence name: to overcome the PAML error: Make sure to separate the sequence from its name by 2 or more spaces
#	for i in *.phy; do sed -i 's/_bp/_bp  /g' $i;done 
    mv *.phy ../results_Phylip_Files
	cd ..
#After editing these phylip files, copy these into a new directory under main called phylip_dir and there you can change the extensions to .phylip and change the | in sequence names to _ (if not changed)
	cp results_Phylip_Files/* phylip_dir/
	cd phylip_dir
	for file in *; do mv "${file}" "${file/.fa_phy.phy/}"; done #change the extensions to .phylip
	for file in *; do mv $file "$file".phylip; done
#	for i in *.phylip; do sed -i 's/|/_/g' $i; done # change the | in species names to _ to avoid possible problems that can come up
    cd ..
}

run_phyml(){

if [ -e "MLtree_dir" ];then rm -rf "MLtree_dir" ; fi  

perl run_PhyML.pl -i phylip_dir -o MLtree_dir

}

cleanup(){
	mkdir temp_files
	mv nohup.out gblocks_nohup_log.txt
	mv *.txt temp_files/
	mv *.out temp_files/
	##remomve sybolic links
	rm goodProteins.fasta goodProteins_CDS.fasta $CLUSTER_TO_PROCESS run_gblocks.pl run_tcoffee.pl pal2nal.pl run_pal2nal.pl run_PhyML.pl faSomeRecords Gblocks
}

## To submit it in batch to cluster, all alignments and their gene trees should be separate like Phylip_dir1 MLtree_dir1 directions for each orthogroup. Furthermore, MLtrees should be marked and branch support values should be removed if present using this script

prepare_MLtreedir(){
#For this step, copy the MLTree outputs from PHYML to mltree_submit_temp directory you created.

cd $SOURCE_FOLDER/
    dir_name=mltree_submit_temp
    if [ -d "$dir_name" ]; then
        echo "Removing $dir_name"
        rm -rf "$dir_name"
    elif [ -f "$dir_name" ]; then
        echo "File with this name already exists, not a directory."
        exit
    fi
    if mkdir "$dir_name"; then
        echo "Clean directory created: $dir_name"
		cd $SOURCE_FOLDER/mltree_submit_temp
				cp $SOURCE_FOLDER/MLtree_dir/*.phylip.tree .
				for i in *.phylip.tree; do sed -i 's@\()-[^:]*\)\([:]\)@):@g' $i; done 
				for i in *.phylip.tree; do sed -i 's@\()[^:]*\)\([:]\)@):@g' $i; done 
				#remove branch support values if existent
				
				###ADJUST Speceies name abbreviation - see line 6-8
				for i in *.phylip.tree; do echo -e $(echo $i | sed 's/.phylip.tree//g')'\t'$(grep -o "8L" $i | grep -c .); done  > nameandinstance.txt
				mv nameandinstance.txt temp
				awk '{s=gensub("[^0-9]*([0-9][0-9]*).*","\\1",$0);print s,$0}' temp |sort -n|awk '{print $2"\t"1}' > nameandinstance.txt
								#first column OG name, second column instance number. Then use the number to copy the phylipdir files as many times as the instance.
				#This step marks each Hsap branch for each tree one by one and makes a new file for every marking event: e.g. marking first Hsap branch is called OG00*_1.phylip.tree_marked and second OG00*_2.phylip.tree_marked etc. THen move the marked tree files in a new directory called mltreemarked and change the extension to *.tree.
				
				###ADJUST Speceies name abbreviation - see line 6-8
				for i in *.phylip.tree; do echo $i; for j in $(seq 1 $(grep -o "8L" $i | grep -c .)); do sed "s@\(8L|[^:]*\)\([:]\)@\1#1\2@$j" $i > "$(echo $i | sed 's/.phylip.tree//g')"_"$j"."phylip.tree"_marked; echo $j;done; done   #to mark different hsap branches one by one and naming them OG*_1 OG*_2...
				mkdir mltreemarked
				mv *_marked mltreemarked
				cd mltreemarked
				for file in *; do mv "${file}" "${file/.phylip.tree_marked/}"; done
				for file in *; do mv $file "$file".tree; done
				#Here first step is to make a list of file names and separate them in different folders called MLtree_dirOG00*.
				for i in *.tree; do echo MLtree_dir"$i"; done > names.txt
				mv names.txt temp
				awk '{s=gensub("[^0-9]*([0-9][0-9]*).*","\\1",$0);print s,$0}' temp |sort -n|awk '{print $2}' > names.txt
				sed -i 's/clus_grp_seq//g' names.txt
				mkdir $(<names.txt cut -d. -f1)
				for i in *.tree; do name=$(echo $i|cut -d. -f1); dir=$(echo $name | sed 's/clus_grp_seq//g') ;mv "$name".* MLtree_dir"$dir"; done
				#And then, you make a list of these new folder names and match them to numbers in a second column; e.g. MLtree_dirOG0004141 to MLtree_dir1... What is important is to sort the matching table so that the order of orthogroups are same in both phylip and mltrees so they are named with the same number e.g. Phylip_dir1 and MLTree_dir1 for OG0004142.
				sed 's/\.tree//g' names.txt > foldernames.txt
				mv foldernames.txt temp
				awk '{s=gensub("[^0-9]*([0-9][0-9]*).*","\\1",$0);print s,$0}' temp |sort -n|awk '{print $2}' > sorted_foldernames.txt
				# sort foldernames.txt > sorted_foldernames.txt
				awk '{printf "%s\t%s\n", $0,"MLtree_dir"NR}' sorted_foldernames.txt > foldernames2.txt
				(while read i; do echo "mv" $(echo $i |  sed -e 's/\t/ /g'); done < foldernames2.txt) | bash
        return 0
    else
        echo "Creating directory failed: $dir_name"
        return 1
    fi 
}

prepare_phylipdir () {
#For this, copy the .phylip files in phylip_submit_temp directory you create. The namesanddistance.txt file comes from preparation of MLTree dirs for hsap_dups and it gives the number of Hsap entries in each tree, a.k.a number of duplications. Using this information, you can copy the phylip file as many as the instance for that orthogroup so at the and you'll have the same number of marked trees (one marked tree for each Hsap branch) and corresponding phylip dirs. The copied phylip files are basically the same but just the names are now OG00_1 and OG00*_2 ... OG00*_n (n=number of hsap branches) where OG00*_1 corresponds to the MLtree with first Hsap branch marked and OG00*_2 to MLtree with second Hsap branch marked etc.
    cd $SOURCE_FOLDER
    dir_name=phylip_submit_temp
    if [ -d "$dir_name" ]; then
        echo "Removing $dir_name"
        rm -rf "$dir_name"
    elif [ -f "$dir_name" ]; then
        echo "File with this name already exists, not a directory."
        exit
    fi
    if mkdir "$dir_name"; then
        echo "Clean directory created: $dir_name"
		cd $SOURCE_FOLDER/phylip_submit_temp
				cp $SOURCE_FOLDER/phylip_dir/*.phylip $SOURCE_FOLDER/phylip_submit_temp
				cp $SOURCE_FOLDER/mltree_submit_temp/nameandinstance.txt $SOURCE_FOLDER/phylip_submit_temp
				mv nameandinstance.txt temp
				awk '{s=gensub("[^0-9]*([0-9][0-9]*).*","\\1",$0);print s,$0}' temp |sort -n|awk '{print $2"\t"1}' > nameandinstance.txt
				for i in *.phylip; do og=$(echo $i|cut -d. -f1); for j in $(grep $og nameandinstance.txt | awk '{print $2}'); do cp $i "$og"_"$j".phylip; done; done
				#to copy the phylip files as many as the instance number / hsap branch number / duplication number in the orthogroup.
				mkdir $SOURCE_FOLDER/phylip_submit_temp/phylipcopied
				mv $SOURCE_FOLDER/phylip_submit_temp/*seq*_*.phylip $SOURCE_FOLDER/phylip_submit_temp/phylipcopied
				cd $SOURCE_FOLDER/phylip_submit_temp/phylipcopied
				for i in *.phylip; do echo Phylip_dir"$i"; done > names.txt
				sed -i 's/clus_grp_seq//g' names.txt
				mv names.txt temp
				awk '{s=gensub("[^0-9]*([0-9][0-9]*).*","\\1",$0);print s,$0}' temp |sort -n|awk '{print $2}' > names.txt
				mkdir $(<names.txt cut -d. -f1)
				for i in *.phylip; do name=$(echo $i|cut -d. -f1);dir=$(echo $name | sed 's/clus_grp_seq//g'); mv "$name".* Phylip_dir"$dir"; done
				#And then, you make a list of these new folder names and match them to numbers in a second column; e.g. Phylip_dirOG0004141 to Phylip_dir1... What is important is to sort the matching table so that the order of orthogroups are same in both phylip and mltrees so they are named with the same number e.g. Phylip_dir1 and MLTree_dir1 for OG0004142.
				sed 's/\.phylip//g' names.txt > foldernames.txt
				mv foldernames.txt temp
				awk '{s=gensub("[^0-9]*([0-9][0-9]*).*","\\1",$0);print s,$0}' temp |sort -n|awk '{print $2}' > sorted_foldernames.txt
				# sort foldernames.txt > sorted_foldernames.txt
				#be careful that the order is the same in MLtree too!
				awk '{printf "%s\t%s\n", $0,"Phylip_dir"NR}' sorted_foldernames.txt > foldernames2.txt
				(while read i; do echo "mv" $(echo $i |  sed -e 's/\t/ /g'); done < foldernames2.txt) | bash
        return 0
    else
        echo "Creating directory failed: $dir_name"
        return 1
    fi 
}


move_final_together () {
#here, move the phylip_dir1 etc and MLtree_dir1 in one folder to use it as input for PAML

    cd $SOURCE_FOLDER
    dir_name=submit
    if [ -d "$dir_name" ]; then
        echo "Removing $dir_name"
        rm -rf "$dir_name"
    elif [ -f "$dir_name" ]; then
        echo "File with this name already exists, not a directory."
        exit
    fi
    if mkdir "$dir_name"; then
        echo "Clean directory created: $dir_name"
		cp -r $SOURCE_FOLDER/phylip_submit_temp/phylipcopied/Phylip_dir* $SOURCE_FOLDER/submit
		cp -r $SOURCE_FOLDER/mltree_submit_temp/mltreemarked/MLtree_dir* $SOURCE_FOLDER/submit
		cp $SOURCE_FOLDER/github_uploaded_scripts/codeml_site_model_scripts/* $SOURCE_FOLDER/submit
		cd $SOURCE_FOLDER/submit
		chmod +x *.sh
		ls
        return 0
    else
        echo "Creating directory failed: $dir_name"
        return 1
    fi 
	
	# cd /scratch-rz/users/shg29ny/PROJECTS_21/Fungi_pr3/PAML_Automate_Updated_v3_reanalysis/submit
	# seq 1 107 | xargs -i sbatch ./codeml_automate_sm_sbatch_gene_tree_optimized.sh {}
	# bash combineM1-M2_v1.sh
}

main

















