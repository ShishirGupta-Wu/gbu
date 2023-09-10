#!/bin/bash

#### Download the latest HPIDB database from https://hpidb.igbb.msstate.edu/about.html#stats locally in the same folder
#### Define the download date 23.08.2020
#### Interactome_scripts folder is needed in the work directory

main(){
#	 set_variables
	 get_scripts
     edit_downloaded_data
	 extract_host_sequences
	 extract_pathogen_sequences
	 extract_interactions
	 clean_up
	 see_stats
	 }


set_variables(){
#     eval "$(conda shell.bash hook)"
#     conda activate /opt/anaconda/envs/shared_env
#     countlines="perl -lne 'END { print $. }'"
#     countseq="grep -c '>'"
      printf '%s\n' "Noothing to set"
}


get_scripts(){
    mkdir bin
    cp Interactome_scripts/orthomclAdjustFasta bin
    wait
    ln -s bin/* .
}

edit_downloaded_data(){
 	unzip hpidb2.mitab.zip
	cp hpidb2.mitab_plus.txt hpidb2.mitab_plus_edit.txt
    sed -i 's/\ /_/g' hpidb2.mitab_plus_edit.txt 
	head hpidb2.mitab_plus_edit.txt
}	 
	 
extract_host_sequences(){
	printf '%s\n' "############# Processing host sequences #############"
	# Column P has host id and Column V has its corresponding sequence; P = 16th column; V = 22nd column
	awk {'print ">", $16, "@", $22'} hpidb2.mitab_plus_edit.txt > tmp1
	sed -i 's/\ //g' tmp1
	sed -i 's/UNIPROT_AC://g' tmp1
	printf '%s\n' "Number of host ids:"
	perl -lne 'END { print $. }' tmp1
	sort tmp1 | uniq > tmp2
	printf '%s\n' "Number of unique host ids:"
	perl -lne 'END { print $. }' tmp2
	tr @ '\n' < tmp2 > host.fasta
	perl orthomclAdjustFasta hosthpidb host.fasta 1
	printf '%s\n' "Number of unique host sequences"
	grep -c '>' hosthpidb.fasta 
	sed -i 's/\./_/g' hosthpidb.fasta 
    sed -i '/^[^>]/s/[U]/X/g' hosthpidb.fasta 
}

extract_pathogen_sequences(){
	printf '%s\n' "############# Processing host sequences #############"
	# Column P has pathogen id and Column V has its corresponding sequence; P = 17th column; V = 23nd column
	awk {'print ">", $17, "@", $23'} hpidb2.mitab_plus_edit.txt > tmp1_p
	sed -i 's/\ //g' tmp1_p
	sed -i 's/UNIPROT_AC://g' tmp1_p
	printf '%s\n' "Number of pathogen ids:"
	perl -lne 'END { print $. }' tmp1_p
	sort tmp1_p | uniq > tmp2_p
	printf '%s\n' "Number of unique pathogen ids:"
	perl -lne 'END { print $. }' tmp2_p
	tr @ '\n' < tmp2_p > pathogen.fasta
	perl orthomclAdjustFasta pathogenhpidb pathogen.fasta 1
	printf '%s\n' "Number of unique pathogen sequences"
	grep -c '>' pathogenhpidb.fasta 
	sed -i 's/\./_/g' pathogenhpidb.fasta 
    sed -i '/^[^>]/s/[U]/X/g' pathogenhpidb.fasta 
}
	 
	 
extract_interactions(){	 
	printf '%s\n' "############# Processing interactions #############"
	awk {'print "hosthpidb|", $16, "@", "pathogenhpidb|", $17'} hpidb2.mitab_plus_edit.txt > tmp_hp
	sed -i 's/UNIPROT_AC://g' tmp_hp
	sed -i 's/\ //g' tmp_hp
	printf '%s\n' "Number of interactions:"
	perl -lne 'END { print $. }' tmp_hp
	uniq -c tmp_hp | sort
	sort tmp_hp | uniq > tmp_hp2
	sed -i 's/@/ /g' tmp_hp2
	perl -p -i -e 's/ /\t/g' tmp_hp2
	printf '%s\n' "Number of unique interactions:"
	perl -lne 'END { print $. }' tmp_hp2
	cp tmp_hp2 final_hpidb_interaction_file
	sed -i 's/\./_/g' final_hpidb_interaction_file
	sed -e "s/[[:space:]]\+/ /g" final_hpidb_interaction_file > HPIDB_interaction_file_nrd.txt_notab
	sed -i 's/\./_/g' HPIDB_interaction_file_nrd.txt_notab
}

clean_up(){
	mkdir temperory_files
	mv tmp* temperory_files/	
	mkdir database_files
	mv hpidb2.mitab.txt hpidb2.mitab_plus_edit.txt hpidb2.mitab_plus.txt final_hpidb_interaction_file host.fasta pathogen.fasta database_files/
	rm -rf orthomclAdjustFasta
	rm -rf bin/
	chmod -R 700 Interactome_scripts
	mkdir inparanoid_host
	mv hosthpidb.fasta inparanoid_host/
	mkdir inparanoid_pathogen
	mv pathogenhpidb.fasta inparanoid_pathogen/
}

see_stats(){
    cd temperory_files
	printf '%s\n' "Number of host ids:"
	perl -lne 'END { print $. }' tmp1
	printf '%s\n' "Number of unique host ids:"
	perl -lne 'END { print $. }' tmp2
	printf '%s\n' "Number of unique host sequences"
	grep -c '>' hosthpidb.fasta 
	
	printf '%s\n' "Number of pathogen ids:"
	perl -lne 'END { print $. }' tmp1_p
	printf '%s\n' "Number of unique pathogen ids:"
	perl -lne 'END { print $. }' tmp2_p
	printf '%s\n' "Number of unique pathogen sequences"
	grep -c '>' pathogenhpidb.fasta
	
	printf '%s\n' "Number of interactions:"
	perl -lne 'END { print $. }' tmp_hp
	printf '%s\n' "Number of unique interactions:"
	perl -lne 'END { print $. }' tmp_hp2
}

# export -f see_stats &>> stats.txt
# result=$(see_stats)
#echo $result &>> stats1.txt

main
