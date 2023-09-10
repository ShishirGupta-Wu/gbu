#!/bin/bash

## Program needs folders Input_files with all the canonical transcripts and the Scripts folder

main(){
	 set_variables
#	 get_scripts
#	 edit_sequences
#	 create_output
	 }

set_variables(){
		# Set variables
		Project_DIR=alicia_p
		Source_Folder=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/2_Sequence_editor/
		Sequence_Folder=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/1_Get_canonical_transcripts/Output/primary_transcripts_NCBI_PROTEOMES
	  
		## declare an array variable (do not write _primary_transcript.fa here)
		declare -ga speciesarray=('GCF_003719485.1' 'GCF_003719475.1' 'GCF_000209065.1' 'GCF_000002445.2' 'GCA_021307395.1' 'GCA_015033655.1' 'GCA_015033625.1' 'GCA_003543875.1' 'GCA_003177105.1' 'GCA_002087225.1' 'GCA_000691245.1' 'GCA_000300495.1' 'GCA_000227395.2' 'GCA_000210295.1' 'GCA_000188675.2');
		
    	## declare abbreviation array variable
		declare -ga speciesabbrivation=('tconor' 'trange' 'tcruCLb' 'tbrubru' 'tvivax' 'tcruYc' 'tcruBc' 'tbruequ' 'tcruDm' 'tthei' 'tgrayi' 'tcrumar' 'tcongol' 'tbrugam' 'tcruSx');
}

get_scripts(){
    ## Making program executable for Unix
	chmod +x Scripts/*
	mkdir bin 
    cp Scripts/orthomclAdjustFasta bin
	wait
	cp Scripts/dos2unix bin
    wait
    ln -s Scripts/dos2unix .
	## 
	mkdir -p analysis_folder Output_edited_sequences_for_orthology
}

edit_sequences(){

    ## Coverting windows opened file into Unix format to avoid ^M error
	./dos2unix $Sequence_Folder/*.fa
	
	cd analysis_folder
    printf '\n%s\n' && printf '\e[1;34m%-6s\e[m' "############# Creating Symbolic links #############" && printf "\n" 
	ln -s $Sequence_Folder/* .
	ls -ltrh
	
    printf '\n%s\n' && printf '\e[1;34m%-6s\e[m' "############# Replacing dots into underscore in all Sequence files  #############" && printf "\n"	
	for FILE in *.fa; do echo 'processing... '$FILE; sed -i 's/\./_/g' $FILE; done
	
	printf '\n%s\n' && printf '\e[1;34m%-6s\e[m' "############# Replacing U into X to avoid BLAST warnings #############" && printf "\n"		
	for FILE in *.fa; do echo 'processing... '$FILE; sed -i '/^[^>]/s/[U]/X/g' $FILE; done
	
	printf '\n%s\n' && printf '\e[1;34m%-6s\e[m' "############# Editing sequence headers: Trimming and adding abbreviation #############" && printf "\n"	
	## loop through the above arrays defined in set_variables()

	j=0
	for i in "${speciesarray[@]}"
	do
		  echo 'processing... '"$i"
		  perl ../Scripts/orthomclAdjustFasta "${speciesabbrivation[$j]}" "$i"_primary_transcript.fa 1
		  head -n5 "${speciesabbrivation[$j]}".fasta
		  printf "\n"
	let j++
	done 
    ls
}

create_output(){
	cd $Source_Folder && rm -rf dos2unix orthomclAdjustFasta
	# mkdir Output_edited_sequences_for_orthology && cd "$_"
	mv analysis_folder/*.fasta Output_edited_sequences_for_orthology/
    printf '%s\n' && printf '\e[1;34m%-6s\e[m' "############# Final Proteomes for Orthology analysis #############" && printf "\n" 
	cd Output_edited_sequences_for_orthology &&  printf '%s\n' *.fasta
	
	#converting fasta extension to .fa for using in orthofinder
	for f in *.fasta; do mv -- "$f" "${f%.fasta}.fa"; 
	done

}

 
main
