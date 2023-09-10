#!/bin/bash

## Program needs Installation of seqkit
## conda activate ncbi_datasets

## chmod +x preprocessing_protein_sequences_v3.sh
## print logfile to verify if seqkit worked fine
## ./preprocessing_protein_sequences_v3.sh  | tee -a logfile.txt

## ** Carefully check log file and fix bad gff related genomes ** quick fix mv file to required fa

## If you are not using translated CDS or mixed refseq and genebank you can run skip_process_NCBI_proteome

main(){
	 set_variables
##	 process_NCBI_proteome
	 skip_process_NCBI_proteome
	 get_all_processed_proteomes
	 }


set_variables(){
	  Set variables
	  Project_DIR=monoc_p
      Source_Folder=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/1_Get_canonical_transcripts
	  NCBI_PROTEOMES=/scratch-rz/users/shg29ny/PROJECT_22/$Project_DIR/0_data_download/input/data_files/ncbi_dataset/data
}

 
process_NCBI_proteome(){	
 
    printf '\n%s\n' "####### Downloaded Datasets: #######"
	cd $NCBI_PROTEOMES
	cat list.txt
	#tree

	
	while read -r line;
		do
		   printf '\n%s\n' "####### Working on: $line #######"
     	   cd "$line"		   
		     
		   seqkit fx2tab -n -l protein.faa | awk 'BEGIN{FS="\t"}{gsub(".+ref\\|", "", $1); gsub("\\|.*$", "", $1); print}' | sort > protein_len.tsv
		   		   
		   cat genomic.gff | awk '$3=="CDS"{gsub(".+;Dbxref=GeneID:", "", $9); gsub(";Name=.+", "", $9); gsub(",Genbank:", "\t", $9); print $9}' | sort -k2,2 | uniq  > gene_protein.tsv
	
		   cat gene_protein.tsv | sed 's/\(.*\):/\1@/' | perl -p -e 's/@/\t/g' > tab_gene_protein.tsv
	  
		   head -n3 tab_gene_protein.tsv protein_len.tsv
		   
		   join -1 2 tab_gene_protein.tsv -2 1 protein_len.tsv | awk 'len[$2]<$3 {len[$2]=$3; id[$2]=$1} END {for (i in id) {print id[i]}}' > protein_longest.txt
		   
		   seqkit grep -r -n -f protein_longest.txt protein.faa  > "$line"_primary_transcript.fa

		   grep -c '>' protein.faa
		   grep -c '>' *_primary_transcript.fa
		   grep -c '>' cds_from_genomic.fna
		   
		   cd ..
		   
		done < list.txt
		wait
		
}

skip_process_NCBI_proteome(){	
 
    printf '\n%s\n' "####### Downloaded Datasets: #######"
	cd $NCBI_PROTEOMES
	cat list.txt
	#tree

	
	while read -r line;
		do
		   printf '\n%s\n' "####### Working on: $line #######"
     	   cd "$line"		   
		     
		   mv protein.faa  "$line"_primary_transcript.fa 	 	   

		   grep -c '>' "$line"_primary_transcript.fa
		   
		   cd ..
		   
		done < list.txt
		wait
		
}


get_all_processed_proteomes(){
      printf '\n%s\n' "####### Preparing Folder for Next Steps #######"
	  cd $Source_Folder
      mkdir -p Output && cd Output
	  mkdir primary_transcripts_NCBI_PROTEOMES
	  
	  cd $NCBI_PROTEOMES	

	while read -r line;
		do
		   printf '\n%s\n' "####### Working on: $line #######"
     	   cd "$line"		   
		     
		   cp *_primary_transcript.fa $Source_Folder/Output/primary_transcripts_NCBI_PROTEOMES/
		   
		   cd ..
		   
		done < list.txt
		wait	  
	  
      printf '\n%s\n' "%%%####### Count in primary transcripts #######%%%"	  
	  
	  cd $Source_Folder/Output/primary_transcripts_NCBI_PROTEOMES/
	  grep -c '>' *.fa 
 	  
}


main
