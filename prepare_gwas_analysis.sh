#!/bin/bash

main(){
	set_variables
    create_folders
 	get_data
	align_reads
	convert_sam_to_bam
	extract_genomic_differences
}

set_variables(){
	  SOURCE_DIR=/storage/evolgen/shg29ny/gwas_test
	  DATA_FOLDER=test_data
	  REFERENCE_FASTA=Aleph1.ref.fa
	  READS=Aleph_new_read_data.fq
	  Reads_file=$(echo $READS | cut -f1 -d".")
}

create_folders(){
	mkdir -p Sequenzdaten
	mkdir -p Sequenzdaten/Index 
}
   
get_data(){
	cp $DATA_FOLDER/* Sequenzdaten
	cd Sequenzdaten
	gunzip *.gz
	
	printf '%s\n' "##### data should have reference_fasta; gff_annotation; and fastq_file. See the contents #####"	
	tree
}
	
align_reads(){
    cd $SOURCE_DIR/Sequenzdaten
	##ADJUST
	Index_file=$(echo $REFERENCE_FASTA | cut -f1 -d".")
		
	printf '%s\n' "##### building index with bowtie2-build #####"
	bowtie2-build -f $REFERENCE_FASTA Index/$Index_file.index -p 16

	printf '%s\n' "##### aligning reads with bowtie2 #####"
	bowtie2 -x Index/$Index_file.index -U $READS -S $Reads_file.sam -p 16 &>> align_stats.txt
	cat align_stats.txt

}

convert_sam_to_bam(){
    cd $SOURCE_DIR/Sequenzdaten
	printf '%s\n' "##### converting sam to bam #####"
	samtools view -bS -o $Reads_file.bam $Reads_file.sam

	printf '%s\n' "##### soring bam file #####"
	samtools sort $Reads_file.bam -o $Reads_file.sorted.bam

	printf '%s\n' "##### indexing the sorted bam file #####"
	samtools index $Reads_file.sorted.bam

	printf '%s\n' "##### look at the files #####"
	samtools tview $Reads_file.sorted.bam $REFERENCE_FASTA
}


extract_genomic_differences(){
    cd $SOURCE_DIR/Sequenzdaten
	
	printf '%s\n' "##### extracting the genomic differences #####"
	
   #bcftools mpileup -t DP -f $REFERENCE_FASTA $Reads_file.sorted.bam | bcftools call -c -> $Reads_file.vcf
	samtools mpileup -u -t DP -f $REFERENCE_FASTA $Reads_file.sorted.bam | bcftools call -c -> $Reads_file.vcf
	cat $Reads_file.vcf
	
	## Compare the .vcf file and graphic output (samtools tview command)
	## samtools tview $Reads_file.sorted.bam $REFERENCE_FASTA
}

main












