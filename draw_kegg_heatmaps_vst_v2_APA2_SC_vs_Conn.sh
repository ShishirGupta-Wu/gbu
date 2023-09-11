#!/bin/bash

### Added function edit_vst_file

## script folder is required

main(){
    set_variables
    extract_genes
	order_vst_file
	select_sample_from_vst_file
    filter_gene_list
    draw_heatmaps
	draw_custom_heatmaps
	clean_up
}

set_variables(){
#    eval "$(conda shell.bash hook)"
#    conda activate /opt/anaconda/envs/shared_env
#    conda list > conda_env.txt
	
	SOURCE_FOLDER=/home/shishir/PROJECTS/Heatmap_Kegg
    PATHWAY_LIST_FILE=$SOURCE_FOLDER/gsea_results_APA2_SC_vs_Conn.txt
	VST_FILE=$SOURCE_FOLDER/vst-transformed_count_table.csv
	
}


extract_genes(){
    
	input=$PATHWAY_LIST_FILE
	while IFS= read -r line
		do
			echo "$line"
			
				Rscript scripts/get_pathway_genes.R ${line} ${line}_geneid.txt &			
		done < "$input"
       wait
	   
}


order_vst_file(){
	# Define input and output filenames
	#input_file="vst-transformed_count_table.csv"
	Ordered_VST_file="ordered_vst-transformed_count_table.csv"

	# Call the R script with input and output filenames as arguments
	#Rscript reorder_vst_rlog_file_bash.R "$input_file" "$output_file"
	Rscript scripts/reorder_vst_rlog_file_bash.R $VST_FILE $Ordered_VST_file
    head -n3 $Ordered_VST_file
}

select_sample_from_vst_file(){
	##ADJUST sample_to_use1 and sample_to_use2

	# Define input and output file names
	 input_file=$SOURCE_FOLDER/ordered_vst-transformed_count_table.csv # here $Ordered_VST_file	
	 output_file="Edited_vst.csv"

	# Define the sample_to_use1 and sample_to_use2 values
	sample_to_use1="APA2_SC"
	sample_to_use2="Conn"

	# Run the R script with arguments

	Rscript scripts/select_vst_rlog_samples_bash.R $Ordered_VST_file $output_file $sample_to_use1 $sample_to_use2
	
	# Rscript select_vst_rlog_samples_bash.R ordered_vst-transformed_count_table.csv Edited_vst.csv APA2_SC Conn2

    head -n5 $output_file

}


filter_gene_list(){

	input=$PATHWAY_LIST_FILE
	while IFS= read -r line
		do
			echo "$line"
			
				awk {'print $2'} ${line}_geneid.txt | sed -e 's/"//g' | sed -e 1d > ${line}_geneid_final
				perl -lpe 's/"/""/g; s/^|$/"/g; s/\t/","/g' < $SOURCE_FOLDER/$output_file > Temp1
				awk -F ","  '{print $7,",",$9}' Temp1 | sed -e 's/ //g' > Temp2
				awk -F ","  '{for (i=1; i<=NF-1; i++) if(i>10) {printf $i FS};{print $NF}}' Temp1 > Temp3
				paste -d ","  Temp2 Temp3 > Temp4
				head -n1 Temp4 > Temp5
				sed -i 's/"//g' Temp4

				awk -F "," 'NR==FNR{a[$1]; next} $2 in a' ${line}_geneid_final Temp4 > Temp6
				
				cat Temp5 Temp6 > Temp7
				awk -F ","  '{for (i=1; i<=NF-1; i++) if(i<2 || i>2) {printf $i FS};{print $NF}}' Temp7 > Temp8
				cat Temp8 | tr '[,]' '[\t]' | sed 's/"//g' | sed 's/ /\t/g' | awk '!a[$1]++' > ${line}_vst.csv


			#	rm -rf Temp*
				wc -l  ${line}_vst.csv ${line}_geneid_final
				
		done < "$input"
       wait

}

draw_heatmaps(){

	input=$PATHWAY_LIST_FILE
	while IFS= read -r line
		do
			echo "$line"
			
				Rscript scripts/create_heatmaps.R ${line}_vst.csv ${line}_${sample_to_use1}_${sample_to_use2}_heatmap.pdf	${line}_${sample_to_use1}_${sample_to_use2}_heatmap_v1.pdf &
				
		done < "$input"
       wait
}





draw_custom_heatmaps(){

	### Please edit the file create_custom_heatmaps_v2.R before running this function
    # check head -n1 of any ${line}_vst.csv for editing samples and colors in create_custom_heatmaps_v2.R
	
	input=$PATHWAY_LIST_FILE
	while IFS= read -r line
		do
			echo "$line"
			
				Rscript scripts/create_custom_heatmaps_v2_APA2_SC_vs_Conn.R ${line}_vst.csv ${line}_${sample_to_use1}_${sample_to_use2}_heatmap_v2.pdf &
				
		done < "$input"
       wait
}


clean_up(){

	mkdir -p intermediate_files pathway_vst_files ${sample_to_use1}_${sample_to_use2}_heatmaps ${sample_to_use1}_${sample_to_use2}_group_heatmaps
	mv *_v2.pdf ${sample_to_use1}_${sample_to_use2}_group_heatmaps
	mv *_v1.pdf ${sample_to_use1}_${sample_to_use2}_heatmaps
	mv *.pdf intermediate_files
	mv *_vst.csv pathway_vst_files
	mv *_geneid.txt *_final intermediate_files
	rm Temp*
}



main
