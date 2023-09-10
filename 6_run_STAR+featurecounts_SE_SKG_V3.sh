#!/bin/bash

main(){
    set_variables
#    log_conda_env
#    create_folders
 #   get_scripts
 #  link_and_rename_read_files
    #concat_annotations
 #   create_index
#    run_read_alignment
#    index_bam_files
#    compile_mapping_stats
 #  count_reads_featureCounts_exon
 #  count_reads_featureCounts_gene
 #  count_reads_featurecounts_rrna
 #  clean_featureCounts_table
 #  generate_sequencing_stats
#   run_deseq2
#   run_deseq2_subset
#    run_deseq2_subset_shishir
### run_deseq2_remove_batch_effect
#    generate_Bokeh_plots
#    generate_custom_plots
    download_KEGG_annotation
   clusterProfiler_KEGG
    download_GO_annotation
    clusterProfiler_GO
    generate_Bokeh_Network_plots
  #map_entrez_to_ensembl_ids
#    generate_package_to_send
}

set_variables(){
#    eval "$(conda shell.bash hook)"
#    conda activate /opt/anaconda/envs/shared_env2
    STAR_BIN=STAR
    READ_SOURCE_FOLDER=../data/RNASeq_data_files/input
    GENOME_FASTA=../data/genome_files/input/GCF_000001405.39_GRCh38.p13_genomic_filtered.fna 
    ANNOTATION_GFF_UNIQ=../data/genome_files/input/GCF_000001405.39_GRCh38.p13_genomic_filtered_uniq.gff
    INPUT_FOLDER=input
    OUTPUT_FOLDER=output
    READ_FOLDER=$INPUT_FOLDER/reads
    GENOME_INDEX_DIR=$OUTPUT_FOLDER/index
    MAPPING_FOLDER=$OUTPUT_FOLDER/mappings
    DESEQ_RESULTS_FOLDER=$OUTPUT_FOLDER/DESeq2
    PATHWAY_RESULTS_FOLDER=$OUTPUT_FOLDER/clusterProfiler
    KEGG_ANNOTATION_FOLDER=$INPUT_FOLDER/KEGG_annotation
    #KEGG internal identifiers
    KEGG_ORGANISM=hsa
    GO_ANNOTATION_FOLDER=$INPUT_FOLDER/GO_annotation
    #GO identifier, whitespaces " " need to be replaced with "_" !
    GO_ORGANISM="Homo_sapiens"
}

log_conda_env(){
    conda list > conda_env.txt
}

create_folders(){
    mkdir -p \
        bin \
        $INPUT_FOLDER \
        $KEGG_ANNOTATION_FOLDER \
        $GO_ANNOTATION_FOLDER \
        $OUTPUT_FOLDER \
        $GENOME_INDEX_DIR \
        $MAPPING_FOLDER \
        $READ_FOLDER \
        $DESEQ_RESULTS_FOLDER
}

get_scripts(){
    cp /scratch-rz/users/shg29ny/tools_git_repo/clean_featureCounts_table.py \
       /scratch-rz/users/shg29ny/tools_git_repo/Mapping_Stat_report.py \
       /scratch-rz/users/shg29ny/tools_git_repo/deseq2_subset.R \
       /scratch-rz/users/shg29ny/tools_git_repo/*.R \
       /scratch-rz/users/shg29ny/tools_git_repo/Bokeh*.py \
	   /scratch-rz/users/shg29ny/tools_git_repo/map_entrez_to_ensembl_ids.R \
       /scratch-rz/users/shg29ny/tools_git_repo/csv2xlsx.py \
	   /scratch-rz/users/shg29ny/tools_git_repo/plot_gene_expression_per_cond.py \
        bin/
}

link_and_rename_read_files(){
    for FILE in $(ls $READ_SOURCE_FOLDER)
    do
#         MOD_FILE_NAME=$(echo $FILE | cut -c 6-)
 	 MOD_FILE_NAME=$(echo $FILE)
        echo $MOD_FILE_NAME
        ln -s ../../$READ_SOURCE_FOLDER/$FILE \
            $READ_FOLDER/$MOD_FILE_NAME
    done
}

# concat_annotations(){
#     cat $ANNOTATION_GFF_UNIQ $SOD1_ANNOTATION_GFF \
#     > $COMBINED_ANNOTATION_GFF
# }

create_index(){
    $STAR_BIN \
        --runThreadN 50 \
        --runMode genomeGenerate \
        --genomeDir $GENOME_INDEX_DIR \
        --genomeFastaFiles $GENOME_FASTA
    mv Log.out $GENOME_INDEX_DIR
}

run_read_alignment(){
    # Increase open files limit for sorting
#   ulimit -n 10000
    ulimit -n 65000

    for LIB in $(ls $READ_FOLDER)
    do
        OUTPUT_FILE_PREFIX=$(echo $LIB | sed -e "s/.fq.bz2/./")
        echo $OUTPUT_FILE_PREFIX
        $STAR_BIN \
            --runMode alignReads \
            --runThreadN 40 \
            --genomeDir $GENOME_INDEX_DIR \
            --readFilesIn $READ_FOLDER/$LIB \
            --readFilesCommand bzcat \
            --sjdbGTFfile $ANNOTATION_GFF_UNIQ \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbGTFtagExonParentGene uniq \
			--outFilterScoreMinOverLread 0.3 \
			--outFilterMatchNminOverLread 0.3 \
            --outFileNamePrefix $MAPPING_FOLDER/$OUTPUT_FILE_PREFIX \
            --outSAMtype BAM SortedByCoordinate
    done
	
#	outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3
}

index_bam_files(){
    for BAM in $(ls $MAPPING_FOLDER | grep bam$)
    do
        /artemis/software/bin/samtools index $MAPPING_FOLDER/$BAM &
    done
    wait
}

compile_mapping_stats(){
    grep "Uniquely mapped reads %" output/mappings/*Log.final.out \
    | sed -e "s/output\/mappings\///" -e s/:.*\|// -e "s/.Log.final.out//"\
    > $OUTPUT_FOLDER/Percentage_uniquely_mapped_reads.csv
}

# IMPORTANT for featureCounts
# -s parameter might need adjustment: 0 count both strands, 1 count reads in sense, 2 count reads in antisense to features
# Parameter values for commonly used kits:
# 0: Takara SMART-Seq v4 Ultra Low Input + Nextera
# 1: NEBNext Multiplex Small RNA
# 2: Illumina TruSeq stranded, Takara SMARTer and SMART-Seq Stranded, NEBNext Ultra II Directional RNA

count_reads_featureCounts_exon(){
    mkdir -p $OUTPUT_FOLDER/featureCounts
    featureCounts \
        -T 42 \
        -s 2 \
        -t exon \
        -g uniq \
	    -M \
        -O \
        --fraction \
        --extraAttributes "ID,gene,gene_biotype,Dbxref,description" \
        -a $ANNOTATION_GFF_UNIQ \
        -o $OUTPUT_FOLDER/featureCounts/counts.csv \
        $(ls $MAPPING_FOLDER/*.Aligned.sortedByCoord.out.bam)
}

count_reads_featureCounts_gene(){
    mkdir -p $OUTPUT_FOLDER/featureCounts
    featureCounts \
        -T 42 \
        -s 2 \
        -f \
        -t gene \
        -g uniq \
        -M \
        -O \
        --fraction \
        --extraAttributes "ID,gene,gene_biotype,Dbxref,description" \
        -a $ANNOTATION_GFF_UNIQ \
        -o $OUTPUT_FOLDER/featureCounts/counts_gene.csv \
        $(ls $MAPPING_FOLDER/*.Aligned.sortedByCoord.out.bam)
}

count_reads_featurecounts_rrna(){
    mkdir -p $OUTPUT_FOLDER/featureCounts
    featureCounts \
        -T 42 \
        -s 2 \
        -f \
        -t rRNA \
        -g uniq \
        -M \
        -O \
        --fraction \
        --extraAttributes "ID,gene,gene_biotype,Dbxref,description" \
        -a $ANNOTATION_GFF_UNIQ \
        -o $OUTPUT_FOLDER/featureCounts/counts_rrna.csv \
        $(ls $MAPPING_FOLDER/*.Aligned.sortedByCoord.out.bam)
}

clean_featureCounts_table(){
    bin/clean_featureCounts_table.py \
        --annotations $ANNOTATION_GFF_UNIQ \
        --input $OUTPUT_FOLDER/featureCounts/counts.csv \
        --output $OUTPUT_FOLDER/featureCounts/counts_cleaned.csv
    bin/clean_featureCounts_table.py \
        --annotations $ANNOTATION_GFF_UNIQ \
        --input $OUTPUT_FOLDER/featureCounts/counts_gene.csv \
        --output $OUTPUT_FOLDER/featureCounts/counts_gene_cleaned.csv
    bin/clean_featureCounts_table.py \
        --annotations $ANNOTATION_GFF_UNIQ \
        --input $OUTPUT_FOLDER/featureCounts/counts_rrna.csv \
        --output $OUTPUT_FOLDER/featureCounts/counts_rrna_cleaned.csv

    ls $OUTPUT_FOLDER/featureCounts/counts*_cleaned.csv | parallel \
        sed -i \
            -e "s:$OUTPUT_FOLDER/mappings/::g" \
            -e "s:.Aligned.sortedByCoord.out.bam::g" \
            $OUTPUT_FOLDER/featureCounts/counts_cleaned.csv

     sed -e "1d" \
        $OUTPUT_FOLDER/featureCounts/counts_cleaned.csv \
    > $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv
}

# samplelist can be used to reorder samples, it does NOT remove them
# Samples given in samplelist are taken out of the BAM order and reinserted at the end of the list in the order provided in samplelist
# Use this to move e.g. no rRNA depletion internal tests to the end of the list
generate_sequencing_stats(){
  python3 bin/Mapping_Stat_report.py \
        --path $PWD \
        --type STAR \
        --output output/sequencing_stats.html \
        --projectnumber Project000013_PRJNA528433 \
        --runids Macrophages \
        --reference GCF_000001405.39_GRCh38.p13  

  killall screen
}

run_deseq2(){
    mkdir -p ${DESEQ_RESULTS_FOLDER}/default
    Rscript bin/deseq2.R \
        --geneCountingFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv \
        --countTansformMethod vst \
        --log2FcCutoff 1.0 \
        --PadjCutoff 0.05 \
        --outputFolder ${DESEQ_RESULTS_FOLDER}/default
}

##         --condRegex '_[1-3]_[^ ]*' \

run_deseq2_subset(){
    mkdir -p ${DESEQ_RESULTS_FOLDER}/default
    Rscript bin/deseq2_subset.R \
        --geneCountingFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv \
        --annoCols 10 \
	    --libSelectRegex '^(Mature|Immature)' \
        --condRegex '_[1-16]' \
		--comparisons "Mature:Immature" \
         --countTansformMethod rlog \
        --log2FcCutoff 1.0 \
        --PadjCutoff 0.05 \
        --outputFolder ${DESEQ_RESULTS_FOLDER}/default
}


run_deseq2_subset_shishir(){
    mkdir -p ${DESEQ_RESULTS_FOLDER}/default
    Rscript bin/deseq2_shi.R \
        --geneCountingFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv \
        --countTansformMethod rlog \
        --libSelectRegex '^(SRR)' \
        --comparisons \
                      "Af2h:Ctrl" \
					  "Af6h:Ctrl" \
        --log2FcCutoff 0.58 \
        --PadjCutoff 0.05 \
        --outputFolder ${DESEQ_RESULTS_FOLDER}/default
		
}

run_deseq2_remove_batch_effect(){
#  mv ${DESEQ_RESULTS_FOLDER}/default ${DESEQ_RESULTS_FOLDER}/default_without_batch

  mkdir -p ${DESEQ_RESULTS_FOLDER}/default
    Rscript bin/deseq2_batch_correction_shi.R \
        --geneCountingFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv \
        --countTansformMethod rlog \
        --libSelectRegex '^(Healthy|ICU|NonICU)' \
        --comparisons \
                      "ICU:Healthy" \
					  "Non-ICU:Healthy" \
					  "ICU:Healthy" \
        --log2FcCutoff 0.58 \
        --PadjCutoff 0.05 \
        --outputFolder ${DESEQ_RESULTS_FOLDER}/default
		
	mkdir -p ${DESEQ_RESULTS_FOLDER}/default_2
    Rscript bin/deseq2_batch_correction_shi2.R \
        --geneCountingFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv \
        --countTansformMethod rlog \
        --libSelectRegex '^(Healthy|ICU|NonICU)' \
        --comparisons \
                      "Non-ICU+ICU:Healthy" \
        --log2FcCutoff 0.58 \
        --PadjCutoff 0.05 \
		--outputFolder ${DESEQ_RESULTS_FOLDER}/default_2		

}


generate_Bokeh_plots(){
    for FILE in $(find ${DESEQ_RESULTS_FOLDER}* -name DESeq2_comparison*.csv)
    do
        COMP=$(basename $FILE | sed -e "s/DESeq2_comparison_//" -e "s/_table.csv//")
        COMP_COND=$(echo $COMP | awk 'BEGIN {FS="_vs_"} {print $1}')
        COMP_REF=$(echo $COMP | awk 'BEGIN {FS="_vs_"} {print $2}')
        python3 bin/Bokeh_plots.py \
            --input_file $FILE \
            --out_file ${FILE%/*}/DESeq2_comparison_${COMP}_interactive_plots.html \
            --condComp $COMP_COND \
            --condRef $COMP_REF \
            --altsym gene
    done
}

generate_custom_plots(){
    bin/plot_gene_expression_per_cond.py \
        --deseq_result_file ${DESEQ_RESULTS_FOLDER}/default/DESeq2_comparison_Af2h_vs_Ctrl_table.csv\
        --plot_file ${DESEQ_RESULTS_FOLDER}/default/beeswarm_Af2h_vs_Ctrl.pdf \
        --gene_symbols CLEC1B CLEC1A AHR SFTPD GSK3B GAA FKBP1A METAP2 METAP1 SENP8 APEX1 TXN1 FASN MTOR PARP1 SOD1 SDS MGLL FLOT1 FLOT2
		
	bin/plot_gene_expression_per_cond.py \
        --deseq_result_file ${DESEQ_RESULTS_FOLDER}/default/DESeq2_comparison_Af6h_vs_Ctrl_table.csv\
        --plot_file ${DESEQ_RESULTS_FOLDER}/default/beeswarm_Af6h_vs_Ctrl.pdf \
        --gene_symbols CLEC1B CLEC1A AHR SFTPD GSK3B GAA FKBP1A METAP2 METAP1 SENP8 APEX1 TXN1 FASN MTOR PARP1 SOD1 SDS MGLL FLOT1 FLOT2
}
		
download_KEGG_annotation(){
    Rscript bin/download_KEGG_annotation.R \
        --geneCountingFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv \
        --annoColNumber 10 \
        --organism ${KEGG_ORGANISM} \
        --annotationFolder $KEGG_ANNOTATION_FOLDER
}

clusterProfiler_KEGG(){
    for FILE in  $(find $DESEQ_RESULTS_FOLDER* -name DESeq2_comparison*.csv)
    do
        COMP=$(basename $FILE | sed -e "s/DESeq2_comparison_//" -e "s/_table.csv//")
        mkdir -p $PATHWAY_RESULTS_FOLDER/KEGG/$(basename ${FILE%/*})/${COMP}
        echo $PATHWAY_RESULTS_FOLDER/KEGG/$(basename ${FILE%/*})/${COMP}
        Rscript bin/clusterProfiler_KEGG_enricher.R \
            --deseqResultsFile $FILE \
            --log2FcCutoff 1.0 \
            --PadjCutoff 0.05 \
            --keggPathID2ExtIDFile $KEGG_ANNOTATION_FOLDER/kegg_path_id_to_entrez_id_${KEGG_ORGANISM}.csv \
            --keggPathID2NameFile $KEGG_ANNOTATION_FOLDER/kegg_path_id_to_name_${KEGG_ORGANISM}.csv \
            --qvalueCutoffs 0.1 0.05 0.01 \
            --outputFolder $PATHWAY_RESULTS_FOLDER/KEGG/$(basename ${FILE%/*})/${COMP}

        Rscript bin/clusterProfiler_KEGG_GSEA.R \
            --deseqResultsFile $FILE \
            --keggPathID2ExtIDFile $KEGG_ANNOTATION_FOLDER/kegg_path_id_to_entrez_id_${KEGG_ORGANISM}.csv \
            --keggPathID2NameFile $KEGG_ANNOTATION_FOLDER/kegg_path_id_to_name_${KEGG_ORGANISM}.csv \
            --qvalueCutoffs 0.1 0.05 0.01 \
            --outputFolder $PATHWAY_RESULTS_FOLDER/KEGG/$(basename ${FILE%/*})/${COMP}
    done
}

download_GO_annotation(){
    Rscript bin/download_GO_annotation.R \
        --geneCountingFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv \
        --annoColNumber 10 \
        --organism ${GO_ORGANISM} \
        --annotationFolder $GO_ANNOTATION_FOLDER
}

clusterProfiler_GO(){
    for FILE in $(find $DESEQ_RESULTS_FOLDER* -name DESeq2_comparison*.csv)
    do
        COMP=$(basename $FILE | sed -e "s/DESeq2_comparison_//" -e "s/_table.csv//")
        mkdir -p $PATHWAY_RESULTS_FOLDER/GO/$(basename ${FILE%/*})/${COMP}
        echo $PATHWAY_RESULTS_FOLDER/GO/$(basename ${FILE%/*})/${COMP}
        Rscript bin/clusterProfiler_GO_enricher.R \
            --deseqResultsFile $FILE \
            --log2FcCutoff 1.0 \
            --PadjCutoff 0.05 \
            --gotranstable $GO_ANNOTATION_FOLDER/go_path_id_to_entrez_id_${GO_ORGANISM}.csv \
            --qvalueCutoffs 0.1 0.05 0.01 \
            --outputFolder $PATHWAY_RESULTS_FOLDER/GO/$(basename ${FILE%/*})/${COMP}

        Rscript bin/clusterProfiler_GO_GSEA.R \
            --deseqResultsFile $FILE \
            --gotranstable $GO_ANNOTATION_FOLDER/go_path_id_to_entrez_id_${GO_ORGANISM}.csv \
            --qvalueCutoffs 0.1 0.05 0.01 \
            --outputFolder $PATHWAY_RESULTS_FOLDER/GO/$(basename ${FILE%/*})/${COMP}
    done
}

generate_Bokeh_Network_plots(){
    for FILE in $(find $PATHWAY_RESULTS_FOLDER -name "enricher_results_*.csv")
    do
        COMP_COND=$(echo $(basename ${FILE%/*}) | awk 'BEGIN {FS="_vs_"} {print $1}')
        COMP_REF=$(echo $(basename ${FILE%/*}) | awk 'BEGIN {FS="_vs_"} {print $2}')
        echo $FILE
        python3 bin/Bokeh_NetworkPlot.py \
            --input_file $FILE \
            --out_file ${FILE%.csv}_Networkplot.html \
            --condComp $COMP_COND \
            --condRef $COMP_REF \
            --plottype Enricher
    done

    for FILE in $(find $PATHWAY_RESULTS_FOLDER -name "gsea_results.csv")
    do
        COMP_COND=$(echo $(basename ${FILE%/*}) | awk 'BEGIN {FS="_vs_"} {print $1}')
        COMP_REF=$(echo $(basename ${FILE%/*}) | awk 'BEGIN {FS="_vs_"} {print $2}')
        echo $FILE
        python3 bin/Bokeh_NetworkPlot.py \
            --input_file $FILE \
            --out_file ${FILE%.csv}_Networkplot.html \
            --condComp $COMP_COND \
            --condRef $COMP_REF \
            --plottype GSEA
    done
	
	killall screen
}

map_entrez_to_ensembl_ids(){
    Rscript bin/map_entrez_to_ensembl_ids.R \
        --annoFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene.csv \
        --outputFile $OUTPUT_FOLDER/featureCounts/Reads_per_gene_ensembl.csv
}

generate_package_to_send(){
    SEND_FOLDER=2020-09-02_Project000011
    mkdir -p $SEND_FOLDER
    cp $OUTPUT_FOLDER/sequencing_stats.html \
        $SEND_FOLDER/
    mkdir -p $SEND_FOLDER/DESeq2
    cp -r $OUTPUT_FOLDER/DESeq2/default/* \
       $SEND_FOLDER/DESeq2
    mkdir -p $SEND_FOLDER/clusterProfiler
    cp -r $OUTPUT_FOLDER/clusterProfiler/KEGG \
        $OUTPUT_FOLDER/clusterProfiler/GO \
        $SEND_FOLDER/clusterProfiler
    find $SEND_FOLDER -name "*.txt" -exec rm {} \;
    find $SEND_FOLDER -name "*csv" -print0 | xargs -0 -n 1 -P24 python2 bin/csv2xlsx.py
    find $SEND_FOLDER -name "*csv" -exec rm {} \;
    zip -r $SEND_FOLDER.zip $SEND_FOLDER
}

main
