library('argparse')
library('dplyr')
library('DESeq2')
library('ggplot2')
library('ggrepel')
library('ggpubr')
library('EnhancedVolcano')
library('RColorBrewer')
library('pheatmap')

main <- function(){
    # Read command line parameters
    parser <- ArgumentParser()
    parser$add_argument("--geneCountingFile")
    parser$add_argument("--countTansformMethod", choices=c('normTransform', 'vst', 'rlog'))
    parser$add_argument("--log2FcCutoff", type="double", default=1.0)
    parser$add_argument("--PadjCutoff", type="double", default=0.05)
    parser$add_argument("--outputFolder")
    args <- parser$parse_args()
    debug <- TRUE

    # Select the pure counting matrix
    
    rawCountTable <- read.table(args$geneCountingFile, sep='\t', header=TRUE, quote='', comment.char='')
 
 #  rawCountTable <- read.csv("output/featureCounts/Reads_per_gene.csv", header=T)
 #   countTable <- round(select(rawCountTable, matches("SMPL_[^ ]*")))
    
     countTable <- round(select(rawCountTable, matches("^SMPL")))
    
    libs <- colnames(countTable)
    annotationTable <- select(rawCountTable, -libs)
    if(debug){print(colnames(countTable))}
    if(debug){print(colnames(annotationTable))}

    # calculate TPM values
    tpmTable <- apply(countTable, 2, tpm, rawCountTable$Combined_exon_length)
    write.table(cbind(annotationTable, tpmTable), file=sprintf("%s/TPM_count_table.csv", args$outputFolder),
        quote=FALSE, sep='\t', row.names=FALSE)

    # Derive the condition and replicate names from lib names
#    conds <- gsub("_RNA[0-9]+$", "", libs) # ADJUST
    conds <- c('Control', 'Control', 'Control', 'KO', 'KO', 'KO')

    samples <- data.frame(row.names=libs, condition=conds, lib=libs)
    if(debug){print(samples)}

    # generate DESeq2 data set
    dds <- DESeqDataSetFromMatrix(
        countData=countTable, colData=samples, design=~condition)
    dds <- DESeq(dds, betaPrior=TRUE)

    # calculated data transformation
    transformed <- get(args$countTansformMethod)(dds)

    # plot PCA and sample distance heatmap
    createPcaAndHeatmap(transformed, args$countTansformMethod, args$outputFolder, lab=libs)

    # generate transformed count table
    transformedTable <- cbind(annotationTable, assay(transformed))
    write.table(transformedTable, file=sprintf(
        "%s/%s-transformed_count_table.csv",
        args$outputFolder, args$countTansformMethod),
        quote=FALSE, sep='\t', row.names=FALSE)

    # calculate normalized counts
    print(sizeFactors(dds))
    normCountTable = data.frame(t(t(countTable)/sizeFactors(dds)))
    colnames(normCountTable) = paste(libs, "normalized", sep = "_")

    # pairwise comparisons	#  c("untreated","treated")
	
	
    comparisons <- list(
        c("SMPL_ctrl_1_218", "SMPL_Npr3_KO_1_217"),
        c("SMPL_ctrl_2_219", "SMPL_Npr3_KO_2_225"),
        c("SMPL_ctrl_3_239", "SMPL_Npr3_KO_3_227"),
        c("SMPL_ctrl_1_218", "SMPL_ctrl_2_219"),
        c("SMPL_ctrl_1_218", "SMPL_ctrl_3_239"),
        c("SMPL_ctrl_2_219", "SMPL_ctrl_3_239"),
        c("SMPL_Npr3_KO_1_217", "SMPL_Npr3_KO_2_225"),
        c("SMPL_Npr3_KO_1_217", "SMPL_Npr3_KO_3_227"),
        c("SMPL_Npr3_KO_2_225", "SMPL_Npr3_KO_3_227")
    ) # ADJUST
	
    for (comp in comparisons)
	{
        compareConditions(comp[1], comp[2], dds, rawCountTable,
            normCountTable, args$log2FcCutoff, args$PadjCutoff, args$outputFolder)
    }

    writeLines(capture.output(sessionInfo()),
    sprintf("%s/session_info.txt", args$outputFolder))
}

tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
}

createPcaAndHeatmap <- function(transformed, countTansformMethod, outputFolder, lab=NULL){
    pdf(sprintf("%s/DESeq2_PCA_and_heatmap_%s.pdf", outputFolder, countTansformMethod))
    # PCA
    pcaData <- plotPCA(transformed, intgroup="condition", returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pca <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
          geom_point(size=3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) +
          coord_fixed() +
          ggtitle(sprintf("PCA (%s-transformed counts)", countTansformMethod)) +
          theme_bw() +
          theme(plot.title = element_text(color="black", size=18, face="bold"))
    print(pca)
    if (!is.null(lab)){
        pca_lab <- pca + geom_text_repel(aes(label=lab), size=2, color="black")
        print(pca_lab)
    }

    # Heatmap
    sampleDists <- dist(t(assay(transformed)))
    sampleDistMatrix <- as.matrix(sampleDists)
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors,
             main="Sample distance heatmap")
    dev.off()
}


compareConditions <- function(Control, KO, dds, rawCountTable,
                              normCountTable, log2FcCutoff, PadjCutoff,
                              outputFolder)
{
    comp <- results(dds, contrast=c('condition', "Control", "KO"), addMLE=FALSE)
    countingAndComparisonResults <- cbind(rawCountTable, normCountTable, comp)
    write.table(countingAndComparisonResults, file=sprintf(
                          "%s/DESeq2_comparison_%s_vs_%s_table.csv",
                          outputFolder, "Control", "KO"),
                quote=FALSE, sep='\t', row.names=FALSE)

    # Generate MA and volcano plot for contrast
    pdf(sprintf("%s/DESeq2_comparison_%s_vs_%s_MA_and_volcano_plot_padj_%s_log2fc_%s.pdf",
                outputFolder, condComp, condRef, PadjCutoff, log2FcCutoff))
    drawMAPlot(comp, PadjCutoff, 2**log2FcCutoff, rawCountTable$gene)
    drawVolcanoPlot(comp, PadjCutoff, log2FcCutoff, rawCountTable$gene)
    dev.off()
}

drawMAPlot <- function(res, alpha, fcCutoff, geneIDs, labelNo=10){
    ylims <- c(-ceiling(max(abs(res$log2FoldChange), na.rm=TRUE)),
               ceiling(max(abs(res$log2FoldChange), na.rm=TRUE)))
    MAPlot <- ggmaplot(data = res, main = "MA plot",
                       fdr = alpha, fc = fcCutoff, size = 1.5,
                       palette = c("red", "blue", "lightgrey"),
                       genenames = geneIDs,
                       legend = "bottom",
                       top = labelNo,
                       font.label = c(10, "plain", "black"),
                       font.x = c(18, "plain", "black"),
                       font.y = c(18, "plain", "black"),
                       font.legend = c(12, "plain", "black"),
                       font.main = c(18, "bold", "black"),
                       xlab = bquote(~Log[2]~ 'base mean'),
                       ylab = bquote(~Log[2]~ 'fold change'),
                       ylim = ylims,
                       ggtheme = theme_bw())
    print(MAPlot)
}

drawVolcanoPlot <- function(res, alpha, fcCutoff, geneIDs, labelNo=10){
    xlims <- c(-ceiling(max(abs(res$log2FoldChange), na.rm=TRUE)),
               ceiling(max(abs(res$log2FoldChange), na.rm=TRUE)))
    resdf <- as.data.frame(res)
    resdf$color <- "lightgrey"
    resdf$class <- "NS"
    resdf$color[resdf$log2FoldChange >= fcCutoff & resdf$padj < alpha] <- "red"
    resdf$class[resdf$log2FoldChange >= fcCutoff & resdf$padj < alpha] <- "Up"
    resdf$color[resdf$log2FoldChange <= -fcCutoff & resdf$padj < alpha] <- "blue"
    resdf$class[resdf$log2FoldChange <= -fcCutoff & resdf$padj < alpha] <- "Down"

    cols <- resdf$color
    names(cols) <- resdf$class

    if(labelNo > 0){
        topgenes <- resdf$padj
        names(topgenes) <- geneIDs
        topgenes <- topgenes[(resdf$padj < alpha) & (
            (resdf$log2FoldChange >= fcCutoff) | (resdf$log2FoldChange <= -fcCutoff))]
        topgenes <- sort(topgenes, decreasing = FALSE, na.last = NA)[1:min(c(labelNo, nrow(topgenes)))]
        selectedLabels <- names(topgenes)
    }
    else{
        selectedLabels <- as.character(c())
    }

    volcanoPlot <- EnhancedVolcano(
        res,
        lab = as.vector(geneIDs),
        x = 'log2FoldChange',
        y = 'padj',
        selectLab = selectedLabels,
        xlim = xlims,
        title = 'Volcano plot',
        subtitle = NULL,
        xlab = bquote(~Log[2]~ 'fold change'),
        ylab = bquote(~-Log[10]~adjusted~italic(P)),
        colCustom = cols,
        pCutoff = alpha,
        FCcutoff = fcCutoff,
        transcriptPointSize = 1.5,
        transcriptLabSize = 3.0,
        legend=c('NS','Log2 FC','Adjusted p-value',
                 'Adjusted p-value & Log2 FC'),
        legendPosition = 'bottom',
        legendLabSize = 12,
        legendIconSize = 2.0,
        colAlpha = 1.0,
        drawConnectors = TRUE)
    print(volcanoPlot)
}

main()
