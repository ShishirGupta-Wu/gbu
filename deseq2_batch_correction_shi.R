library('argparse')
library('dplyr')
library('stringr')
library('limma')
library('DESeq2')
library('ggplot2')
library('ggrepel')
library('ggpubr')
library('EnhancedVolcano')
library('RColorBrewer')
library('pheatmap')

main <- function(){
    options(warn=1)
    # Read command line parameters
    parser <- ArgumentParser()
    parser$add_argument("--geneCountingFile")
    parser$add_argument("--annoCols", type="integer", default=10)
    parser$add_argument("--libSelectRegex", required=TRUE)
    parser$add_argument("--condRegex", required=FALSE)
    parser$add_argument("--comparisons", nargs="+", required=TRUE)
    parser$add_argument("--countTansformMethod", choices=c('normTransform', 'vst', 'rlog'))
    parser$add_argument("--log2FcCutoff", type="double", default=1.0)
    parser$add_argument("--PadjCutoff", type="double", default=0.05)
    parser$add_argument("--outputFolder")
    args <- parser$parse_args()
    debug <- TRUE

    rawCountTable <- read.table(args$geneCountingFile, sep='\t', header=TRUE,
                                quote='', comment.char='')

    # Select the pure counting matrix
    countTable <- round(select(rawCountTable, matches(args$libSelectRegex)))
    libs <- colnames(countTable)
	libs
    annotationTable <- rawCountTable[, 1:args$annoCols]
    if(debug){print(colnames(countTable))}
    if(debug){print(colnames(annotationTable))}

    # calculate TPM values
    tpmTable <- apply(countTable, 2, tpm, rawCountTable$Combined_exon_length)
    write.table(cbind(annotationTable, tpmTable), file=sprintf("%s/TPM_count_table.csv", args$outputFolder),
        quote=FALSE, sep='\t', row.names=FALSE)

    # Derive the condition and replicate names from lib names
    # conds <- gsub("_RNA[0-9]+$", "", libs)
  
    conds <- c('Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature', 'Immature',  'Immature', 'Immature', 'Mature', 'Mature', 'Mature', 'Mature', 'Mature','Mature', 'Mature', 'Mature', 'Mature', 'Mature', 'Mature', 'Mature', 'Mature', 'Mature', 'Mature', 'Mature')
		
    #  batch <- str_extract(libs, "_[1-16]$")
	
	####IMPORTANT##### - This is an example of 4 batches
	######Assigning batch based on batch affected PCA and from sample details. Please check header by     head -n 1 output/featureCounts/Reads_per_gene.csv  and assign batch number accordingly in order of Reads_per_gene.csv
    #	Immature_10     Immature_11     Immature_12     Immature_13     Immature_14 Immature_15      Immature_16     Immature_1      Immature_2      Immature_3      Immature_4      Immature_5      Immature_6      Immature_7      Immature_8      Immature_9      Mature_10   Mature_11        Mature_12       Mature_13       Mature_14       Mature_15       Mature_16       Mature_1        Mature_2        Mature_3        Mature_4        Mature_5        Mature_6    Mature_7 Mature_8        Mature_9

    batch <- c('3', '3', '3', '4', '4', '4', '4', '1', '1', '1', '1', '2', '2', '2', '2', '3', '3', '3', '3', '4', '4', '4', '4', '1', '1', '1', '1', '2', '2', '2', '2', '3') 
	
		
    samples <- data.frame(row.names=libs, condition=conds, lib=libs, batch=batch)
    if(debug){print(samples)}

    # generate DESeq2 data set
    dds <- DESeqDataSetFromMatrix(
        countData=countTable, colData=samples, design=~batch+condition)
    dds <- DESeq(dds, betaPrior=TRUE)

    # plot genes of interest
    # plotGenesOfInterest(dds, rawCountTable, genesOfInterest, args$outputFolder)

    # calculated data transformation
    transformed <- get(args$countTansformMethod)(dds)

    # plot genes of interest heatmap
    # plotGenesOfInterestHeatmap(transformed, profiles, annotationTable,
    #                            args$outputFolder)

    # plot PCA and sample distance heatmap
    createPcaAndHeatmap(transformed, args$countTansformMethod, batch, args$outputFolder, lab=libs)

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

    # pairwise comparisons
    for (comp in strsplit(args$comparison, ":", fixed=TRUE)){
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

createPcaAndHeatmap <- function(transformed, countTansformMethod, batch, outputFolder, lab=NULL){
    pdf(sprintf("%s/DESeq2_PCA_and_heatmap_%s.pdf", outputFolder, countTansformMethod))
    # PCA
    # assay(transformed) <- removeBatchEffect(assay(transformed), batch)
    transformed_batch_corrected <- transformed
    assay(transformed_batch_corrected) <- removeBatchEffect(assay(transformed_batch_corrected), batch)
    pcaData <- plotPCA(transformed, intgroup=c("condition", "batch"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
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
    pcaData <- plotPCA(transformed_batch_corrected, intgroup=c("condition", "batch"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=batch)) +
          geom_point(size=3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance")) +
          coord_fixed() +
          ggtitle(sprintf("PCA (%s-transformed counts\nbatch corrected)", countTansformMethod)) +
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
    sampleDists <- dist(t(assay(transformed_batch_corrected)))
    sampleDistMatrix <- as.matrix(sampleDists)
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors,
             main="Sample distance heatmap batch corrected")
    dev.off()
}

compareConditions <- function(condComp, condRef, dds, rawCountTable,
                              normCountTable, log2FcCutoff, PadjCutoff,
                              outputFolder){
    comp <- results(dds, contrast=c('condition', condComp, condRef), addMLE=TRUE)
    countingAndComparisonResults <- cbind(rawCountTable, normCountTable, comp)
    write.table(countingAndComparisonResults, file=sprintf(
                          "%s/DESeq2_comparison_%s_vs_%s_table.csv",
                          outputFolder, condComp, condRef),
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

plotGenesOfInterestHeatmap <- function(transformed, profiles, annotationTable,
                                       outputFolder){
    pdf(sprintf("%s/genes_of_interest_heatmaps.pdf", outputFolder), width=12)
    transformedDF <- assay(transformed)
    transformedDF <- as.data.frame(transformedDF)
    transformedDF$Geneid <- annotationTable$gene
    transformedDF$id <- annotationTable$ID
    transformedDF <- na.omit(transformedDF)
    profileCols <- colnames(profiles)
    profileList <- c()
    values <- matrix(ncol = length(colnames(transformedDF)),
                     dimnames=list(c(), colnames(transformedDF)), byrow=TRUE)
    # for (curProfile in profileCols){
        curProfile <- "Megakaryocyte"
        geneList <- profiles[curProfile][,1]
        curVals <- transformedDF[transformedDF$Geneid %in% geneList,]
        profileList <- c(profileList, rep(curProfile, length(curVals$Geneid)))
        values <- rbind(values, curVals)

        values_df <- data.frame(curVals)
        values_df <- values_df[-c(1),]

        rownames(values_df)<-values_df$Geneid
        values_df$Geneid <- NULL
        values_df$id <- NULL
        pheatmap(as.matrix(values_df), main = curProfile, cluster_rows = FALSE,
                 cluster_cols = FALSE, scale = "row")
    # }
    # values <- data.frame(values)
    # values <- values[-c(1),]
    #
    # rownames(values)<-values$id
    # profileList <- data.frame(profileList)
    # rownames(profileList) <- values$id
    # values$Geneid <- NULL
    # values$id <- NULL
    #
    # colorsList <- c("#e6194B", "#f58231","#ffe119", "#bfef45", "#3cb44b",
    #                 "#42d4f4", "#4363d8", "#911eb4", "#f032e6", "#800000",
    #                 "#9A6324", "#808000", "#469990")
    # names(colorsList) <- profileCols
    # colPal <- list(profileList = colorsList)
    #
    # pheatmap(as.matrix(values), annotation_colors=colPal, show_rownames = FALSE,
    #          cellwidth = 20, annotation_row=profileList, cluster_rows = FALSE,
    #          cluster_cols=FALSE, scale = "row")
    dev.off()
}

plotGenesOfInterest <- function(dds, rawCountTable, genesOfInterest, outputFolder){
    pdf(sprintf("%s/genes_of_interest.pdf", outputFolder))
    for (gene in genesOfInterest){
        plotData <- plotCounts(dds, which(rawCountTable$gene==gene),
                         intgroup = "condition", returnData = TRUE)
        g <- ggplot(plotData, aes(x = condition, y = count)) +
             geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
        g <- g + stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1),
                              geom="crossbar", color="red", width=0.5)
        g <- g + coord_flip() +
                 ggtitle(sprintf("Gene %s", gene)) +
                 theme(plot.title = element_text(color="black", size=18, face="bold"))
        print(g)
    }
    dev.off()
}

main()
