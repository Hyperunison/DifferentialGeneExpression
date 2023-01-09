###############################################################################
## Differential gene expression module                                       ##

# This module reads 
#   * a design file
#   * a model file 
#   * a RSEM read count matrix (rnaseq pipeline output)

# and produces the differential gene expression analyses specified in the 
# design/model files. 

# Outputs: DEseqw differential gene epxression output, written out as a text 
# file. 


## load relevant files
# Design file
FNdesign <- "design.table.tsv"

dfDesign <- read.delim(
    file = FNdesign,
    sep = "\t",
    stringsAsFactors = F
)

# Model file
FNmodel <- "model.table.tsv"

dfModel <- read.delim(
    file = FNmodel,
    sep = "\t",
    stringsAsFactors = F
)

# RSEM count file
FNrsem <- "rsem.merged.gene_counts.tsv"

dfRSEM <- read.delim(
    file = FNrsem,
    sep = "\t",
    stringsAsFactors = F
)

# the count file needs to be made into an integer count matrix
mRSEM <- dfRSEM
row.names(mRSEM) <- mRSEM$gene_id
mRSEM$gene_id <- NULL
mRSEM$transcript_id.s. <- NULL

# make sure it's integer, as required by DESeq2
mRSEM <- data.matrix(
    round(
        mRSEM
    )
)

# TPM file
FNtpm <- "rsem.merged.gene_tpm.tsv"

dfTPM <- read.delim(
    file = FNtpm,
    sep = "\t",
    stringsAsFactors = F
)

# setting parameters
primaryAlignmentGeneID <- "gene_id"

DEseq2contrastTable <- data.frame(NULL)

## Obtain gene annotation
# The variable below need to be set in a project-specific fashion
# Required if geneIDs, e.g. Ensembl have been used for alignment
# In this case, the NCBI -GRCh38 default genome was used and annotation is in
# gene names already.

## Do DGE analysis
dfDGE <- dfModel

dfDGE <- dfDGE[dfDGE$test == "Wald",]

for (i in 1:nrow(dfDGE)){
    designFormula <- as.formula(as.vector(dfDGE[i, "model"]))
    addCols <- as.vector(dfDGE[i, "model"])
    addCols <- gsub("~", "", addCols)
    addCols <- unlist(strsplit(addCols, "[+]"))
    addCols <- gsub(" ", "", addCols)
    addCols <- sapply(addCols, function(x) unlist(strsplit(x, ":")))
    addCols <- unique(Reduce(c, addCols))

    addCols <- addCols[addCols != "condition"]

    selCols <- c("sample.id", as.vector(dfDGE[i,"comparisonID"]), addCols,"sample.group")

    colData = unique(dfDesign[, selCols])
    names(colData) <- gsub(as.vector(paste0("^",dfDGE[i,"comparisonID"], "$")), "condition", names(colData))
    rownames(colData) = as.vector(colData$sample.id)
    colData$sample.id <- NULL



    if (!as.logical(dfDGE[i,"normalizeAllSamplesTogether"])) {
        colData = droplevels(data.frame(colData[colData$condition != "",]))
        colData <- colData[order(colData$condition),]
    } else {
        colData[colData$condition == "", "condition"] <- "rest"
    }


    contrasts = sort(as.vector(unique(colData[,"condition"])), decreasing = FALSE)
    contrasts = contrasts[contrasts != "rest"]
    contrasts = contrasts[contrasts != ""]

    colData$condition <- gsub("^1_", "", colData$condition)
    colData$condition <- gsub("^2_", "", colData$condition)


    fCols <- c("condition", addCols)
    for (j in 1:length(fCols)){
        if (!is.numeric(colData[,fCols[j]])){
            colData[,fCols[j]] <- as.factor(colData[,fCols[j]])
        }

    }

    colData[,"condition"] = as.factor(colData[,"condition"])
    colData$sample.group <- as.factor(colData$sample.group)




    contrasts <- gsub("^1_", "", contrasts)
    contrasts <- gsub("^2_", "", contrasts)

    #Create contrast vector
    #contrast.vector = c([condition],[1_diff.gene set, e.g. mt],[2_baseline, e.g. wt])
    #if (contrasts[2] != "scr"){
    #  contrasts = rev(contrasts)
    #}
    sel.col = contrasts

    contrast.vector = append("condition", contrasts)
    colName = paste(contrasts, collapse = "_vs_")
    #colName

    if (as.logical(dfDGE[i,"normalizeAllSamplesTogether"])) {
        raw.counts.temp = mRSEM
    } else {
        raw.counts.temp = mRSEM[,rownames(colData)]
    }

    colData$condition <- as.factor(colData$condition)


    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = raw.counts.temp[,row.names(colData)],
        colData   = colData,
        design    = designFormula
    )

    #dds$condition <- factor(dds$condition, levels=contrasts)
    if (as.vector(dfDGE[i, "betaPrior"]) == "TRUE"){
        betaPrior <- TRUE
    } else {
        betaPrior <- FALSE
    }

    #temporary fix
    parallelProcessing <- F

    # deseq2 processing
    dds <- DESeq2::DESeq(
        dds,
        test = as.vector(dfDGE[i, "test"]),
        parallel = parallelProcessing,
        betaPrior = betaPrior
    )

    res <- DESeq2::results(dds, contrast = contrast.vector)
    #https://support.bioconductor.org/p/83773/
    #res <- results(dds, contrast=list("conditioncell_type_A","conditioncell_type_B"))



    res = data.frame(res)
    names(res) = paste(names(res), colName, sep="_")
    res[[primaryAlignmentGeneID]] = rownames(res)


    names(res) = gsub("log2FoldChange", "logFC", names(res))
    names(res) = gsub(
        "logFC",
        paste("contrast_", i, "_logFC", sep=""),
        names(res)
    )

    names(res) = gsub(
        "padj",
        paste("contrast_", i, "_padj", sep=""),
        names(res)
    )

    names(res) = gsub(
        "stat",
        paste("contrast_", i, "_stat", sep=""),
        names(res)
    )

    res$baseMean <- log2(res$baseMean)
    names(res) = gsub(
        "baseMean",
        paste("contrast_", i, "_lg2BaseMean", sep=""),
        names(res)
    )

    #Remove all rows without a padj
    padj.col = grep("padj", names(res))[1]
    res[,padj.col][is.na(res[,padj.col])] = ""
    res = res[res[,padj.col] != "", ]
    res[,padj.col] <- as.numeric(res[,padj.col])

    ## Add log10p column ##
    padj  <- names(res)[grep("_padj_", names(res))]
    lg10p <- gsub("padj", "lg10p", padj)

    for (z in 1:length(padj)){
        preprocess <- as.numeric(res[,padj[z]])
        minNum <- min(preprocess[preprocess != 0])
        preprocess[preprocess == 0] <- minNum

        # if (length(grep("padj_LRT", padj[i])) > 0){
        #     preprocess <- as.numeric(res[,padj[z]])
        #     minNum <- min(preprocess[preprocess != 0])
        #     preprocess[preprocess == 0] <- minNum
        # } else {
        #     preprocess <- as.numeric(res[,padj[z]])
        # }

        temp <- -1*log10(preprocess)
        #temp[temp >= 50] = 50
        res[,lg10p[z]] <- temp
    }

    col.vector = c(
        primaryAlignmentGeneID,
        names(res)[grep("contrast", names(res))]
    )

    res = res[,col.vector]

    ## Make all numeric columns numeric ##
    res[,grep("contrast_", names(res))] <- apply(res[,grep("contrast_", names(res))], 2, as.numeric)



    ## Add to result array ##
    if (nrow(DEseq2contrastTable) == 0){
        DEseq2contrastTable <- res
    } else {
        DEseq2contrastTable <- merge(
            DEseq2contrastTable,
            res,
            by.x = primaryAlignmentGeneID,
            by.y = primaryAlignmentGeneID,
            all = TRUE
        )
        DEseq2contrastTable[is.na(DEseq2contrastTable)] <- 0
    }

} ## End of DGE loop

## Write differential gene expression result table to file
FNdgeOut <- "DGEResulTtable.tsv"
write.table(
    DEseq2contrastTable,
    FNdgeOut,
    sep = "\t",
    row.names = F
)




