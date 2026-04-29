biomex.readFile <- function(filePath, sep = "auto", colClasses = NA, nrows = -1) {
  # Read data
  data <- tryCatch({
    if (grepl(".zip$",filePath)) {
      data <- fread(cmd = paste("unzip -p",filePath), sep = sep, colClasses = colClasses, nrows = nrows)
    } else {
      data <- fread(filePath, sep = sep, colClasses = colClasses, nrows = nrows)
    }
  }, error = function(err) {
    print(err)
    return(NULL)
  })
  
  return(data)
}

biomex.loadCustomSets <- function(setFilePath, experimentInformation) {
  setTable <- biomex.readFile(filePath = setFilePath)
  
  # Invalid file
  if (is.null(setTable)) return(NULL)
  if (ncol(setTable) < 2) return(NULL)
  
  setNames <- unique(setTable[[1]])
  
  sets <- list()
  for (setName in setNames)
  {
    set <- setTable[setTable[[1]] == setName]
    set <- set[[2]]
    
    sets[setName] <- list(set)
  }
  
  return(sets)
}

biomex.performCompetitiveSetEnrichment <- function(differentialExpression, featureSets, minimumSetSize, randomSeed) {
  differentialExpression$entity <- as.character(differentialExpression$entity)
  
  # Define sets and genes
  sets <- c()
  features <- c()
  for (i in 1:length(featureSets)) {
    setName <- names(featureSets[i])
    setContent <- featureSets[[i]]
    
    sets <- c(sets, rep(setName, length(setContent)))
    features <- c(features, setContent)
  }
  
  # Creating the terms...
  TERM2NAME <- data.frame(term = names(featureSets), name = names(featureSets), stringsAsFactors = FALSE)
  TERM2FEATURE <- data.frame(term = sets, gene = features, stringsAsFactors = FALSE)
  
  features <- bitr(differentialExpression$entity, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  idx <- features %>%
    dplyr::select(SYMBOL)%>%
    duplicated() %>%
    which()
  
  if(length(idx > 0)) features <- features[-idx, ]
  
  for(i in 1:nrow(differentialExpression)) {
    differentialExpression$entity[i] <- features[features$SYMBOL == differentialExpression$entity[i], "ENTREZID"]
  }
  
  logFC <- differentialExpression$n
  names(logFC) <- differentialExpression$entity
  
  logFC <- sort(logFC, decreasing = TRUE)
  
  # Evaluate whether more than two sets are included for GSEA
  index <- which(TERM2FEATURE$gene %in% names(logFC))
  includedPerSet <- table(as.factor(TERM2FEATURE[index, 1]))
  setsIncluded <- length(which(includedPerSet >= minimumSetSize))
  
  if(setsIncluded <= 1) return (list(errorMessage = "Not enough sets"))
  
  # Competitive set enrichment analysis
  set.seed(randomSeed)
  cseaResult <- GSEA(geneList = logFC, exponent = 1, minGSSize = minimumSetSize, 
                     pvalueCutoff = 1, pAdjustMethod = 'fdr', TERM2GENE = TERM2FEATURE,
                     TERM2NAME = TERM2NAME, verbose = TRUE, seed = TRUE)
  
  # Refine results
  cseaTable <- as.data.table(cseaResult@result)
  cseaTable$Direction <- ifelse(cseaTable$NES <= 0, "Down", "Up")
  cseaTable <- cseaTable[, c("ID", "enrichmentScore", "NES", "Direction", "pvalue", "p.adjust", "core_enrichment")] 
  colnames(cseaTable) <- c("Set", "Enrichment score", "NES", "Direction", "Pvalue", "Adjusted.pvalue", "Enriched_genes")
  
  return(cseaTable)
}

# Function to convert ENTREZ IDs to Gene Symbols
convert_entrez_to_symbols <- function(entrez_str) {
  # Split the string by "/"
  entrez_ids <- unlist(strsplit(entrez_str, "/"))
  
  # Convert ENTREZ IDs to Gene Symbols
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = entrez_ids,
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")
  
  # Collapse into a single string with "/" separator
  gene_symbols <- paste(gene_symbols, collapse = "/")
  return(gene_symbols)
}
