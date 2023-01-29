#' Calculate the distance for enhancers to their nearest TSS,
#' which was borrowed from Signac
#'
#' @keywords internal
#'
DistanceToTSS <- function(peaks, genes, distance = 2e+05,
                          sep = c("-", "-")) {

  tss <- GenomicRanges::resize(x = genes, width = 1, fix = 'start') # find the TSS location
  genes.extended <- suppressWarnings(
    expr = Signac::Extend(
      x = tss, upstream = distance, downstream = distance
    )
  ) # extand the genomic range from the TSS till downstream/upstream 200000 bp
  overlaps <- GenomicAlignments::findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  ) # find the peaks overlapped with the extended genomic ranges of genes
  hit_matrix <- Matrix::sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  ) # build a sparse matrix to record the overlaps between peaks and extended genomic ranges of genes
  rownames(x = hit_matrix) <- Signac::GRangesToString(grange = peaks, sep = sep) # use peak names as the row names
  colnames(x = hit_matrix) <- genes.extended$gene_name # use gene names as the column names
  hit_matrix
}



#' Retain the longest transcripts for all exons
#'
#' @keywords internal
#' @import data.table
#'
CollapseToLongestTranscript <- function(ranges) {

  range.df <- data.table::as.data.table(x = ranges) # transform a GRanges object into a data frame
  range.df$strand <- as.character(x = range.df$strand) # transform GRanges into character strings
  range.df$strand <- ifelse( # check whether the strand information is available
    test = range.df$strand == "*", # ambiguous
    yes = "+", # treat all sequences by only using their positive strands
    no = range.df$strand # use the provided strands
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]]),
    "gene_name"
  ] # merge exons into the longest transcripts as a representative of the gene
  colnames(x = collapsed) <- c(
    "gene_name", "seqnames", "start", "end", "strand", "gene_biotype"
  )
  gene.ranges <- GenomicRanges::makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = T # the information not used to build the data frame
    # will be retained in meta data
  )
  gene.ranges
}



#' Retain genes or enhancers which are within a nearby range of at least one enhancer or gene,
#' respectively
#'
#' @importFrom dplyr %>% filter
#'
#' @keywords internal
#'
filter_nearby_genes <- function(obj, distance = 1e+6, peak.assay = "ATAC") {

  # calculate the nearby genes
  gene.coords <- CollapseToLongestTranscript(ranges = Signac::Annotation(object = obj[[peak.assay]]))
  peaks <- Signac::StringToGRanges(rownames(obj[[peak.assay]])) # get peak coordinates
  distance.df <- Matrix::summary(DistanceToTSS(peaks = peaks, genes = gene.coords,
                                       distance = distance)) # distance matrix
  peak.names <- rownames(obj[[peak.assay]]) # get peak names
  gene.names <- gene.coords$gene_name # get gene names


  data.table::rbindlist(pbmcapply::pbmclapply(1:nrow(distance.df), function(i) {
    return(list(peak = peak.names[distance.df[i, 1]], gene = gene.names[distance.df[i, 2]]))
  }, mc.cores = min(parallel::detectCores(), nrow(distance.df))), fill = T) %>%
    dplyr::filter(gene %in% rownames(obj[["RNA"]]))
}



#' Link matrix to granges
#' 
#' @keywords internal
#' 
#' @param linkmat A sparse matrix with genes in the rows and peaks in the
#' columns
#' @param gene.coords Genomic coordinates for each gene
#' @return Returns a GRanges object
#' @importFrom GenomicRanges resize start width GRanges makeGRangesFromDataFrame
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics sort
#' 
LinksToGRanges <- function(linkmat, gene.coords, sep = c("-", "-")) {
 
  # get TSS for each gene
  tss <- resize(gene.coords, width = 1, fix = 'start')
  gene.idx <- sapply(
    X = rownames(x = linkmat),
    FUN = function(x) {
      which(x = x == tss$gene_name)[[1]]
    }
  )
  tss <- tss[gene.idx]
  
  # get midpoint of each peak
  peak.ranges <- Signac::StringToGRanges(
    regions = colnames(x = linkmat),
    sep = sep
  )
  midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2)
  
  # convert to triplet form
  dgtm <- as(object = linkmat, Class = "dgTMatrix")
  
  # create dataframe
  df <- data.frame(
    chromosome = as.character(x = seqnames(x = peak.ranges)[dgtm@j + 1]),
    tss = start(x = tss)[dgtm@i + 1],
    pk = midpoints[dgtm@j + 1],
    score = dgtm@x,
    gene = rownames(x = linkmat)[dgtm@i + 1],
    peak = colnames(x = linkmat)[dgtm@j + 1]
  )
  
  # work out start and end coords
  df$start <- ifelse(test = df$tss < df$pk, yes = df$tss, no = df$pk)
  df$end <- ifelse(test = df$tss < df$pk, yes = df$pk, no = df$tss)
  df$tss <- NULL
  df$pk <- NULL
  
  # convert to granges
  gr.use <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  return(sort(x = gr.use))
}



#' Load the genomic annotations of an organism
#'
#' @keywords internal
#'
load_database <- function(org = "hg38") {

  ifelse (grepl("^mm", org) | org == "mouse", enddb <- "org.Mm.eg.db",
          enddb <- "org.Hs.eg.db")
  return(enddb)
}



#' Convert ENSEMBL IDs to gene symbols
#'
#' @keywords internal
#'
ensembl_to_symbol <- function(ensembl.ll, org = org) {

  org.db <- load_database(org = org)


  # Remove the version numbers
  key.ll <- ensembl.ll
  if (grepl("^mm", org)) {
    key.ll <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1", key.ll)
  } else {
    key.ll <- gsub("(ENSG[0-9]+)\\.[0-9]+", "\\1", key.ll)
  }


  # Map Ensembl IDs to gene symbols
  quiet(require(org.db, character.only = T))
  symbol.ll <- AnnotationDbi::mapIds(x = get(org.db),
                      keys = key.ll,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")


  return(symbol.ll)
}



#' Load the annotation database of an organism
#'
#' @keywords internal
#'
load_annotation <- function(org = "hg38") {

  ifelse (grepl("^mm", org), enddb <- "EnsDb.Mmusculus.v75",
          enddb <- "EnsDb.Hsapiens.v75")


  enddb
}



#' Generate a \code{GRanges} object composed of gene annotations
#'
#' @keywords internal
#'
build_gene_GRanges <- function(org = "hg38") {

  enddb <- load_annotation(org = org)
  require(enddb, character.only = T) # if get the error of "cli", please run "update.packages("cli")"
  annotations <- Signac::GetGRangesFromEnsDb(get(enddb))
  ensembldb::seqlevelsStyle(annotations) <- 'UCSC'
  gene.gr <- CollapseToLongestTranscript(annotations)

  return(gene.gr)
}



#' Link enhancers to genes within a distance cutoff
#'
#' @importFrom Matrix summary
#' @importFrom dplyr %>%
#'
#' @keywords internal
#'
link_peaks_to_genes <- function(peak.obj = c("chrX-192989-220023", "chr2-178095031-178129859"),
                                gene.obj = c("PLCXD1", "NFE2L2"),
                                org = "hg38", distance = 5e+05) {

  peak.gr <- Signac::StringToGRanges(peak.obj)
  gene.annotation <- build_gene_GRanges(org = org)
  if (grepl("^ENSG|ENSMUG", gene.annotation$gene_name[1])) {
    symbol.ll <- ensembl_to_symbol(ensembl.ll = gene.annotation$gene_name,
                                   org = org)
    gene.gr <- gene.annotation[which(gene.annotation$gene_name %in% symbol.ll)]
  } else {
    gene.gr <- gene.annotation[which(gene.annotation$gene_name %in% gene.obj)]
  }


  # Link peaks to genes
  if (is.numeric(distance)) {
    message ("Finding nearby genes for each peak within ", distance, " bp ...")
    summ <- DistanceToTSS(peaks = peak.gr, genes = gene.gr,
                          distance = distance, sep = c("-", "-")) %>% summary
  } else if (distance == "gene") { # to-do : find the closest gene for each peak using Signac
    message ("Finding the closest peak for each gene ...")
    summ <- Signac::nearest(gene.gr, peak.gr)
    summ <- summ[which(!is.na(summ))]
    summ <- data.frame(i = summ, j = seq_along(summ), x = rep(1, length(summ)))
  }
  cis.gr <- peak.gr[summ$i]
  GenomicRanges::mcols(cis.gr)$gene <- gene.gr$gene_name[summ$j]


  return(cis.gr)
}



#' Link peaks to genes using \code{Signac} functions
#'
#' @keywords internal
#' 
link_peaks <- function(object, peak.assay = "ATAC", expression.assay = "RNA", 
                        peak.slot = "counts", 
                        expression.slot = "data", method = "pearson", gene.coords = NULL, 
                        distance = 5e+05, min.distance = NULL, min.cells = 10, genes.use = NULL, 
                        n_sample = 200, pvalue_cutoff = 0.05, score_cutoff = 0.05, 
                        gene.id = FALSE, verbose = TRUE) {
  
  if (!requireNamespace(package = "qlcMatrix", quietly = TRUE)) {
    stop("Please install qlcMatrix: install.packages('qlcMatrix')")
  }
  if (!inherits(x = object[[peak.assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  if (!is.null(x = min.distance)) {
    if (!is.numeric(x = min.distance)) {
      stop("min.distance should be a numeric value")
    }
    if (min.distance < 0) {
      warning("Requested a negative min.distance value, setting min.distance to zero")
      min.distance <- NULL
    } else if (min.distance == 0) {
      min.distance <- NULL
    }
  }
  features.match <- c("GC.percent", "count", "sequence.length")
  if (method == "pearson") {
    cor_method <- qlcMatrix::corSparse
  } else if (method == "spearman") {
    cor_method <- SparseSpearmanCor
  } else {
    stop("method can be one of 'pearson' or 'spearman'.")
  }
  if (is.null(x = gene.coords)) {
    annot <- Signac::Annotation(object = object[[peak.assay]])
    if (is.null(x = annot)) {
      stop("Gene annotations not found")
    }
    gene.coords <- CollapseToLongestTranscript(ranges = annot)
  }
  meta.features <- Seurat::GetAssayData(object = object, assay = peak.assay, 
                                slot = "meta.features")
  if (!(all(c("GC.percent", "sequence.length") %in% colnames(x = meta.features)))) {
    stop("DNA sequence information for each peak has not been computed.\n", 
         "Run RegionsStats before calling this function.")
  }
  if (!("count" %in% colnames(x = meta.features))) {
    data.use <- Seurat::GetAssayData(object = object[[peak.assay]], 
                             slot = "counts")
    hvf.info <- Signac::FindTopFeatures(object = data.use, verbose = FALSE)
    hvf.info <- hvf.info[rownames(meta.features), , drop = FALSE]
    meta.features <- cbind(meta.features, hvf.info)
  }
  peak.data <- Seurat::GetAssayData(object = object, assay = peak.assay, 
                            slot = peak.slot)
  expression.data <- Seurat::GetAssayData(object = object, assay = expression.assay, 
                                  slot = expression.slot)
  peakcounts <- rowSums(x = peak.data > 0)
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, , drop = FALSE]
  if (!is.null(x = genes.use)) {
    genes.keep <- intersect(x = names(x = genes.keep[genes.keep]), 
                            y = genes.use)
  }
  expression.data <- expression.data[genes.keep, , drop = FALSE]
  if (verbose) {
    message("Testing ", nrow(x = expression.data), " genes and ", 
            sum(peaks.keep), " peaks")
  }
  genes <- rownames(x = expression.data)
  if (gene.id) {
    gene.coords.use <- gene.coords[gene.coords$gene_id %in% 
                                     genes, ]
    gene.coords.use$gene_name <- gene.coords.use$gene_id
  } else {
    gene.coords.use <- gene.coords[gene.coords$gene_name %in% 
                                     genes, ]
  }
  if (length(x = gene.coords.use) == 0) {
    stop("Could not find gene coordinates for requested genes")
  }
  if (length(x = gene.coords.use) < nrow(x = expression.data)) {
    message("Found gene coordinates for ", length(x = gene.coords.use), 
            " genes")
  }
  peaks <- Signac::StringToGRanges(rownames(object[[peak.assay]])) # Author added this line
  # peaks <- Signac::granges(x = object[[peak.assay]])
  peaks <- peaks[peaks.keep]
  peak_distance_matrix <- DistanceToTSS(peaks = peaks, genes = gene.coords.use, 
                                        distance = distance)
  if (!is.null(x = min.distance)) {
    peak_distance_matrix_min <- DistanceToTSS(peaks = peaks, 
                                              genes = gene.coords.use, distance = min.distance)
    peak_distance_matrix <- peak_distance_matrix - peak_distance_matrix_min
  }
  if (sum(peak_distance_matrix) == 0) {
    stop("No peaks fall within distance threshold\n", "Have you set the proper genome and seqlevelsStyle for ", 
         peak.assay, " assay?")
  }
  genes.use <- colnames(x = peak_distance_matrix)
  all.peaks <- rownames(x = peak.data)
  message ("Dimensions of peak.data:", nrow(peak.data), " x ", ncol(peak.data))
  peak.data <- t(x = peak.data)
  coef.vec <- c()
  gene.vec <- c()
  zscore.vec <- c()
  # future::nbrOfWorkers() <- future::availableCores() # Author added this line
  if (T) {
  # if (future::nbrOfWorkers() > 1) {
    # invisible(require(future.apply))
    mylapply <- pbmcapply::pbmclapply
  } else {
    mylapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
  }
  res <- mylapply(X = seq_along(along.with = genes.use), mc.cores = max(parallel::detectCores() / 2, 
                                                                        1),
                  FUN = function(i) {
    peak.use <- as.logical(x = peak_distance_matrix[, genes.use[[i]]])
    gene.expression <- t(x = expression.data[genes.use[[i]], 
                                             , drop = FALSE])
    gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
    if (sum(peak.use) < 2) {
      return(list(gene = NULL, coef = NULL, zscore = NULL))
    }
    else {
      peak.access <- peak.data[, peak.use, drop = FALSE]
      coef.result <- cor_method(X = peak.access, Y = gene.expression)
      rownames(x = coef.result) <- colnames(x = peak.access)
      coef.result <- coef.result[abs(x = coef.result) > 
                                   score_cutoff, , drop = FALSE]
      if (nrow(x = coef.result) == 0) {
        return(list(gene = NULL, coef = NULL, zscore = NULL))
      }
      else {
        peaks.test <- rownames(x = coef.result)
        trans.peaks <- all.peaks[!grepl(pattern = paste0("^", 
                                                         gene.chrom), x = all.peaks)]
        meta.use <- meta.features[trans.peaks, ]
        pk.use <- meta.features[peaks.test, ]
        bg.peaks <- lapply(X = seq_len(length.out = nrow(x = pk.use)), 
                           FUN = function(x) {
                             Signac::MatchRegionStats(meta.feature = meta.use, 
                                              query.feature = pk.use[x, , drop = FALSE], 
                                              features.match = features.match, n = n_sample, 
                                              verbose = FALSE)
                           })
        bg.access <- peak.data[, unlist(x = bg.peaks), 
                               drop = FALSE]
        bg.coef <- cor_method(X = bg.access, Y = gene.expression)
        rownames(bg.coef) <- colnames(bg.access)
        zscores <- vector(mode = "numeric", length = length(x = peaks.test))
        for (j in seq_along(along.with = peaks.test)) {
          coef.use <- bg.coef[(((j - 1) * n_sample) + 
                                 1):(j * n_sample), ]
          z <- (coef.result[j] - mean(x = coef.use))/sd(x = coef.use)
          zscores[[j]] <- z
        }
        names(x = coef.result) <- peaks.test
        names(x = zscores) <- peaks.test
        zscore.vec <- c(zscore.vec, zscores)
        gene.vec <- c(gene.vec, rep(i, length(x = coef.result)))
        coef.vec <- c(coef.vec, coef.result)
      }
      gc(verbose = FALSE)
      pval.vec <- pnorm(q = -abs(x = zscore.vec))
      links.keep <- pval.vec < pvalue_cutoff
      if (sum(x = links.keep) == 0) {
        return(list(gene = NULL, coef = NULL, zscore = NULL))
      }
      else {
        gene.vec <- gene.vec[links.keep]
        coef.vec <- coef.vec[links.keep]
        zscore.vec <- zscore.vec[links.keep]
        return(list(gene = gene.vec, coef = coef.vec, 
                    zscore = zscore.vec))
      }
    }
  })
  gene.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                              1))
  coef.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                              2))
  zscore.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                                3))
  if (length(x = coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    return(object)
  }
  peak.key <- seq_along(along.with = unique(x = names(x = coef.vec)))
  names(x = peak.key) <- unique(x = names(x = coef.vec))
  coef.matrix <- Matrix::sparseMatrix(i = gene.vec, j = peak.key[names(x = coef.vec)], 
                              x = coef.vec, dims = c(length(x = genes.use), max(peak.key)))
  rownames(x = coef.matrix) <- genes.use
  colnames(x = coef.matrix) <- names(x = peak.key)
  links <- LinksToGRanges(linkmat = coef.matrix, gene.coords = gene.coords.use)
  z.matrix <- Matrix::sparseMatrix(i = gene.vec, j = peak.key[names(x = zscore.vec)], 
                           x = zscore.vec, dims = c(length(x = genes.use), max(peak.key)))
  rownames(x = z.matrix) <- genes.use
  colnames(x = z.matrix) <- names(x = peak.key)
  z.lnk <- LinksToGRanges(linkmat = z.matrix, gene.coords = gene.coords.use)
  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(q = -abs(x = links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]
  Signac::Links(object = object[[peak.assay]]) <- links
  return(object)
}
