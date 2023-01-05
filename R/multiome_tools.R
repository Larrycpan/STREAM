#' Subset a \code{Seurat} object
#'
#' @import dplyr
#'
#' @keywords internal
#'
subset_object <- function(LTMG.obj, object, peak.assay = 'ATAC',
                          atac.dis, max.peaks = 3000, links.df = NULL,
                          n.blocks = 100) {

  atac.dis <- atac.dis[intersect(unique(links.df[links.df$gene %in%
                                                   rownames(object[["RNA"]])]$peak),
                                 rownames(atac.dis)),]
  block.list <- split(LTMG.obj@BiCluster@CoCond_cell,
                      f = LTMG.obj@BiCluster@CoCond_cell$Condition) %>%
    sapply(., "[[", "cell_name") %>% head(n = n.blocks)


  return( parallel::mclapply(block.list, function(block) {
    n.cutoff <- 0
    r.sum <- Matrix::rowSums(atac.dis[, block]) # rowwise sums
    r.sum <- r.sum[which(r.sum > n.cutoff)] # remove the zero sums
    r.sum <- r.sum[order(r.sum, decreasing = T)] # rank the row sums
    block.peaks <- names(r.sum) # get the consistent peaks
    if (length(block.peaks) < 1) {
      return(NULL)
    }
    if (length(block.peaks) > max.peaks) {
      block.peaks <- block.peaks[1:max.peaks]
    }

    return(list(peaks = block.peaks, cells = block))
    # return(subset(x = object, features = c(rownames(object[['RNA']]), block.peaks),
    #               cells = block))
  }, mc.cores = parallel::detectCores()) )
}



#' Build a list of heterogeneous graphs based on the input list of \code{Seurat} objects
#'
#' Function org_to_DB used to get the database based on organism
#'
#' @keywords internal
#'
org_to_DB <- function(org = "hg38") {
  
  message ("The organism version is ", org, ".")
  
  
  # Return the database
  if (grepl("mm", org)) {
    message ("Loading the database EnsDb.Mmusculus.v79 ...")
    require(EnsDb.Mmusculus.v79)
    return(EnsDb.Mmusculus.v79)
  } else {
    message ("Loading the database EnsDb.Hsapiens.v86 ...")
    require(EnsDb.Hsapiens.v86)
    return(EnsDb.Hsapiens.v86)
  }
}



generate_windows <- function(window, genomic_coords) {
  
  if (!is(genomic_coords, "data.frame")) {
    chr_maxes <- read.table(genomic_coords)
  }
  else {
    chr_maxes <- genomic_coords
  }
  names(chr_maxes) <- c("V1", "V2")
  win_ranges <- plyr::ddply(chr_maxes, plyr::.(V1), function(x) {
    r <- seq(from = 1, to = x$V2[1], by = window/2)
    l <- r + window - 1
    data.frame(start = r, end = l)
  })
  gr <- GenomicRanges::GRanges(win_ranges$V1,
                               ranges = IRanges::IRanges(win_ranges$start,
                                                         win_ranges$end))
  
  
  return(gr)
}



get_genomic_range <- function(grs, cds, win) {
  
  end1 <- as.numeric(as.character(GenomicRanges::end(grs[win])))
  end2 <- as.numeric(as.character(GenomicRanges::start(grs[win])))
  win_range <- cds[(monocle3::fData(cds)$bp1 < end1 & monocle3::fData(cds)$bp1 >
                      end2) | (monocle3::fData(cds)$bp2 < end1 & monocle3::fData(cds)$bp2 > end2),
  ]
  win_range <- win_range[as.character(monocle3::fData(win_range)$chr) ==
                           gsub("chr", "", as.character(GenomicRanges::seqnames(grs[win]))),
  ]
  monocle3::fData(win_range)$mean_bp <- (as.numeric(as.character(monocle3::fData(win_range)$bp1)) +
                                           as.numeric(as.character(monocle3::fData(win_range)$bp2)))/2
  
  
  return(win_range)
}



find_distance_parameter <- function(dist_mat,
                                    gene_range,
                                    maxit,
                                    null_rho,
                                    s,
                                    distance_constraint,
                                    distance_parameter_convergence) {
  
  if (sum(dist_mat > distance_constraint) / 2 < 1) {
    return("No long edges")
  }
  
  found <- FALSE
  starting_max <- 2
  distance_parameter <- 2
  distance_parameter_max <- 2
  distance_parameter_min <- 0
  it <- 0
  while(found != TRUE & it < maxit) {
    vals <- monocle3::exprs(gene_range)
    # message (vals, "\n")
    cov_mat <- cov(base::as.data.frame(t(vals)))
    diag(cov_mat) <- Matrix::diag(cov_mat) + 1e-4
    
    rho <- get_rho_mat(dist_mat, distance_parameter, s)
    
    GL <- glasso::glasso(cov_mat, rho)
    big_entries <- sum(dist_mat > distance_constraint)
    
    if (((sum(GL$wi[dist_mat > distance_constraint] != 0) / big_entries) > 0.05) |
        (sum(GL$wi == 0) / (nrow(GL$wi) ^ 2) < 0.2 ) ) {
      longs_zero <- FALSE
    } else {
      longs_zero <- TRUE
    }
    
    if (longs_zero != TRUE | (distance_parameter == 0)) {
      distance_parameter_min <- distance_parameter
    } else {
      distance_parameter_max <- distance_parameter
    }
    new_distance_parameter <- (distance_parameter_min +
                                 distance_parameter_max)/2
    
    if(new_distance_parameter == starting_max) {
      new_distance_parameter <- 2 * starting_max
      starting_max <- new_distance_parameter
    }
    
    if (distance_parameter_convergence > abs(distance_parameter -
                                             new_distance_parameter)) {
      found <- TRUE
    } else {
      distance_parameter <- new_distance_parameter
    }
    it <- it + 1
  }
  if (maxit == it) warning ("maximum iterations hit")
  
  
  return(distance_parameter)
}



calc_dist_matrix <- function(gene_range) {
  
  dist_mat <- Matrix::as.matrix(stats::dist(monocle3::fData(gene_range)$mean_bp))
  row.names(dist_mat) <- colnames(dist_mat) <- row.names(monocle3::fData(gene_range))
  
  
  return(dist_mat)
}



estimate_distance_parameter <- function(cds,
                                        window=500000,
                                        maxit=100,
                                        s=0.75,
                                        sample_num = 100,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        max_elements = 200,
                                        genomic_coords = NULL,
                                        max_sample_windows = 500) {
  
  # require(monocle3)
  # require(SummarizedExperiment)
  # require(SingleCellExperiment)
  
  
  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(assertthat::is.number(window))
  assertthat::assert_that(assertthat::is.count(maxit))
  assertthat::assert_that(assertthat::is.number(s), s < 1, s > 0)
  assertthat::assert_that(assertthat::is.count(sample_num))
  assertthat::assert_that(assertthat::is.count(distance_constraint))
  assertthat::assert_that(distance_constraint < window)
  assertthat::assert_that(assertthat::is.number(distance_parameter_convergence))
  if (!is.data.frame(genomic_coords)) {
    assertthat::is.readable(genomic_coords)
  }
  assertthat::assert_that(assertthat::is.count(max_sample_windows))
  
  grs <- generate_windows(window, genomic_coords)
  
  monocle3::fData(cds)$chr <- gsub("chr", "", monocle3::fData(cds)$chr)
  monocle3::fData(cds)$bp1 <- as.numeric(as.character(monocle3::fData(cds)$bp1))
  monocle3::fData(cds)$bp2 <- as.numeric(as.character(monocle3::fData(cds)$bp2))
  
  distance_parameters <- list()
  distance_parameters_calced <- 0
  it <- 0
  
  while(sample_num > distance_parameters_calced & it < max_sample_windows) {
    # message (it, "\n")
    it <- it + 1
    win <- sample(seq_len(length(grs)), 1)
    # message (win, "\n")
    GL <- "Error"
    win_range <- get_genomic_range(grs, cds, win)
    
    if (nrow(monocle3::exprs(win_range))<=1) {
      next()
    }
    if (nrow(monocle3::exprs(win_range)) > max_elements) {
      next()
    }
    
    dist_matrix <- calc_dist_matrix(win_range)
    
    distance_parameter <- find_distance_parameter(dist_matrix,
                                                  win_range,
                                                  maxit = maxit,
                                                  null_rho = 0,
                                                  s,
                                                  distance_constraint = distance_constraint,
                                                  distance_parameter_convergence =
                                                    distance_parameter_convergence)
    
    if (!is(distance_parameter, "numeric")) next()
    distance_parameters = c(distance_parameters, distance_parameter)
    distance_parameters_calced <- distance_parameters_calced + 1
  }
  
  if(length(distance_parameters) < sample_num)
    warning(paste0("Could not calculate sample_num distance_parameters (",
                   length(distance_parameters), " were calculated) - see ",
                   "documentation details"))
  if(length(distance_parameters) == 0)
    stop("No distance_parameters calculated")
  
  unlist(distance_parameters)
}



get_rho_mat <- function(dist_matrix, distance_parameter, s) {
  xmin <- 1000
  out <- (1-(xmin/dist_matrix)^s) * distance_parameter
  out[!is.finite(out)] <- 0
  out[out < 0] <- 0
  
  
  out
}



generate_cicero_models <- function(cds, distance_parameter, s = 0.75, window = 5e+05,
                                   max_elements = 200, genomic_coords = NULL) {
  
  # require(monocle3)
  # require(SummarizedExperiment)
  # require(SingleCellExperiment)
  
  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(assertthat::is.number(distance_parameter))
  assertthat::assert_that(assertthat::is.number(s), s < 1,
                          s > 0)
  assertthat::assert_that(assertthat::is.number(window))
  assertthat::assert_that(assertthat::is.count(max_elements))
  if (!is.data.frame(genomic_coords)) {
    assertthat::is.readable(genomic_coords)
  }
  grs <- generate_windows(window, genomic_coords)
  monocle3::fData(cds)$chr <- gsub("chr", "", monocle3::fData(cds)$chr)
  monocle3::fData(cds)$bp1 <- as.numeric(as.character(monocle3::fData(cds)$bp1))
  monocle3::fData(cds)$bp2 <- as.numeric(as.character(monocle3::fData(cds)$bp2))
  outlist <- pbmcapply::pbmclapply(seq_len(length(grs)), 
                                   mc.cores = max(parallel::detectCores() / 2, 1),
                                   function(win) {
                                     GL <- "Error"
                                     win_range <- get_genomic_range(grs, cds, win)
                                     if (nrow(monocle3::exprs(win_range)) <= 1) {
                                       return("Zero or one element in range")
                                     }
                                     if (nrow(monocle3::exprs(win_range)) > max_elements) {
                                       return("Too many elements in range")
                                     }
                                     dist_matrix <- calc_dist_matrix(win_range)
                                     rho_mat <- get_rho_mat(dist_matrix, distance_parameter,
                                                            s)
                                     vals <- Matrix::as.matrix(monocle3::exprs(win_range))
                                     cov_mat <- cov(t(vals))
                                     Matrix::diag(cov_mat) <- diag(cov_mat) + 1e-04
                                     GL <- glasso::glasso(cov_mat, rho_mat)
                                     colnames(GL$w) <- row.names(GL$w) <- row.names(vals)
                                     colnames(GL$wi) <- row.names(GL$wi) <- row.names(vals)
                                     return(GL)
                                   })
  names_df <- as.data.frame(grs)
  names(outlist) <- paste(names_df$seqnames, names_df$start,
                          names_df$end, sep = "_")
  
  
  outlist
}



assemble_connections <- function (cicero_model_list, silent = FALSE) {
  
  require(monocle3)
  require(SummarizedExperiment)
  require(SingleCellExperiment)
  
  types <- vapply(cicero_model_list, FUN = class, FUN.VALUE = "character")
  char_hbn <- cicero_model_list[types == "character"]
  gl_only <- cicero_model_list[types == "list"]
  if (!silent) {
    print(paste("Successful cicero models: ", length(gl_only)))
    print("Other models: ")
    print(table(unlist(char_hbn)))
    print(paste("Models with errors: ", sum(is.null(cicero_model_list))))
  }
  cors <- lapply(gl_only, function(gl) {
    cors <- stats::cov2cor(gl$w)
    data.table::melt(as.data.table(cors, keep.rownames = TRUE),
                     measure = patterns("[0-9]"))
  })
  cors <- data.table::rbindlist(cors)
  names(cors) <- c("Var1", "Var2", "value")
  data.table::setkey(cors, "Var1", "Var2")
  cors_rec <- as.data.frame(cors[, list(mean_coaccess = reconcile(value)),
                                 by = "Var1,Var2"])
  names(cors_rec) <- c("Peak1", "Peak2", "coaccess")
  cors_rec <- cors_rec[cors_rec$Peak1 != cors_rec$Peak2, ]
  
  
  cors_rec
}



#' 
#' @keywords internal
#' 
make_atac_cds <- function(input, binarize = FALSE ) {
  
  if (is(input, "character")) {
    assertthat::is.readable(input)
    intersect_lean <- as.data.frame(data.table::fread(input,
                                                      header = FALSE))
  }
  else if (class(input) %in% c("matrix", "data.frame")) {
    intersect_lean <- input
  }
  else {
    stop("Input must be file path, matrix, or data.frame")
  }
  assertthat::assert_that(assertthat::are_equal(ncol(intersect_lean),
                                                3))
  assertthat::assert_that(is.logical(binarize))
  names(intersect_lean) <- c("site_name", "cell_name", "read_count")
  assertthat::assert_that(is.numeric(intersect_lean$read_count))
  intersect_lean$site_name <- as.factor(intersect_lean$site_name)
  intersect_lean$cell_name <- as.factor(intersect_lean$cell_name)
  cellinfo <- data.frame(cells = levels(intersect_lean$cell_name))
  row.names(cellinfo) <- cellinfo$cells
  cellinfo$temp <- seq_len(nrow(cellinfo))
  dhsinfo <- data.frame(site_name = levels(intersect_lean$site_name))
  dhsinfo <- cbind(dhsinfo, split_peak_names(dhsinfo$site_name))
  row.names(dhsinfo) <- dhsinfo$site_name
  names(dhsinfo) <- c("site_name", "chr", "bp1", "bp2")
  dhsinfo$chr <- gsub("chr", "", dhsinfo$chr)
  dhsinfo <- dhsinfo[order(as.character(dhsinfo$chr),
                           as.numeric(as.character(dhsinfo$bp2))),
  ]
  intersect_lean_ord <- intersect_lean[order(intersect_lean$site_name,
                                             intersect_lean$cell_name), ]
  dhsinfo <- dhsinfo[order(dhsinfo$site_name), ]
  cellinfo <- cellinfo[order(cellinfo$cells), ]
  intersect_lean_ord$site_name <- factor(intersect_lean_ord$site_name)
  intersect_lean_ord$cell_name <- factor(intersect_lean_ord$cell_name)
  intersect_lean_ord$site_name_num <- as.numeric(intersect_lean_ord$site_name)
  intersect_lean_ord$cell_name_num <- as.numeric(intersect_lean_ord$cell_name)
  if (binarize)
    intersect_lean_ord$read_count <- as.numeric(intersect_lean_ord$read_count >
                                                  0)
  sparse_intersect <- Matrix::sparseMatrix(i = intersect_lean_ord$site_name_num,
                                           j = intersect_lean_ord$cell_name_num,
                                           x = intersect_lean_ord$read_count)
  invisible(require(SingleCellExperiment))
  atac_cds <- suppressWarnings(monocle3::new_cell_data_set(methods::as(sparse_intersect,
                                                                       "sparseMatrix"),
                                                           cell_metadata = cellinfo,
                                                           gene_metadata = dhsinfo))
  monocle3::pData(atac_cds)$temp <- NULL
  monocle3::fData(atac_cds)$chr <- as.character(monocle3::fData(atac_cds)$chr)
  monocle3::fData(atac_cds)$bp1 <- as.numeric(as.character(monocle3::fData(atac_cds)$bp1))
  monocle3::fData(atac_cds)$bp2 <- as.numeric(as.character(monocle3::fData(atac_cds)$bp2))
  atac_cds <- atac_cds[order(monocle3::fData(atac_cds)$chr, monocle3::fData(atac_cds)$bp1),
  ]
  atac_cds <- monocle3::detect_genes(atac_cds)
  
  
  atac_cds
}



run_cicero <- function(cds,
                       genomic_coords,
                       window = 500000,
                       silent=FALSE,
                       sample_num = 100 ) {
  # Check input
  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(is.logical(silent))
  assertthat::assert_that(assertthat::is.number(window))
  assertthat::assert_that(assertthat::is.count(sample_num))
  if (!is.data.frame(genomic_coords)) {
    assertthat::is.readable(genomic_coords)
  }
  
  if (!silent) print("Starting Cicero")
  if (!silent) print("Calculating distance_parameter value")
  # require(monocle3)
  # require(SummarizedExperiment)
  # require(SingleCellExperiment)
  distance_parameters <- estimate_distance_parameter(cds, window=window,
                                                     maxit=100, sample_num = sample_num,
                                                     distance_constraint = 250000,
                                                     distance_parameter_convergence = 1e-22,
                                                     genomic_coords = genomic_coords)
  
  mean_distance_parameter <- mean(unlist(distance_parameters))
  
  if (!silent) print("Running models")
  cicero_out <-
    generate_cicero_models(cds,
                           distance_parameter = mean_distance_parameter,
                           window = window,
                           genomic_coords = genomic_coords)
  
  if (!silent) print("Assembling connections")
  all_cons <- assemble_connections(cicero_out, silent=silent)
  
  if (!silent) print("Done")
  all_cons
}



df_for_coords <- function (coord_strings) {
  
  coord_strings <- gsub(",", "", coord_strings)
  coord_cols <- as.data.frame(split_peak_names(coord_strings), 
                              stringsAsFactors = FALSE)
  names(coord_cols) <- c("chr", "bp1", "bp2")
  coord_cols$Peak <- coord_strings
  coord_cols$bp1 <- as.numeric(coord_cols$bp1)
  coord_cols$bp2 <- as.numeric(coord_cols$bp2)
  coord_cols
}



split_peak_names <- function(inp) {
  
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}



make_cicero_cds <- function(cds, reduced_coordinates, k = 50, summary_stats = NULL,
                            size_factor_normalize = TRUE, silent = FALSE,
                            return_agg_info = FALSE) {
  
  assertthat::assert_that(is(cds, "cell_data_set"))
  assertthat::assert_that(is.data.frame(reduced_coordinates) |
                            is.matrix(reduced_coordinates))
  assertthat::assert_that(assertthat::are_equal(nrow(reduced_coordinates),
                                                nrow(monocle3::pData(cds))))
  assertthat::assert_that(setequal(row.names(reduced_coordinates),
                                   colnames(cds)))
  assertthat::assert_that(assertthat::is.count(k) & k > 1)
  assertthat::assert_that(is.character(summary_stats) | is.null(summary_stats))
  if (!is.null(summary_stats)) {
    assertthat::assert_that(all(summary_stats %in% names(monocle3::pData(cds))),
                            msg = paste("One of your summary_stats is missing",
                                        "from your pData table. Either add a", "column with the name in",
                                        "summary_stats, or remove the name", "from the summary_stats parameter.",
                                        collapse = " "))
    assertthat::assert_that(sum(vapply(summary_stats, function(x) {
      !(is(monocle3::pData(cds)[, x], "numeric") | is(monocle3::pData(cds)[,
                                                                           x], "integer"))
    }, 1)) == 0, msg = paste("All columns in summary_stats must be",
                             "of class numeric or integer.", collapse = " "))
  }
  assertthat::assert_that(is.logical(size_factor_normalize))
  assertthat::assert_that(is.logical(silent))
  assertthat::assert_that(is.logical(return_agg_info))
  reduced_coordinates <- as.data.frame(reduced_coordinates)
  reduced_coordinates <- reduced_coordinates[colnames(cds),
  ]
  nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
  row.names(nn_map) <- row.names(reduced_coordinates)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  good_choices <- seq_len(nrow(nn_map))
  choice <- sample(seq_len(length(good_choices)), size = 1,
                   replace = FALSE)
  chosen <- good_choices[choice]
  good_choices <- good_choices[good_choices != good_choices[choice]]
  it <- 0
  k2 <- k * 2
  get_shared <- function(other, this_choice) {
    k2 - length(union(cell_sample[other, ], this_choice))
  }
  while (length(good_choices) > 0 & it < 5000) {
    it <- it + 1
    choice <- sample(seq_len(length(good_choices)), size = 1,
                     replace = FALSE)
    new_chosen <- c(chosen, good_choices[choice])
    good_choices <- good_choices[good_choices != good_choices[choice]]
    cell_sample <- nn_map[new_chosen, ]
    others <- seq_len(nrow(cell_sample) - 1)
    this_choice <- cell_sample[nrow(cell_sample), ]
    shared <- sapply(others, get_shared, this_choice = this_choice)
    if (max(shared) < 0.9 * k) {
      chosen <- new_chosen
    }
  }
  cell_sample <- nn_map[chosen, ]
  if (!silent) {
    combs <- utils::combn(nrow(cell_sample), 2)
    shared <- apply(combs, 2, function(x) {
      k2 - length(unique(as.vector(cell_sample[x, ])))
    })
    message(paste0("Overlap QC metrics:\nCells per bin: ",
                   k, "\nMaximum shared cells bin-bin: ", max(shared),
                   "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ",
                   median(shared)))
    if (mean(shared)/k > 0.1)
      warning("On average, more than 10% of cells are shared between paired bins.")
  }
  exprs_old <- monocle3::exprs(cds)
  mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in%
                   cell_sample[x, , drop = FALSE])
  if (return_agg_info) {
    row.names(mask) <- colnames(exprs_old)
    colnames(mask) <- paste0("agg_", chosen)
    agg_map <- reshape2::melt(mask)
    agg_map <- agg_map[agg_map$value, ]
    agg_map$value <- NULL
    names(agg_map) <- c("cell", "agg_cell")
  }
  mask <- Matrix::Matrix(mask)
  new_exprs <- exprs_old %*% mask
  pdata <- monocle3::pData(cds)
  new_pcols <- "agg_cell"
  if (!is.null(summary_stats)) {
    new_pcols <- c(new_pcols, paste0("mean_", summary_stats))
  }
  new_pdata <- plyr::adply(cell_sample, 1, function(x) {
    sub <- pdata[x, ]
    df_l <- list()
    df_l["temp"] <- 1
    for (att in summary_stats) {
      df_l[paste0("mean_", att)] <- mean(sub[, att])
    }
    data.frame(df_l)
  })
  new_pdata$agg_cell <- paste("agg_", chosen, sep = "")
  new_pdata <- new_pdata[, new_pcols, drop = FALSE]
  row.names(new_pdata) <- new_pdata$agg_cell
  colnames(new_exprs) <- new_pdata$agg_cell
  fdf <- monocle3::fData(cds)
  new_pdata$temp <- NULL
  cicero_cds <- suppressWarnings(monocle3::new_cell_data_set(new_exprs,
                                                             cell_metadata = new_pdata, gene_metadata = fdf))
  cicero_cds <- monocle3::detect_genes(cicero_cds, min_expr = 0.1)
  cicero_cds <- monocle3::estimate_size_factors(cicero_cds)
  if (any(!c("chr", "bp1", "bp2") %in% names(monocle3::fData(cicero_cds)))) {
    monocle3::fData(cicero_cds)$chr <- NULL
    monocle3::fData(cicero_cds)$bp1 <- NULL
    monocle3::fData(cicero_cds)$bp2 <- NULL
    monocle3::fData(cicero_cds) <- cbind(monocle3::fData(cicero_cds),
                                         df_for_coords(row.names(monocle3::fData(cicero_cds))))
  }
  if (size_factor_normalize) {
    cicero_cds <- suppressWarnings(monocle3::new_cell_data_set(Matrix::t(Matrix::t(monocle3::exprs(cicero_cds)) /
                                                                           monocle3::pData(cicero_cds)$Size_Factor),
                                                               cell_metadata = monocle3::pData(cicero_cds),
                                                               gene_metadata = monocle3::fData(cicero_cds)))
  }
  if (return_agg_info) {
    return(list(cicero_cds, agg_map))
  }
  
  
  cicero_cds
}



reconcile <- function(values) {
  if (length(values) == 1) return(values)
  if (sum(values >= 0) == length(values)) return(mean(values))
  if (sum(values <= 0) == length(values)) return(mean(values))
  if (sum(values == 0) == length(values)) return(0)
  return(NA_real_)
}



#' @keywords internal
#'
#' @importFrom dplyr %>%
#' @importFrom monocle3 estimate_size_factors preprocess_cds reduce_dimension detect_genes
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#'
build_graph <- function(obj.list, obj = NULL, rna.dis, atac.dis,
                        distance = 250000, cicero.covar = -Inf,
                        org.gs = BSgenome.Hsapiens.UCSC.hg38,
                        signac.score = -Inf, signac.pval = Inf,
                        min.cells = 10, ifWeighted = F,
                        peak.assay = 'ATAC') {

  # Link enhancers to genes
  Seurat::DefaultAssay(obj) <- peak.assay
  obj <- Signac::LinkPeaks(object = obj, peak.assay = peak.assay, expression.assay = "RNA",
                   expression.slot = "data", distance = distance,
                   pvalue_cutoff = signac.pval, score_cutoff = signac.score)
  signac.links <- stats::setNames(data.frame(Signac::Links(obj)$peak,
                             Signac::Links(obj)$gene), c("node1", "node2"))
  message ("Built ", nrow(signac.links), " enhancer-gene relations from the total dataset.")


  # Link enhancers to enhancers
  x <- Seurat::GetAssayData(object = obj, slot = "data", assay = peak.assay)
  summ <- Matrix::summary(x)
  # require(monocle3)
  # require(SummarizedExperiment)
  # require(SingleCellExperiment)
  # convert the matrix into a sparse matrix

  cicero.data <- data.frame(Origin = rownames(x)[summ$i],
                            Destination = colnames(x)[summ$j],
                            Weight      = summ$x) # transform the sparse matrix into a data frame
  input.cds <- make_atac_cds(cicero.data, binarize = TRUE ) %>% detect_genes
  # input.cds <- monocle3::detect_genes(input.cds)
  input.cds <- input.cds[Matrix::rowSums(monocle3::exprs(input.cds)) != 0, ] %>% estimate_size_factors %>%
    preprocess_cds(method = "LSI", verbose = FALSE) %>%
    reduce_dimension(reduction_method = 'UMAP', preprocess_method = "LSI")
  umap.coords <- SingleCellExperiment::reducedDims(input.cds)$UMAP # obtain the UMAP coordinates
  cicero.cds <- make_cicero_cds(input.cds, reduced_coordinates = umap.coords )
  genome.info <- data.frame(org.gs@seqinfo@seqnames,
                            org.gs@seqinfo@seqlengths) # genome sequence lengths
  colnames(genome.info) <- c("seqnames", "seqlengths") # rename the columns
  cicero.links <- run_cicero(cds = cicero.cds, genomic_coords = genome.info,
                             window = distance * 2 ) # build peak-peak linkages using cicero
  colnames(cicero.links) <- c('node1', 'node2', 'weight')
  cicero.links$node2 <- as.character(cicero.links$node2) # convert factors into characters
  coaccess.links <- data.table::rbindlist(apply(cicero.links, 1, function(r) {
    ifelse(r[1] >= r[2], return(list(node1 = r[1], node2 = r[2], weight = r[3])),
           return(list(node1 = r[2], node2 = r[1], weight = r[3])))
  }), fill = T) %>% dplyr::distinct()
  coaccess.links <- coaccess.links[, c("node1", "node2")]


  # Make graphs
  return( pbmcapply::pbmclapply(obj.list, function(x) {
    x.signac <- signac.links[signac.links$node1 %in% x$peaks,]
    if (nrow(x.signac) < 1) {
      return(NULL)
    }
    if (ifWeighted) {
      x.signac <- cbind(x.signac,
                        weight = sapply(1:nrow(x.signac), function(i) {
                          xx <- x.signac[i,]
                          sum(atac.dis[xx$node1, x$cells] > 0 &
                            rna.dis[xx$node2, x$cells] > 0)
                        })
                          ) %>% dplyr::filter(weight > 0)
    }
    x.cicero <- coaccess.links[coaccess.links$node1 %in% x$peaks &
                                 coaccess.links$node2 %in% x$peaks,]
    if (ifWeighted) {
      x.cicero <- cbind(x.cicero,
                        weight = sapply(1:nrow(x.cicero), function(i) {
                          xx <- x.cicero[i,]
                          sum(atac.dis[xx$node1, x$cells] > 0 &
                                atac.dis[xx$node2, x$cells] > 0)
                        })
                          ) %>% dplyr::filter(weight > 0)
    }


    ig <- igraph::graph_from_data_frame(rbind(x.signac, x.cicero), directed = F)
    if (ifWeighted) {
      igraph::E(ig)$weight <- length(x$cells) - c(x.signac$weight, x.cicero$weight)
    }

    ig
  }, mc.cores = max(1, parallel::detectCores() / 2)) )
}



# build_graph <- function(obj.list, obj = NULL,
#                         distance = 500000, cicero.covar = 0,
#                         org.gs = BSgenome.Hsapiens.UCSC.hg38,
#                         signac.score = 0, signac.pval = 0.5,
#                         min.cells = 10,
#                         peak.assay = 'ATAC') {
#
#   # # Libraries
#   # library(Signac)
#   # library(cicero)
#   # # library(org.gs)
#   # library(pbapply)
#   # library(igraph)
#
#
#   `%!in%` <- Negate(`%in%`) # define the negation of %in%
#   return(pbmclapply(obj.list, function(i) {
#     x <- subset(x = obj, features = c(rownames(rna.dis), i$peaks),
#                 cells = i$cells)
#     signac.links <- link_signac(x, distance = distance,
#                                 signac.score = signac.score,
#                                 signac.pval = signac.pval,
#                                 min.cells = min.cells,
#                                 peak.assay = peak.assay)
#     if (is.null(signac.links) | nrow(signac.links) < 1) {
#       return(NULL)
#     }
#     signac.links
#
#
#     # # Link a pair of enhancers
#     # x.atac <- GetAssayData(object = x, slot = 'data',
#     #                        assay = peak.assay) # extract the ATAC assay
#     # ifelse (nrow(x.atac) <= 1000 | ncol(x.atac) <= 1000,
#     #         cicero.links <- link_cor(x.atac, distance = distance * 2,
#     #                                  cicero.covar = cicero.covar),
#     #         cicero.links <- tryCatch(link_cicero(x.atac,
#     #                                              distance = distance * 2, org.gs = org.gs,
#     #                                              cicero.covar = cicero.covar),
#     #                                  error = function(e) {
#     #                                    message ('Error: something is wrong with Cicero.\n')
#     #                                    0
#     #                                  }))
#     # if ((is.null(cicero.links) | "data.frame" %!in%
#     #     class(cicero.links)) &
#     #     nrow(x.atac) > 1000 &
#     #     ncol(x.atac) > 1000) { # error exists
#     #   cicero.links <- link_cor(x.atac, distance = distance * 2,
#     #                            cicero.covar = cicero.covar)
#     # }
#     #
#     #
#     # quiet(ifelse ("data.frame" %in% class(cicero.links),
#     #               links <- rbind(signac.links, cicero.links),
#     #               links <- signac.links)) # bind two data frames
#     # G <- graph_from_data_frame(links, directed = F) # build graph
#     # E(G)$weight <- links$weight # assign the edge weights
#     #
#     #
#     # G
#   }, mc.cores = detectCores()))
# }



#' Construct Steiner forest based on minimum spanning trees
#'
#' @importFrom dplyr %>%
#'
#' @keywords internal
#'
MST_to_SFP <- function(G, ter) {

  comp.info <- igraph::components(graph = G, mode = "weak") # calculate the connected components
  splited.ter <- as.factor(comp.info$membership[ter]) %>% split(., f = factor(.)) # split terminals into sets
  Reduce(rbind, lapply(seq_along(splited.ter), function(i) {
    if (length(splited.ter[[i]]) < 2) { # only one gene
      v.edges <- igraph::incident(graph = G, v = names(splited.ter[[i]]), mode = "all") # get all neighbors
      return( as.data.frame(igraph::ends(graph = G, es = v.edges[which.min(v.edges$weight)])) )
    } else if (length(splited.ter[[i]]) == 2) { # multiple genes
      return(as.data.frame( igraph::ends(graph = G,
                                es = igraph::shortest_paths(graph = G, from = names(splited.ter[[i]])[[1]],
                                                    to = names(splited.ter[[i]])[[2]], output = "epath",
                                                    weights = igraph::E(G)$weight) %>% .$epath %>%
                                  `[[` (1)) )) # the shortest path
    } else {
      path <- igraph::shortest_paths(graph = G, from = names(splited.ter[[i]])[[1]],
                             to = names(splited.ter[[i]])[[2]], output = "both",
                             weights = igraph::E(G)$weight)
      path.vs <- igraph::as_ids(path$vpath[[1]]) # nodes
      for (j in 3 : length(splited.ter[[i]])) {
        path <- igraph::shortest_paths(graph = G, from = path.vs,
                               to = names(splited.ter[[i]])[[j]], output = "both",
                               weights = igraph::E(G)$weight)
        path.vs <- union(path.vs, igraph::as_ids(path$vpath[[1]])) # add the new nodes
      }
      sub.G <- igraph::induced_subgraph(graph = G, vids = path.vs) # build subgraph
      sub.mst <- igraph::mst(graph = sub.G,
                     weights = igraph::E(sub.G)$weight,
                     algorithm = "prim") # get the minimum spanning tree


      return( as.data.frame(igraph::ends(graph = sub.mst, es = igraph::E(sub.mst))) )
    }
  }))
}



#' Remove edges connecting two Steiner nodes, i.e., enhancers
#'
#' @keywords internal
#'
remove_steiner <- function(df) {
  Reduce(rbind, lapply(1:nrow(df), function(i) {
    if (grepl("^chr", df[i, 2])) {
      return(NULL)
    }

    return(df[i,])
  }))
}



#' Score hybrid biclusters (HBCs)
#'
#' @import dplyr
#' @import igraph
#'
#' @keywords internal
#'
score_HBC <- function(HBC, m = NULL, KL = "min.exp", Q = NULL,
                      P = NULL, G = NULL, VERY_SMALL = 0.00001) {

  if (KL == "KL") {
    # Score HBCs using KL divergence
    R <- Matrix::t(apply(m[HBC$genes, HBC$cells], 1, function(r) {
      n <- sum(r > 0)
      return(c(n, ncol(m[HBC$genes, HBC$cells]) - n))
    })) %>% `+` (VERY_SMALL) %>% `/` (ncol(m[HBC$genes, HBC$cells]) +
                                        2 * VERY_SMALL)
    C <- Matrix::t(apply(m[HBC$genes, HBC$cells], 2, function(cc) {
      n <- sum(cc > 0)
      return(c(n, nrow(m[HBC$genes, HBC$cells]) - n))
    })) %>% `+` (VERY_SMALL) %>% `/` (nrow(m[HBC$genes, HBC$cells]) +
                                        2 * VERY_SMALL)
    return(mean(sapply(rownames(R), function(j) {
      return(R[j, 1] * log(R[j, 1] / Q[j, 1]) + R[j, 2] * log(R[j, 2] / Q[j, 2]))
    })) +
      mean(sapply(rownames(C), function(j) {
        return(C[j, 1] * log(C[j, 1] / P[j, 1]) + C[j, 2] * log(C[j, 2] / P[j, 2]))
      })))
  } else if (KL == "sum.weight") {
    # Score HBCs using the total weights of enhancer-gene relations
    return( length(E(G)[HBC$genes %--% HBC$peaks]) * length(HBC$cells) -
      sum(E(G)[HBC$genes %--% HBC$peaks]$weight) )
  } else if (KL == "mean.weight") {
    # Score HBCs using the mean weights of enhancer-gene relations
    return(length(HBC$cells) -
             mean(E(G)[HBC$genes %--% HBC$peaks]$weight))
  } else if (KL == "density.weight") {
    # Score HBCs using the mean weights of enhancer-gene relations normalized
    # by the numbers of genes and enhancers
    return( (length(E(G)[HBC$genes %--% HBC$peaks]) * length(HBC$cells) -
             sum(E(G)[HBC$genes %--% HBC$peaks]$weight)) /
              length(HBC$genes) / length(HBC$peaks) )
  } else if (KL == "square") {
    return( length(HBC$genes) * length(HBC$cells) )
  } else {
    return ( min(length(HBC$genes), length(HBC$cells)) )
  }
}



#' Rearrange the columns of a data frame
#'
#' @keywords internal
#'
rearrange_cols <- function(df) {

  Reduce(rbind, lapply(1 : nrow(df), function(i) {
    if (!grepl("^chr", df[i, 1])) {
      return(list(df[i, 2:1]))
    }

    return(df[i, ])
  }))
}



#' Seeding based on Steiner forest problem (SFP) model
#'
#' @importFrom dplyr %>%
#'
#' @keywords internal
#'
SFP_seeding <- function(G.list, obj.list, bound.TFs, binding.CREs, block.list,
                        TFGene.pairs, rna.dis, atac.dis, KL = "min.exp", P = NULL,
                        Q = NULL, score.cutoff = 1, TOP_TFS = Inf, ifWeighted = FALSE,
                        quantile.cutoff = 4) {

  seed.es <- Reduce(rbind, pbmcapply::pbmclapply(seq_along(obj.list), function(i) {
    G <- G.list[[i]] # the graph
    G.vs <- igraph::as_ids(V(G)) # all nodes
    peak.vs <- G.vs[grep("^chr", G.vs)] # peak nodes
    terminal <- intersect(block.list[[i]], G.vs) # gene nodes
    if (length(terminal) == 0) {
      return(NULL)
    }
    G <- igraph::induced_subgraph(graph = G, vids = c(peak.vs, terminal))
    es.df <- MST_to_SFP(G, terminal) %>% rearrange_cols # solve the SFP model using MST
    x.df <- remove_steiner(es.df) # remove the edges whose two nodes are peaks
    colnames(x.df) <- c('steiner_node', 'terminal_node') # name the columns
    x.df <- cbind(terminal = rep(i, nrow(x.df)), x.df)


    return(x.df)
  }, mc.cores = max(1, parallel::detectCores()) / 2) )


  seeds <- Reduce(c, pbmcapply::pbmclapply(unique(seed.es$terminal), function(i) {
    G <- G.list[[i]]
    block.cells <- obj.list[[i]]$cells # the cells where genes are coexpressed
    df <- seed.es[seed.es$terminal == i,] # get the data frame
    gene.TFs <- TFGene.pairs[[i]]$gene.TFs # gene-TF relations of this cell cluster
    TF.genes <- TFGene.pairs[[i]]$TF.genes # TF-gene relations of this cell cluster
    TFs <- Reduce(c, gene.TFs[unique(df$terminal_node)]) # get the regulating TFs
    names(TFs) <- NULL
    TF.freq <- table(TFs) # the occurring frequency of TFs
    TF.freq <- TF.freq[which(TF.freq > score.cutoff)] # remove the TFs regulating no gene
    if (length(TF.freq) < 2) {
      return(NULL)
    }
    TF.freq <- TF.freq[order(TF.freq, decreasing = T)] # filter TFs
    TF.top <- names(TF.freq[1 : min(length(TF.freq), TOP_TFS)]) # select the top-ranked TFs
    ter.seeds <- lapply(TF.top, function(tt) {
      # message (tt)
      overlap.genes <- intersect(TF.genes[[tt]], unique(df$terminal_node)) # genes regulated by the TF
      if (length(overlap.genes) <= score.cutoff) {
        return(NULL)
      }
      overlap.peaks <- intersect(unique(df$steiner_node[df$terminal_node %in% overlap.genes]),
                                 binding.CREs[[tt]]) # peaks bound the the TF and meanwhile linked to
      if (length(overlap.peaks) < 1) {
        return(NULL)
      }
      overlap.genes <- df$terminal_node[df$terminal_node %in% overlap.genes &
                                          df$steiner_node %in% overlap.peaks]
      atac.ratio <- ifelse (length(overlap.peaks) > 1,
                            quantile((rowSums(atac.dis[overlap.peaks, block.cells, drop = FALSE])) /
                              length(block.cells))[quantile.cutoff],
                            # the cutoff of consistency between accessibility and expression

                            sum(atac.dis[overlap.peaks, block.cells]) / length(block.cells)
      )
      order.cells <- apply(rna.dis[overlap.genes, block.cells, drop = F], 2, sum) %>%
        sort(., decreasing = T) %>% names() # sort cells according to the sum of expression values
      HBC <- list(terminal = i, TF = tt, genes = overlap.genes, peaks = overlap.peaks,
                  cells = order.cells, atac.ratio = atac.ratio, score = 0, weight = 1)
      HBC$score <- score_HBC(HBC = HBC, KL = KL, m = rna.dis, P = P, Q = Q, G = G)
      if (ifWeighted) {
        HBC$weight <- sum(atac.dis[HBC$peaks, HBC$cells]) /
                            sum(rna.dis[HBC$genes, HBC$cells])
      }
      return(HBC)
    })
  }, mc.cores = max(1, parallel::detectCores()) / 2) )
  # }) )
  seeds <- seeds[sapply(seeds, is.not.null)]
  seeds <- seeds[sapply(seeds, function(x) {
    length(x$genes) > 1 & length(x$cells) > 1
  })]
  score.ll <- sapply(seeds, "[[", ("score")) # obtain the scores
  seeds <- seeds[order(score.ll, decreasing = TRUE )] # sort the seeds in decreasing order of scores


  return(seeds[sapply(seeds, "[[", ("score")) > score.cutoff])
}



#' Get the list of RNA and ATAC matrices
#'
#' @keywords internal
#'
get_matrix_list <- function(m, obj.list, assay = "RNA") {

  if (assay == "RNA") {
    return( parallel::mclapply(obj.list, function(x) {
      return(m[, x$cells])
    }, mc.cores = max(1, parallel::detectCores() / 2)) )
  } else {
    return( parallel::mclapply(obj.list, function(x) {
      return(m[x$peaks, x$cells])
    }, mc.cores = max(1, parallel::detectCores() / 2)) )
  }
}



#' Whether the seed is eligible based on the similarity with other
#' seeds related to the same TF
#'
#' @keywords internal
#'
intra_eligible_seed <- function(s, h, gene.cutoff = 1.0) {

  cutoff <- gene.cutoff * min(length(s$genes), length(h$genes))
  overlap.genes <- intersect(s$genes, h$genes)
  if (length(overlap.genes) > cutoff) {
    return(F)
  }


  return(T)
}



#' Get the top-ranked genes, 100 by default
#'
#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom igraph as_ids neighbors ends E %--% V
#' @importFrom Matrix rowSums
#'
#'
get_top_genes_peaks <- function(HBC, cand.genes, cand.peaks,
                                top.ngenes = 100, G,
                                rna.m, atac.m) {

  cand.genes <- setdiff(cand.genes, HBC$genes) %>%
    intersect(., rownames(rna.m)) # the genes to choose from
  if (length(cand.genes) < 1) {
    return(NULL)
  }
  cand.genes <- (rna.m[cand.genes, , drop = F] > 0) %>% rowSums %>%
    sort(., decreasing = T) %>%
    head(n = top.ngenes) %>% names # sort genes in decreasing
  cand.peaks <- setdiff(cand.peaks, HBC$peaks) %>% intersect(., as_ids(V(G)))
  # cand.peaks <- intersect(Reduce(union, parallel::mclapply(cand.genes, function(v) {
  #   as_ids(neighbors(G, v, "all"))
  # }, mc.cores = max(1, parallel::detectCores() / 2)) ), cand.peaks)
  end.mat <- ends(G, E(G)[cand.peaks %--% cand.genes])
  cand.peaks <- unique(end.mat[, 1])
  if (length(cand.peaks) < 1) {
    return(NULL)
  }
  cand.peaks <- (atac.m[cand.peaks, , drop = F] > 0) %>% rowSums %>%
    sort(., decreasing = T) %>% names
  top.peaks.genes <- lapply(1:ncol(end.mat), function(i) {
    unique(end.mat[, i])
  })
  if ("matrix" %in% class(top.peaks.genes)) {
    top.peaks.genes <- list(top.peaks.genes[, 1], top.peaks.genes[, 2])
  }


  top.peaks.genes
}



#' Expand the core part of the hybrid bicluster (HBC)
#'
#' @importFrom dplyr %>%
#' @importFrom igraph E %--% ego as_ids
#' @importFrom Matrix rowSums
#'
#' @keywords internal
#'
expand_core <- function(HBC, top.genes.peaks, rna.m, atac.m, G,
                        KL = "min.exp", Q = NULL, P = NULL,
                        DELTA = 0.0015, closure = FALSE, ego.order = 1,
                        min.cells = 10) {

  choice.genes <- top.genes.peaks[[2]] # genes to choose
  choice.peaks <- top.genes.peaks[[1]] # peaks to choose
  end.mat <- as.data.frame(ends(G, E(G)[c(HBC$genes, choice.genes) %--% c(HBC$peaks, choice.peaks)]))
  choice.gene.peaks <- split(end.mat, f = end.mat[, 2]) %>% sapply(., "[[", 1)
  iter <- 0
  score.end <- F
  while (T) {
    if (length(choice.genes) < 1) {
      break
    }
    choice.info <- parallel::mclapply(choice.genes, function(g) {
      g.prof <- rna.m[g, HBC$cells] %>% sort(., decreasing = T)
      g.cells <- g.prof[which(g.prof > 0)] %>% names
      g.lcs <- qualV::LCS(g.cells, HBC$cells) # the longest common path
      g.cells <- g.lcs$LCS
      cell.cutoff <- g.lcs$LLCS * HBC$atac.ratio # the cutoff of cell numbers
      if (length(intersect(HBC$peaks, choice.gene.peaks[[g]])) > 0) {
        return(list(gene = g, peak = NA, cells = g.cells))
      }
      g.peaks <- choice.gene.peaks[[g]]
      # g.peaks <- as.data.frame(ends(G, E(G)[g %--% choice.peaks])) %>% `[[` (1) %>% unique

      if (length(g.peaks) < 1) {
        return(NULL)
      }
      if (closure) {
        closure.peaks <- ego(G, order = ego.order, nodes = g.peaks, mode = "all") %>%
          sapply(., function(pp) {
          ppp <- as_ids(pp)
          ppp[grepl("^chr", ppp)]
        })
        atac.tmp <- atac.m[Reduce("union", closure.peaks), g.cells]
        p.cells <- lapply(closure.peaks, function(pp) {
          sum(apply(atac.m[pp, g.cells, drop = F], 2, function(ppp) {
            sum(ppp) > 0
          }))
        })
        names(p.cells) <- g.peaks
        p.cells
      } else {
        p.cells <- atac.m[g.peaks, g.cells, drop = F] %>% rowSums # calculate the qualified peaks
      }
      qual.peaks <- p.cells[which(p.cells >= cell.cutoff)] %>%
        names %>% head(n = 1) # select peaks

      if (length(qual.peaks) < 1) {
        return(NULL)
      }


      return(list(gene = g, peak = qual.peaks, cells = g.cells))
    }, mc.cores = max(1, parallel::detectCores() / 2) )


    choice.info <- choice.info[sapply(choice.info, is.not.null)]
    qual.info <- choice.info[!is.na(choice.info %>% lapply(., "[[", ("gene")) %>% unlist)]
    if (length(qual.info) < 1) {
      break
    }
    add.id <- qual.info %>% lapply(., "[[", ("cells")) %>% sapply(., length) %>% which.max
    if (length(qual.info[[add.id]]$cells) < min.cells) {
      break
    }
    add.HBC <- HBC
    add.HBC$genes <- c(add.HBC$genes, qual.info[[add.id]]$gene)
    if (!is.na(qual.info[[add.id]]$peak)) {
      add.HBC$peaks <- c(add.HBC$peaks, qual.info[[add.id]]$peak)
    }
    add.HBC$cells <- qual.info[[add.id]]$cells
    add.HBC$score <- score_HBC(add.HBC, KL = KL, P = P, Q = Q, m = m, G = G)
    iter <- iter + 1

    if (add.HBC$score < HBC$score & length(HBC$genes) > 2) {
      break
    }
    choice.genes <- setdiff(choice.genes, qual.info[[add.id]]$gene)
    HBC <- add.HBC
  }


  HBC
}



#' Add other genes or peaks that are consistent with the hybrid bicluster (HBC)
#'
#' @importFrom dplyr %>%
#' @importFrom igraph as_ids V ends %--% ego
#' @importFrom Matrix rowSums colSums
#'
#' @keywords internal
#'
expand_fuzzy <- function(HBC, G, cand.genes, cand.peaks, rna.m, atac.m,
                         g.cutoff = 1.0, m, mm, ego.order = 1, closure = F,
                         P = NULL, Q = NULL, KL = "min.exp") {

  cand1.genes <- setdiff(cand.genes, HBC$genes) %>% intersect(., as_ids(V(G))) %>%
    intersect(., rownames(m)) # remove the included genes
  cand1.peaks <- setdiff(cand.peaks, HBC$peaks) %>% intersect(., as_ids(V(G))) %>%
    intersect(., rownames(mm)) # remove the included peaks
  if (length(cand1.genes) < 1) {
    return(HBC)
  }
  gene.cutoff <- g.cutoff * length(HBC$cells)
  cand2.genes <- cand1.genes[(rna.m[cand1.genes, HBC$cells,
                                    drop = FALSE] > 0) %>%
                               rowSums >= gene.cutoff]
  if (length(cand2.genes) < 1) {
    return(HBC)
  }
  retain.genes <- cand2.genes
  new.gene.links <- ends(G, E(G)[retain.genes %--% HBC$peaks]) %>%
    apply(., 2, unique) # new genes linked to old peaks
  if ("array" %in% class(new.gene.links)) { # in cases we got a matrix
    new.gene.links <- list(new.gene.links[, 1], new.gene.links[, 2])
  }
  if (length(new.gene.links) > 1) {
    add.genes <- c(HBC$genes, new.gene.links[[1]]) # add genes directly linked to old peaks
    retain.genes <- setdiff(retain.genes, new.gene.links[[1]]) # other genes to be added
  }
  extreme.cutoff <- HBC$atac.ratio
  peak.cutoff <- length(HBC$cells) * extreme.cutoff
  if (closure) {
    closure.peaks <- ego(G, order = ego.order, nodes = cand1.peaks, mode = "all") %>%
        sapply(., function(pp) {
        ppp <- as_ids(pp)
        ppp[grepl("^chr", ppp)]
    })
    atac.tmp <- atac.m[Reduce("union", closure.peaks), HBC$cells]
    p.cells <- parallel::mclapply(closure.peaks, function(pp) {
        sum(apply(atac.m[pp, HBC$cells, drop = F], 2, function(ppp) {
        sum(ppp) > 0
        }))
    }, mc.cores = max(1, parallel::detectCores() / 2))
    names(p.cells) <- cand1.peaks
    retain.peaks <- names(p.cells)[p.cells >= peak.cutoff]
  } else {
     retain.peaks <- cand1.peaks[atac.m[cand1.peaks, HBC$cells, drop = FALSE] %>%
                                rowSums >= peak.cutoff]
  }
  cand.es <- E(G)[retain.genes %--% retain.peaks]
  if (length(cand.es) < 1) {
    return(HBC)
  }
  cand.es.df <- cbind(as.data.frame(ends(graph = G, es = cand.es)),
                      cand.es$weight) # convert to a data frame
  colnames(cand.es.df) <- c("peak", "gene", "weight")
  add.genes <- c(HBC$genes, unique(cand.es.df$gene)) # add new genes
  add.peaks <- c(HBC$peaks, unique(cand.es.df$peak[tapply(seq_along(cand.es.df$weight),
                                                          cand.es.df$gene, min)]))
  HBC$genes <- add.genes
  HBC$peaks <- add.peaks
  # HBC$cells <- add.cells
  HBC$score <- score_HBC(HBC, KL = KL, Q = Q, P = P, m = m, G = G) # score the HBC


  HBC
}



#' Expand HBC by adding cells
#'
#' @keywords internal
#'
expand_cells <- function(HBC = NULL, m = NULL, mm = NULL, quantile.cutoff = 4,
                         P = NULL, Q = NULL, G = NULL, KL = "min.exp") {

  gene.cells <- names(which(apply(m[HBC$genes, setdiff(colnames(m),
                                                       HBC$cells)], 2, sum) >=
                length(HBC$genes)))
  peak.cells <- names(which(apply(mm[HBC$peaks, gene.cells], 2, sum) >=
    quantile(apply(mm[HBC$peaks, HBC$cells], 2, sum))[quantile.cutoff]))
  HBC$cells <- c(HBC$cells, peak.cells)
  HBC$score <- score_HBC(HBC, KL = KL, Q = Q, P = P, m = m, G = G)


  HBC
}



#' Expand a hybrid bicluster (HBC)
#'
#' @keywords internal
#'
#'
#' @importFrom igraph %--% ends E
#'
expand_HBC <- function(HBC, cand.genes, cand.peaks, quantile.cutoff = 4,
                       rna.m, atac.m, dual = F, G, ego.order = 1,
                       top.ngenes = 5, c.cutoff = 1.0, closure = F,
                       KL = "min.exp", rna.dis = NULL, atac.dis = NULL,
                       Q = NULL, P = NULL, min.cells = 10) {

  top.genes.peaks <- get_top_genes_peaks(HBC = HBC, cand.peaks = cand.peaks,
                                         cand.genes = cand.genes,
                                         G = G, rna.m = rna.m, atac.m = atac.m,
                                         top.ngenes = top.ngenes)
  if (is.null(top.genes.peaks)) {
    return(HBC)
  }
  HBC <- expand_core(HBC = HBC, top.genes.peaks = top.genes.peaks,
                     rna.m = rna.m, atac.m = atac.m, closure = closure,
                     G = G, KL = KL, Q = Q, P = P, ego.order = ego.order,
                     min.cells = min.cells) # expand the core part of the HBC
  HBC <- expand_fuzzy(HBC = HBC, G = G, cand.genes = cand.genes, cand.peaks = cand.peaks,
                      m = rna.dis > 0, mm = atac.dis, rna.m = rna.m, atac.m = atac.m,
                      g.cutoff = c.cutoff, KL = KL, closure = closure,
                      ego.order = ego.order,
                      P = P, Q = Q)
  HBC <- expand_cells(HBC = HBC, m = rna.dis > 0, mm = atac.dis,
                      quantile.cutoff = quantile.cutoff)
  HBC.ends <- ends(G, E(G)[HBC$peaks %--% HBC$genes])
  HBC$genes <- unique(HBC.ends[, 2])
  HBC$peaks <- unique(HBC.ends[, 1])
  HBC.links <- Signac::StringToGRanges(HBC.ends[, 1])
  GenomicRanges::mcols(HBC.links)$gene <- HBC.ends[, 2]
  HBC$links <- HBC.links
  HBC$weight <- sum(atac.dis[HBC$peaks, HBC$cells]) / sum(rna.dis[HBC$genes, HBC$cells])


  HBC
}



#' Check whether a seed is eligible
#'
#' @keywords internal
#'
exist_seed <- function(s, HBCs, same.terminal = T,
                       intra.cutoff = 1.0, inter.cutoff = 0.80,
                       peak.cutoff = 0.80 ) {

  for (h in HBCs) {

    # Filter seeds based on terminal
    if (same.terminal & s$terminal == h$terminal) {
      return(T)
    }


    # Filter seeds based on similarity
    if (s$TF == h$TF) {
      if (!intra_eligible_seed(s = s, h = h,
                               gene.cutoff = intra.cutoff)) {
        return(F)
      }
    } else {
      if (!inter_eligible_seed(s = s, h = h,
                               gene.cutoff = inter.cutoff,
                               peak.cutoff = peak.cutoff)) {
        return(F)
      }
    }
  }


  return(T)
}



#' Check whether a seed is eligible based on the similarity with other
#' seeds related to different TFs
#'
#' @keywords internal
#'
inter_eligible_seed <- function(s, h,
                                gene.cutoff = 0.80,
                                peak.cutoff = 0.80) {

  # Check overlapped genes
  cutoff <- gene.cutoff * min(length(s$genes), length(h$genes)) # the cutoff of gene number
  overlap.genes <- intersect(s$genes, h$genes)
  if (length(overlap.genes) <= cutoff) {
    return(T)
  }
  cutoff <- peak.cutoff * min(length(s$peaks), length(h$peaks)) # the cutoff of gene number
  overlap.peaks <- intersect(s$peaks, h$peaks)
  if (length(overlap.peaks) <= cutoff) {
    return(T)
  }


  return(F)
}



#' Hybrid biclustering
#'
#' @import dplyr
#'
#' @keywords internal
#'
hybrid_biclust <- function(seeds = NULL, rna.list = NULL, atac.list = NULL,
               top.ngenes = 5, bound.TFs = NULL, same.terminal = F,
               binding.CREs = NULL, G.list = NULL, TFGene.pairs = NULL,
               c.cutoff = 1.0, KL = "min.exp", org.gs = org.gs, ego.order = 1, closure = F,
               intra.cutoff = 1.0, inter.cutoff = 0.80,
               peak.cutoff = 0.80, quantile.cutoff,
               rna.dis = NULL, atac.dis = NULL,
               min.cells = 10 ) {

  # Hybrid biclustering
  HBCs <- list() # the set of HBCs in nested list, which is the same as that of seeds
  id <- 0
  idd <- 0
  Q <- NULL
  P <- NULL
  if (KL == "KL") {

    # proportion of nonzero and zero values in each row
    Q <- Matrix::t(apply(rna.dis, 1, function(r) {
      n <- sum(r > 0)
      return(c(n, ncol(rna.dis) - n))
    })) %>% `+` (VERY_SMALL) %>% `/` (ncol(rna.dis) + 2 * VERY_SMALL)


    # proportion of nonzero and zero values in each column
    P <- Matrix::t(apply(rna.dis, 2, function(cc) {
      n <- sum(cc > 0)
      return(c(n, nrow(rna.dis) - n))
    })) %>% `+` (VERY_SMALL) %>% `/` (nrow(rna.dis) + 2 * VERY_SMALL)
  }


  # Initializes the progress bar
  message ("Hybrid biclustering on the transcriptome and chromatin accessibility matrices ...")
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(seeds), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar


  # for (i in 1:5) {
  for (i in seq_along(seeds)) {
    HBC <- seeds[[i]]
    Sys.sleep(0.1) # Remove this line and add your code
    utils::setTxtProgressBar(pb, i)


    idd <- idd + 1
    # message ("\nProcessing the putative seed : ", idd, " ...")
    if (!exist_seed(s = HBC, HBCs = HBCs, same.terminal = same.terminal,
                    intra.cutoff = intra.cutoff,
                    inter.cutoff = inter.cutoff,
                    peak.cutoff = peak.cutoff)) { # Ineligible
      next
    }
    id <- id + 1


    # message("Processing the ", id, " hybrid bicluster ...")
    HBC <- expand_HBC(HBC = HBC, cand.peaks = binding.CREs[[HBC$TF]] ,
                      rna.m = rna.list[[HBC$terminal]] , atac.m = atac.list[[HBC$terminal]] ,
                      cand.genes = TFGene.pairs[[HBC$terminal]]$TF.genes[[HBC$TF]] ,
                      G = G.list[[HBC$terminal]] , rna.dis = rna.dis , atac.dis = atac.dis ,
                      top.ngenes = top.ngenes , c.cutoff = c.cutoff,
                      ego.order = ego.order, closure = closure,
                      quantile.cutoff = quantile.cutoff,
                      KL = KL, Q = Q, P = P, min.cells = min.cells)
    if (length(HBC$genes) < 3) {
      next
    }
    HBC$seed <- idd
    HBCs[[id]] <- HBC
   }


  return(HBCs)
}



# merge_HBCs <- function(HBCs, stat = T, phyper.cutoff = 0.05,
#                        rna.dis, atac.dis) {
#
#   if (length(HBCs) < 2) {
#     return(HBCs)
#   }
#   HBCs <- HBCs[HBCs %>% sapply(., "[[", ("genes")) %>% sapply(., length) > 2]
#   HBCs <- HBCs[HBCs %>% sapply(., "[[", ("cells")) %>% sapply(., length) > 2]
#   sorted.HBCs <- HBCs[order(HBCs %>% sapply(., "[[", ("score")), decreasing = T)]
#   HBC.TFs <- unique(sapply(sorted.HBCs, "[[", "TF"))
#   new.HBCs <- do.call("c", pbmclapply(HBC.TFs, function(tf) {
#     TF.HBCs <- sorted.HBCs[sapply(sorted.HBCs, "[[", "TF") == tf]
#     if (length(TF.HBCs) < 2) {
#       return(TF.HBCs)
#     }
#     if (stat) {
#       TF.cells <- Reduce(union, unlist(sapply(TF.HBCs, "[[", "cells")))
#       N <- length(TF.cells)
#       comb.pairs <- combn(seq_along(TF.HBCs), 2)
#       adj.p.cutoff <- phyper.cutoff * length(TF.HBCs) * length(TF.HBCs)
#       tri.df <- rbindlist(lapply(1:ncol(comb.pairs), function(j) {
#         l1 <- comb.pairs[1, j]
#         l2 <- comb.pairs[2, j]
#         h1 <- TF.HBCs[[l1]]$cells
#         h2 <- TF.HBCs[[l2]]$cells
#         m <- length(h1)
#         k <- length(h2)
#         q <- length(intersect(h1, h2))
#         weight <- phyper(q - 1, m, N - m, k, lower.tail = F)
#         return(list(node1 = l1, node2 = l2, weight = weight))
#       }), fill = T) %>% dplyr::filter(weight <= adj.p.cutoff) %>%
#         dplyr::select(node1, node2) # build an igraph object
#       if (nrow(tri.df) < 1) { # no need to merge
#         return(TF.HBCs)
#       }
#       TF.cliques <- max_cliques(graph_from_data_frame(d = tri.df, directed = F)) # find maximal cliques
#       merged.TF.HBCs <- lapply(TF.cliques, function(k) {
#         HBC.list <- as.numeric(as_ids(k)) # the terminals to obtain the merged HBC
#         clique.HBCs <- TF.HBCs[HBC.list] # the HBCs belonging to this clique
#         terminal <- as.vector(sapply(clique.HBCs, "[[", "terminal"))
#         genes <- unique(unlist(lapply(clique.HBCs, "[[", "genes")))
#         peaks <- unique(unlist(lapply(clique.HBCs, "[[", "peaks")))
#         cells <- unique(unlist(lapply(clique.HBCs, "[[", "cells")))
#         atac.ratio <- mean(apply(atac.dis[peaks, cells], 1, sum)) / length(cells)
#         new.HBC <- list(terminal = terminal, TF = tf, genes = genes, peaks = peaks,
#                         cells = cells, atac.ratio = atac.ratio, score = 0)
#         new.HBC$score <- score_HBC(new.HBC, KL = 6)
#
#
#         new.HBC
#       })
#       ifelse(length(TF.cliques) < 2, covered.HBCs <- as.numeric(as_ids(TF.cliques[[1]])),
#              covered.HBCs <- unique(unlist(sapply(TF.cliques, function(y) {
#                as.numeric(as_ids(y))
#              })))) # the set of merged HBCs
#       isolated.TF.HBCs <- TF.HBCs[setdiff(seq_along(TF.HBCs), covered.HBCs)] # the isolated HBC
#       return(c(merged.TF.HBCs, isolated.TF.HBCs)) # merge the new HBCs
#     }
#   }, mc.cores = detectCores()))
#
#
#   return(new.HBCs[order(sapply(new.HBCs, "[[", "score"), decreasing = T)])
# }



# patch_HBCs <- function(merged.HBCs, binding.CREs, x, peak.ratio = NULL,
#                        peak.assay = "ATAC", distance = Inf) {
#
#   # # Libraries
#   # library(Seurat)
#   # library(pbmcapply)
#   # library(dplyr)
#   # library(Signac)
#
#
#   `%!in%` <- Negate(`%in%`) # define the negation of %in%
#
#
#   rna.m <- binarize(GetAssayData(x, slot = "data", assay = "RNA"))
#   gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = x[[peak.assay]]))
#   atac.m <- binarize(GetAssayData(x, slot = "data", assay = peak.assay))
#   patched.HBCs <- pbmclapply(seq_along(merged.HBCs), function(i) {
#     HBC <- merged.HBCs[[i]]
#     HBC.rna <- rna.m[setdiff(rownames(rna.m), HBC$genes), HBC$cells, drop = F]
#     genes.keep <- names(which(rowSums(HBC.rna) == ncol(HBC.rna))) %>%
#       intersect(gene.coords$gene_name) # genes to keep
#     if (length(genes.keep) < 1) { # No gene has coordinate annotation
#       add.HBC <- HBC
#     } else if (!is.finite(distance)) { # Without constraints of peaks
#       HBC.seqs <- unique(seqnames(gene.coords[gene.coords$gene_name %in% HBC$genes]))
#       genes.keep.coords <- gene.coords[gene.coords$gene_name %in% genes.keep] # coordinates of genes to add
#       add.HBC <- HBC
#       add.HBC$genes <- c(add.HBC$genes, genes.keep.coords[seqnames(genes.keep.coords) %in% HBC.seqs]$gene_name)
#       add.HBC$score <- score_HBC(HBC = add.HBC, KL = 6)
#     } else {
#       HBC.gene.coords <- gene.coords[which(gene.coords$gene_name %in% genes.keep)] # coordinates to keep
#       peak.HBC.distance <- DistanceToTSS(peaks = StringToGRanges(HBC$peaks), genes = HBC.gene.coords,
#                                          distance = distance) # get the distance between genes and included peaks
#       add.genes <- c(names(which(colSums(peak.HBC.distance) > 0))) # directly add the genes
#       if (length(add.genes) < 1) {
#         add.HBC <- HBC
#       } else {
#         HBC.gene.coords <- HBC.gene.coords[HBC.gene.coords$gene_name %!in% add.genes]
#         HBC.atac <- atac.m[setdiff(binding.CREs[[HBC$TF]], HBC$peaks), HBC$cells, drop = F]
#         peak.ratio <- quantile(rowSums(atac.m[HBC$peaks, HBC$cells,
#                                               drop = F]))[[3]] / ncol(HBC.atac) # 25% quantile
#         peaks.keep <- names(which(rowSums(HBC.atac) >= peak.ratio * ncol(HBC.atac))) # peaks to keep
#         if (length(peaks.keep) < 1) { # Only add the genes linked to the peaks incorporated in the HBC
#           add.HBC <- HBC
#           add.HBC$genes <- c(add.HBC$genes, add.genes)
#           add.HBC$score <- score_HBC(HBC = add.HBC, KL = 6)
#         } else { # Have qualified peaks
#           peaks.keep.GR <- StringToGRanges(peaks.keep) # get GRange objects
#           peak_distance_matrix <- DistanceToTSS(peaks = peaks.keep.GR, genes = HBC.gene.coords,
#                                                 distance = distance) # get the distance between peaks and genes
#           if (sum(peak_distance_matrix) == 0) {
#             add.HBC <- HBC
#             add.HBC$genes <- c(add.HBC$genes, add.genes)
#             add.HBC$score <- score_HBC(HBC = add.HBC, KL = 6)
#           } else {
#             colnames(peak_distance_matrix) <- HBC.gene.coords$gene_name # rename the columns
#             row.col <- get_coherent_peak_gene_pairs(peak_distance_matrix = peak_distance_matrix,
#                                                     HBC.rna = HBC.rna, HBC.atac = HBC.atac)
#             add.HBC <- HBC
#             add.HBC$genes <- c(HBC$genes, names(row.col)) # update the genes
#             add.HBC$peaks <- union(HBC$peaks, row.col) # update the peaks
#             add.HBC$score <- score_HBC(HBC = add.HBC, KL = 6)
#           }
#         }
#       }
#     }
#     add.gene.cells <- names(which(colSums(rna.m[HBC$genes,
#                                                 setdiff(colnames(rna.m),
#                                                         HBC$cells),
#                                                 drop = F]) >= length(HBC$genes)))
#     if (is.finite(distance)) {
#       peak.ratio <- quantile(colSums(atac.m[HBC$peaks, HBC$cells,
#                                             drop = F]))[[3]] / length(HBC$peaks)
#     } else {
#       peak.ratio <- quantile(colSums(atac.m[HBC$peaks, HBC$cells,
#                                             drop = F]))[[4]] / length(HBC$peaks)
#     }
#     peak.cell.cutoff <- peak.ratio * length(HBC$peaks)
#     add.peak.cells <- names(which(colSums(atac.m[HBC$peaks, add.gene.cells,
#                                                  drop = F]) >= peak.cell.cutoff))
#     add.HBC$cells <- c(add.HBC$cells, add.peak.cells)
#     add.HBC$score <- score_HBC(add.HBC, KL = 6) # score the HBC
#
#
#     return(add.HBC)
#   }, mc.cores = min(detectCores(), length(merged.HBCs)))
#
#
#   return(patched.HBCs[!sapply(patched.HBCs, is.null)])
# }



#' Calculate the pairwise similarity between hybrid biclusters (HBCs)
#'
#' @importFrom dplyr %>% filter
#'
#' @keywords internal
#'
compute_original_sim <- function(HBCs, features = "genes") {

  return(Reduce(rbind, parallel::mclapply(seq_along(HBCs), function(i) {
    sq1 <- length(HBCs[[i]][[features]]) * length(HBCs[[i]]$cells)

    sq1.sim <- lapply(seq_along(HBCs), function(j) {
      if(i < j) {
        return(NA)
      } else if (i == j) {
        return(list(HBC1 = i, HBC2 = i, Sim = 1, Same.TF = T))
      }

      sq2 <- length(HBCs[[j]][[features]]) * length(HBCs[[j]]$cells)
      sq <- length(intersect(HBCs[[i]][[features]], HBCs[[j]][[features]])) *
        length(intersect(HBCs[[i]]$cells, HBCs[[j]]$cells))
      sim <- sq / min(sq1, sq2) # original similarity

      ifelse(HBCs[[i]]$TF == HBCs[[j]]$TF,
             return(list(HBC1 = i, HBC2 = j, Sim = sim,
                         Same.TF = T)),
             return(list(HBC1 = i, HBC2 = j, Sim = sim,
                         Same.TF = F)))
    })

    return(data.table::rbindlist(sq1.sim[!is.na(sq1.sim)], fill = T))
  }, mc.cores = min(parallel::detectCores(), length(HBCs)))) )
  # %>%
  #   dplyr::filter(HBC1 != HBC2))
}



# normalize_sim <- function(sim.df, HBCs) {
#
#   same.df <- sim.df[sim.df$Same.TF, ] # pairs regulated by the same TF
#   diff.df <- sim.df[!sim.df$Same.TF, ] # pairs regulated by different TFs
#   norm.same.df <- same.df
#   norm.same.df$Sim <- (norm.same.df$Sim - mean(norm.same.df$Sim)) / sd(norm.same.df$Sim)
#   norm.diff.df <- diff.df
#   norm.diff.df$Sim <- (norm.diff.df$Sim - mean(norm.diff.df$Sim)) / sd(norm.diff.df$Sim)
#   return(rbind(norm.same.df, norm.diff.df, data.frame(HBC1 = seq_along(HBCs),
#                                                       HBC2 = seq_along(HBCs),
#                                                       Sim = rep(1, length(HBCs)),
#                                                       Same.TF = rep(T, length(HBCs)))))
#
# }



#' Calculate pairwise similarity between hybrid biclusters (HBCs)
#'
#' @keywords internal
#'
compute_sim <- function(HBCs) {

  message ("Calculating the pairwise similarity between hybrid biclusters (HBCs) ...")
  gene.sim.df <- compute_original_sim(HBCs = HBCs) # compute the original similarity overlaps on expression level
  peak.sim.df <- compute_original_sim(HBCs = HBCs, features = "peaks")
  # norm.gene.df <- normalize_sim(sim.df = gene.sim.df, HBCs = HBCs) # normalize the similarity matrix
  # norm.peak.df <- normalize_sim(sim.df = peak.sim.df, HBCs = HBCs) # normalize the similarity matrix


  norm.gene.df <- gene.sim.df
  norm.gene.df$Sim <- scale(norm.gene.df$Sim)
  norm.peak.df <- peak.sim.df
  norm.peak.df$Sim <- scale(norm.peak.df$Sim)


  norm.sim.df <- norm.gene.df
  norm.sim.df$Sim <- pmin(norm.gene.df$Sim, norm.peak.df$Sim)
  # norm.sim.df$Sim <- parallel::mclapply(1:nrow(norm.gene.df), function(i) {
  #   min(norm.gene.df$Sim[i], norm.peak.df$Sim[i])
  # }, mc.cores = parallel::detectCores() ) %>% unlist
  sim.m <- Matrix::sparseMatrix(i = norm.sim.df$HBC1, j = norm.sim.df$HBC2,
                        x = norm.sim.df$Sim) # generate a similarity matrix


  sim.m
}



#' Calculate the submodular optimization value after adding a hybrid bicluster (HBC)
#'
#' @keywords internal
#'
# compute_add <- function(query, hit, sim.m, HBC.max.sim) {
#
#   # return( unlist(sapply(query, function(i) {
#   #   add.sim <- max(sim.m[i, tail(hit, n = 1)],
#   #                  sim.m[tail(hit, n = 1), i])
#   #   ifelse (HBC.max.sim[[i]] < add.sim, return(add.sim),
#   #           return(HBC.max.sim[[i]]))
#   # })) )
#   return(apply(sim.m[query, hit, drop = FALSE], 1, max))
# }

compute_add <- function(query, hit, sim.m, HBC.max.sim) {

  return(apply(sim.m[query, hit, drop = FALSE], 1, max))
}



#' Select the next hybrid bicluster (HBC) to add
#'
#' @keywords internal
#'
# select_HBC <- function(HBCs, sim.m, HBC.flag, HBC.max.sim,
#                        regulon.ids) {
#
#   max.sim.ll <- parallel::mclapply(seq_along(HBCs), function(x) {
#     if (!HBC.flag[x]) {
#       return(-1)
#     } # do not consider the HBCs that have been already selected
#     return(compute_add(query = seq_along(HBCs), hit = c(regulon.ids, x),
#                        sim.m = sim.m,
#                        HBC.max.sim = HBC.max.sim))
#     # compute the objective value
#   }, mc.cores = max(1, parallel::detectCores() / 2) )
#   which.HBC <- which.max(unlist(sapply(max.sim.ll, sum)))
#
#
#   return(list(which.HBC, max.sim.ll[[which.HBC]]))
# }

select_HBC <- function(HBCs, sim.m, HBC.flag, HBC.max.sim,
                       step = 30,
                       regulon.ids) {

  max.sim.ll <- parallel::mclapply(seq_along(HBCs), function(x) {
    if (!HBC.flag[x]) {
      return(-1)
    } # do not consider the HBCs that have been already selected
    return( compute_add(query = seq_along(HBCs) , hit = c(regulon.ids, x),
                       sim.m = sim.m,
                       HBC.max.sim = HBC.max.sim) )
    # compute the objective value
  }, mc.cores = max(1, parallel::detectCores() / 2) )
  # which.HBC <- which.max(unlist(sapply(max.sim.ll, sum)))
  which.HBC <- head(order(unlist(sapply(max.sim.ll, sum)), decreasing = TRUE),
                    n = step)
  which.val <- compute_add(query = seq_along(HBCs) , hit = c(regulon.ids, which.HBC),
              sim.m = sim.m,
              HBC.max.sim = HBC.max.sim)


  return(list(which.HBC, which.val))
}



#' Calculate the optimization function value
#'
#' @importFrom dplyr %>%
#'
#' @keywords internal
#'
calculate_sil <- function(selected.HBCs, links.df) {

  use.genes <- Reduce(union, lapply(selected.HBCs, "[[", "genes")) # all genes included in all regulons
  use.peaks <- Reduce(union, lapply(selected.HBCs, "[[", "peaks")) # all peaks included in all regulons
  use.links <- links.df[links.df$gene %in% use.genes &
                          links.df$peak %in% use.peaks, , drop = F] # useful links
  inner.links <- sum(unlist(parallel::mclapply(1:nrow(use.links), function(i) {
    whether.include <- F
    for (j in seq_along(selected.HBCs)) {
      if (use.links[i, 1] %in% selected.HBCs[[j]]$peaks &
          use.links[i, 2] %in% selected.HBCs[[j]]$genes) {
        whether.include <- T # include
        break
      }
    }
    return(whether.include)
  }, mc.cores = max(1, parallel::detectCores() / 2))) )
  outer.links <- nrow(use.links) - inner.links # links not included
  sil.score <- inner.links - outer.links # Silhouette score
  message ("\neRegulon objective value (n = ", length(selected.HBCs), ") : ", sil.score, ".")


  sil.score # silhouette score
}



#' Perform submodular optimization
#'
#' @keywords internal
#'
# sub_mod <- function(HBCs, sim.m, G.list, n.cells, rna.list, block.list,
#                     obj, peak.assay = "ATAC", links.df,
#                     min.eRegs = 100, submod.mode = 1,
#                     distance = 500000) {
#
#   regulon.ids <- c() # the list of HBC ids
#   obj.list <- c() # the list of objective values
#   HBC.flag <- rep(T, length(HBCs))
#   HBC.max.sim <- rep(-1, length(HBCs))
#   max.n <- 0 # the number of HBCs to select which yields the maximum objective value
#   max.obj <- -2 # the current maximum objective value
#   # links.df <- filter_nearby_genes(obj, distance = distance,
#   #                                 peak.assay = peak.assay)
#
#
#   message ("Optimizing the set of enhancer regulons (eRegulons) via submodular function ...")
#   pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
#                               max = length(HBCs), # Maximum value of the progress bar
#                               style = 3,    # Progress bar style (also available style = 1 and style = 2)
#                               width = 50,   # Progress bar width. Defaults to getOption("width")
#                               char = "=")   # Character used to create the bar
#   while (1) {
#     utils::setTxtProgressBar(pb, length(regulon.ids))
#     double.output <- select_HBC(HBCs = HBCs, sim.m = sim.m,
#                                 HBC.flag = HBC.flag,
#                                 HBC.max.sim = HBC.max.sim,
#                                 regulon.ids = regulon.ids) # select the next HBC
#     add.id <- double.output[[1]]
#     HBC.max.sim <- double.output[[2]]
#     regulon.ids <- c(regulon.ids, add.id) # add the new HBC
#     HBC.flag[add.id] <- F # mask this selected HBC
#     obj.list <- c(obj.list, calculate_sil(selected.HBCs = HBCs[regulon.ids],
#                                           links.df = links.df))
#     if (length(regulon.ids) >= length(HBCs)) {
#       break
#     }
#   }
#
#   max.n <- tail(which(obj.list == max(obj.list[(min.eRegs + 1) :
#                                                  length(HBCs)])), n = 1)
#   message (max.n, " enhancer regulons (eRegulons) yield the maximum score.")
#   eRegs <- HBCs[regulon.ids[1 : max.n]] # select regulons
#
#
#   return(list(eRegs = eRegs, obj = obj.list))
# }
#'
sub_mod <- function(HBCs, sim.m, G.list, n.cells, rna.list, block.list,
                    obj, peak.assay = "ATAC", links.df,
                    min.eRegs = 100, step = 50) {

  regulon.ids <- c() # the list of HBC ids
  obj.list <- c() # the list of objective values
  HBC.flag <- rep(TRUE, length(HBCs))
  HBC.max.sim <- rep(-1, length(HBCs))
  max.n <- 0 # the number of HBCs to select which yields the maximum objective value
  max.obj <- -2 # the current maximum objective value
  # links.df <- filter_nearby_genes(obj, distance = distance,
  #                                 peak.assay = peak.assay)


  message ("Optimizing the set of enhancer regulons (eRegulons) via submodular function ...")
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                              max = length(HBCs), # Maximum value of the progress bar
                              style = 3,    # Progress bar style (also available style = 1 and style = 2)
                              width = 50,   # Progress bar width. Defaults to getOption("width")
                              char = "=")   # Character used to create the bar
  while (1) {
    utils::setTxtProgressBar(pb, length(regulon.ids))
    double.output <- select_HBC(HBCs = HBCs, sim.m = sim.m,
                                HBC.flag = HBC.flag, step = step,
                                HBC.max.sim = HBC.max.sim,
                                regulon.ids = regulon.ids) # select the next HBC
    add.id <- double.output[[1]]
    HBC.max.sim <- double.output[[2]]
    regulon.ids <- c(regulon.ids, add.id) # add the new HBC
    HBC.flag[add.id] <- FALSE # mask this selected HBC
    obj.list <- c(obj.list, calculate_sil(selected.HBCs = HBCs[regulon.ids],
                                          links.df = links.df) )
    if (length(regulon.ids) > length(HBCs)) {
      break
    }
  }


  obj.ids <- seq_along(obj.list) * step
  obj.ids[length(obj.ids)] <- min(obj.ids[length(obj.ids)], length(HBCs))
  names(obj.list) <- obj.ids
  # max.n <- tail(which(obj.list == max(obj.list[(min.eRegs + 1) :
  #                                                length(HBCs) / step ])), n = 1)
  eRegs <- HBCs[regulon.ids[1 : max(min.eRegs, as.numeric(names(which.max(obj.list))))]] # select regulons
  message (length(eRegs), " enhancer regulons (eRegulons) yield the maximum score.")


  return(list(eRegs = eRegs, obj = obj.list))
}



# expand_eGRNs <- function(obj, submod.HBCs, peak.assay = 'ATAC',
#                          distance = 250000, expand.cutoff = 0.70) {
#
#   # # Libraries
#   # library(Seurat)
#   # library(Signac)
#   # library(pbapply)
#
#
#   # Expansion
#   rna.m <- GetAssayData(object = obj, slot = "data", assay = "RNA") > 0 # get the RNA expression matrix
#   atac.m <- GetAssayData(object = obj, slot = "data", assay = peak.assay) > 0
#   gene.coords <- CollapseToLongestTranscript(ranges = Annotation(object = obj[[peak.assay]])) # coordinates
#   distance.df <- DistanceToTSS(peaks = StringToGRanges(rownames(atac.m)),
#                                genes = gene.coords,
#                                distance = distance) # distance matrix
#   expanded.eGRNs <- pblapply(submod.HBCs, function(x) {
#     add.genes <- setdiff(names(which(apply(rna.m[, x$cells], 1, sum) >=
#                                        length(x$cells) * expand.cutoff)),
#                          x$genes) # select additional genes
#     add.coords <- gene.coords[gene.coords$gene_name %in% add.genes] # subset the genes to add
#     overlap.genes <- intersect(add.genes, colnames(distance.df)) # some genes have not coordinates
#     add.dist <- summary(distance.df[, overlap.genes]) # get the distance matrix
#     x$gene.status <- c(rep(T, length(x$genes)), rep(F, length(add.genes))) # record the status of genes
#     x$genes <- c(x$genes, add.genes) # extend genes
#     x
#   })
#
#
#   expanded.eGRNs
# }



#' Convert a pair of scRNA-seq and scATAC-seq mastrices into a \code{Seurat} object
#'
#'
#' @keywords internal
#'
rna_atac_matrices_to_Seurat <- function(rna_counts = NULL, atac_counts = NULL,
                                        org = "hg38", frag.file = NULL,
                                        sep = c("-", "-")) {

  # Parameters
  message ("Loaded scRNA-seq and scATAC-seq matrices of sizes: ", nrow(rna_counts),
           " x ", ncol(rna_counts), " and ", nrow(atac_counts), " x ", ncol(atac_counts),
           ".")


  # Create Seurat object
  pbmc <- Seurat::CreateSeuratObject(counts = rna_counts)
  pbmc[["percent.mt"]] <- Seurat::PercentageFeatureSet(pbmc, pattern = "^MT-")


  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <- Signac::StringToGRanges(rownames(atac_counts), sep = sep)
  grange.use <- BSgenome::seqnames(grange.counts) %in% GenomeInfoDb::standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- Signac::GetGRangesFromEnsDb(ensdb = org_to_DB(org = org))
  ensembledb::seqlevelsStyle(annotations) <- 'UCSC'
  GenomeInfoDb::genome(annotations) <- org


  # Add fragments
  if (!is.null(frag.file)) {
    chrom_assay <- Signac::CreateChromatinAssay(
      counts = atac_counts,
      sep = sep,
      genome = org,
      fragments = frag.file,
      min.cells = 10,
      annotation = annotations
    )
  } else {
    chrom_assay <- GenomeInfoDb::CreateChromatinAssay(
      counts = atac_counts,
      sep = sep,
      genome = org,
      min.cells = 10,
      annotation = annotations
    )
  }
  pbmc[["ATAC"]] <- chrom_assay


  pbmc
}
