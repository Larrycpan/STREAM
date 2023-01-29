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



#'
#' @keywords internal
#' 
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



#'
#' @keywords internal
#' 
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


#'
#' @keyword internal
#' 
find_distance_parameter <- function(dist_mat,
                                    gene_range,
                                    maxit,
                                    null_rho,
                                    s,
                                    distance_constraint,
                                    distance_parameter_convergence) {
  
  if (sum(dist_mat > distance_constraint)/2 < 1) {
    return("No long edges")
  }
  
  found <- FALSE
  starting_max <- 2
  distance_parameter <- 2
  distance_parameter_max <- 2
  distance_parameter_min <- 0
  it <- 0
  while(found != TRUE & it < maxit) {
    vals <- as.matrix(monocle3::exprs(gene_range))
    cov_mat <- cov(t(vals))
    diag(cov_mat) <- diag(cov_mat) + 1e-4
    
    rho <- get_rho_mat(dist_mat, distance_parameter, s)
    
    GL <- glasso::glasso(cov_mat, rho)
    big_entries <- sum(dist_mat > distance_constraint)
    
    if (((sum(GL$wi[dist_mat > distance_constraint] != 0)/big_entries) > 0.05) |
        (sum(GL$wi == 0)/(nrow(GL$wi)^2) < 0.2 ) ) {
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
  if (maxit == it) warning("maximum iterations hit")
  return(distance_parameter)
}



#'
#' @keyword internal
#' 
calc_dist_matrix <- function(gene_range) {
  
  dist_mat <- Matrix::as.matrix(stats::dist(monocle3::fData(gene_range)$mean_bp))
  row.names(dist_mat) <- colnames(dist_mat) <- row.names(monocle3::fData(gene_range))
  
  
  return(dist_mat)
}


#'
#' @keyword internal
#' 
estimate_distance_parameter <- function(cds,
                                        window=500000,
                                        maxit=100,
                                        s=0.75,
                                        sample_num = 100,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        max_elements = 200,
                                        genomic_coords = NULL,
                                        max_sample_windows = 500 ) {
  
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
    it <- it + 1
    win <- sample(seq_len(length(grs)), 1)
    GL <- "Error"
    win_range <- get_genomic_range(grs, cds, win)
    
    if (nrow(monocle3::exprs(win_range))<=1) {
      next()
    }
    if (nrow(monocle3::exprs(win_range)) > max_elements) {
      next()
    }
    
    dist_matrix <- calc_dist_matrix(win_range)
    
    
    distance_parameter <- find_distance_parameter(dist_mat = dist_matrix,
                                                  gene_range = win_range,
                                                  maxit = maxit,
                                                  null_rho = 0,
                                                  s = s,
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



#'
#' @keyword internal
#' 
get_rho_mat <- function(dist_matrix, distance_parameter, s) {
  xmin <- 1000
  out <- (1-(xmin/dist_matrix)^s) * distance_parameter
  out[!is.finite(out)] <- 0
  out[out < 0] <- 0
  
  
  out
}



#'
#' @keyword internal
#' 
generate_cicero_models <- function(cds, distance_parameter, s = 0.75, window = 5e+05,
                                   max_elements = 200, genomic_coords = NULL) {

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



#'
#' @keyword internal
#' 
assemble_connections <- function (cicero_model_list, silent = FALSE) {
  
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
    data.table::melt(data.table::as.data.table(cors, keep.rownames = TRUE),
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


#'
#' @keyword
#' 
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



#' 
#' @keyword internal
#' 
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



#' 
#' @keyword internal
#' 
split_peak_names <- function(inp) {
  
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}



#' 
#' @keyword
#' 
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
    cicero_cds <- suppressWarnings(monocle3::new_cell_data_set(Matrix::t(Matrix::t(monocle3::exprs(cicero_cds)) / monocle3::pData(cicero_cds)$Size_Factor), cell_metadata = monocle3::pData(cicero_cds),
                                                               gene_metadata = monocle3::fData(cicero_cds)))
  }
  if (return_agg_info) {
    return(list(cicero_cds, agg_map))
  }
  
  
  cicero_cds
}



#' 
#' @keyword internal
#' 
reconcile <- function(values) {
  if (length(values) == 1) return(values)
  if (sum(values >= 0) == length(values)) return(mean(values))
  if (sum(values <= 0) == length(values)) return(mean(values))
  if (sum(values == 0) == length(values)) return(0)
  return(NA_real_)
}



#' Overlap two lists of \code{GRanges} objects
#'
#' @keywords internal
#' 
#' @importFrom S4Vectors queryHits subjectHits
#' 
overlap_peak_lst <- function(lst1, lst2) {
  
  overlaps <- GenomicRanges::findOverlaps(
    query = lst1,
    subject = lst2,
    type = 'any',
    select = 'all'
  ) # find the peaks overlapped with the extended genomic ranges of genes
  Matrix::sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = lst1), length(x = lst2))
  ) # build a sparse matrix to record the overlaps between peaks and extended genomic ranges of genes
}
