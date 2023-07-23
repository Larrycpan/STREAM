#' Subset a \code{Seurat} object based on \code{IRIS-FGM} blocks
#'
#' @import dplyr
#' @keywords internal
#'
subset_object <- function(CoCond_cell, object, peak.assay = 'ATAC', 
                          atac.dis, max.peaks = 3000, links.df = NULL,
                          n.blocks = 100) {
  
  atac.dis <- atac.dis[intersect(unique(links.df[links.df$gene %in%
                                                   rownames(object[["RNA"]])]$peak),
                                 rownames(atac.dis)),]

  if (length(unique(CoCond_cell$Condition)) == 0) {
    stop ("No bicluster was discovered by IRIS-FGM\n", 
          "Please perform quality control or denoising before analysis")
  } else if (length(unique(CoCond_cell$Condition)) < 2) {
    block.cells <- list(unique(CoCond_cell$cell_name))
  } else {
    block.cells <- split(CoCond_cell,
                        f = CoCond_cell$Condition) %>%
      sapply(., "[[", "cell_name") %>% head(n = n.blocks)
  }
  
  
  return( parallel::mclapply(block.cells, function(block) {
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
  }, mc.cores = parallel::detectCores()) )
}



#' 
#' @importFrom S4Vectors queryHits subjectHits
#' 
#' @keywords internal
#' 
subset_GR_intervals <- function(GR = NULL, distance = 5e+05, 
                    peaks = NULL) {
  
  return(peaks[unique(queryHits(GenomicRanges::findOverlaps(
    query = Signac::StringToGRanges(peaks), 
    subject = Signac::Extend(
      x = GR, 
      upstream = distance, 
      downstream = distance
    ),
    type = 'any',
    select = 'all'
  )))])
}



#' @keywords internal
#'
#' @importFrom dplyr %>%
#' @importFrom monocle3 estimate_size_factors preprocess_cds reduce_dimension detect_genes
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#'
build_graph <- function(obj.list, obj = NULL, rna.dis, atac.dis,
                        distance = 5e+5, cicero.covar = -Inf,
                        org.gs = BSgenome.Hsapiens.UCSC.hg38,
                        signac.score = -Inf, signac.pval = Inf,
                        min.cells = 10, ifWeighted = FALSE, 
                        filter_peaks_for_cicero = FALSE,
                        peak.assay = 'ATAC', out.dir = "./") {

  # Link enhancers to genes
  Seurat::DefaultAssay(obj) <- peak.assay
  obj <- link_peaks(object = obj, peak.assay = peak.assay, expression.assay = "RNA",
                   expression.slot = "data", distance = distance,
                   pvalue_cutoff = signac.pval, score_cutoff = signac.score)
  signac.links <- stats::setNames(data.frame(obj[[peak.assay]]@links$peak,
                                             obj[[peak.assay]]@links$gene), 
                                  c("node1", "node2"))
  message ("Built ", nrow(signac.links), " enhancer-gene relations from the total dataset.\n", 
           "Saved Signac enhancer-gene links to file: ", out.dir, "Signac_links.qsave")
  qs::qsave(signac.links, paste0(out.dir, "Signac_links.qsave"))


  # Link enhancers to enhancers
  x <- Seurat::GetAssayData(object = obj, slot = "data", assay = peak.assay)
  if (filter_peaks_for_cicero) {
    suppressWarnings(x <- x[subset_GR_intervals(GR = obj[[peak.assay]]@links, distance = distance, 
                        peaks = rownames(x)), , drop = FALSE])
    message ("Filtered enhancers before running Cicero, \n", 
             "Obtaining ", nrow(x), " enhancers.")
  }
  summ <- Matrix::summary(x)
  invisible(require(monocle3))
  cicero.data <- data.frame(Origin = rownames(x)[summ$i],
                            Destination = colnames(x)[summ$j],
                            Weight      = summ$x) # transform the sparse matrix into a data frame
  input.cds <- make_atac_cds(cicero.data, binarize = TRUE ) %>% detect_genes
  input.cds <- input.cds[Matrix::rowSums(monocle3::exprs(input.cds)) != 0, ] %>% estimate_size_factors %>%
    preprocess_cds(method = "LSI", verbose = FALSE) %>%
    reduce_dimension(reduction_method = 'UMAP', preprocess_method = "LSI")
  umap.coords <- SingleCellExperiment::reducedDims(input.cds)$UMAP # obtain the UMAP coordinates
  cicero.cds <- make_cicero_cds(cds = input.cds, reduced_coordinates = umap.coords )
  genome.info <- data.frame(org.gs@seqinfo@seqnames,
                            org.gs@seqinfo@seqlengths) # genome sequence lengths
  colnames(genome.info) <- c("seqnames", "seqlengths") # rename the columns
  cicero.links <- run_cicero(cds = cicero.cds, genomic_coords = genome.info,
                             window = distance ) # build peak-peak linkages using Cicero
  colnames(cicero.links) <- c('node1', 'node2', 'weight')
  cicero.links$node2 <- as.character(cicero.links$node2) # convert factors into characters
  coaccess.links <- data.table::rbindlist(apply(cicero.links, 1, function(r) {
    ifelse(r[1] >= r[2], return(list(node1 = r[1], node2 = r[2], weight = r[3])),
           return(list(node1 = r[2], node2 = r[1], weight = r[3])))
  }), fill = T) %>% dplyr::distinct()
  coaccess.links <- coaccess.links[, c("node1", "node2")]
  message ("Built ", nrow(coaccess.links), " co-accessible enhancer-enhancer links.\n", 
           "Saved the Cicero enhancer-enhancer links to file: ", out.dir, "Cicero_links.qsave")
  qs::qsave(coaccess.links, paste0(out.dir, "Cicero_links.qsave"))


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
  }, mc.cores = max(1, parallel::detectCores() / 4)) )
}



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
#' @importFrom dplyr %>% arrange
#'
#' @keywords internal
#'
SFP_seeding <- function(motif.obj = NULL, G.list, obj.list, bound.TFs, binding.CREs, block.list,
                        TFGene.pairs, rna.dis, atac.dis, KL = "min.exp", P = NULL,
                        Q = NULL, score.cutoff = 1, TOP_TFS = Inf, ifWeighted = TRUE,
                        quantile.cutoff = 4, peak.assay = "ATAC",
                        ifPutativeTFs = FALSE, padj.cutoff = 0.05 ) {

  # Identify enhancer-gene relations to construct seeds
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
  message ("Identified ", nrow(seed.es), " enhancer-gene relations using SFP model.")


  # Built seeds for TFs of which sites are incorporated in JASPAR2022
  seeds <- Reduce("c", pbmcapply::pbmclapply(unique(seed.es$terminal), function(i) {
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
      overlap.genes <- intersect(TF.genes[[tt]], unique(df$terminal_node)) # genes regulated by the TF
      if (length(overlap.genes) <= score.cutoff) {
        return(NULL)
      }
      overlap.peaks <- intersect(unique(df$steiner_node[df$terminal_node %in% overlap.genes]),
                                 binding.CREs[[tt]]) # peaks bound the the TF and meanwhile linked to
      if (length(overlap.peaks) < 1) {
        return(NULL)
      }
      overlap.genes <- unique(df$terminal_node[df$terminal_node %in% overlap.genes &
                                          df$steiner_node %in% overlap.peaks])
      atac.ratio <- ifelse (length(overlap.peaks) > 1,
                            quantile((rowSums(atac.dis[overlap.peaks, 
                                                       block.cells, drop = FALSE])) /
                              length(block.cells))[quantile.cutoff],
                            # the cutoff of consistency between accessibility and expression

                            sum(atac.dis[overlap.peaks, block.cells]) / length(block.cells)
      )
      order.cells <- apply(rna.dis[overlap.genes, block.cells, drop = F], 2, sum) %>%
        sort(., decreasing = T) %>% names() # sort cells according to the sum of expression values
      HBC <- list(terminal = i, Tier = 1, TF = tt, genes = overlap.genes,
                  peaks = overlap.peaks,
                  cells = order.cells, atac.ratio = atac.ratio, score = 0, weight = 1)
      HBC$score <- score_HBC(HBC = HBC, KL = KL, m = rna.dis, P = P, Q = Q, G = G)
      if (ifWeighted) {
        HBC$weight <- sum(atac.dis[HBC$peaks, HBC$cells]) /
                            sum(rna.dis[HBC$genes, HBC$cells])
      }
      return(HBC)
    })
    
    
    return(ter.seeds)
  }, mc.cores = max(1, parallel::detectCores()) / 2) )
  
  
  # Putative TFs not included in JASPAR 2022
  if (ifPutativeTFs) {
    puta.seeds <- pbmcapply::pbmclapply(
      unique(seed.es$terminal), 
      mc.cores = max(1, parallel::detectCores()) / 2, 
      function(i) {
        
        # Test enrichment
        block.cells <- obj.list[[i]]$cells # the cells where genes are coexpressed
        features <- 
          intersect(rownames(Signac::GetMotifData(object = motif.obj, 
                                                  assay = peak.assay, 
                                                  slot = "data")), 
                    seed.es[seed.es$terminal == i, "steiner_node"])
        suppressWarnings(enriched.motifs <- tryCatch(Signac::FindMotifs(
          object = motif.obj,
          features = features, 
          assay = peak.assay
        ), error = function (e) { 0 }))
        if (is.numeric(enriched.motifs)) {
          return(NULL)
        }
        tt <- enriched.motifs[enriched.motifs$p.adjust < padj.cutoff, , 
                              drop = FALSE] %>% arrange(pvalue) %>% 
          head(n = 1) %>% pull("motif.name")
        if (identical(tt, character(0))) {
          return(NULL)
        }
        puta.genes <- unique(seed.es[seed.es$terminal == i,
                                     "terminal_node"])
        HBC <- list(terminal = i, Tier = 2, TF = tt, 
                    genes = puta.genes,
                    peaks = features,
                    cells = apply(rna.dis[puta.genes, block.cells, 
                                          drop = F], 2, sum) %>%
                      sort(., decreasing = T) %>% names(), 
                    # sort cells according to the sum of 
                    # expression values
                    score = 0, weight = 1)
        HBC$atac.ratio <- 
          ifelse (length(features) > 1,
                  quantile((rowSums(atac.dis[features, 
                                             HBC$cells,
                                             drop = FALSE])) / 
                             length(HBC$cells))[quantile.cutoff],
                  sum(atac.dis[features, HBC$cells]) / 
                    length(HBC$cells)
          )
        HBC$score <- score_HBC(HBC = HBC, KL = KL, m = rna.dis, 
                               P = P, Q = Q, G = G)
        if (ifWeighted) {
          HBC$weight <- sum(atac.dis[HBC$peaks, HBC$cells]) /
            sum(rna.dis[HBC$genes, HBC$cells])
        }
        
        return(HBC)
      }
    )
    seeds <- c(seeds, puta.seeds)
  }
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
  end.mat <- ends(G, E(G)[cand.peaks %--% cand.genes])
  cand.peaks <- unique(end.mat[, 1])
  if (length(cand.peaks) < 1) {
    return(NULL)
  }
  cand.peaks <- (atac.m[cand.peaks, , drop = FALSE] > 0) %>% rowSums %>%
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
  score.end <- FALSE
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
                                                       HBC$cells), 
                                    drop = FALSE], 2, sum) >=
                length(HBC$genes)))
  peak.cells <- names(which(apply(mm[HBC$peaks, gene.cells, drop = FALSE], 2, sum) >=
    quantile(apply(mm[HBC$peaks, HBC$cells, drop = FALSE], 2, sum))[quantile.cutoff]))
  HBC$cells <- c(HBC$cells, peak.cells)
  HBC$score <- score_HBC(HBC, KL = KL, Q = Q, P = P, m = m, G = G)
  HBC
}



#' Expand a hybrid bicluster (HBC)
#'
#' @keywords internal
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
  if (is.not.null(top.genes.peaks)) {
    HBC <- expand_core(HBC = HBC, top.genes.peaks = top.genes.peaks,
                       rna.m = rna.m, atac.m = atac.m, closure = closure,
                       G = G, KL = KL, Q = Q, P = P, ego.order = ego.order,
                       min.cells = min.cells) # expand the core part of the HBC
    if (is.null(HBC)) {
      return(NULL)
    }
    HBC <- expand_fuzzy(HBC = HBC, G = G, cand.genes = cand.genes, cand.peaks = cand.peaks,
                        m = rna.dis > 0, mm = atac.dis, rna.m = rna.m, atac.m = atac.m,
                        g.cutoff = c.cutoff, KL = KL, closure = closure,
                        ego.order = ego.order,
                        P = P, Q = Q)
  }
  HBC <- expand_cells(HBC = HBC, m = rna.dis > 0, mm = atac.dis,
                      quantile.cutoff = quantile.cutoff)
  HBC.ends <- ends(G, E(G)[HBC$peaks %--% HBC$genes])
  HBC$genes <- unique(HBC.ends[, 2])
  HBC$peaks <- unique(HBC.ends[, 1])
  HBC.links <- Signac::StringToGRanges(HBC.ends[, 1])
  if (is.null(HBC) | length(HBC.links) < 1) {
    return(NULL)
  }
  GenomicRanges::mcols(HBC.links)$gene <- HBC.ends[, 2]
  HBC$links <- HBC.links
  HBC$weight <- sum(atac.dis[HBC$peaks, HBC$cells]) / sum(rna.dis[HBC$genes, HBC$cells])
  return(HBC)
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
               top.ngenes = 5, bound.TFs = NULL, same.terminal = FALSE,
               binding.CREs = NULL, G.list = NULL, TFGene.pairs = NULL,
               c.cutoff = 1.0, KL = "min.exp", org.gs = org.gs, ego.order = 1, closure = FALSE,
               intra.cutoff = 1.0, inter.cutoff = 0.80,
               peak.cutoff = 0.80, quantile.cutoff = 4,
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

  for (i in seq_along(seeds)) {
    HBC <- seeds[[i]]
    Sys.sleep(0.1) # Remove this line and add your code
    utils::setTxtProgressBar(pb, i)


    idd <- idd + 1
    if (!exist_seed(s = HBC, HBCs = HBCs, same.terminal = same.terminal,
                    intra.cutoff = intra.cutoff,
                    inter.cutoff = inter.cutoff,
                    peak.cutoff = peak.cutoff)) { # Ineligible
      next
    }


    HBC <- expand_HBC(HBC = HBC, cand.peaks = binding.CREs[[HBC$TF]] ,
                      rna.m = rna.list[[HBC$terminal]] , atac.m = atac.list[[HBC$terminal]] ,
                      cand.genes = TFGene.pairs[[HBC$terminal]]$TF.genes[[HBC$TF]] ,
                      G = G.list[[HBC$terminal]] , rna.dis = rna.dis , atac.dis = atac.dis ,
                      top.ngenes = top.ngenes , c.cutoff = c.cutoff,
                      ego.order = ego.order, closure = closure,
                      quantile.cutoff = quantile.cutoff,
                      KL = KL, Q = Q, P = P, min.cells = min.cells)
    if (is.null(HBC) | length(HBC$genes) < 3) {
      next
    }
    id <- id + 1
    HBC$seed <- idd
    HBCs[[id]] <- HBC
   }


  return(HBCs)
}



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
}



#' Calculate pairwise similarity between hybrid biclusters (HBCs)
#'
#' @keywords internal
#'
compute_sim <- function(HBCs) {

  message ("Calculating the pairwise similarity between hybrid biclusters (HBCs) ...")
  gene.sim.df <- compute_original_sim(HBCs = HBCs) # compute the original similarity overlaps on expression level
  peak.sim.df <- compute_original_sim(HBCs = HBCs, features = "peaks")
  norm.gene.df <- gene.sim.df
  norm.gene.df$Sim <- scale(norm.gene.df$Sim)
  norm.peak.df <- peak.sim.df
  norm.peak.df$Sim <- scale(norm.peak.df$Sim)


  norm.sim.df <- norm.gene.df
  norm.sim.df$Sim <- pmin(norm.gene.df$Sim, norm.peak.df$Sim)
  sim.m <- Matrix::sparseMatrix(i = norm.sim.df$HBC1, j = norm.sim.df$HBC2,
                        x = norm.sim.df$Sim) # generate a similarity matrix
  sim.m
}



#' Calculate the submodular optimization value after adding a hybrid bicluster (HBC)
#'
#' @keywords internal
#'
compute_add <- function(query, hit, sim.m, HBC.max.sim) {

  return(apply(sim.m[query, hit, drop = FALSE], 1, max))
}



#' Select the next hybrid bicluster (HBC) to add
#'
#' @keywords internal
#'
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
sub_mod <- function(HBCs, sim.m, G.list, n.cells, rna.list, block.list,
                    obj, peak.assay = "ATAC", links.df,
                    min.eRegs = 100, step = 50) {

  regulon.ids <- c() # the list of HBC ids
  obj.list <- c() # the list of objective values
  HBC.flag <- rep(TRUE, length(HBCs))
  HBC.max.sim <- rep(-1, length(HBCs))
  max.n <- 0 # the number of HBCs to select which yields the maximum objective value
  max.obj <- -2 # the current maximum objective value


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



#' Convert a pair of scRNA-seq and scATAC-seq mastrices into a \code{Seurat} object
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
  ensembldb::seqlevelsStyle(annotations) <- 'UCSC'
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
    chrom_assay <- Signac::CreateChromatinAssay(
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
