#' Identify enhancer regulons (eRegulons) from jointly profiled scRNA-seq and scATAC-seq data
#' 
#' @description Most transcriptional regulation analyses, including prediction of transcription factor 
#' (TF) binding sites and identification of enhancer-gene relations require the discovery of enhancer 
#' regulons (eRegulons), which are active in a set of cells. A TF with its set of target enhancers 
#' and target genes is called an eRegulon. This function takes as input a \code{Seurat} 
#' object (containing both scRNA-seq and scATAC-seq) and then simultaneously identifies eRegulons as well 
#' as the cell sets where eRegulons are active. In addition, this function determines the number of eRegulons 
#' in a dataset by submodular optimization.
#' 
#' @param obj A \code{Seurat} object composed of both scRNA-seq and scATAC-seq assays.
#' @param top.peaks The number of top-ranked peaks to identify the core part of hybrid biclusters (HBCs),
#' 3000 by default.
#' @param min.cells The cutoff of minimum number of cells for quality control (QC), 10 by default.
#' @param out.dir The directory to save the intermediate or final results, "./" by default.
#' @param org The organism version, hg38 by default.
#' @param top.ngenes The number of genes composing the core part of an HBC, 5 by default.
#' @param c.cutoff The cutoff of consistency during hybrid biclustering process, 1.0 by default.
#' @param n.blocks The cutoff of the maximum number of blocks output by IRIS-FGM, 100 by default.
#' @param min.eRegs The minimum number of enhancer regulons (eRegulons) to output, 100 by default.
#' @param peak.assay The scATAC-seq assay, "ATAC" by default.
#' @param distance The distance cutoff to build enhancer-enhancer relations, 5e+05 by default.
#' @param BlockOverlap The cutoff of maximum overlap between blocks output by \code{IRIS-FGM}, 0.50 by default.
#' @param Extension The consistency level to expand a block by \code{IRIS-FGM}, 0.70 by default.
#' @param intra.cutoff The cutoff to calculate pairwise similarity among HBCs associated with the same TFs,
#' 1.0 by default.
#' @param inter.cutoff The cutoff to compute pairwise similarity among genes in HBCs associated with different TFs,
#' 0.80 by default.
#' @param peak.cutoff The cutoff to quantify pairwise similarity among enhancers in HBCs associated with 
#' different TFs, 0.80 by default.
#' @param var.genes The number of highly variable genes to predict used to identify the core part of HBCs,
#' 3000 by default.
#' @param KL Which method to use for measuring the score of HBCs, "min.exp" by default, i.e.,
#' the smaller one between the numbers of genes and cells in a HBC.
#' @param quantile.cutoff The quantile cutoff of the ratio of HBC cells, where enhancers are accessible,
#' 4 by default, indicating the top-25% among ranks.
#' @param submod.step The step size of, i.e., the number of HBCs to add in each step during iteration,
#'  for submodular optimization, 30 by default.
#' @param filter_peaks_for_cicero Whether filter the peaks in a neighborhood of the \code{Signac} results 
#' before running \code{cicero}, FALSE by default.
#' @param candidate.TFs The list of candidate TFs used to identify eRegulons and eGRNs, NULL by default.
#'
#' @rdname run_stream
#' @importFrom dplyr %>%
#' @export
#' 
#' @return When running on a \code{Seurat} object,
#' returns a list of eRegulons saved in a nested list, each of which contains the following attributes:
#'    \item{terminal} {The \code{IRIS-FGM} block used to predict the eRegulon.}
#'    \item{Tier} {The tier of the TF-enhancer relations: 1 represents JASPAR annotations; 
#'    2 denotes motif scanning.}
#'    \item{TF} {The TF of the eRegulon.}
#'    \item{genes} {Genes of the eRegulon.}
#'    \item{peaks} {Enhancers of the eRegulon.}
#'    \item{cells} {Cells where the eRegulon is active.}
#'    \item{atac.ratio} {The ratio of cells where the eRegulon enhancers are 
#'    accessible against cells in which the eRegulon genes are expressed.}
#'    \item{score} {The eRegulon score.}
#'    \item{weight} {The eRegulon weight.}
#'    \item{links} {The enhancer-gene relations saved in \code{GRanges} object.}
#'    \item{seed} {The seed to obtain the eRegulon.}
#' 
#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' @references Chang, Y., Allen, C., Wan, C., Chung, D., Zhang, C., Li, Z. and Ma, Q., 2021. 
#' IRIS-FGM: an integrative single-cell RNA-Seq interpretation system for functional gene module analysis. 
#' Bioinformatics, 37(18), pp.3045-3047.
#' @references Xie, J., Ma, A., Zhang, Y., Liu, B., Cao, S., Wang, C., ... & Ma, Q. (2020). 
#' QUBIC2: a novel and robust biclustering algorithm for analyses and interpretation of 
#' large-scale RNA-Seq data. Bioinformatics, 36(4), 1143-1149.
#' @references Stuart, T., Srivastava, A., Madad, S., Lareau, C. A., & Satija, R. (2021). 
#' Single-cell chromatin state analysis with Signac. Nature methods, 18(11), 1333-1341.
#' @references Pliner, H. A., Packer, J. S., McFaline-Figueroa, J. L., Cusanovich, 
#' D. A., Daza, R. M., Aghamirzaie, D., ... & Trapnell, C. (2018). 
#' Cicero predicts cis-regulatory DNA interactions from single-cell chromatin accessibility data. 
#' Molecular cell, 71(5), 858-871.
#' @references Castro-Mondragon, J. A., Riudavets-Puig, R., Rauluseviciute, I., Berhanu Lemma, R., 
#' Turchi, L., Blanc-Mathieu, R., ... & Mathelier, A. (2022). 
#' JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles. 
#' Nucleic acids research, 50(D1), D165-D173.
#' 
run_stream <- function(obj = NULL,
                       candidate.TFs = NULL,
                       peak.assay = "ATAC",
                       var.genes = 3000,
                       top.peaks = 3000,
                       ifPutativeTFs = FALSE,
                       min.cells = 10,
                       out.dir = "./",
                       org = "hg38",
                       top.ngenes = 5,
                       c.cutoff = 1.0,
                       n.blocks = 100,
                       distance = 5e+05,
                       filter_peaks_for_cicero = FALSE,
                       ifWeighted = TRUE,
                       cicero.covar = -Inf,
                       signac.score = -Inf,
                       signac.pval = Inf,
                       intra.cutoff = 1.0,
                       inter.cutoff = 0.80,
                       peak.cutoff = 0.80,
                       score.cutoff = 1,
                       KL = "min.exp",
                       quantile.cutoff = 4,
                       BlockOverlap = 0.50,
                       Extension = 0.90,
                       submod.step = 30,
                       min.eRegs = 100, 
                       url.link = "https://figshare.com/ndownloader/files/38794185"
                       ) {

  # Set start time
  start.time <- Sys.time()
  message ("Began running STREAM2 at time: \n", start.time)
  
  
  # Check parameters
  if (!dir.exists(out.dir)) {
    message ("Creating the directory: ", out.dir, " to save the intermediate or final results ...")
    dir.create(out.dir)
  }


  # Quality control
  invisible(peaks.use <- rownames(obj[[peak.assay]])[as.vector(as.character(BSgenome::seqnames(obj[[peak.assay]]@ranges)) %in% GenomeInfoDb::standardChromosomes(obj[[peak.assay]]@ranges))])
  obj <- subset(obj,
                features = c(
                  rownames(obj[["RNA"]])[!grepl("\\.", rownames(obj[["RNA"]]))],
                  peaks.use
                ))
  message ("Filtered out features not a gene symbol or belonging to non-standard chromosomes,\n",
           "leading to ", nrow(obj[["RNA"]]), " genes and ", nrow(obj[[peak.assay]]), " enhancers.")
  qs::qsave(obj, paste0(out.dir, "Obj_filtered.qsave"))
  message ("Saved the Seurat object after filtering to file: ", paste0(out.dir, "Obj_filtered.qsave."))


  # Libraries
  set.seed(1234)
  quiet(if (org == "mm10") {
    require(BSgenome.Mmusculus.UCSC.mm10)
    org.gs <- BSgenome.Mmusculus.UCSC.mm10
  } else if (org == "mm9") {
    require(BSgenome.Mmusculus.UCSC.mm9)
    org.gs <- BSgenome.Mmusculus.UCSC.mm9
  } else if (org == "hg19") {
    require(BSgenome.Hsapiens.UCSC.hg19)
    org.gs <- BSgenome.Hsapiens.UCSC.hg19
  } else {
    require(BSgenome.Hsapiens.UCSC.hg38)
    org.gs <- BSgenome.Hsapiens.UCSC.hg38
  })
  message ("Loaded full genome sequences for ", org, ".")


  # Filter the genes that have annotations
  Seurat::DefaultAssay(obj) <- peak.assay
  obj <- Signac::RegionStats(object = obj, assay = peak.assay,
                     genome = org.gs)
  message ("Computed the GC content, region lengths, and dinucleotide base frequencies",
           " for regions in the ", peak.assay, " assay and add them to the feature metadata.")
  qs::qsave(obj, paste0(out.dir, "Obj_quality_ctr.qsave"))
  links.df <- filter_nearby_genes(obj = obj, peak.assay = peak.assay)
  obj <- subset(x = obj, features = c(unique(links.df$peak),
                                      unique(links.df$gene)))
  obj <- obj[, Matrix::colSums(obj[['RNA']]@counts) > 0 &
               Matrix::colSums(obj[[peak.assay]]@counts) > 0]
  message ("Removed the cells with zero count of RNA or ATAC, \n",
           "leading to ", ncol(obj), " cells.")
  message ("Filtered the enhancers that are beyond 500 kb from any genes in the Seurat object.\n",
           length(unique(links.df$peak)), " enhancers and ", length(unique(links.df$gene)),
           " genes were retained.")


  # Find highly variable genes and top-ranked enhancers
  Seurat::DefaultAssay(obj) <- "RNA"
  suppressWarnings(obj <- Seurat::SCTransform(obj, verbose = FALSE, variable.features.n = var.genes))
  Seurat::DefaultAssay(obj) <- peak.assay
  obj <- Signac::RunTFIDF(obj) # frequency-inverse document frequency (TF-IDF) normalization
  obj <- Signac::FindTopFeatures(obj, min.cutoff = "q5")
  obj <- subset(x = obj, features = c(obj[[peak.assay]]@var.features,
                                        obj[['SCT']]@var.features))
  message ("Saved the Seurat object composed of ", nrow(obj[['RNA']]), " genes (",
           length(obj[["SCT"]]@var.features), " variable genes) and ",
           nrow(obj[[peak.assay]]), " enhancers (", length(obj[[peak.assay]]@var.features),
           " top-ranked enhancers) to file: ", out.dir,
           "Obj_var_genes_top_enhs.qsave.")
  qs::qsave(obj, paste0(out.dir, "Obj_var_genes_top_enhs.qsave"))


  options(timeout = 2000)
  load(url(url.link))
  message ("Loaded TF binding sites from JASPAR 2022.")
  TF.CRE.pairs <- find_TFBS(peaks = rownames(obj[[peak.assay]]),
                            TFBS.list = TFBS.list, org = org, 
                            candidate.TFs = candidate.TFs)
  qs::qsave(TF.CRE.pairs, paste0(out.dir, "TF_binding_sites_on_enhs.qsave"))
  bound.TFs <- TF.CRE.pairs$CRE
  binding.CREs <- TF.CRE.pairs$TF
  rm(TF.CRE.pairs)
  message ("Identified ", length(binding.CREs), " TFs binding ",
           length(bound.TFs), " enhancers,\n",
           "which were saved to file: ", out.dir, "TF_binding_sites_on_enhs.qsave.")
  
  
  # Identify TFs tentatively binding enhancers based on the position weight matrices (PWMs) in JASPAR 2022
  if (ifPutativeTFs) {
    
    # Get a list of motif position frequency matrices from the JASPAR database
    quiet(require(JASPAR2022))
    pfm <- TFBSTools::getMatrixSet(
      x = JASPAR2022,
      opts = list(collection = "CORE", tax_group = 'vertebrates', 
                  all_versions = FALSE)
    )
    puta.TF.peak.pairs <- run_motifmatchr(pfm = pfm, 
                                          peaks = rownames(obj[[peak.assay]]), 
                                          org.gs = org.gs, 
                                          candidate.TFs = candidate.TFs)
    qs::qsave(puta.TF.peak.pairs, paste0(out.dir, "Putative_TF_binding_sites_on_enhs.qsave"))
    puta.bound.TFs <- puta.TF.peak.pairs$puta.bound.TFs
    puta.binding.peaks <- puta.TF.peak.pairs$puta.binding.peaks
    rm(puta.TF.peak.pairs)
    message ("Identified ", length(puta.binding.peaks), " TFs putatively binding ",
             length(puta.bound.TFs), " enhancers,\n",
             "which were saved to file: ", out.dir, "Putative_TF_binding_sites_on_enhs.qsave.")
  }


  # LTMG modeling
  LTMG.obj <- call_LTMG(obj = obj)
  message ("Finished LTMG modeling for ", nrow(LTMG.obj@LTMG@LTMG_discrete),
           " genes across ", ncol(LTMG.obj@LTMG@LTMG_discrete), " cells.\n",
           "There are ", length(unique(range(LTMG.obj@LTMG@LTMG_discrete))),
           " transcriptional regulatory states (TRSs) identified.")


  # Please uncomment the following command and make sure to set a correct working directory
  # so that the following command will generate intermediate files.
  LTMG.obj <- IRISFGM::CalBinaryMultiSignal(LTMG.obj) # partition each gene into TRSs
  LTMG.obj <- IRISFGM::RunBicluster(LTMG.obj, DiscretizationModel = "LTMG", OpenDual = FALSE,
                         NumBlockOutput = n.blocks, BlockOverlap = BlockOverlap,
                         Extension = Extension,
                         BlockCellMin = min.cells)
  message ("Identified ", max(LTMG.obj@BiCluster@CoReg_gene$Condition),
           " IRIS-FGM biclusters from the LTMG matrix.\n",
           "Saved the LTMG object to file: ", out.dir, "LTMG.qsave.")
  qs::qsave(LTMG.obj, paste0(out.dir, "LTMG.qsave"))


  # Get the list of Seurat objects
  rna.dis <- subset(x = LTMG.obj@LTMG@LTMG_discrete,
                    rownames(LTMG.obj@LTMG@LTMG_discrete) %in%
                      rownames(Seurat::GetAssayData(obj, assay = "RNA", slot = "data")))
  atac.dis <- binarize(Seurat::GetAssayData(object = obj, slot = 'data',
                                    assay = peak.assay))
  link.ratio <- sapply(apply(links.df, 2, unique), length)
  obj.list <- subset_object(LTMG.obj = LTMG.obj, object = obj, links.df = links.df,
                            atac.dis = atac.dis, n.blocks = n.blocks,
                            max.peaks = round(nrow(rna.dis) * link.ratio[1] / link.ratio[2]),
                            peak.assay = peak.assay)
  block.list <- split(LTMG.obj@BiCluster@CoReg_gene,
                      f = LTMG.obj@BiCluster@CoReg_gene$Condition) %>%
    sapply(., "[[", "Gene")
  flags <- sapply(obj.list, is.not.null)
  obj.list <- obj.list[flags]
  block.list <- block.list[flags]


  # Construct heterogeneous graphs
  G.list <- build_graph(obj.list = obj.list, obj = obj, rna.dis = rna.dis,
                        atac.dis = atac.dis, ifWeighted = ifWeighted,
                      distance = distance, cicero.covar = cicero.covar,
                      org.gs = org.gs, peak.assay = peak.assay,
                      signac.score = signac.score, signac.pval = signac.pval,
                      filter_peaks_for_cicero = filter_peaks_for_cicero,
                      min.cells = min.cells, out.dir = out.dir)
  flags <- sapply(G.list, is.not.null) # whether the graph is empty
  G.list <- G.list[flags] # filter the graphs
  obj.list <- obj.list[flags] # filter the Seurat objects
  block.list <- block.list[flags]
  qs::qsave(obj.list, paste0(out.dir, "Obj_list.qsave"))
  qs::qsave(G.list, paste0(out.dir, "Graph_list.qsave"))
  if (length(G.list) < 1) {
    stop ("No heterogeneous graph was constructed")
  }
  message ("List containing ", length(obj.list), " Seurat objects were saved to file: ",
           out.dir, "Obj_list.qsave.")
  message ("Constructed ", length(G.list), " heterogeneous graphs,\n",
           "which were saved to file: ", out.dir, "Graph_list.qsave.")


  # Discover cell-subpopulation-active TF-target pairs
  TFGene.pairs <- find_TF_gene(G.list = G.list, bound.TFs = bound.TFs,
                               binding.CREs = binding.CREs)
  qs::qsave(TFGene.pairs, paste0(out.dir, "TF_target_pairs.qsave"))
  message ("Saved the TF-target pairs to file: ", out.dir, "TF_target_pairs.qsave.")

  
  # Discover putative cell-subpopulation-active TF-target pairs (based on motif enrichment)
  if (ifPutativeTFs) {
    puta.TF.gene.pairs <- find_TF_gene(G.list = G.list, bound.TFs = puta.bound.TFs,
                                       binding.CREs = puta.binding.peaks)
    qs::qsave(puta.TF.gene.pairs, paste0(out.dir, "Putative_TF_target_pairs.qsave"))
    message ("Saved the putative TF-target pairs to file: ", out.dir, "Putative_TF_target_pairs.qsave.")
  }
  
  
  # Add motif annotations (if demanded)
  if (ifPutativeTFs) {
    if (grepl("^mm", org)) {
      tax.id <- 10090
    } else {
      tax.id <- 9606
    }
    quiet(require(JASPAR2022))
    pfm <- TFBSTools::getMatrixSet(
      x = JASPAR2022,
      opts = list(collection = "CORE", 
                  tax_group = 'vertebrates', 
                  species = tax.id, 
                  all_versions = FALSE)
    )
    quiet(require(motifmatchr))
    quiet( obj <- Signac::AddMotifs(
      object = obj,
      genome = org.gs,
      pfm = pfm
    ) )
    motif.obj <- subset(obj, features = 
                          intersect(rownames(obj[[peak.assay]]@motifs), 
                                                  rownames(obj[[peak.assay]])))
  } else {
    motif.obj <- NULL
  }
  

  # Seeding based upon Steiner forest problem (SFP) model
  seeds <- SFP_seeding(motif.obj = motif.obj, G.list = G.list, obj.list = obj.list, 
                       block.list = block.list,
                       bound.TFs = bound.TFs, binding.CREs = binding.CREs,
                       TFGene.pairs = TFGene.pairs, score.cutoff = score.cutoff,
                       quantile.cutoff = quantile.cutoff, ifPutativeTFs = ifPutativeTFs,
                       rna.dis = rna.dis, atac.dis = atac.dis, KL = KL, 
                       peak.assay = peak.assay)
  if (length(seeds) < 1) {
    stop ("No seeds was identified")
  }
  message (length(seeds), " seeds are identified for hybrid biclustering,\n",
           "which were saved to file: ", out.dir, "Seeds.qsave.")
  qs::qsave(seeds, paste0(out.dir, "Seeds.qsave"))


  # Get the list of RNA and ATAC matrices
  rna.list <- get_matrix_list(m = rna.dis, obj.list = obj.list, assay = "RNA") # RNA matrix
  atac.list <- get_matrix_list(m = atac.dis, obj.list = obj.list, assay = peak.assay) # ATAC matrix
  rm(obj.list)


  # Hybrid biclustering for Tier-1 seeds (TF binding sites from JASPAR 2022 database)
  message ("Biclustering based on Tier-1 TF-enhancer binding annotations ...")
  HBCs <- hybrid_biclust(seeds = seeds[(sapply(seeds, "[[", "Tier") == 1)], 
                         rna.list = rna.list, atac.list = atac.list,
                         top.ngenes = top.ngenes, bound.TFs = bound.TFs,
                         binding.CREs = binding.CREs, G.list = G.list,
                         TFGene.pairs = TFGene.pairs, peak.cutoff = peak.cutoff,
                         c.cutoff = c.cutoff, KL = KL, org.gs = org.gs,
                         intra.cutoff = intra.cutoff, inter.cutoff = inter.cutoff,
                         quantile.cutoff = quantile.cutoff,
                         rna.dis = rna.dis, atac.dis = atac.dis, min.cells = min.cells)
  
  
  # Hybrid biclustering for Tier-2 seeds (TF binding sites from motif enrichment)
  if (ifPutativeTFs) {
    message ("Biclustering based on Tier-2 motif enrichment ...")
    puta.HBCs <- hybrid_biclust(seeds = seeds[(sapply(seeds, "[[", "Tier") == 2)], 
                           rna.list = rna.list, atac.list = atac.list,
                           top.ngenes = top.ngenes, bound.TFs = puta.bound.TFs,
                           binding.CREs = puta.binding.peaks, G.list = G.list,
                           TFGene.pairs = puta.TF.gene.pairs, peak.cutoff = peak.cutoff,
                           c.cutoff = c.cutoff, KL = KL, org.gs = org.gs,
                           intra.cutoff = intra.cutoff, inter.cutoff = inter.cutoff,
                           quantile.cutoff = quantile.cutoff,
                           rna.dis = rna.dis, atac.dis = atac.dis, min.cells = min.cells)
    HBCs <- c(HBCs, puta.HBCs)
  }
  
  
  # Basic information of the identified HBCs
  gene.range <- sapply(HBCs, "[[", "genes") %>% sapply(., length) %>% range
  peak.range <- sapply(HBCs, "[[", "peaks") %>% sapply(., length) %>% range
  cell.range <- sapply(HBCs, "[[", "cells") %>% sapply(., length) %>% range
  message ("Identified ", length(HBCs), " hybrid biclusters (HBCs),\n",
           "which were saved to file: ", out.dir, "HBCs.qsave.\n",
           "HBCs contain ", gene.range[1], "-", gene.range[2], " genes, ",
           peak.range[1], "-", peak.range[2], " enhancers, ",
           cell.range[1], "-", cell.range[2], " cells, and ",
           length(unique(sapply(HBCs, "[[", "TF"))), " TFs.")
  qs::qsave(HBCs, paste0(out.dir, "HBCs.qsave"))


  # Free some memory
  rm(TFGene.pairs)
  rm(atac.list)
  rm(atac.dis)
  rm(bound.TFs)


  # Submodular optimization
  if (length(HBCs) <= min.eRegs) {
    submod.HBCs <- HBCs
  } else {
    message("Performing submodular optimization ...")
    sim.m <- compute_sim(HBCs = HBCs) # calculate the pairwise similarity between HBCs
    submod.obj <- sub_mod(HBCs = HBCs, sim.m = sim.m, rna.list = rna.list,
                          G.list = G.list, links.df = links.df,
                          block.list = block.list, n.cells = ncol(rna.dis),
                          obj = obj, min.eRegs = min.eRegs, step = submod.step,
                          peak.assay = peak.assay) # submodular optimization
    rm(sim.m)
    submod.HBCs <- submod.obj$eRegs
    message ("HBC scores after submodular optimization were written to file: ", 
             out.dir, "Submodular_scores.qsave.")
    qs::qsave(submod.obj$obj, paste0(out.dir, "Submodular_scores.qsave"))
    
    
    # Set end time
    end.time <- Sys.time()
    message ("Finished running STREAM2 at time: \n", end.time, "\n\n", 
             "Total running time: ", end.time - start.time)


    return(submod.obj$eRegs)
  }


  return(submod.HBCs)
}



#' Simulate a jointly profiled scRNA-seq and scATAC-seq dataset in which
#' several enhancer regulons (eRegulons) are contained. 
#' 
#' @importFrom dplyr %>%
#'
#' @export
#' @rdname create_rna_atac
#'
#' @param obj A \code{Seurat} object used as the prototype to generate simulated dataset.
#' @param ntfs The number of eRegulons (TFs) to include in the simulated dataset.
#' @param ngenes The average number of genes in an eRegulon.
#' @param ncells The average number of cells in an eRegulon.
#' @param org The organism, hg38 by default.
#' @param atac.assay The scATAC-seq assay.
#' @param gene.links The average number of enhancers linked to each gene in an eRegulon.
#' @param distance The maximum distance between a gene and its linked enhancers.
#' @param all.genes the number of genes in the simulated \code{Seurat} object.
#' @param all.enhs The number of enhancers in the simulated \code{Seurat} object.
#' @param all.cells The number of cells in the simulated \code{Seurat} object.
#'
#' @return Returns a list composed of an eRegulon list and a \code{Seurat} object.
#' each item of the simulated eRegulon contains the following attributes:
#'    \item{TF} {The TF of the eRegulon.}
#'    \item{genes} {Genes of the eRegulon.}
#'    \item{peaks} {Enhancers of the eRegulon.}
#'    \item{cells} {Cells where the eRegulon is active.}
#' 
#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' @references Castro-Mondragon, J. A., Riudavets-Puig, R., Rauluseviciute, I., Berhanu Lemma, R., 
#' Turchi, L., Blanc-Mathieu, R., ... & Mathelier, A. (2022). 
#' JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles. 
#' Nucleic acids research, 50(D1), D165-D173.
#' @references Garcia-Alonso, L., Holland, C. H., Ibrahim, M. M., Turei, D., & Saez-Rodriguez, J. (2019). 
#' Benchmark and integration of resources for the estimation of human transcription factor activities. 
#' Genome research, 29(8), 1363-1375.
#' 
create_rna_atac <- function(obj = NULL, ntfs = 5, ngenes = 100,
                            ncells = 100, all.genes = 1000, all.enhs = 3000, all.cells = 1000,
                            org = "hg38", atac.assay = "ATAC", gene.links = 2,
                            distance = 5e+05, url.link = "https://figshare.com/ndownloader/files/38794185"
                            ) {

  # Parameters
  message ("Loading TF binding sites from JASPAR 2022 ...")
  options(timeout = 2000)
  load(url(url.link))
  sites <- TFBS.list[[org]]
  rm(TFBS.list)
  message ("There are ", length(unique(sites$TF)), " TFs from the JASPAR 2022 database.\n",
           ntfs, " eRegulons regulated by different TFs will be generated.\n",
           "The Seurat object contains ", nrow(obj[["RNA"]]),
           " genes, ",
           nrow(obj[[atac.assay]]), " enhancers, and ", ncol(obj), " cells.")


  # Select TFs
  set.seed(123)
  cand.tfs <- unique(sites$TF)
  cand.tfs <- cand.tfs[!grepl("::", cand.tfs)]
  tf.lst <- unique(sites$TF)[sample(seq_along(cand.tfs), size = ntfs)]


  # Select target genes
  if (grepl("^mm", org)) {
    tf.target <- dorothea::dorothea_mm[, c("tf", "target")]
  } else {
    tf.target <- dorothea::dorothea_hs[, c("tf", "target")]
  }
  message ("Loaded ", length(unique(tf.target$tf)), " TFs collected in DoRothEA database.")
  tf.target.overlap <- tf.target[tf.target$tf %in% tf.lst,] %>% split(., .$tf)
  message ("There are ", length(names(tf.target.overlap)),
           " common TFs between JASPAR and DoRothEA.")
  ngene.lst <- stats::rpois(ntfs, ngenes)
  gene.lst <- pbmcapply::pbmclapply(seq_along(tf.target.overlap), function(i) {
    x <- tf.target.overlap[[i]]
    x$target[sample(1:nrow(x), min(ngene.lst[i], nrow(x)))]
  }, mc.cores = parallel::detectCores())
  names(gene.lst) <- names(tf.target.overlap)
  message ("Generated regulons composed of ", range(sapply(gene.lst, length))[1],
           "-", range(sapply(gene.lst, length))[2], " genes.")


  # Select enhancers bound by the selected TFs
  cis.links <- suppressWarnings(link_peaks_to_genes(peak.obj = rownames(obj[[atac.assay]]),
                                   gene.obj = Reduce("union", gene.lst), distance = distance,
                                   org = org))
  message ("Identified ", length(unique(Signac::GRangesToString(cis.links))), " enhancers within ",
           distance, " from the selected genes.")
  enhs <- Signac::GRangesToString(cis.links)
  tf.sites <- data.frame(tf = sites$TF[sites$TF %in% tf.lst],
                         peak = Signac::GRangesToString(sites$peak[sites$TF %in% tf.lst])) %>%
    split(., .$tf) %>% sapply(., "[[", 2)


  # Convert TF binding sites to enhancers
  tf.enhs <- pbmcapply::pbmclapply(tf.sites, function(x) {
    hit.m <- overlap_peak_lst(lst1 = Signac::StringToGRanges(x), lst2 = cis.links)
    summ <- Matrix::summary(hit.m)
    enhs[unique(summ$j)]
  }, mc.cores = parallel::detectCores())


  # Ensure each gene was linked to at least one enhancer
  en.regs <- pbmcapply::pbmclapply(names(tf.enhs), function(x) {
    Reduce("c", lapply(gene.lst[[x]], function(y) {
      cis.links[cis.links$gene == y & enhs %in% tf.enhs[[x]]] %>%
        head(n = max(1, rpois(1, gene.links)))
    }))
  }, mc.cores = parallel::detectCores())
  names(en.regs) <- names(tf.enhs)
  message ("Generated ", length(en.regs), " eRegulons composed of ",
           range(sapply(en.regs, length))[1], "-", range(sapply(en.regs, length))[2],
           " enhancer-gene relations.")


  # Build hybrid bicluster (HBCs) based on eRegulons
  rna.m <- Seurat::GetAssayData(object = obj, slot = "counts", assay = "RNA")
  atac.m <- Seurat::GetAssayData(object = obj, slot = "counts", assay = atac.assay)
  hbc.lst <- c()
  ncell.lst <- stats::rpois(length(en.regs), ncells)


  # Initializes the progress bar
  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(en.regs), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar


  for (i in seq_along(en.regs)) {
    Sys.sleep(0.1) # Remove this line and add your code


    x <- list(
      TF = names(en.regs)[i],
      genes = unique(en.regs[[i]]$gene),
      enhancers = unique(Signac::GRangesToString(en.regs[[i]]))
    )


    xcells <- sample(1:ncol(obj), ncell.lst[[i]], replace = FALSE)
    rna.m[x$genes, xcells] <- 1 + Reduce("rbind", parallel::mclapply(apply(rna.m[x$genes, xcells], 1, max),
                                                           function(x) rep(x, length(xcells)),
                                                           mc.cores = parallel::detectCores()))
    atac.m[x$enhancers, xcells] <- 1 + Reduce("rbind", parallel::mclapply(apply(atac.m[x$enhancers,
                                                                                       xcells], 1, max),
                                                                function(x) rep(x, length(xcells)),
                                                                mc.cores = parallel::detectCores()))
    x$cells <- colnames(obj)[xcells]
    hbc.lst[[i]] <- x


    utils::setTxtProgressBar(pb, i)
  }


  # Add genes, enhancers, and cells not associated with eRegulons
  hbc.rna <- rna.m[Reduce("union", sapply(hbc.lst, "[[", "genes")),
                   Reduce("union", sapply(hbc.lst, "[[", "cells"))]
  hbc.atac <- atac.m[Reduce("union", sapply(hbc.lst, "[[", "enhancers")),
                   Reduce("union", sapply(hbc.lst, "[[", "cells"))]
  message ("Simulated hybrid biclusters (HBCs) cover ", nrow(hbc.rna),
           " genes, ", nrow(hbc.atac), " enhancers, and ", ncol(hbc.rna),
           " cells.")


  # Add background genes, enhancers, and cells
  choice.genes <- setdiff(rownames(rna.m), rownames(hbc.rna))
  choice.genes <- choice.genes[!grepl("^A", choice.genes) & !grepl("\\.", choice.genes)]
  bg.genes <- sample(choice.genes,
                     all.genes - nrow(hbc.rna))
  bg.enhs <- sample(setdiff(rownames(atac.m), rownames(hbc.atac)),
                    all.enhs - nrow(hbc.atac))
  bg.cells <- sample(setdiff(colnames(atac.m), colnames(hbc.rna)),
                     all.cells - ncol(hbc.rna))
  message ("We have ", length(bg.genes), " genes, ", length(bg.enhs), " enhancers, and ",
           length(bg.cells), " cells as background.")


  # Generate the simulated matrices
  simul.rna <- rna.m[Reduce("union", list(rownames(hbc.rna), bg.genes, tf.lst)), 
                     c(colnames(hbc.rna), bg.cells)]
  simul.atac <- atac.m[c(rownames(hbc.atac), bg.enhs), c(colnames(hbc.atac), bg.cells)]
  rand.rna <- simul.rna[sample(nrow(simul.rna)), sample(ncol(simul.rna))]
  rand.atac <- simul.atac[sample(nrow(simul.atac)), sample(ncol(simul.atac))]
  message ("Generated simulated scRNA-seq and scATAC-seq matrices of sizes: ",
           nrow(rand.rna), " x ", ncol(rand.rna), " and ",
           nrow(rand.atac), " x ", ncol(rand.atac), ".")


  # Create the Seurat object
  simul.obj <- rna_atac_matrices_to_Seurat(rna_counts = rand.rna,
                                           atac_counts = rand.atac,
                                           org = org)


  message ("The simulated scRNA-seq and scATAC-seq data contains:\n",
           "1. Hybrid biclusters (HBCs) in nested list,\n",
           "2. Seurat object containing the scRNA-seq and scATAC-seq data.")
  return( list(HBCs = hbc.lst, Seurat = simul.obj) )
}



#' Predict cell-type-specific enhancer-driven gene regulatory networks (eGRNs)
#' 
#' @description Transcription factors (TFs) interact with chromatin regions (herein we call them enhancers), 
#' to regulate the downstream target genes. Joint profiling of scRNA-seq and scATAC-seq shed light upon 
#' gene expression and chromatin accessibility in each cell, providing a great opportunity to inspect 
#' the underlying gene regulatory mechanisms in different cell types. Enhancer-driven 
#' gene regulatory networks (eGRNs) are gene regulatory networks (GRNs) in which TFs regulate their target genes 
#' via binding accessible enhancers in a cell type/state/subpopulation. This function takes as input a 
#' \code{Seurat} object and a list of eRegulons, and then outputs cell-type-specific eGRNs.
#' 
#' @importFrom dplyr %>%
#' @importFrom Seurat DefaultAssay SCTransform RunPCA RunUMAP FindMultiModalNeighbors FindClusters
#' @importFrom Signac RunTFIDF FindTopFeatures RunSVD
#'
#' @export
#' @rdname get_cts_en_GRNs
#'
#' @param obj An \code{Seurat} object.
#' @param celltype The metadata column indicating the cell types or clusters, "seurat_clusters" by default.
#' @param en.regs A list of cell-type-specific eRegulons, each of which contains the following attributes:
#'    \item{terminal} {The \code{IRIS-FGM} block used to predict the eRegulon.}
#'    \item{Tier} {The tier of the TF-enhancer relations: 1 represents JASPAR annotations; 
#'    2 denotes motif scanning.}
#'    \item{TF} {The TF of the eRegulon.}
#'    \item{genes} {Genes of the eRegulon.}
#'    \item{peaks} {Enhancers of the eRegulon.}
#'    \item{cells} {Cells where the eRegulon is active.}
#'    \item{atac.ratio} {The ratio of cells where the eRegulon enhancers are 
#'    accessible against cells in which the eRegulon genes are expressed.}
#'    \item{score} {The eRegulon score.}
#'    \item{weight} {The eRegulon weight.}
#'    \item{links} {The enhancer-gene relations saved in \code{GRanges} object.}
#'    \item{seed} {The seed to obtain the eRegulon.}
#' @param peak.assay The chromatin accessibility assay, "ATAC" by default.
#' @param rna.dims The number of dimensions for RNA dimension reduction, 50 by default.
#' @param atac.dims The number of dimensions for ATAC dimension reduction, 50 by default.
#' @param out.dir The directory to save the intermediate results or final results, "./" by default.
#' @param padj.cutoff The cutoff of adjusted p-value of hyper-geometric test, 0.05 by default.
#' 
#' @return Returns a list of eGRNs in each cell type saved in \code{GRanges} object, the name of each of which 
#' is cell type, and the \code{GRanges} object contains metadata columns "gene" and "TF".
#'
#'#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' 
get_cts_en_GRNs <- function(obj = NULL, celltype = "seurat_clusters",
                            en.regs = NULL, peak.assay = "ATAC",
                            rna.dims = 50, atac.dims = 50,
                            padj.cutoff = 0.05,
                            out.dir = "./") {

  # Information
  message ("The Seurat object contains ", nrow(obj[['RNA']]), " genes, ", nrow(obj[[peak.assay]]),
           " enhancers, and ", ncol(obj), " cells.")


  if (!celltype %in% colnames(obj@meta.data)) {
    # RNA analysis
    if (!"SCT" %in% names(obj@assays)) {
      DefaultAssay(obj) <- "RNA"
      obj <- SCTransform(obj, verbose = FALSE)
    } else {
      DefaultAssay(obj) <- "SCT"
    }
    obj <- obj %>% RunPCA() %>% RunUMAP(dims = 1:rna.dims, reduction.name = 'umap.rna',
                                        reduction.key = 'rnaUMAP_')
    message ("Reduced dimensions for transcriptome.")


    # ATAC analysis
    DefaultAssay(obj) <- peak.assay
    if (is.null(obj[[peak.assay]]@var.features)) {
      obj <- RunTFIDF(obj)
      obj <- FindTopFeatures(obj, min.cutoff = 'q0')
    }
    obj <- RunSVD(obj)
    obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:atac.dims,
                   reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    message ("Reduced dimensions for chromatin accessibility.")


    # WNN analysis
    obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                                   dims.list = list(1:rna.dims, 2:atac.dims))
    obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                   reduction.key = "wnnUMAP_")
    obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
    message ("Identified ", length(levels(obj$seurat_clusters)), " cell clusters.\n",
             "Cell clusters are saved in obj$seurat_clusters.\n",
             "The Seurat object after cell clustering is saved in: ", out.dir,
             "Obj_clustered.qsave.")
    celltype <- "seurat_clusters"
  }


  # Hyper-geometric test
  type.cells <- split(obj[[celltype]], f = obj[[celltype]]) %>% sapply(., rownames)
  pval.m <- Reduce("rbind", pbmcapply::pbmclapply(seq_along(en.regs), function(i) {
    group1 <- en.regs[[i]]$cells
    sapply(seq_along(type.cells), function(j) {
      return(phyper(length(intersect(group1, type.cells[[j]])),
             length(type.cells[[j]]),
             ncol(obj) - length(type.cells[[j]]),
             length(group1),
             lower.tail = TRUE))
    })
  }, mc.cores = parallel::detectCores()))
  rownames(pval.m) <- seq_along(en.regs)
  colnames(pval.m) <- names(type.cells)
  message ("Finished hyper-geometric test of ", length(en.regs),
           " eRegulons against ", length(type.cells),
           ' cell types, obtaining p-values ranging from ',
           range(pval.m)[1], " to ", range(pval.m)[2])


  # Bonferroni correction
  padj.m <- pval.m * length(en.regs) * length(type.cells)
  padj.m[padj.m > 1] <- 1.0
  message ("Finished Bonferroni correction and obtained adjusted p-values ranging from ",
           range(padj.m)[1], " to ", range(padj.m)[2])


  # Identify cell-type-specific eRegulons
  cts.reg.ids <- pbmcapply::pbmclapply(1:ncol(padj.m), function(j) {
    return(which(padj.m[, j] < padj.cutoff))
  }, mc.cores = parallel::detectCores())
  names(cts.reg.ids) <- colnames(padj.m)


  # Construct cell-type-specific eGRNs
  cts.reg.ids <- cts.reg.ids[sapply(cts.reg.ids, length) > 0]
  cts.en.grns <- pbmcapply::pbmclapply(names(cts.reg.ids), function(i) {
    links <- unlist(as(lapply(en.regs[cts.reg.ids[[i]]], function(y) {
      GenomicRanges::mcols(y$links)$TF <- y$TF
      y$links
    }), "GRangesList"))
    return(list(links = links,
                cell.type = i,
                cells = type.cells[[i]]))
  }, mc.cores = parallel::detectCores())
  names(cts.en.grns) <- names(cts.reg.ids)
  message ("Identified enhancer gene regulatory networks (eGRNs) for cell types (n = ", 
           length(cts.reg.ids), ") : ",
           paste(names(cts.reg.ids), collapse = ", ") )


  return(cts.en.grns)
}



#' Identify cell-type-specific enhancer regulons (eRegulons)
#' 
#' @description Enhancer regulons (eRegulons) are active in various cell subpopulations. Some 
#' of them may be enriched in one or more cell types. The successful identification of eRegulons 
#' at the single-cell level can improve the detection of heterogeneous transcriptional regulatory 
#' mechanisms across various cell types and allows for reliable constructions of global gene regulatory 
#' networks encoded in complex diseases. Hence, it is critical to study cell-type-specific eRegulons. 
#' eRegulons specific to a cell type are called cell-type-specific eRegulons. This function takes as 
#' input a \code{Seurat} object (composed of scRNA-seq and scATAC-seq) and cell-type-specific enhancer-drive 
#' gene regulatory networks (eGRNs) and then identify the cell-type-specific eRegulons in each cell type/cluster.
#' 
#' @importFrom dplyr %>%
#'
#' @export
#' @rdname get_cts_en_regs
#'
#' @param obj An \code{Seurat} object.
#' @param celltype The metadata column indicating the cell types or clusters, "seurat_clusters" by default.
#' @param cts.en.grns Cell-type-specific eGRNs saved in \code{GRanges} object, the name of each of which 
#' is cell type, and the \code{GRanges} object contains metadata columns "gene" and "TF".
#' @param peak.assay The chromatin accessibility assay, "ATAC" by default.
#' @param de.genes A list of differentially expressed genes (DEGs).
#' @param out.dir The directory to save the intermediate results or final results, "./" by default.
#' @param accessibility Whether perform differential accessibility analysis, FALSE by default.
#' @param padj.cutoff The cutoff of adjusted p-values of differential expression, 0.05 by default.
#'
#' @return Returns a list of cell-type-specific eRegulons,  each of which contains the following attributes:
#'    \item{TF} {The TF of the cell-type-specific eRegulons.}
#'    \item{genes} {Genes of the cell-type-specific eRegulons.}
#'    \item{enhancers} {Enhancers of the cell-type-specific eRegulons.}
#'    \item{cells} {Cells where the celltype.}
#'    \item{links} {The enhancer-gene relations saved in \code{GRanges} object.}
#'    \item{celltype} {The celltype of the cell-type-specific eRegulons.}
#'
#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' 
get_cts_en_regs <- function(obj = NULL, peak.assay = "ATAC", de.genes = NULL,
                            cts.en.grns = NULL, accessibility = FALSE, out.dir = "./",
                            min.pct = 0.25, logfc.threshold = 0.25, padj.cutoff = 0.05) {

  # Information
  message ("The Seurat object contains ", nrow(obj[['RNA']]), " genes, ", nrow(obj[[peak.assay]]),
           " enhancers, and ", ncol(obj), " cells.")


  # Differential analyses
  if (is.null(de.genes)) {
    message ("Need to filter genes based on differential expression.\n", 
             "Predicting differentially expressed genes (DEGs) ...")
    Seurat::DefaultAssay(obj) <- "RNA"
    Idents(obj) <- stats::setNames(obj@meta.data[, celltype], colnames(obj))
    future::plan("multiprocess", workers = 4)
    de.genes <- Seurat::FindAllMarkers(obj, only.pos = TRUE,
                                       min.pct = min.pct, logfc.threshold = logfc.threshold)
  }
  de.genes <- de.genes[de.genes$p_val_adj < padj.cutoff,]
  message ("We have ", nrow(de.genes), " differentially expressed genes (DEGs).")
  if (accessibility) {
    message ("Need to filter genes based on differential accessibility.\n", 
             "Predicting differentially accessible regions (DARs) ...")
    Seurat::DefaultAssay(obj) <- peak.assay
    Idents(obj) <- stats::setNames(obj@meta.data[, celltype], colnames(obj))
    future::plan("multiprocess", workers = 4)
    da.peaks <- Seurat::FindAllMarkers(obj, only.pos = TRUE,
                                    min.pct = min.pct, logfc.threshold = logfc.threshold)
    da.peaks <- da.peaks[da.peaks$p_val_adj < padj.cutoff,]
    message ("We have ", nrow(da.peaks), " differentially accessible regions (DARs).")
  }


  # Build cell-type-specific eRegulons
  cts.en.regs <- do.call("c", pbmcapply::pbmclapply(cts.en.grns, function(x) {
    links <- x$links
    splitted <- split(links, f = links$TF)
    lapply(seq_along(splitted), function(i) {
      cts.en.reg <- list(
        TF = names(splitted[i]),
        celltype = x$cell.type,
        cells = x$cells,
        genes = intersect(unique(splitted[[i]]$gene),
                          de.genes[de.genes$cluster == x$cell.type, "gene"])
        )
      if (accessibility) {
        cts.en.reg$enhancers <- intersect(unique(Signac::GRangesToString(splitted[[i]])),
                                          da.peaks[da.peaks$cluster == x$cell.type, "gene"])
      } else {
        cts.en.reg$enhancers <- unique(Signac::GRangesToString(splitted[[i]]))
      }
      cts.en.reg$links <- splitted[[i]][Signac::GRangesToString(splitted[[i]]) %in%
                                          cts.en.reg$enhancers &
                                          splitted[[i]]$gene %in% cts.en.reg$genes]
      cts.en.reg$enhancers <- Signac::GRangesToString(cts.en.reg$links)


      return(cts.en.reg)
    })
  }, mc.cores = parallel::detectCores()))
  names(cts.en.regs) <- NULL
  message ("Identified ", length(cts.en.regs), " cell-type-specific eRegulons.")
  
  
  return(cts.en.regs)
}



#' Calculate the precision, recall, and f-scores of overlaps between 
#' two \code{GRanges} objects indicating enhancer-gene relations.
#' 
#' @description Given two \code{GRanges} objects, each of which has the meta column named 
#' "gene", this function calculates the overlaps between them. Based on the calculated overlaps, 
#' this function computes precision, recall, and f-score. This function aims to assess the enhancer-gene 
#' relations in eRegulons or eGRNs. We may use the enhancer-target pair databases, e.g., EnhancerAtlas or 
#' scEnhancer of the same tissues or cell lines.
#' 
#' @import dplyr
#' @export
#' @rdname intersect_enhancer_gene_relations
#' 
#' @param x The first \code{GRanges} object saving enhancer-gene relations.
#' @param y The second \code{GRanges} object saving enhancer-gene relations.
#' @return Return a \code{data.frame} indicating overlapped \code{GRanges} objects, 
#' containing the following columns:
#'    \item{x.peak} {The enhancer in the first \code{GRanges} object for each pair of 
#'    overlapped \code{GRanges} objects.}
#'    \item{y.peak} {The enhancer in the second \code{GRanges} object for each pair of 
#'    overlapped \code{GRanges} objects.}
#'    \item{gene} {The gene for each pair of 
#'    overlapped \code{GRanges} objects.}
#'
#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' @references Gao, T., & Qian, J. (2020). EnhancerAtlas 2.0: an updated resource with enhancer 
#' annotation in 586 tissue/cell types across nine species. Nucleic acids research, 48(D1), D58-D64.
#' @references Gao, T., Zheng, Z., Pan, Y., Zhu, C., Wei, F., Yuan, J., ... & Qian, J. (2022). 
#' scEnhancer: a single-cell enhancer resource with annotation across hundreds of tissue/cell 
#' types in three species. Nucleic acids research, 50(D1), D371-D379.
#' 
intersect_enhancer_gene_relations <- function(x, y) {
  
  # Calculate intersections
  overlap.genes <- intersect(unique(x$gene), unique(y$gene))
  x.overlap <- x[x$gene %in% overlap.genes]
  y.overlap <- y[y$gene %in% overlap.genes]
  query.subject <- GenomicAlignments::findOverlaps(query = x.overlap,
                                                   subject = y.overlap)
  # library(Repitools)
  x.peaks <- Signac::GRangesToString(x)
  x.genes <- x$gene
  y.peaks <- Signac::GRangesToString(y)
  overlap.df <- data.table::rbindlist(lapply(seq_along(query.subject), function(i) {
    list(x.peak = x.peaks[queryHits(query.subject)[i]], 
         y.peak = y.peaks[subjectHits(query.subject)[i]], 
         gene = x.genes[queryHits(query.subject)[i]])
  }))
  return(overlap.df)
}



#' Calculate the precision, recall, and f-scores of the overlaps between 
#' two lists of \code{GRanges} objects indicating enhancer-gene relations
#' 
#' @description This function has the same functionality as the function "intersect_enhancer_gene_relations". 
#' The only difference is that this function aims to compare two lists of \code{GRanges} objects.
#' 
#' @export
#' @rdname intersect_enhancer_gene_relations_in_batch
#' 
#' @param link.pairs The first list of \code{GRanges} objects saving enhancer-gene relations.
#' @param ep.ll The second list of \code{GRanges} objects saving enhancer-gene relations.
#' @param only.overlap Only consider the \code{GRanges} objects of which genes were overlapped against databases, 
#' TRUE by default.
#' @param max.score Which score will be used to select the best query-hit pairs of \code{GRanges} objects,
#' "precision" by default.
#' 
#' @return Returns a \code{data.frame} to indicate the query-hit pairs as well as precision, recall, and f-score. 
#' The \code{data.frame} contains the following columns:
#'    \item{EP} {The ID of overlapped enhancer-gene pairs in databases.}
#'    \item{precision} {The precision of the overlaps between enhancer-gene relations between
#'     enhancer regulons (eRegulons) and that in databases.}
#'    \item{recall} {The recall of the overlaps between enhancer-gene relations between
#'     eRegulons and that in databases.}
#'    \item{fscore} {The f-score of the overlaps between enhancer-gene relations between
#'     eRegulons and that in databases.}
#' 
#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' @references Gao, T., & Qian, J. (2020). EnhancerAtlas 2.0: an updated resource with enhancer 
#' annotation in 586 tissue/cell types across nine species. Nucleic acids research, 48(D1), D58-D64.
#' @references Gao, T., Zheng, Z., Pan, Y., Zhu, C., Wei, F., Yuan, J., ... & Qian, J. (2022). 
#' scEnhancer: a single-cell enhancer resource with annotation across hundreds of tissue/cell 
#' types in three species. Nucleic acids research, 50(D1), D371-D379.
#' 
intersect_enhancer_gene_relations_in_batch <- function(link.pairs, ep.ll, 
                                                       only.overlap = FALSE, 
                                                       max.score = "precision") {
  
  # Load the peak-gene linkages
  if (only.overlap) {
    message ("Evaluating enhancer-gene relations only for the ones overlapped with enhancers ...\n", 
             "There are in total ", length(link.pairs), " enhancer-gene relations before overlapping.")
    link.pairs <- link.pairs[unique(queryHits(findOverlaps(link.pairs, do.call("c", ep.ll))))]
    message ("There are in total ", length(link.pairs), " enhancer-gene relations after overlapping.")
  }
  n.pairs <- length(link.pairs)
  if (n.pairs < 1) {
    warning ("There is no enhancer-gene linkages input!")
    if (only.overlap) {
      warning ("Please note: we only consider the enhancers overlapped with enhancers!\n", 
               "Enhancers not overlapped with enhancers may still exist.")
    }
    return(data.frame(EP = NA, precision = NA, recall = NA, fscore = NA))
  }
  message ("There are in total ", n.pairs, " enhancer-gene pairs input.")
  
  
  # Calculate the overlaps
  message ("Performing enrichment analysis against the enhancer-target pairs in databases ...")
  score.dt <- data.table::rbindlist(pbapply::pblapply(seq_along(ep.ll), function(i) {
    ee <- ep.ll[[i]]
    overlap.genes <- intersect(unique(link.pairs$gene), unique(ee$gene))
    overlap.query <- link.pairs[link.pairs$gene %in% overlap.genes]
    overlap.subject <- ee[ee$gene %in% overlap.genes]
    overlap.hit <- GenomicAlignments::findOverlaps(query = overlap.query, subject = overlap.subject)
    hit.query <- length(unique(queryHits(overlap.hit)))
    hit.subject <- length(unique(subjectHits(overlap.hit)))
    precision <- hit.query / n.pairs
    recall <- hit.subject / length(ee)
    fscore <- 2 * precision * recall / (precision + recall)
    list(EP = i, precision = precision, recall = recall, fscore = fscore)
  }))
  score.dt <- na.omit(score.dt)
  max.scores <- score.dt[which.max(score.dt[[max.score]]),]
  return(max.scores)
}



#' Calculate the p-values of overlaps between two \code{GRanges} objects
#' 
#' @description Given two \code{GRanges} objects, this function calculates the overlaps between them. 
#' Based on the calculated overlaps, 
#' this function relies upon \code{regionR} to perform permutation test to assess the significance of 
#' overlaps between the two \code{GRanges} objects. Finally, a p-value will be calculated. Usually, we 
#' perform comparison for the enhancer set of an eRegulon against a series of ChIP-seq peaks in the same tissues 
#' or cell lines.
#' 
#' @export
#' @rdname intersect_peaks
#' 
#' @param x The first \code{GRanges} object or \code{data.frame}.
#' @param y The first \code{GRanges} object or \code{data.frame}.
#' @param n.times The number of times of permutation, 1000 by default.
#' @param alternative The direction of alternative hypothesis, "greater" by default.
#' 
#' @return Return a numeric p-value.
#'
#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' @references Gel, B., Díez-Villanueva, A., Serra, E., Buschbeck, M., Peinado, M. A., & Malinverni, R. (2016). 
#' regioneR: an R/Bioconductor package for the association analysis of genomic regions based on permutation tests. 
#' Bioinformatics, 32(2), 289-291.
#' 
intersect_peaks <- function(x, y, n.times = 100, alternative = "greater") {
  
  # Calculate intersections
  message ("Perform significance test between two GRange objects composed of ", 
           length(x), " and ", length(y), " peaks.")
  regioneR::overlapPermTest(A = x, B = y, ntimes = n.times, alternative = alternative)
}



#' Calculate the p-values of overlaps between two lists of \code{GRanges} objects
#' 
#' @description This function has the same functionality as the function "intersect_peaks". The 
#' only difference is that this function compares two lists of \code{GRanges} objects and returns 
#' a \code{data.frame} composed of the IDs in each \code{GRanges} object as well as the p-values.
#' 
#' @export
#' @rdname intersect_peaks_in_batch
#' 
#' @param x.ll The first list of \code{GRanges} objects.
#' @param y.ll The first list of \code{GRanges} objects.
#' @param n.times The number of times of permutation, 1000 by default.
#' 
#' @return Returns a \code{data.frame} composed of the IDs of significantly overlapped 
#' \code{GRanges} objects and p-values.
#'
#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' @references Gel, B., Díez-Villanueva, A., Serra, E., Buschbeck, M., Peinado, M. A., & Malinverni, R. (2016). 
#' regioneR: an R/Bioconductor package for the association analysis of genomic regions based on permutation tests. 
#' Bioinformatics, 32(2), 289-291.
#' 
intersect_peaks_in_batch <- function(x.ll, y.ll, n.times = 100) {
  
  message ("Perform significance test between two lists composed of ", 
           length(x.ll), " and ", length(y.ll), " GRange objects.")
  pval.df <- data.frame(numeric(), numeric(), 
                        numeric())
  for (i in seq_along(x.ll)) {
    for (j in seq_along(y.ll)) {
      pval <- suppressMessages(intersect_peaks(x = x.ll[[i]], y = y.ll[[j]], n.times = n.times))
      pval.df <- rbind(pval.df, c(i, j, pval))
    }
  }
  message ("Finished performing enrichment analysis between two lists.\n", 
           nrow(pval.df), " pairs of peak lists have been compared.")
  colnames(pval.df) <- c("query", "subject", "pval")
  return(pval.df)
}



#' Perform enrichment analysis for eRegulon genes against gene ontology (GO) terms 
#' and KEGG pathways
#' 
#' @description This function aims to calculate the enrichment of genes in each eRegulon against 
#' gene ontology (GO) terms or KEGG patthways. Finally, this function returns a list of \code{data.frame} 
#' objects, each of which represents one functional genomics databases, e.g., Biological Process, 
#' Molecular Function, Cellular Component, and KEGG pathways (human or mouse).
#' 
#' @export
#' @rdname enrich_genes
#' 
#' @param regs The list of enhancer regulons (eRegulons) or cell-type-specific eRegulons.
#' @param dbs The list of databases to run enrichment analysis, c("GO", "KEGG") by default.
#' @param org The organism, "human" by default.
#' 
#' @references Li, Y., Ma, A., Wang, Y., Wang, C., Chen, S., Fu, H., Liu, B. and Ma, Q., 2022. 
#' Enhancer-driven gene regulatory networks inference from single-cell RNA-seq and ATAC-seq data. 
#' bioRxiv, pp.2022-12.
#' @references Kuleshov, M. V., Jones, M. R., Rouillard, A. D., Fernandez, N. F., Duan, Q., 
#' Wang, Z., ... & Ma'ayan, A. (2016). Enrichr: a comprehensive gene set enrichment analysis web 
#' server 2016 update. Nucleic acids research, 44(W1), W90-W97.
#' 
enrich_genes <- function(regs = NULL, dbs = c("GO", "KEGG"), 
                         org = "human" ) {
  
  # Determine databases
  genes.ll <- sapply(regs, "[[", "genes")
  databases <- c()
  if ("GO" %in% dbs) {
    message ("Loading GO databases ...")
    databases <- c("GO_Molecular_Function_2018",
                   "GO_Cellular_Component_2018",
                   "GO_Biological_Process_2018")
  }
  if ("KEGG" %in% dbs) {
    if (org == "mouse") {
      databases <- c(databases, "KEGG_2019_Mouse")
    } else {
      databases <- c(databases, "KEGG_2019_Human")
    }
    message ("Loading KEGG database of ", org, " ...")
  }
  
  
  # Run enrichment analyses
  message ("Running enrichment analyses for ", length(genes.ll), " gene lists ...")
  quiet(require(enrichR))
  enriched <- pbmcapply::pbmclapply(genes.ll, 
                                    mc.cores = 4, 
                                    function(x) {
                                      enrichr(genes = x, databases = databases)
                                    })
  names(enriched) <- seq_along(genes.ll)
  
  
  # Parse the results
  message ("Merging enrichment analyses results ...")
  processed <- pbmcapply::pbmclapply(databases, 
                                     mc.cores = parallel::detectCores(), 
                                     function(x) {
                                       message ("Preparing enrichment analysis results for ", x, " ...")
                                       x.prosessed <- Reduce("rbind", lapply(seq_along(enriched), function(i) {
                                         cbind(Id = rep(names(enriched)[i], nrow(enriched[[i]][[x]])), enriched[[i]][[x]])
                                       }) )
                                     })
  names(processed) <- databases
  message ("There are ", paste(sapply(processed, nrow), collapse = ", "), 
           " terms/pathways for ", paste(names(processed), collapse = ", "), 
           ".")
  return(processed)
}
