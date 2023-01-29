#' Given a list of TFs, output associated proteins from BioGrid interacting with these TFs.
#' 
#' @description Collaborations among proteins, e.g., transcription factors (TFs), contribute to the regulation 
#' of target genes via binding to the same set of enhancers. To narrow down the scope of 
#' enhancer regulons (eRegulons), we encourage users to select a TF list when running \code{stream2}. This function 
#' aims to enlarge the TF list by incorporating the proteins interacting with TFs in the user-provided list.
#' 
#' @keywords internal
#' @export
#' @rdname get_associated_prot
#' 
#' @import dplyr
#' 
#' @param key.TFs The list of interested TFs.
#' @param db The database of protein relations, including BioGrid, IntAct, STRING, and CORUM, 
#' BioGrid by default. Other databases have been under development and will be added soon.
#' @param org The organism, human by default.
#' 
#' @return This function returns a list of TF names.
#' 
get_associated_prot <- function(key.TFs = NULL, db = "BioGrid", org = "human") {
  
  # Parameters
  message ("Identifying proteins interacting with the ", 
           length(key.TFs), " TFs provided by user.")
  
  
  # Load database
  if (db == "STRING") {
    message ("Loaded the STRING database.")
    if (org == "mouse") {
      species <- 10090
    } else {
      species <- 9606
    }
    tmp.db <- quiet(STRINGdb::STRINGdb$new(species = species))
    tmp.graph <- quiet(tmp.db$get_graph())
    ensembl.lst <- igraph::as_adj_list(tmp.graph, mode = "all") %>% lapply(., function(x) {
      strsplit(igraph::as_ids(seq = x), split = "\\.") %>% sapply(., "[[", 2)
    })
    names(ensembl.lst) <- strsplit(names(ensembl.lst), split = "\\.") %>% 
      sapply(., "[[", 2)
    if (org == "mouse") {
      mart <- biomaRt::useMart(host = 'grcm39.ensembl.org',
                      biomart = 'ENSEMBL_MART_ENSEMBL',
                      dataset = 'mmusculus_gene_ensembl')
    } else {
      mart <- biomaRt::useMart(host = 'grch38.ensembl.org',
                               biomart = 'ENSEMBL_MART_ENSEMBL',
                               dataset = 'hsapiens_gene_symbol')
    }
    mart_results <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                  "ensembl_peptide_id"),
                                   filters = "ensembl_peptide_id", 
                                   values = protein_ids,
                                   mart = mart)
    ix <- match(protein_ids, mart_results$ensembl_peptide_id)
    ix <- ix[!is.na(ix)]
    relation.lst <- pbmcapply::pbmclapply(ensembl.lst, mc.cores = parallel::detectCores(), 
                                          function(protein_ids) {
                                            newnames <- protein_ids
                                            newnames[match(mart_results[ix,'ensembl_peptide_id'], 
                                                           newnames)] <-
                                              mart_results[ix, 'ensembl_gene_symbol']
                                            newnames
                                          })
    names(relation.lst)[match(mart_results[ix,'ensembl_peptide_id'], 
                              names(relation.lst))] <- mart_results[ix, 'ensembl_gene_symbol']
  } else {
    message ("Loaded the BioGrid database.")
    if (org == "mouse") {
      data(MouseBioGRIDInteractionOfficial)
      tmp.lst <- MouseBioGRIDInteractionOfficial
    } else {
      data(HumanBioGRIDInteractionOfficial)
      tmp.lst <- HumanBioGRIDInteractionOfficial
    }
    relation.lst <- lapply(tmp.lst, "[[", "interactors")
    names(relation.lst) <- sapply(tmp.lst, "[[", "name")
  }
  if (org == "mouse") {
    taxid <- 10090
  } else {
    taxid <- 9606
  }
  # tmp.lst <- PItools::fullInteractome(taxid = taxid, database = database,
  #                            format = "tab27", clean = TRUE,
  #                            protein_only = TRUE,
  #                            directory = NULL) # keep data files inside R library - default
  # relation.lst <- tmp.lst # To-do: add the conversion steps!
  message ("Loaded ", length(relation.lst), " protein relations from ", 
           db, " database of ", org, ".")
  
  
  # Obtain proteins
  TFs <- c(
    unique(unlist(relation.lst[intersect(names(relation.lst), key.TFs)])),
    key.TFs
  )
  message ("Identified ", length(TFs), " TFs in total.")
  return(TFs)
}
