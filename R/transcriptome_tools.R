#' \code{IRIS-FGM} modeling
#'
#' @keywords internal
#'
call_LTMG <- function(obj) {

  ltmg.df <- as.data.frame(obj@assays$RNA@data) # get the expression matrix saved in a data frame
  ltmg.obj <- IRISFGM::ProcessData(IRISFGM::CreateIRISFGMObject(ltmg.df),
                          normalization = "cpm", IsImputation = F)
  IRISFGM::RunLTMG(ltmg.obj, Gene_use = "all")
}
