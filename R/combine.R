#' Title
#'
#' @param PD_list
#'
#' @return
#' @export
#'
#' @examples
combine <- function(PD_list){
  #make sure each element in PD_list is of class ProtData
  # Ensure each element in PD_list is of class ProtData
  if (!all(sapply(PD_list, inherits, "ProtData"))) {
    stop("All elements in PD_list must be of class ProtData.")
  }

  # Make sure the metadata colnames are the same for all PDs, warning/error if not
  meta_colnames <- lapply(PD_list, function(pd) colnames(pd@prot_meta))
  if (!all(sapply(meta_colnames, function(cols) identical(cols, meta_colnames[[1]])))) {
    warning("Metadata column names are not identical across all ProtData objects.")
  }

  # Make sure there is at least 1 common colname in the condition file between all PDs
  common_condition_cols <- Reduce(intersect, lapply(PD_list, function(pd) colnames(pd@condition)))
  if (length(common_condition_cols) == 0) {
    stop("There is no common column in the condition data between the ProtData objects.")
  }

  datas <- data.frame()
  conds <- data.frame()

  for (pd in PD_list){
    cond <- pd@condition
    dat <- pd@data
    meta <- pd@prot_meta

    dat$Protein.Group <- meta[,1]

    # create df with all experiments
    if (nrow(datas) == 0) {
      datas <- dat
    } else {
      datas <- datas %>%
        merge(dat, by = "Protein.Group", all = TRUE)
    }


    #create metadata for all experiments, keeping only the common columns
    if (nrow(conds) == 0) {
      conds <- cond
    } else {
      # Identify the matching columns between current meta and accumulated metas
      common_cols <- intersect(colnames(conds), colnames(cond))
      # Subset both 'metas' and 'meta' to keep only the common columns
      conds <- dplyr::bind_rows(conds[, common_cols, drop = FALSE], cond[, common_cols, drop = FALSE])
    }
  }
  return(create_protdata(datas, conds))
}
