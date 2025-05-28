#' @importFrom magrittr %>%

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
  metas <- data.frame()

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


    # merge protein metadata
    if (nrow(metas) == 0) {
      metas <- meta
    } else {
      common_cols <- intersect(colnames(metas), colnames(meta))
      metas <- metas[,common_cols]
      meta <- meta[,common_cols]
      meta_filtered <- dplyr::anti_join(meta, metas, by = names(meta)[1])
      if (nrow(meta_filtered) > 0) {
        metas <- dplyr::bind_rows(metas, meta_filtered)
      }
    }

    # merge conditions
    if (nrow(conds) == 0) {
      conds <- cond
    } else {
      # Identify the matching columns between current cond and accumulated cond
      common_cols <- intersect(colnames(conds), colnames(cond))
      # Subset both 'conds' and 'cond' to keep only the common columns
      conds <- dplyr::bind_rows(conds[, common_cols, drop = FALSE], cond[, common_cols, drop = FALSE])
    }
  }
  tasd<<-datas
  dsfd<<-metas
  datas <- datas[, c("Protein.Group", setdiff(names(datas), "Protein.Group"))]
  #merge metas and datas
  num_metas_cols <- length(names(metas))
  names(metas)[1] <- "Protein.Group"
  merged_df <- merge(metas, datas, by = "Protein.Group")
  print("hi")
  names(merged_df)[1] <- names(metas)[1]

  return(create_protdata(merged_df, conds, intensity_cols = (num_metas_cols+1):length(merged_df)))
}
