#' @importFrom magrittr %>%

#' Title
#'
#' @param adat
#' @param condition
#'
#' @return
#' @export
#'
#' @examples
create_protdata_from_olink <- function(npx, condition = NULL, filter = TRUE) {
  npx <- as.data.frame(npx)
  dat <- olink_sample_out(npx, filter)
  return(create_protdata(dat, intensity_cols = c(4:length(colnames(dat))), condition, method = "Olink"))
}


olink_sample_out=function(my_npx, filter){
  if(filter){
    npx_wide <- my_npx |>
      dplyr::mutate(NPX = ifelse(NPX < LOD, NA, NPX)) |>
      # dplyr::filter(SampleType == "SAMPLE") |>
      # dplyr::filter(AssayType == "assay") |>
      dplyr::select(SampleID, UniProt,Assay, OlinkID, NPX) |>
      tidyr::pivot_wider(names_from = SampleID, values_from = NPX, values_fn = mean) |>
      dplyr::rename(Protein_Group = UniProt, Genes = Assay)
  }else{
    npx_wide <- my_npx |>
      # dplyr::filter(AssayType == "assay") |>
      dplyr::select(SampleID, UniProt,Assay, OlinkID, NPX) |>
      tidyr::pivot_wider(names_from = SampleID, values_from = NPX,values_fn = mean) |>
      dplyr::rename(Protein_Group = UniProt, Genes = Assay)
  }
  return(npx_wide)
}
