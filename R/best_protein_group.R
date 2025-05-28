#' @importFrom magrittr %>%

bestProteins <- function(PD, gene_col, prot_col){
  prot_col = "Protein_Group"
  gene_col = "Genes"
  PD = dat_pro
  dat_original <- getData(PD)
  dat <- as.data.frame(getData(PD)) %>%
    # Calculte missing values
    dplyr::mutate(missing_value = rowSums(is.na(.)))%>%
    # Calculate median values
    dplyr::mutate(median = matrixStats::rowMedians(as.matrix(dplyr::select(., -contains("missing_value")) %>%
                                           dplyr::select_if(is.numeric)), na.rm = TRUE)) %>%
    merge(dplyr::select(getProtMeta(PD),c(gene_col, prot_col)),by='row.names') %>%
    #Identify unique Protein_Group with the least missing values and highest median intensity
    dplyr::group_by(..gene_col) %>%
    dplyr::filter(missing_value == min(missing_value)) %>%
    dplyr::slice(which.max(median)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(missing_value, median))
}
