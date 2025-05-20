

#' Title
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
setGeneric("do_t_test", function(object, treatment_samples, control_samples, meta_col) standardGeneric("do_t_test"))


#' Title
#'
#' @param ProtData
#'
#' @return
#' @export
#'
#' @examples
setMethod("do_t_test", "ProtData", function(object, treatment_samples, control_samples) {

  col <- names(object@prot_meta)
  DT <- cbind(object@prot_meta, object@data)
  meta <- object@condition

  n_treatment <- length(treatment_samples)
  n_control <- length(control_samples)
  #ttest table
  # Initial data transformation
  DT_ttest <- DT[,c(col, treatment_samples, control_samples)]

  # Convert NA to 0 and add missing value calculations
  DT_ttest <- DT_ttest %>%
    # Convert NA to 0
    dplyr::mutate(across(where(is.numeric), ~ dplyr::coalesce(., 0))) %>%
    # Calculate missing values
    dplyr::mutate(
      missing_value = rowSums(dplyr::select(., -col) == 0),  # Total missing values per row, excluding PTM column
      missing_value_c = rowSums(dplyr::select(., all_of(control_samples)) == 0),  # Missing values in control samples per row
      missing_value_t = rowSums(dplyr::select(., all_of(treatment_samples)) == 0)  # Missing values in treatment samples per row
    ) %>%
    # Filter features with > 50% missing values
    dplyr::filter(
      !(missing_value_t > (n_treatment / 2) & missing_value_t < n_treatment),
      !(missing_value_c > (n_control / 2) & missing_value_c < n_control),
      missing_value != (n_treatment + n_control)
    ) %>%
    # Select relevant columns
    dplyr::select(-missing_value, -missing_value_c, -missing_value_t)%>%
    as.data.frame()
  # Restore row names from PTM column and remove the PTM column
  # rownames(DT_ttest) <- DT_ttest[,col]
  # DT_ttest[,col]=NULL

  prot_meta <- DT_ttest[,col]
  DT_ttest[,col] <- NULL

  # Perform t-test  on treatment and control columns
  t_test <- apply(DT_ttest, 1, function(x){
    a =factor(c(rep('treatment',n_treatment),
                rep("control",n_control)),
              levels = c('treatment',"control"))
    fvalue=var.test(x~a)
    if (!is.na(fvalue$p.value)){
      if (fvalue$p.value > 0.05){
        result <- t.test(x~a, var.equal = T)
      }else{
        result <- t.test(x~a, var.equal = F)
      }
    }
    treatment_estimate <- as.numeric(unlist(result$estimate[1]))
    control_estimate <- as.numeric(unlist(result$estimate[2]))
    return(data.table::data.table('P.Value'=result$p.value,
                      'treatment_estimate'=treatment_estimate,
                      'control_estimate'=control_estimate)
    )
  })
  t_test <- data.table::rbindlist(t_test)
  result_ttest <- cbind(prot_meta, DT_ttest, t_test)
  # Merge DT and result_ttest based on PTM column
  #result_ttest <- merge(DT[, name, drop=FALSE], result_ttest, by.x=col, by.y='row.names')
  # Calculate logFC and adjust P-values using dplyr
  result_ttest <- result_ttest %>%
    dplyr::mutate(logFC = log2((treatment_estimate+1) / (control_estimate + 1)),
           adj.P.Val = p.adjust(P.Value, method='BH'))
  return(result_ttest)
})

#' Helper function to filter out sparse proteins
filter_features <- function(DT_limma, control_samples, treatment_samples, alpha){
  n_treatment <- length(treatment_samples)
  n_control <- length(control_samples)
  if (n_treatment>3&n_control>3) {
    # Drop rows (protein groups) with > 50% missingness in samples
    DT_limma$missing_value= apply(DT_limma, 1, function(x) sum(x==0))
    DT_limma$missing_value_c= apply(DT_limma[,control_samples], 1, function(x) sum(x==0))
    DT_limma$missing_value_t= apply(DT_limma[,treatment_samples], 1, function(x) sum(x==0))
    DT_limma <- DT_limma[!(DT_limma$missing_value_t >(n_treatment * alpha)&DT_limma$missing_value_t < n_treatment),]
    DT_limma <- DT_limma[!(DT_limma$missing_value_c >(n_control * alpha)&DT_limma$missing_value_c <n_control),]
    DT_limma <- DT_limma[DT_limma$missing_value != (n_treatment+n_control),]
    DT_limma[,grep('missing_value',colnames(DT_limma))]=NULL
  } else{
    # Drop rows (protein groups) with missingness in samples
    DT_limma$missing_value= apply(DT_limma, 1, function(x) sum(x==0))
    DT_limma$missing_value_c= apply(DT_limma[,control_samples], 1, function(x) sum(x==0))
    DT_limma$missing_value_t= apply(DT_limma[,treatment_samples], 1, function(x) sum(x==0))
    DT_limma <- DT_limma[DT_limma$missing_value_t %in% c(0,n_treatment),]
    DT_limma <- DT_limma[DT_limma$missing_value_c %in% c(0,n_control),]
    DT_limma <- DT_limma[DT_limma$missing_value != (n_treatment+n_control),]
    DT_limma[,grep('missing_value',colnames(DT_limma))]=NULL
  }
  return(DT_limma)

}


#' Perform limma differential expression on a ProtData object
#'
#' This function takes a `ProtData` object and two vectors containing the column names
#' of the treatment and control samples. It filters out proteins found in fewer than 50%
#' of samples, performs a limma-based differential expression analysis, and returns
#' a data frame with metadata, intensities, log fold changes, and p-values.
#'
#' @param object A `ProtData` object containing protein intensities, metadata, and condition info.
#' @param treatment_samples Character vector of column names representing treatment samples.
#' @param control_samples Character vector of column names representing control samples.
#'
#' @return A data frame containing filtered proteins with metadata, intensity values, log fold change,
#' p-value, and adjusted p-value.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- do_limma(my_protdata, treatment_samples = c("T1", "T2"), control_samples = c("C1", "C2"))
#' }
setGeneric("do_limma", function(object, treatment_samples, control_samples) standardGeneric("do_limma"))

#' @describeIn do_limma Method for ProtData objects
setMethod("do_limma", "ProtData", function(object, treatment_samples, control_samples) {
  #data
  meta_cols <- names(object@prot_meta)
  DT <- cbind(object@prot_meta, object@data)
  meta <- object@condition

  # treatment_samples=grep(treatment,colnames(Log2_DT),value = T)
  # control_samples=grep(control,colnames(Log2_DT),value = T)
  DT_limma <- DT[,c(treatment_samples, control_samples)]
  n_treatment <- length(treatment_samples)
  n_control <- length(control_samples)
  # Convert NA to 0
  DT_limma[is.na(DT_limma)] <- 0

  # Filter out sparse proteins
  DT_limma <- filter_features(DT_limma, control_samples, treatment_samples, alpha = 0.5)

  #design
  group_list <- factor(c(rep('treatment',n_treatment),
                         rep("control",n_control)),
                       levels = c('treatment',"control"))
  limma_design <- model.matrix(~0+group_list)
  colnames(limma_design) <- levels(group_list)
  rownames(limma_design) <- colnames(DT_limma)
  cont.matrix <- limma::makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = limma_design)

  #limma
  fit <- limma::lmFit(DT_limma, limma_design)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2, trend=TRUE)

  result_limma <- limma::topTable(fit2, coef=1,n=Inf)
  result_limma=merge(DT[,c(meta_cols,treatment_samples, control_samples)],result_limma,by.x=0,by.y=0)
  result_limma$Row.names <- NULL
  return(result_limma)
  }
)










#' Plot a Volcano Plot for Differential Expression Results
#'
#' This function generates a volcano plot based on log fold changes (logFC) and
#' adjusted p-values from differential expression analysis. It highlights genes
#' that pass the specified thresholds for logFC and FDR (false discovery rate).
#' Optionally, it can label genes of interest or the top up/downregulated genes.
#'
#' @param DT.original A data frame containing the differential expression results.
#'        It should have columns `logFC` (log fold change) and `adj.P.Val` (adjusted p-value).
#' @param label_col The column name (as a string) used for labeling genes in the plot. Default is `NULL`.
#'        If `NULL`, the function will select the first column of `DT.original` for labeling.
#' @param lfc_threshold A numeric value representing the threshold for log fold change (default is 1).
#'        Genes with logFC greater than or equal to this value are labeled as "UP",
#'        and genes with logFC less than or equal to the negative of this value are labeled as "DOWN".
#' @param fdr_threshold A numeric value representing the false discovery rate (FDR) threshold (default is 0.01).
#'        Genes with an adjusted p-value greater than or equal to this threshold will be labeled as "Others".
#' @param labelgene A character vector of gene names to be labeled in the plot (default is `NULL`).
#'        If provided, only these genes will be labeled in the plot.
#'
#' @return A `ggplot2` object representing the volcano plot.
#' @export
#'
#' @examples
#' # Example data
#' dt <- data.frame(
#'   Gene = paste("Gene", 1:100),
#'   logFC = rnorm(100),
#'   adj.P.Val = runif(100)
#' )
#'
#' # Plot volcano plot for the top genes
#' plot_volcano(dt, label_col = "Gene", lfc_threshold = 1, fdr_threshold = 0.05)
#'
#' # Plot with specific gene labels
#' plot_volcano(dt, label_col = "Gene", lfc_threshold = 1, fdr_threshold = 0.05, labelgene = c("Gene1", "Gene2"))
plot_volcano <- function(DT.original, label_col = NULL, lfc_threshold=1, fdr_threshold=0.01, labelgene=NULL) {
  if(is.null(label_col)){
    label_col = names(DT.original)[1]
  }
  options(ggrepel.max.overlaps = Inf)
  DT <- DT.original

  # Set initial group to 'Others' and update based on thresholds
  DT <- DT %>%
    dplyr::mutate(Group = 'Others',
           Group = dplyr::if_else(logFC >= lfc_threshold, 'UP', Group),
           Group = dplyr::if_else(logFC <= -lfc_threshold, 'DOWN', Group),
           Group = dplyr::if_else(adj.P.Val >= fdr_threshold, 'Others', Group),
           labeltext = '')

  # If labelgene is provided, update labeltext accordingly
  if (!is.null(labelgene)) {
    DT <- DT %>%
      dplyr::mutate(labeltext = dplyr::if_else(!!rlang::sym(label_col) %in% label_gene, !!rlang::sym(label_col), labeltext))
  } else{
    # Select top 5 genes for UP and DOWN groups
    top5_gene <- DT %>%
      dplyr::filter(Group == 'UP') %>%
      dplyr::arrange(desc(logFC)) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::bind_rows(
        DT %>%
          dplyr::filter(Group == 'DOWN') %>%
          dplyr::arrange(logFC) %>%
          dplyr::slice_head(n = 5)
      )

    # Label top 5 genes using the specified label column
    DT <- DT %>%
      dplyr::mutate(labeltext = dplyr::if_else(!!rlang::sym(label_col) %in% top5_gene[[label_col]], !!rlang::sym(label_col), ''))
  }




  g <- ggplot2::ggplot(DT, ggplot2::aes(x = logFC, y = -log10(adj.P.Val))) +
    ggplot2::geom_point(ggplot2::aes(color = Group)) +
    ggplot2::scale_color_manual(breaks = c("DOWN", "Others", "UP"),
                       values = c("#67a9cf", "#969696", "#ef8a62")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom") +
    ggrepel::geom_label_repel(
      data = subset(DT),
      ggplot2::aes(label = labeltext),
      size = 5,
      box.padding = grid::unit(0.35, "lines"),
      point.padding = grid::unit(0.3, "lines")) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = lfc_threshold, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = -lfc_threshold, linetype = "dashed") +
    ggplot2::theme_classic()

  return(g)
}



