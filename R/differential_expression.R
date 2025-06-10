#' @importFrom magrittr %>%

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
  # sort by adjusted o value and log fold change
  result_limma_sorted <- result_limma[order(result_limma$adj.P.Val, -abs(result_limma$logFC)), ]
  return(result_limma_sorted)
  }
)


#' Perform DEA using the condition labels of the protdata object
#'
#' @param object
#' @param condition
#' @param control_grouop
#' @param treatment_group
#'
#' @return
#' @export
#'
#' @examples
setGeneric("do_limma_by_condition", function(object, condition, control_group, treatment_group) standardGeneric("do_limma_by_condition"))

#' @describeIn do_limma Method for ProtData objects
setMethod("do_limma_by_condition", "ProtData", function(object, condition, control_group, treatment_group) {
  meta <- object@condition
  conditions <- names(meta)

  if (!(condition %in% conditions)) {
    stop("The 'condition' must be a column name of the ProtData condition metadata.")
  }

  groups <- unique(meta[[condition]])
  if (!(control_group %in% groups && treatment_group %in% groups)) {
    stop("Both control and treatment groups must be valid entries in the condition metadata.")
  }

  control_samples <- rownames(meta %>% dplyr::filter(.data[[condition]] == control_group))
  treatment_samples <- rownames(meta %>% dplyr::filter(.data[[condition]] == treatment_group))

  if (length(control_samples) < 2 || length(treatment_samples) < 2) {
    stop("Each of the control and treatment groups must contain at least 2 samples.")
  }

  return(ProtPipe::do_limma(object, treatment_samples, control_samples))
})








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
      dplyr::mutate(labeltext = dplyr::if_else(!!rlang::sym(label_col) %in% labelgene, !!rlang::sym(label_col), labeltext))
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
      data = DT,
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

# Enrichment Functions
#' Title
#'
#' @param DE
#' @param org
#'
#' @return
#' @export
#'
#' @examples
# add_entrez <- function(DE, org = org.Hs.eg.db::org.Hs.eg.db, gene_col = "Genes"){
#   DE[[gene_col]] <- sapply(strsplit(DE[[gene_col]], ";"), `[`, 1)
#   dplyr::left_join(
#     DE,
#     clusterProfiler::bitr(
#       DE$Genes,
#       fromType = "SYMBOL",
#       toType = "ENTREZID",
#       OrgDb = org
#     ),
#     by = c("Genes" = "SYMBOL")
#   )
# }
add_entrez <- function(DE, org = org.Hs.eg.db::org.Hs.eg.db, gene_col = "Genes") {
  DE[[gene_col]] <- sapply(strsplit(DE[[gene_col]], ";"), `[`, 1)

  # Remove NA or empty strings before mapping
  genes_to_map <- DE[[gene_col]]
  genes_to_map <- genes_to_map[!is.na(genes_to_map) & genes_to_map != ""]

  if (length(genes_to_map) == 0) {
    warning("No valid gene symbols found to map.")
    return(NULL)
  }

  # Use tryCatch to prevent app crash
  mapping <- tryCatch({
    clusterProfiler::bitr(
      genes_to_map,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org
    )
  }, error = function(e) {
    warning("Gene mapping failed: ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(mapping) || nrow(mapping) == 0) {
    warning("No Entrez IDs mapped.")
    return(NULL)
  }

  dplyr::inner_join(
    DE,
    mapping,
    by = setNames("SYMBOL", gene_col)
  )
}


#' Title
#'
#' @param gene_id
#' @param all_gene_vector
#' @param enrich_pvalue
#'
#' @return
#' @export
#'
#' @examples
enrich_go <- function(gene_id, all_gene_vector, enrich_pvalue = 1, org = org.Hs.eg.db::org.Hs.eg.db) {
  cat("Processing GO\n")

  GO <- tryCatch({
    clusterProfiler::enrichGO(
      gene          = gene_id,
      universe      = names(all_gene_vector),
      OrgDb         = org,
      ont           = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = enrich_pvalue,
      qvalueCutoff  = enrich_pvalue,
      readable      = TRUE
    )
  }, error = function(e) {
    warning("GO enrichment failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(GO) && !is.null(GO@result) && nrow(GO@result) > 0) {
    return(GO)
  } else {
    return(NULL)
  }
}


#' Title
#'
#' @param gene_id
#' @param all_gene_vector
#' @param enrich_pvalue
#'
#' @return
#' @export
#'
#' @examples
enrich_kegg <- function(gene_id, all_gene_vector, enrich_pvalue = 1, org = org.Hs.eg.db::org.Hs.eg.db, organism = 'hsa') {
  cat("Processing KEGG\n")

  KEGG <- tryCatch({
    clusterProfiler::enrichKEGG(
      gene         = gene_id,
      organism     = organism,
      universe     = all_gene_vector,
      pvalueCutoff = enrich_pvalue,
      qvalueCutoff = enrich_pvalue
    )
  }, error = function(e) {
    warning("KEGG enrichment failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(KEGG)) {
    KEGG <- DOSE::setReadable(KEGG, OrgDb = org, keyType = "ENTREZID")
    if (!is.null(KEGG@result) && nrow(KEGG@result) > 0) {
      return(KEGG)
    }
  }

  return(NULL)
}

#' Title
#'
#' @param gene_list
#' @param enrich_pvalue
#' @param org
#'
#' @return
#' @export
#'
#' @examples
gse_go <- function(gene_list, enrich_pvalue = 1, org = org.Hs.eg.db::org.Hs.eg.db) {
  cat("Processing GSEA GO\n")

  GO <- tryCatch({
    clusterProfiler::gseGO(
      geneList     = gene_list,
      OrgDb        = org,
      ont          = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff  = enrich_pvalue,
      verbose       = FALSE
    )
  }, error = function(e) {
    warning("GSEA GO failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(GO)) {
    GO <- DOSE::setReadable(GO, OrgDb = org, keyType = "ENTREZID")
    if (!is.null(GO@result) && nrow(GO@result) > 0) {
      return(GO)
    }
  }

  return(NULL)
}

#' Title
#'
#' @param gene_list
#' @param enrich_pvalue
#' @param org
#' @param organism
#'
#' @return
#' @export
#'
#' @examples
gse_kegg <- function(gene_list, enrich_pvalue = 1, org = org.Hs.eg.db::org.Hs.eg.db, organism = 'hsa') {
  cat("Processing GSEA KEGG\n")

  KEGG <- tryCatch({
    clusterProfiler::gseKEGG(
      geneList     = gene_list,
      organism     = organism,
      pvalueCutoff = enrich_pvalue,
      verbose      = FALSE
    )
  }, error = function(e) {
    warning("GSEA KEGG failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(KEGG)) {
    KEGG <- DOSE::setReadable(KEGG, OrgDb = org, keyType = "ENTREZID")
    if (!is.null(KEGG@result) && nrow(KEGG@result) > 0) {
      return(KEGG)
    }
  }

  return(NULL)
}

do_enrichment <- function(DE, lfc_threshold, fdr_threshold, enrich_pvalue){

}
###pathway analysis
#' Title
#'
#' @param DE
#' @param lfc_threshold
#' @param fdr_threshold
#' @param enrich_pvalue
#'
#' @return
#' @export
#'
#' @examples
enrich_pathways = function(DE, lfc_threshold=1, fdr_threshold=0.01, enrich_pvalue=0.05, go_org = org.Hs.eg.db, kegg_org = 'hsa', gene_col = "Genes"){
  datas <- list()
  plots <- list()

  DT <- ProtPipe::add_entrez(DE, org = go_org, gene_col = gene_col)

  # Check if any genes were successfully mapped
  if (is.null(DT)) {
    warning("No genes were mapped to Entrez IDs. Skipping enrichment.")
    return(NULL)
  }

  #initialize all plots and dataframes
  datas$go_up <- datas$go_down <- datas$kegg_up <- datas$kegg_down <- datas$gse_go <- datas$gse_kegg <- NULL
  plots$go_up_dotplot <- plots$go_up_barplot <- plots$go_down_dotplot <- plots$go_down_barplot <- NULL
  plots$kegg_up_dotplot <- plots$kegg_up_barplot <- plots$kegg_down_dotplot <- plots$kegg_down_barplot <- NULL
  plots[["gse_go_dotplot"]] <- plots[["gse_go_emapplot"]] <- plots[["gse_kegg_dotplot"]] <- plots[["gse_kegg_emapplot"]] <- NULL
  expected_ontologies <- c("BP", "MF", "CC")
  for (ont in expected_ontologies) {
    plots[[paste0("go_up_dotplot_", tolower(ont))]] <- NULL
    plots[[paste0("go_up_barplot_", tolower(ont))]] <- NULL
    plots[[paste0("go_down_dotplot_", tolower(ont))]] <- NULL
    plots[[paste0("go_down_barplot_", tolower(ont))]] <- NULL
  }


  ## up and down regulated genes
  up_genes=DT[which(DT$logFC>=lfc_threshold&DT$adj.P.Val<=fdr_threshold),]
  down_genes=DT[which(DT$logFC<=(-lfc_threshold)&DT$adj.P.Val<=fdr_threshold),]

  # over representation enrichment for upregulated genes ###############################
  if (length(up_genes)>0){
    go_up <- ProtPipe::enrich_go(up_genes$ENTREZID, DT$ENTREZID, org = go_org, enrich_pvalue = enrich_pvalue)
    kegg_up <- ProtPipe::enrich_kegg(up_genes$ENTREZID, DT$ENTREZID, org = go_org, organism = kegg_org, enrich_pvalue = enrich_pvalue)

    if(!is.null(go_up)){
      datas$go_up <- go_up@result

      p <- enrichplot::dotplot(go_up, showCategory = 10, split = "ONTOLOGY") +
        ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")
      plots$go_up_dotplot <- p

      g <- barplot(go_up, showCategory = 10, split = "ONTOLOGY") +
        ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")
      plots$go_up_barplot <- g

      #individual GO ontolgies
      ontologies <- unique(go_up@result$ONTOLOGY)

      # Loop over ontologies and generate separate plots
      for (ont in ontologies) {
        # Subset the go_up object by ontology
        go_up_subset <- clusterProfiler::filter(go_up, ONTOLOGY == ont)
        if(nrow(go_up_subset) > 0){
          # Create dotplot
          dp <- enrichplot::dotplot(go_up_subset, showCategory = 10) +
            ggplot2::ggtitle(paste("GO Dotplot -", ont))
          plots[[paste0("go_up_dotplot_", tolower(ont))]] <- dp

          # Create barplot
          bp <- barplot(go_up_subset, showCategory = 10) +
            ggplot2::ggtitle(paste("GO Barplot -", ont))
          plots[[paste0("go_up_barplot_", tolower(ont))]] <- bp
        }
      }
    }
    if(!is.null(kegg_up)){
      datas$kegg_up <- kegg_up@result

      p <- enrichplot::dotplot(kegg_up, showCategory = 10)
      plots$kegg_up_dotplot <- p

      g <- barplot(kegg_up, showCategory = 10)
      plots$kegg_up_barplot <- g
    }
  }
    # over representation enrichment for downregulated genes ###############################
    if (length(down_genes)>0){
      go_down <- ProtPipe::enrich_go(down_genes$ENTREZID, DT$ENTREZID, org = go_org, enrich_pvalue = enrich_pvalue)
      kegg_down <- ProtPipe::enrich_kegg(down_genes$ENTREZID, DT$ENTREZID, org = go_org, organism = kegg_org, enrich_pvalue = enrich_pvalue)

      # check if enriched pathways exist
      if(!is.null(go_down)){
        datas$go_down <- go_down@result

        p <- enrichplot::dotplot(go_down, showCategory = 10, split = "ONTOLOGY") +
          ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")
        plots$go_down_dotplot <- p
        g <- barplot(go_down, showCategory = 10, split = "ONTOLOGY") +
          ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")
        plots$go_down_barplot <- g

        #individual GO ontolgies
        ontologies <- unique(go_down@result$ONTOLOGY)

        # Loop over ontologies and generate separate plots
        for (ont in ontologies) {
          # Subset the go_up object by ontology
          go_down_subset <- clusterProfiler::filter(go_down, ONTOLOGY == ont)
          if(nrow(go_down_subset) > 0){
            # Create dotplot
            dp <- enrichplot::dotplot(go_down_subset, showCategory = 10) +
              ggplot2::ggtitle(paste("GO Dotplot -", ont))
            plots[[paste0("go_down_dotplot_", tolower(ont))]] <- dp

            # Create barplot
            bp <- barplot(go_down_subset, showCategory = 10) +
              ggplot2::ggtitle(paste("GO Barplot -", ont))
            plots[[paste0("go_down_barplot_", tolower(ont))]] <- bp
          }

        }
        }
        if(!is.null(kegg_down)){
          datas$kegg_down <- kegg_down@result

          p <- enrichplot::dotplot(kegg_down, showCategory = 10)
          plots$kegg_down_dotplot <- p

          g <- barplot(kegg_down, showCategory = 10)
          plots$kegg_down_barplot <- g
        }
      }
  #Gene Set Enrichment ##############
  df_ordered <- DT[order(DT$logFC, decreasing = TRUE), ]
  ordered_genes <- df_ordered$logFC
  names(ordered_genes) <- df_ordered$ENTREZID
  ordered_genes_unique <- ordered_genes[!duplicated(names(ordered_genes))]

  gse_go <- ProtPipe::gse_go(ordered_genes_unique, org = go_org, enrich_pvalue = enrich_pvalue)
  gse_kegg <- ProtPipe::gse_kegg(ordered_genes_unique, org = go_org, organism = kegg_org, enrich_pvalue = enrich_pvalue)

  if(!is.null(gse_go)){
    datas$gse_go <- gse_go@result
    if (nrow(gse_go) > 0) {
      p=enrichplot::dotplot(gse_go, showCategory=10, split=".sign") + facet_grid(.~.sign)
      plots[["gse_go_dotplot"]] <- p
      x2 <- enrichplot::pairwise_termsim(gse_go)
      p=enrichplot::emapplot(x2)
      plots[["gse_go_emapplot"]] <- p
    }

  }
  if(!is.null(gse_kegg)){
    datas$gse_kegg <- gse_kegg@result
    if (nrow(gse_kegg) > 0) {
      p=enrichplot::dotplot(gse_kegg, showCategory=10, split=".sign") + facet_grid(.~.sign)
      plots[["gse_kegg_dotplot"]] <- p
      x2 <- enrichplot::pairwise_termsim(gse_kegg)
      p=enrichplot::emapplot(x2)
      plots[["gse_kegg_emapplot"]] <- p
    }
  }
  return(list(results = datas, plots = plots))
}







