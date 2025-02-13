do_Limma <- function(protdata, baseCondition, control, treatment, lfc_threshold, fdr_threshold) {
  Log2_DT <- protdata@data
  Log2_DT$Genes <- protdata@prot_meta$PG.Genes
  treatment_samples <- which(protdata@condition$base_condition == treatment)
  control_samples <- which(protdata@condition$base_condition == control)
  
  if (length(treatment_samples) == 0 | length(control_samples) == 0) {
    stop("No samples found for treatment or control. Check baseCondition mapping.")
  }
  
  base::print(paste0("Selected samples - Treatment: ", paste(treatment_samples, collapse = ", ")))
  base::print(paste0("Selected samples - Control: ", paste(control_samples, collapse = ", ")))
  
  DT_limma <- Log2_DT[, c(treatment_samples, control_samples), drop = FALSE]
  DT_limma[is.na(DT_limma)] <- 0
  DT_limma$Genes <- Log2_DT$Genes
  
  DT_limma <- as.data.table(DT_limma)
  DT_limma[, missing_value := rowSums(.SD == 0), .SDcols = -c("Genes")]
  DT_limma <- DT_limma[missing_value < (length(treatment_samples) + length(control_samples)) * 0.5, ][, missing_value := NULL]
  
  group_list <- factor(c(rep("treatment", length(treatment_samples)), rep("control", length(control_samples))))
  limma_design <- model.matrix(~0 + group_list)
  colnames(limma_design) <- levels(group_list)
  
  rownames(limma_design) <- colnames(DT_limma[, -c("Genes")])
  
  cont.matrix <- makeContrasts(contrasts = paste0("treatment-control"), levels = limma_design)
  
  data_matrix <- as.matrix(DT_limma[, -c("Genes")])
  
  fit <- lmFit(data_matrix, limma_design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, trend = TRUE)
  
  result_limma <- topTable(fit2, coef = 1, n = Inf)
  result_limma <- as.data.frame(result_limma)
  
  result_limma$Genes <- DT_limma$Genes
  
  return(result_limma)
}



do_Ttest <- function(protdata, baseCondition, control, treatment, lfc_threshold, fdr_threshold) {
  Log2_DT <- protdata@data
  Log2_DT$Genes <- protdata@prot_meta$PG.Genes
  # Extract sample column names for treatment and control
  treatment_samples <- which(protdata@condition$base_condition == treatment)
  control_samples <- which(protdata@condition$base_condition == control)
  
  if (length(treatment_samples) == 0 | length(control_samples) == 0) {
    stop("No samples found for treatment or control. Check baseCondition mapping.")
  }
  
  base::print(paste0("Selected samples - Treatment: ", paste(treatment_samples, collapse = ", ")))
  base::print(paste0("Selected samples - Control: ", paste(control_samples, collapse = ", ")))
  
  DT_ttest <- Log2_DT[, c(treatment_samples, control_samples), drop = FALSE]
  DT_ttest[is.na(DT_ttest)] <- 0
  DT_ttest$Genes <- Log2_DT$Genes
  
  DT_ttest <- as.data.table(DT_ttest)
  DT_ttest[, missing_value := rowSums(.SD == 0), .SDcols = -c("Genes")]
  DT_ttest <- DT_ttest[missing_value < (length(treatment_samples) + length(control_samples)) * 0.5, ][, missing_value := NULL]
    treatment_data <- DT_ttest[, ..treatment_samples]
  control_data <- DT_ttest[, ..control_samples]
  
  ttest_results <- DT_ttest[, .(
    Genes = Genes,
    p_value = t.test(as.numeric(unlist(.SD[, ..treatment_samples])), 
                     as.numeric(unlist(.SD[, ..control_samples])))$p.value,
    mean_diff = mean(as.numeric(unlist(.SD[, ..treatment_samples]))) - 
      mean(as.numeric(unlist(.SD[, ..control_samples])))
  ), by = .(Genes)]
  
  ttest_results[, adj_p_value := p.adjust(p_value, method = "fdr")]
  
  ttest_results <- ttest_results[
    abs(mean_diff) >= lfc_threshold & adj_p_value <= fdr_threshold
  ]
  
  return(ttest_results)
}


plot_volcano <- function(DT.original, 
                         out_dir = "./", 
                         output_filename = "volcano_plot", 
                         label_col = "Genes", 
                         lfc_threshold = 1, 
                         fdr_threshold = 0.05, 
                         labelgene = NULL) {
  DT <- as.data.table(DT.original)
  
  required_cols <- c("logFC", "adj.P.Val", label_col)
  if (!all(required_cols %in% colnames(DT))) {
    stop("The input data must contain the columns: logFC, adj.P.Val, and the specified label_col.")
  }
  
  DT[, Group := "Others"]
  DT[logFC >= lfc_threshold, Group := "UP"]
  DT[logFC <= -lfc_threshold, Group := "DOWN"]
  DT[adj.P.Val >= fdr_threshold, Group := "Others"]
  
  DT[, labeltext := ""]
  
  if (!is.null(labelgene)) {
    label_gene <- fread(labelgene)  # Assuming labelgene is a file path
    DT[label_col %in% label_gene$gene, labeltext := get(label_col)]
  } else {
    top5_gene <- DT[Group == "UP"][order(-logFC)][1:5]
    top5_gene <- rbind(top5_gene, DT[Group == "DOWN"][order(logFC)][1:5])
    DT[get(label_col) %in% top5_gene[[label_col]], labeltext := get(label_col)]
  }
  
  # Create the volcano plot
  g <- ggplot(DT, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Group), alpha = 0.8) +
    scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                       values = c("#67a9cf", "#969696", "#ef8a62")) +
    theme_minimal(base_size = 14) +
    labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", color = "Group") +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "black") + 
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black") +
    geom_label_repel(
      aes(label = labeltext),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      max.overlaps = Inf
    ) +
    theme(legend.position = "bottom")
  
  # Save the plot
  output_path <- paste0(out_dir, output_filename, "_volcano.pdf")
  ggsave(g, filename = output_path, width = 8, height = 8)
  
  cat(paste0("Volcano plot saved to: ", output_path, "\n"))
  
  return(g)  # ggplot object for add ons and customs
}
