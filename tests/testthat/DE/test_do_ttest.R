test_that("correctly perform ttest", {
  df <- readRDS("EXAMPLES/VIRUS/virus_data.rds")
  meta <- readRDS("EXAMPLES/VIRUS/virus_metadata.rds")

  dat_pro <- create_protdata(df, intensity_cols = c(3:147), meta)

  treatment_samples <- c("EBV_5_5d_1", "EBV_5_5d_1", "EBV_5_5d_1")
  control_samples <- c("EBV_0_5d_1", "EBV_0_5d_2", "EBV_0_5d_3")

  t <- do_t_test(dat_pro, treatment_samples = treatment_samples, control_samples = control_samples)
})


test_that("correctly perform ttest", {
  df <- readRDS("EXAMPLES/VIRUS/virus_data.rds")
  meta <- readRDS("EXAMPLES/VIRUS/virus_metadata.rds")

  dat_pro <- create_protdata(df, intensity_cols = c(4:148), meta)

  treatment_samples <- c("EBV_5_5d_1", "EBV_5_5d_1", "EBV_5_5d_1")
  control_samples <- c("EBV_0_5d_1", "EBV_0_5d_2", "EBV_0_5d_3")


  Log2_DP <- ProtPipe::log2_transform(dat_pro)

  DE <- do_limma(Log2_DP, treatment_samples = treatment_samples, control_samples = control_samples)
  plot_volcano(DE, label_col = "Genes")

  df_with_entrez <- add_entrez(DE)
  up_genes <- df_with_entrez %>%
    dplyr::filter(adj.P.Val < 0.01) %>%
    dplyr::filter(logFC > 2) %>%
    .$ENTREZID

  go_up <- ProtPipe::enrich_go(up_genes, df_with_entrez$ENTREZID)
  kegg_up <- ProtPipe::enrich_kegg(up_genes, df_with_entrez$ENTREZID)

  df_ordered <- df_with_entrez[order(df_with_entrez$logFC, decreasing = TRUE), ]
  ordered_genes <- df_ordered$logFC
  names(ordered_genes) <- df_ordered$ENTREZID
  ordered_genes_unique <- ordered_genes[!duplicated(names(ordered_genes))]


  gse_go <- ProtPipe::gse_go(ordered_genes_unique)
  gse_kegg <- ProtPipe::gse_kegg(ordered_genes_unique)

  en <- enrich_pathways(DE)

})
meta <- meta %>%
  mutate(sample = gsub("24hr", "1d", sample))%>%
  mutate(sample = gsub("48hr", "2d", sample))%>%
  mutate(time = gsub("24hr", "1d", time))%>%
  mutate(time = gsub("48hr", "2d", time))
