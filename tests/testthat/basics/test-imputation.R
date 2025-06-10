test_that("correctly perform ttest", {
  df <- readRDS("EXAMPLES/VIRUS/virus_data.rds")
  meta <- readRDS("EXAMPLES/VIRUS/virus_metadata.rds")

  dat_pro <- create_protdata(df, condition = meta)

  impute_zero <- ProtPipe::impute(dat_pro, 0)
  impute_min <- ProtPipe::impute_min(dat_pro)
  impute_minimal <- ProtPipe::impute_minimal(dat_pro)

  ProtPipe::plot_pg_counts(dat_pro)
  ProtPipe::plot_pg_counts(dat_pro, condition = "viral.exposure")
  ProtPipe::plot_CVs(dat_pro, condition = "viral.exposure")

  ProtPipe::plot_pg_intensities(dat_pro)
  ProtPipe::plot_correlation_heatmap(dat_pro)
  ProtPipe::plot_correlation_heatmap(dat_pro, condition = "time")

  ProtPipe::plot_PCs(dat_pro, condition = "viral.exposure")
  ProtPipe::plot_umap(dat_pro, condition = "concentration")
  ProtPipe::plot_hierarchical_cluster(dat_pro)

  dat_pro <- ProtPipe::log2_transform(dat_pro)
  ProtPipe::plot_proteomics_heatmap(dat_pro)
})
