test_that("correctly perform ttest", {
  df <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")

  dat_pro <- create_protdata(df)

  normalize_median <- ProtPipe::median_normalize(dat_pro)
  normalize_mean <- ProtPipe::mean_normalize(dat_pro)

  ProtPipe::plot_CVs(normalize_mean, condition = "base_condition")


  ProtPipe::plot_pg_counts(normalize_median)
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
