test_that("correctly perform ttest", {
  df <- data.table::fread("EXAMPLES/DIANN/report.pg_matrix.tsv")
  meta <- data.table::fread("EXAMPLES/DIANN/metadata.csv")
  dat_pro <- create_protdata(df, intensity_cols = c(5:138), meta)

  treatment_samples <- c("SLAM_1", "SLAM_2", "SLAM_3", "SLAM_4", "SLAM_5")
  control_samples <- c("SLAM_6", "SLAM_7", "SLAM_8", "SLAM_9", "SLAM_10")

  t <- do_t_test(dat_pro, treatment_samples = treatment_samples, control_samples = control_samples)
})


test_that("correctly perform ttest", {
  df <- data.table::fread("EXAMPLES/DIANN/report.pg_matrix.tsv")
  meta <- data.table::fread("EXAMPLES/DIANN/metadata.csv")
  dat_pro <- create_protdata(df, intensity_cols = c(5:138), meta)

  Log2_DP <- ProtPipe::log2_transform(dat_pro)

  treatment_samples <- c("SLAM_1", "SLAM_2", "SLAM_3", "SLAM_4", "SLAM_5")
  control_samples <- c("SLAM_6", "SLAM_7", "SLAM_8", "SLAM_9", "SLAM_10")
  DE <- do_limma(Log2_DP, treatment_samples = treatment_samples, control_samples = control_samples)
  plot_volcano(DE)
})
