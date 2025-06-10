library(testthat)
setwd("../..")

#DIANN input with condition file
test_that("correctly make a prot_data object from DIANN input and condition file", {
  df <- data.table::fread("~/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/SLAM/Liver/report.pg_matrix.tsv")
  meta <- data.table::fread("~/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/SLAM/Liver/metadata.csv")
  dat_pro <- create_protdata(df, intensity_cols = c(5:138), meta)
  dat_pro_norm <- log2_transform(dat_pro)
  data_pro_scaled <- ProtPipe::scale(dat_pro_norm)
  data_pro_impute <- impute(data_pro_scaled, 0)
  ProtPipe::plot_PCs(data_pro_impute, condition = "Preperation batch")
  data_pro_corrected <- batch_correct(data_pro_impute,"Preparation batch")

  print(plot_PCs(data_pro_impute, condition = "Preparation batch"))
  print(plot_PCs(data_pro_corrected, condition = "Preparation batch"))

  print(ProtPipe::plot_umap(data_pro_corrected, condition = "Preparation batch"))

  print(plot_PCs(data_pro_impute, condition = "age"))
  print(plot_PCs(data_pro_corrected, condition = "age"))

  data_corrected <- limma::removeBatchEffect(data, dat_pro@condition$`Preparation batch`) %>% as.data.frame()
  dat_pro@data <- data_corrected
  print(plot_PCs(dat_pro, condition = "age"))

  expect_s4_class(dat_pro, "ProtData")
})
