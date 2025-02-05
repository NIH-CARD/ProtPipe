library(testthat)
setwd("../../..")

#spectronaut input
test_that("correctly remove outliers from a prot_data object", {
  df1 <- data.table::fread("~/OneDrive - National Institutes of Health/projects/SLAM/Muscle/report.PG_matrix.tsv")
  meta1 <- data.table::fread("~/OneDrive - National Institutes of Health/projects/SLAM/Muscle/metadata.csv")
  dat_pro1 <- create_protdata(df1)
  dat_pro2 <- remove_outliers(dat_pro1)
  expect_s4_class(dat_pro1, "ProtData")
})
