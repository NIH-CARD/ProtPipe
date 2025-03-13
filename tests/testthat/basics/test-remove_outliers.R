library(testthat)
setwd("../../..")

#spectronaut input
test_that("correctly remove outliers from a prot_data object", {
  df <- data.table::fread("EXAMPLES/DIANN/report.pg_matrix.tsv")
  meta <- data.table::fread("EXAMPLES/DIANN/metadata.csv")
  dat_pro1 <- create_protdata(df, intensity_cols = c(5:138), meta)
  dat_pro2 <- remove_outliers(dat_pro1)
  expect_s4_class(dat_pro2, "ProtData")
})
