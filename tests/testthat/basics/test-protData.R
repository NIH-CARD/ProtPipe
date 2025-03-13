library(testthat)
setwd("../../..")

#spectronaut input
test_that("correctly make a prot_data object", {
  # Print the current working directory
  print(getwd())
  dat <- data.table::fread("EXAMPLES/DIFF_ABUNDANCE/iPSC.csv")
  dat_pro <<- create_protdata(dat, intensity_cols = c(3:44))
  expect_s4_class(dat_pro, "ProtData")
})

#DIANN input
test_that("correctly make a prot_data object from DIANN input", {
  df <- data.table::fread("EXAMPLES/DIANN/report.pg_matrix.tsv")
  dat_pro <<- create_protdata(df, intensity_cols = c(5:138))
  expect_s4_class(dat_pro, "ProtData")
})

#DIANN input with condition file
test_that("correctly make a prot_data object from DIANN input and condition file", {
  df <- data.table::fread("EXAMPLES/DIANN/report.pg_matrix.tsv")
  meta <- data.table::fread("EXAMPLES/DIANN/metadata.csv")
  dat_pro <<- create_protdata(df, intensity_cols = c(5:138), meta)
  expect_s4_class(dat_pro, "ProtData")
})
