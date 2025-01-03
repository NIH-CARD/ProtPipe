library(testthat)
setwd("../..")

#spectronaut input
test_that("correctly make a prot_data object", {
  # Print the current working directory
  print(getwd())
  dat <- data.table::fread("EXAMPLES/DIFF_ABUNDANCE/iPSC.csv")
  dat_pro <<- create_protdata(dat)
  expect_s4_class(dat_pro, "ProtData")
})

#DIANN input
test_that("correctly make a prot_data object from DIANN input", {
  dat <- data.table::fread("~/OneDrive - National Institutes of Health/projects/SLAM/Heart/report.PG_matrix.tsv")
  dat_pro <<- create_protdata(dat)
  expect_s4_class(dat_pro, "ProtData")
})
