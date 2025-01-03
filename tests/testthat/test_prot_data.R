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
  df <- data.table::fread("~/OneDrive - National Institutes of Health/projects/SLAM/Heart/report.PG_matrix.tsv")
  dat_pro <<- create_protdata(df)
  expect_s4_class(dat_pro, "ProtData")
})

#DIANN input with condition file
test_that("correctly make a prot_data object from DIANN input and condition file", {
  df <- data.table::fread("~/OneDrive - National Institutes of Health/projects/SLAM/Heart/report.PG_matrix.tsv")
  meta <- data.table::fread("~/OneDrive - National Institutes of Health/projects/SLAM/Heart/metadata.csv")
  dat_pro <<- create_protdata(df)
  expect_s4_class(dat_pro, "ProtData")

  df1 <- data.table::fread("~/OneDrive - National Institutes of Health/projects/SLAM/Muscle/report.PG_matrix.tsv")
  meta1 <- data.table::fread("~/OneDrive - National Institutes of Health/projects/SLAM/Muscle/metadata.csv")
  dat_pro1 <<- create_protdata(df1)
  expect_s4_class(dat_pro1, "ProtData")
})
