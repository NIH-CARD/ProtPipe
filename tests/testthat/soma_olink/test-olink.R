library(testthat)
setwd("../../..")

load("EXAMPLES/olink/npx_data1.rda")
file = "EXAMPLES/olink/npx_data1.xlsx"
metafile <- "EXAMPLES/olink/npx_data1_meta_original.csv"

npx <- OlinkAnalyze::read_NPX(file)

meta <- read.delim(metafile, sep = ";")

test_that("correctly make a prot_data object from Olink without LOD filtering", {
  dat_pro <- create_protdata_from_olink(npx, filter = F)
  expect_s4_class(dat_pro, "ProtData")
  cols = ncol(dat_pro@data)
  expect_equal(cols, 158)
  rows = nrow(dat_pro@data)
  expect_equal(rows, 184)
})

test_that("correctly make a prot_data object from Olink with LOD filtering", {
  dat_pro <- create_protdata_from_olink(npx)
  expect_s4_class(dat_pro, "ProtData")
  cols = ncol(dat_pro@data)
  expect_equal(cols, 158)
  rows = nrow(dat_pro@data)
  expect_equal(rows, 184)
})
