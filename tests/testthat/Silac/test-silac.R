library(testthat)
setwd("../../..")

dat <<- data.table::fread("EXAMPLES/SILAC/pSILAC_CARD_Challenge_Report.tsv")

#test silac constructor
test_that("correctly make a SilacData object from tsv input", {
  dfas <<- create_silac(dat, intensity_cols = c(5:ncol(dat)))
})
