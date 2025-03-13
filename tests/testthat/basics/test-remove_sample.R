library(testthat)
setwd("../../..")

test_that("remove a sample from a prot_data object", {
  dat <- data.table::fread("EXAMPLES/DIFF_ABUNDANCE/iPSC.csv")
  dat_pro <- create_protdata(dat, intensity_cols = c(3:44))
  a <- num_samples(dat_pro)
  dat_pro <- remove_sample(dat_pro, "Day0_1")
  b <- num_samples(dat_pro)
  print(a-b)
  expect_equal(a-b, 1)
})

test_that("remove 3 samples from a prot_data object", {
  dat <- data.table::fread("EXAMPLES/DIFF_ABUNDANCE/iPSC.csv")
  dat_pro <- create_protdata(dat, intensity_cols = c(3:44))
  a <- num_samples(dat_pro)
  dat_pro1 <- remove_sample(dat_pro, c("Day0_1", "Day0_2", "Day0_3"))
  b <- num_samples(dat_pro1)
  print(a-b)
  expect_equal(a-b, 3)
})

test_that("remove 0 samples from a prot_data object", {
  dat <- data.table::fread("EXAMPLES/DIFF_ABUNDANCE/iPSC.csv")
  dat_pro <- create_protdata(dat, intensity_cols = c(3:44))
  a <- num_samples(dat_pro)
  dat_pro1 <- remove_sample(dat_pro, c())
  b <- num_samples(dat_pro1)
  print(a-b)
  expect_equal(a-b, 0)
})
