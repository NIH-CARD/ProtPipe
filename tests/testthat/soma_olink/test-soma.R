library(testthat)
setwd("../../..")

dat <- SomaDataIO::read_adat("EXAMPLES/soma/example_data_v5.0_plasma.adat")

test_that("create protdata object from somascan data without filtering", {
  soma_pro <- create_protdata_from_soma(dat, filter = FALSE)
  expect_s4_class(soma_pro, "ProtData")
})

test_that("create protdata object from somascan data with buffer filtering", {
  soma_pro <- create_protdata_from_soma(dat)
  expect_s4_class(soma_pro, "ProtData")
})

test_that("make qc figures from soma data", {
  dat <- SomaDataIO::read_adat("EXAMPLES/soma/example_data_v5.0_plasma.adat")
  dat_csf <- SomaDataIO::read_adat("/Users/epsteinjm/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/projects/CSF_soma/data/CHI-24-010_v5.0_CSF.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.AddLOD.anmlSMP.20241010.adat")
  soma_pro <- create_protdata_from_soma(dat_csf)
  expect_s4_class(soma_pro, "ProtData")
  print(plot_pg_counts(soma_pro))
})
