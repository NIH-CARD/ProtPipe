test_that("correctly make a prot_data object from DIANN input and condition file", {
  df <- data.table::fread("EXAMPLES/DIANN/report.pg_matrix.tsv")
  meta <- data.table::fread("EXAMPLES/DIANN/metadata.csv")
  dat_pro <- create_protdata(df, intensity_cols = c(5:138), meta)
  dat_pro_unique <- create_protdata(df, intensity_cols = c(5:138), meta) %>%
    unique_data(col = "Genes")
  expect_s4_class(dat_pro, "ProtData")
})

test_that("correctly make a prot_data object from DIANN peptide input and condition file", {
  df <- data.table::fread("EXAMPLES/DIANN/report.pr_matrix.tsv")
  meta <- data.table::fread("EXAMPLES/DIANN/metadata.csv")
  dat_pro <- create_protdata(df, intensity_cols = c(5:138), meta)
  dat_pro_unique <- create_protdata(df, intensity_cols = c(5:138), meta) %>%
    unique_data(col = "Genes")
  expect_s4_class(dat_pro, "ProtData")
})
