remove_outliers <-function(PD, sds = 3){
  dat <- getData(PD)
  cond <- getCondition(PD)
  dat.long <- getDataLong(PD)

  pgcounts <- data.table::as.data.table(table(dat.long$Sample))
  colnames(pgcounts) <- c("Sample", "N")

  stdev <- sd(pgcounts[,N])
  mean_count <- mean(pgcounts[,N])
  min_protein_groups <- floor(mean_count - (sds * stdev))
  max_protein_groups <- ceiling(mean_count + (sds * stdev))
  cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']'))
  low_count_samples <- as.character(pgcounts[N < min_protein_groups, Sample])
  high_count_samples <- as.character(pgcounts[N > max_protein_groups, Sample])
  if(length(low_count_samples)==0) {
    cat('\nINFO: No low group count samples to remove\n')
  } else {
    cat(paste0('\nINFO: runing low-count outlier ', low_count_samples))
    cat('\n\n')
    print(pgcounts[Sample %in% low_count_samples])
    cat('\n')
    dat[, colnames(dat) %in% low_count_samples]= NULL    # remove sample columns from wide table
    dat.long <- dat.long[! (Sample %in% low_count_samples)] # remove rows from long table
    cond[rownames(cond) %in% low_count_samples, ]= NULL
  }
  if(length(high_count_samples)==0) {
    cat('INFO: No high group count samples to remove\n')
  } else {
    cat(paste0('\nINFO: runing high-count outlier ', high_count_samples))
    cat('\n')
    print(pgcounts[Sample %in% high_count_samples])
    dat[, colnames(dat) %in% high_count_samples]= NULL     # remove sample columns from wide table
    dat.long <- dat.long[! (Sample %in% high_count_samples)] # remove rows from long table
    cond[rownames(cond) %in% high_count_samples, ]= NULL
  }
}
