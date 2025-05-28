#' @importFrom magrittr %>%

#' Title
#'
#' @param object
#' @param sds
#'
#' @return
#' @export
#'
#' @examples
setGeneric("remove_outliers", function(object, sds = 3) standardGeneric("remove_outliers"))

#' Title
#'
#' @param ProtData
#'
#' @return
#' @export
#'
#' @examples
setMethod("remove_outliers",
          "ProtData",
          function(object, sds = 3){
            N_values <- colSums(!is.na(getData(object)))
            pgcounts <- data.frame(Sample = names(N_values), N = N_values)

            dat <- object@data
            cond <- object@condition

            stdev <- sd(pgcounts[,'N'])
            mean_count <- mean(pgcounts[,'N'])
            min_protein_groups <- floor(mean_count - (sds * stdev))
            max_protein_groups <- ceiling(mean_count + (sds * stdev))
            cat(paste0('INFO: Tolerating protein group counts in the range [', min_protein_groups,',',max_protein_groups,']\n'))
            low_count_samples <- as.character(pgcounts[pgcounts$N < min_protein_groups, 'Sample'])
            if(length(low_count_samples)>0){
              print("removing the following samples with low protein counts:")
              print(low_count_samples)
            }
            high_count_samples <- as.character(pgcounts[pgcounts$N > max_protein_groups, 'Sample'])
            if(length(high_count_samples)>0){
              print("removing the following samples with high protein counts:")
              print(high_count_samples)
            }

            outliers <- c(low_count_samples, high_count_samples)
            if (length(outliers) >0){
              dat <- dat[,!(colnames(dat) %in% outliers),drop = FALSE]
              cond <- cond[!rownames(cond) %in% outliers, ,drop = FALSE]

              object@data <- dat
              object@condition <- cond
            }
            return(object)
          }
)
