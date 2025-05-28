#' @importFrom magrittr %>%

#' Title
#'
#' @param PD A prot_data object
#' @param by_condition True if you want to group by condition
#'
#' @return
#' @export
#'
#' @examples
plot_pg_counts <- function(PD, condition = NULL) {

  # get the number of protein groups per sample
  N_values <- colSums(!is.na(getData(PD)))
  pgcounts <- data.frame(Sample = names(N_values), N = N_values)

  # Order samples by ascending counts
  n_samples <- nrow(pgcounts)
  if (is.null(condition)){
    if (n_samples > 20) {
      p=ggplot2::ggplot(pgcounts, ggplot2::aes(x=Sample, y=N)) +
        ggplot2::geom_bar(stat="identity", fill="#67a9cf")+
        ggplot2::theme_classic()+
        ggplot2::labs(fill = "",x="",y='Number of Protein Groups')+
        ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))
    }
    if (n_samples <= 20) {
      p=ggplot2::ggplot(pgcounts, ggplot2::aes(x=Sample, y=N,fill=Condition)) +
        ggplot2::geom_bar(stat="identity")+
        ggplot2::theme_classic()+
        ggplot2::labs(fill = "",x="",y='Number of Protein Groups')+
        ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 90))+
        ggplot2::geom_text(ggplot2::aes(label=N, y=N + (0.05*max(pgcounts$N))))
    }
  }else{
  # group by condition
    condition_file <- getCondition(PD)
    if (condition %in% colnames(condition_file)){
      condition_file <- condition_file %>%
        dplyr::mutate(Sample = rownames(condition_file)) %>%
        dplyr::select(c(!!rlang::sym(condition), Sample))
      pgcounts <- pgcounts %>%
        dplyr::left_join(condition_file, by = "Sample") %>%
        dplyr::rename(Condition = !!rlang::sym(condition))  # Rename the dynamic column to 'Condition'
      summary_data <- pgcounts %>%
        dplyr::group_by(Condition) %>%
        dplyr::summarize(mean = mean(N), sd = sd(N)) %>%
        dplyr::arrange(Condition)
      p=ggplot2::ggplot(summary_data, ggplot2::aes(x=as.factor(Condition), y=mean)) +
        ggplot2::geom_bar(stat="identity",fill="#67a9cf", position= ggplot2::position_dodge())+
        ggplot2::theme_classic()+
        ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                               position=ggplot2::position_dodge(.9))+
        ggplot2::labs(fill = "",x="",y='Number of Protein Groups')
    }else{
      stop("the selected condition does not appear in the condition file")
    }
  }
  return(p)
}

#' Title
#'
#' @param PD prot_data object
#'
#' @return ggplot boxplot of intensities for each sample
#' @export
#'
#' @examples
plot_pg_intensities <- function(PD) {
  # Assuming PD@data is a data frame with proteins as rows and samples as columns
  dat <- PD@data  # Or however you access your data frame in the wide format
  dat_long <- reshape2::melt(dat,
                            measure.vars=names(dat)[sapply(dat, function(x) all(is.numeric(x)))],
                            variable.name='Sample',
                            value.name='Intensity')
  dat_long <- dat_long[dat_long$Intensity>0,]
  dat_long <<- dat_long[rowSums(is.na(dat_long)) < ncol(dat_long),]
  # Plot the boxplot
  g <- ggplot2::ggplot(dat_long, ggplot2::aes(x = Sample, y = log10(Intensity))) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#67a9cf") +
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "", x = "", y = "Log10 Protein Intensity") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::geom_boxplot(width = 0.1) +
    ggplot2::geom_hline(color = '#ef8a62', linetype = 'dashed',
                        ggplot2::aes(yintercept = quantile(log10(Intensity), 0.50, na.rm = TRUE)))

  return(g)
}

#' Title
#'
#' @param PD
#' @param condition
#' @param min_samples
#'
#' @return
#' @export
#'
#' @examples
plot_CVs <- function(PD, condition, min_samples = 2){
  condition_file <- getCondition(PD)
  if (!condition %in% colnames(condition_file)){
    stop("the selected condition does not appear in the condition file")
  }

  intensities_t <- as.data.frame(t(PD@data))
  intensities_t[[condition]] <- condition_file[[condition]]

  # Split by condition
  groups <- split(intensities_t, intensities_t[[condition]])

  # Compute CVs per protein per group
  cv_list <- lapply(names(groups), function(cond) {
    group_data <- groups[[cond]][, !(colnames(groups[[cond]]) %in% condition), drop = FALSE]

    if (nrow(group_data) < min_samples) {
      return(NULL)  # Not enough samples
    }

    # CV per protein across samples in this group
    cv_values <- apply(group_data, 2, function(x) {
      if (all(is.na(x))) return(NA)
      stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
    })

    data.frame(
      Protein = names(cv_values),
      CV = cv_values,
      Condition = cond,
      stringsAsFactors = FALSE
    )
  })

  # Combine and filter
  cv_df <- do.call(rbind, cv_list)
  cv_df <- cv_df[!is.na(cv_df$CV), ]

  # Plot
  p <- ggplot2::ggplot(cv_df, ggplot2::aes(x = Condition, y = CV)) +
    ggplot2::geom_violin(fill = "#67a9cf", color = "black", trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.3, size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "", y = "Coefficient of Variation (CV)") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = stats::quantile(cv_df$CV, 0.5, na.rm = TRUE)),
                        color = "#ef8a62", linetype = "dashed")

  return(p)
}

#' Title
#'
#' @param PD
#' @param method
#'
#' @return
#' @export
#'
#' @examples
get_spearman <- function(PD, method = 'spearman') {
  DT <- getData(PD)
  #### Pairwise correlations between sample columns
  #dt.samples <- DT[,-c(1:2)]     # Ignore info columns (subset to only intensity values)
  dt.samples <- DT[, sapply(DT, is.numeric)] #better way of getting just intensity columns

  dt.corrs <- cor(log2(as.matrix(na.omit(dt.samples))+1), method=method)

  # Format correlations as 3 digits
  dt.corrs <- data.table::data.table(reshape2::melt(dt.corrs, measure.vars=dt.corrs[,rn], value.name='Spearman'))
  dt.corrs <- dt.corrs[! is.na('Spearman')]
  data.table::setnames(dt.corrs, c('Var1', 'Var2'), c('SampleA','SampleB'))
  dt.corrs <- dt.corrs %>% dplyr::mutate(Spearman = round(Spearman, 3))

  return(dt.corrs[])
}

## correlation
#' Title
#'
#' @param DT.corrs
#'
#' @return
#' @export
#'
#' @examples
plot_correlation_heatmap <- function(PD, condition = NULL) {
  DT.corrs <- get_spearman(PD)
  n_samples <- length(unique(DT.corrs[,'SampleA']))
  max_limit <- max(DT.corrs$Spearman)
  min_limit <- min(DT.corrs$Spearman)
  mid_limit <- as.numeric(format(((max_limit + min_limit) / 2), digits=3))
  if (is.null(condition)){
    g <- ggplot2::ggplot(DT.corrs, ggplot2::aes(x=SampleA, y=SampleB, fill=Spearman, label = NA)) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(color='gray10') +
      ggplot2::theme_classic() +
      ggplot2::scale_fill_gradient2(low = "skyblue", high = "tomato1", mid = "white",
                                    midpoint = mid_limit, limit = c(min_limit,max_limit),
                                    space = "Lab", breaks=c(min_limit, mid_limit, max_limit),
                                    name="Spearman\nCorrelation\n") +
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank())
  }else{
    condition_file <- getCondition(PD)
    if (condition %in% colnames(condition_file)){
      condition_map <- setNames(condition_file[[condition]], rownames(condition_file))

      #Reorder the levels of SampleA and SampleB so that samples with the same condition appear together
      DT.corrs$SampleA <- factor(DT.corrs$SampleA, levels = names(condition_map)[order(condition_map[DT.corrs$SampleA])])
      DT.corrs$SampleB <- factor(DT.corrs$SampleB, levels = names(condition_map)[order(condition_map[DT.corrs$SampleB])])

      DT.corrs <- DT.corrs
      g <- ggplot2::ggplot(DT.corrs, ggplot2::aes(x=SampleA, y=SampleB, fill=Spearman, label = NA)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(color='gray10') +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_gradient2(low = "skyblue", high = "tomato1", mid = "white",
                                      midpoint = mid_limit, limit = c(min_limit,max_limit),
                                      space = "Lab", breaks=c(min_limit, mid_limit, max_limit),
                                      name="Spearman\nCorrelation\n") +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank())+

        # Update x and y axis labels with conditions
        ggplot2::scale_x_discrete(labels = function(x) {
          # Show the condition only for the first sample in each group of same-condition samples
          labels <- condition_map[x]
          labels[duplicated(labels)] <- ""  # Blank out duplicates
          return(labels)
        }) +
        ggplot2::scale_y_discrete(labels = function(x) {
          # Show the condition only for the first sample in each group of same-condition samples
          labels <- condition_map[x]
          labels[duplicated(labels)] <- ""  # Blank out duplicates
          return(labels)
        })
    }else{
      stop("the selected condition does not appear in the condition file")
    }
  }
  return (g)

}
