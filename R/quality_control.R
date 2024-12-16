#' Title
#'
#' @param PD A prot_data object
#' @param by_condition True if you want to group by condition
#'
#' @return
#' @export
#'
#' @examples
plot_pg_counts <- function(PD, by_condition = F) {

  # get the number of protein groups per sample
  pgcounts <- data.table::as.data.table(table(getDataLong(PD)$Sample))
  colnames(pgcounts) <- c("Sample", "N")
  pgcounts$Condition=as.factor(gsub('_[0-9]+$','',pgcounts$Sample))

  # Order samples by ascending counts
  n_samples <- nrow(pgcounts)
  if (!by_condition){
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
  DT <- getDataLong(PD)
  n_samples <- length(unique(DT$Sample))
  g <- ggplot2::ggplot(DT, ggplot2::aes(x=Sample, y=log10(Intensity))) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill="#67a9cf") +
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "",x="",y='Log10 Protein Intensity') +
    ggplot2::theme(axis.text.x = ggplot2::element_text( angle=90)) +
    ggplot2::geom_boxplot(width=0.1) +
    ggplot2::geom_hline(color='#ef8a62', linetype='dashed',  ggplot2::aes(yintercept=quantile(log10(DT$Intensity), 0.50)))

  return(g)
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
  dt.samples <- DT[,-c(1:2)]     # Ignore info columns (subset to only intensity values)
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
plot_correlation_heatmap <- function(DT.corrs) {
  n_samples <- length(unique(DT.corrs[,'SampleA']))
  max_limit <- max(DT.corrs$Spearman)
  min_limit <- min(DT.corrs$Spearman)
  mid_limit <- as.numeric(format(((max_limit + min_limit) / 2), digits=3))
  g <- ggplot2::ggplot(DT.corrs, ggplot2::aes(x=SampleA, y=SampleB, fill=Spearman, label=Spearman)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(color='gray10') +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_gradient2(low = "skyblue", high = "tomato1", mid = "white",
                         midpoint = mid_limit, limit = c(min_limit,max_limit),
                         space = "Lab", breaks=c(min_limit, mid_limit, max_limit),
                         name="Spearman\nCorrelation\n") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank())
  return (g)
}
