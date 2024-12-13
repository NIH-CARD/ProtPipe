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
  DT=getDataLong(PD)
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
