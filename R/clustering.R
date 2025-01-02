## PCA
#' Title
#'
#' @param PD
#'
#' @return
#' @export
#'
#' @examples
get_PCs <- function(PD, condition = NA) {
  DT <- getData(PD)
  out <- list()
  ##cluster data(na=0)
  cluster_data <- DT %>%
    dplyr::select_if(is.numeric) %>%
    # Replace NA values with 0
    dplyr::mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    dplyr::filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data %>%
    dplyr::mutate(across(everything(), ~ log2(. + 1)))

  ##PCA and plot
  pca_data <- t(log2_cluster_data)
  pca=stats::prcomp(pca_data, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
  out$summary <- data.table::as.data.table(t(summary(pca)$importance), keep.rownames=T)
  setnames(out$summary, c('component','stdv','percent','cumulative'))
  out$summary$percent=round(out$summary$percent*100, digits = 2)
  pca_df = as.data.frame(pca$x)[,1:5]
  pca_df$Sample=rownames(pca_df)
  if(is.na(condition)){
    pca_df$Condition=gsub('_[0-9]+$','',rownames(pca_df))
  }else {
    pca_df$Condition=getCondition(PD)[[condition]]
  }
  out$components <- pca_df
  return(out)
}


#' Title
#'
#' @param PD
#'
#' @return
#' @export
#'
#' @examples
plot_PCs <- function(PD, condition = NA) {
  PCA <- get_PCs(PD, condition)
  p <- ggplot2::ggplot(PCA$components, ggplot2::aes(x = PC1, y = PC2, color = Condition)) +
    ggplot2::geom_point(size=4) +
    ggplot2::xlab(paste0("PC1","(",PCA$summary$percent[1],"%)")) +
    ggplot2::ylab(paste0("PC2","(",PCA$summary$percent[2],"%)")) +
    ggplot2::theme_classic()
  return(p)
}

## HC cluster
#' Title
#'
#' @param PD
#'
#' @return
#' @export
#'
#' @examples
plot_hierarchical_cluster <- function(PD) {
  DT <- getData(PD)
  cluster_data <- DT %>%
    dplyr::select_if(is.numeric) %>%
    # Replace NA values with 0
    dplyr::mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    dplyr::filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ log2(. + 1)))

  dist_mat <- stats::dist(t(log2_cluster_data)) #
  hc_cluster <- stats::hclust(dist_mat,method = "complete")
  g <- ggdendro::ggdendrogram(hc_cluster, rotate=TRUE) + ggplot2::labs(title='Hierarchical clustering')
  return(g)
}

## umap
#' Title
#'
#' @param PD
#' @param neighbors
#'
#' @return
#' @export
#'
#' @examples
get_umap <- function(PD, neighbors = 15, condition = NA) {
  DT <- getData(PD)
  ##cluster data(na=0)
  cluster_data <- DT %>%
    dplyr::select_if(is.numeric)%>%
    # Replace NA values with 0
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    dplyr::filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data %>%
    dplyr::mutate(across(everything(), ~ log2(. + 1)))

  set.seed(100)
  DT.umap <- umap::umap(t(log2_cluster_data), n_neighbors=neighbors)
  DT.out <- data.table::as.data.table(DT.umap$layout, keep.rownames=TRUE)
  data.table::setnames(DT.out, c('Sample', 'UMAP1', 'UMAP2'))
  if(is.na(condition)){
    DT.out$Condition=gsub('_[0-9]+$','',DT.out$Sample)
  }else {
    DT.out$Condition=getCondition(PD)[[condition]]
  }
  return(DT.out[])
}

#' Title
#'
#' @param PD
#'
#' @return
#' @export
#'
#' @examples
plot_umap <- function(PD, neighbors = 15, condition = NA) {
  DT <- get_umap(PD, neighbors = neighbors, condition = condition)
  g <- ggplot2::ggplot(DT, ggplot2::aes(x=UMAP1, y=UMAP2, color=Condition)) +
    ggplot2::geom_point(size=4) +
    ggplot2::theme_classic()
  return(g)
}
