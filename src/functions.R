
## FUNCTIONS #####################################################################################

tryTo <- function(infomsg='', cmd,  errormsg='ERROR: failed!') {
    # tryTo is a simple tryCatch wrapper, taking a command + message.
    # If the try block fails, an error is printed and R quits.
    cat(paste0(infomsg, '\n'))
    tryCatch(cmd,
        error = function(e){
            cat(paste0(errormsg, '\n', e))
            quit(status=1)
        },
        finally = cat('')
    )
}

replace_NAs <- function(DT, sdcols, newvalue) {
    # Within data.table `DT`,
    # for `sdcols` specified columns,
    # replaces all NA with `newvalue`
    DT=data.table(DT)
    DT[, (sdcols) := lapply(.SD, function(x) {ifelse(is.na(x),newvalue,x)}), .SDcols=sdcols]
}

replace_values <- function(DT, sdcols, oldvalue, newvalue) {
    # Within data.table `DT`,
    # for `sdcols` specified columns,
    # replaces all NA with `newvalue`
    DT[, (sdcols) := lapply(.SD, function(x) {ifelse(x==oldvalue,newvalue,x)}), .SDcols=sdcols]
}

ezwrite <- function(x, output_dir, output_filename) {
    # Wrapper for fwrite that uses standard TSV output defaults.
    # Concatenates output directory and filename for final output location.
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    fwrite(x, file=paste0(output_dir, '/', output_filename),
        quote=F,
        row.names=F,
        col.names=T,
        sep='\t')

}

convert_log_to_raw <- function(DT.original, log_base) {
  DT <- DT.original
  samplenames <- colnames(DT[,-c(1:2)])
  DT[, (samplenames) := lapply(.SD, function(x) log_base^x), .SDcols=samplenames]
  return(DT[])
}

#total proteomics############
standardize_format <- function(DT.original) {
    # Accepts an input protein group intensity data.table, whether spectronaut or DIA-NN format,
    # and restructures into one consistent style for downstream processing
    DT <- DT.original
    if("Protein.Ids" %in% colnames(DT)) {
      print("DIAnn input")
        DT[, 'Protein.Ids' := NULL]
        DT[, 'Protein.Names' := NULL]
        DT[, 'First.Protein.Description' := NULL]
        setnames(DT, 'Protein.Group', 'Protein_Group')
    }
    else if('EG.PrecursorId' %in% colnames(DT)) {
      print("Spectronaut input")
      setnames(DT, 'EG.PrecursorId', 'Peptide_Sequence')
      setnames(DT, 'PG.Genes', 'Genes')
        DT=as.data.frame(DT)
        # Use only Protein_Group and Genes
        col_select=c('Peptide_Sequence','Genes',grep('raw',colnames(DT),value = T))
        DT=DT[, col_select]
        #as number
        DT[,grep('raw',colnames(DT))]=as.data.frame(apply(DT[,grep('raw',colnames(DT))],2,as.numeric))
        DT=data.table(DT)
    }
    else if('PG.ProteinGroups' %in% colnames(DT)) {
      print("Spectronaut input")
      setnames(DT, 'PG.ProteinGroups', 'Protein_Group')
      setnames(DT, 'PG.Genes', 'Genes')
      DT=as.data.frame(DT)
      # Use only Protein_Group and Genes
      col_select=c('Protein_Group','Genes',grep('PG.Quantity',colnames(DT),value = T))
      DT=DT[, col_select]
      #as number
      DT[,grep('PG.Quantity',colnames(DT))]=as.data.frame(apply(DT[,grep('PG.Quantity',colnames(DT))],2,as.numeric))
      DT=data.table(DT)
    }
    else if('Peptide Sequence' %in% colnames(DT)) {
      print("FragPipe input")
      setnames(DT, 'Gene', 'Genes')
      colnames(DT)=gsub("\\s", "_",colnames(DT))
      # Use only Protein_Group and Genes
      DT=as.data.frame(DT)
      col_select=c('Peptide_Sequence','Genes',grep('[0-9]_Intensity',colnames(DT),value = T))
      DT=DT[, col_select]
      DT=data.table(DT)
    }

    # Remove leading directories for sample names
    # e.g. /path/to/sample1.mzML -> sample1.mzML
    setnames(DT, basename(colnames(DT)))

    # Remove trailing file extensions
    extensions <- '.mzML$|.mzml$|.RAW$|.raw$|.dia$|.DIA$|_Intensity'
    extension_samplenames <-  colnames(DT)[colnames(DT) %like% extensions]
    trimmed_samplenames <- gsub(extensions, '', extension_samplenames)
    setnames(DT, extension_samplenames, trimmed_samplenames)
    return(DT[])
}

trim_colnames <- function(DT) {
    colnames_out <- gsub(pattern="\\[.*\\] ", replacement='', x=colnames(DT))   # trim leading [N]
    colnames_out <- gsub(pattern="\\..*\\.PG\\.Quantity|\\.PG\\.Quantity|\\..*Quantity.*", replacement='', x=colnames_out)   # remove suffix
    return(colnames_out)
}

unique_data <- function(DT, col,output_dir,output_filename) {
  DT_unique <- DT %>%
    # Calculate missing values for all columns except specified ones
    mutate(missing_value = rowSums(is.na(select(., -contains("Protein_Group|Genes|Precursor|Peptide_Sequence"))))) %>%

    # Calculate median values for numeric columns except specified ones
    mutate(median = rowMedians(as.matrix(select(., -contains("Protein_Group|Genes|Precursor|Peptide_Sequence|missing_value")) %>%
                                           select_if(is.numeric)), na.rm = TRUE)) %>%

    # Group by the specified column dynamically
    group_by(!!sym(col))%>%

    # Filter rows with minimum missing values and maximum median intensity
    filter(missing_value == min(missing_value)) %>%
    slice(which.max(median)) %>%
    ungroup() %>%

    # Remove the temporary columns
    select(-c(missing_value, median))

  # Write the result to a CSV file
  ezwrite(DT_unique, output_dir = output_dir,output_filename = output_filename)
  return(DT_unique)
}


melt_intensity_table <- function(DT) {
    # Converts intensity data.table to long format
    # info_cols <- c('Protein_Group', 'Genes', 'First_Protein_Description')
  DT.long <- melt(DT,
                  measure.vars=names(DT)[sapply(DT, function(x) all(is.numeric(x)))],
                  variable.name='Sample',
                  value.name='Intensity')
  DT.long=data.table(DT.long)
  DT.long=DT.long[DT.long$Intensity>0,]
  return(DT.long)
}

plot_pg_counts <- function(DT.long, output_dir) {
  pgcounts <- DT.long[, .N, by=Sample]
  pgcounts$Condition=as.factor(gsub('_[0-9]+$','',pgcounts$Sample))
  # Order samples by ascending counts
  ezwrite(pgcounts, output_dir, 'protein_group_nonzero_counts.tsv')
  n_samples <- nrow(pgcounts)
  if (n_samples > 20) {
    p=ggplot(pgcounts, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="#67a9cf")+
      theme_classic()+
      labs(fill = "",x="",y='Number of Protein Groups')+
      scale_x_discrete(guide = guide_axis(angle = 90))
  }
  if (n_samples <= 20) {
    p=ggplot(pgcounts, aes(x=Sample, y=N,fill=Condition)) +
      geom_bar(stat="identity")+
      theme_classic()+
      labs(fill = "",x="",y='Number of Protein Groups')+
      scale_x_discrete(guide = guide_axis(angle = 90))+
      geom_text(aes(label=N, y=N + (0.05*max(pgcounts$N))))
  }

  if (n_samples>50){
    ggsave(plot = p,filename = paste0(output_dir, 'protein_group_nonzero_counts.pdf'),width = n_samples/10,height = 6)
    cat(paste0('   -> ', output_dir, 'protein_group_nonzero_counts.pdf', '\n'))
  }else {
    ggsave(plot = p,filename = paste0(output_dir, 'protein_group_nonzero_counts.pdf'),width = 8,height = 6)
    cat(paste0('   -> ', output_dir, 'protein_group_nonzero_counts.pdf', '\n'))
  }
  # group with sd
  summary_data <- pgcounts %>%
    group_by(Condition) %>%
    summarize(mean = mean(N), sd = sd(N)) %>%
    arrange(Condition)
  p=ggplot(summary_data, aes(x=as.factor(Condition), y=mean)) +
    geom_bar(stat="identity",fill="#67a9cf", position=position_dodge())+
    theme_classic()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    labs(fill = "",x="",y='Number of Phospho-sites')
  ggsave(plot = p,filename = paste0(output_dir,'/','protein_group_nonzero_counts_Condition.pdf'),height=5,width = 4)
  cat(paste0('   -> ', output_dir, '/','protein_group_nonzero_counts_Condition.pdf', '\n'))
  return(pgcounts)
  }

plot_pep_counts <- function(DT, output_dir, output_filename) {
  n_samples <- nrow(DT)
  if (n_samples > 20) {
    p=ggplot(DT, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="#67a9cf")+
      theme_classic()+
      labs(fill = "",x="Sample",y='Number of Peptides')+
      scale_x_discrete(guide = guide_axis(angle = 90))
  }
  if (n_samples <= 20) {
    p=ggplot(DT, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="#67a9cf")+
      theme_classic()+
      labs(fill = "",x="Sample",y='Number of Peptides')+
      scale_x_discrete(guide = guide_axis(angle = 90))+
      geom_text(aes(label=N, y=N + (0.05*max(pep_counts$N))))
  }

  if (n_samples>50){
    ggsave(plot = p,filename = paste0(output_dir, output_filename),width = n_samples/10,height = 6)
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
  }else {
    ggsave(plot = p,filename = paste0(output_dir, output_filename),width = 8,height = 6)
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
  }
}

plot_pg_thresholds <- function(DT, output_dir, output_filename) {
    # G
    g <- ggplot(DT, aes(x=Threshold, y=N, color=Sample)) +
            geom_line() +
            geom_point(shape=21, alpha=0.5) +
            theme_classic() +
            labs(x=paste0('Minimum Log[', opt$log_base, '](Intensity) Threshold'),
            y=paste0('N Protein Groups where Log[', opt$log_base, '](Intensity)> 0'))
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g,
        filename=paste0(output_dir, output_filename),
        width=20,
        height=20,
        units='cm'
    )
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
}

plot_density <- function(DT.original, output_dir, output_filename) {
    # Currently UNUSED as the beeswarm function serves the purpose well
    DT <- DT.original
    intensity_median <- median(DT[,Intensity])
    n_samples <- length(unique(DT[,Sample]))
    dat.quantiles <- DT[, list(
                    'q025'=quantile(Intensity, 0.025),
                    'q25'=quantile(Intensity, 0.25),
                    'q50'=quantile(Intensity, 0.50),
                    'q75'=quantile(Intensity, 0.75),
                    'q975'=quantile(Intensity, 0.975)
                    ), by=Sample]

    dat.legend <- melt(dat.quantiles[which.min(as.numeric(Sample))], measure.vars=colnames(dat.quantiles[,-1]))
    dat.legend[, qlabel := tstrsplit(variable, split='q')[2]]
    dat.legend[, qlabel := paste0('0.', qlabel)]
    dat.legend[, qlabel := as.numeric(qlabel)]

    g <- ggplot(DT, linetype='solid', aes(x=Intensity)) +
        geom_density(fill='gray80') +
        theme_few() +
        facet_grid(Sample~., switch='y') +
        geom_vline(xintercept=intensity_median, color='#ef8a62') +
        geom_vline(data=dat.quantiles, linetype='solid',  alpha=0.7, aes(xintercept=q50))+
        geom_vline(data=dat.quantiles, linetype='dashed', alpha=0.7,  aes(xintercept=q25))+
        geom_vline(data=dat.quantiles, linetype='dashed', alpha=0.7,  aes(xintercept=q75))+
        geom_vline(data=dat.quantiles, linetype='dotted', alpha=0.7,  aes(xintercept=q025))+
        geom_vline(data=dat.quantiles, linetype='dotted', alpha=0.7,  aes(xintercept=q975))+
        theme(strip.text.y.left = element_text(angle = 0, hjust=0.5, vjust=0.5)) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) +
        labs(x=paste0('Log[', opt$log_base, '](Intensity)'), title='Intensity Distribution across Samples') +
        geom_label(data=dat.legend, aes(x=value, y=0.285, label=qlabel)) +
        theme(panel.border = element_blank()) +
        ylim(0,0.3)
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g,
        filename=paste0(output_dir, output_filename),
        height=2.5*n_samples,
        width=30,
        units='cm'
    )
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
}

## Global Normalization function by Median with NA handling
median_normalize_intensity <- function(DT.original) {
    DT <- DT.original
    # Get global median of intensity values
    global_median <- median(DT[, Intensity])
    DT[, 'sample_median' := median(Intensity), by=Sample]
    DT[, 'global_median' := global_median]
    DT[, 'NormInt' := Intensity * (global_median / sample_median)]
    DT[, c('Intensity', 'sample_median', 'global_median') := NULL]
    setnames(DT, 'NormInt', 'Intensity')
    return(DT[])
}
## Global Normalization function by Mean with NA handling
mean_normalize_intensity <- function(DT.original) {
  DT <- DT.original
  # Get global mean of intensity values
  global_mean <- mean(DT[, Intensity])
  DT[, 'sample_mean' := mean(Intensity), by=Sample]
  DT[, 'global_mean' := global_mean]
  DT[, 'NormInt' := Intensity * (global_mean / sample_mean)]
  DT[, c('Intensity', 'sample_mean', 'global_mean') := NULL]
  setnames(DT, 'NormInt', 'Intensity')
  return(DT[])
}

plot_pg_intensities <- function(DT, output_dir, output_filename) {
  DT=data.table(DT)
    n_samples <- length(unique(DT$Sample))
    g <- ggplot(DT, aes(x=Sample, y=log10(Intensity))) +
        geom_boxplot(outlier.shape = NA, fill="#67a9cf") +
        theme_classic() +
        labs(fill = "",x="",y='Log10 Protein Intensity') +
        theme(axis.text.x = element_text( angle=90)) +
        geom_boxplot(width=0.1) +
        geom_hline(color='#ef8a62', linetype='dashed',  aes(yintercept=quantile(log10(DT$Intensity), 0.50)))

    if (n_samples>50){
      ggsave(plot = g,filename = paste0(output_dir, output_filename),width = n_samples/10,height =6)
      cat(paste0('   -> ', output_dir, output_filename, '\n'))
    }else{
      ggsave(plot = g,filename = paste0(output_dir, output_filename),width = 8,height = 6)
      cat(paste0('   -> ', output_dir, output_filename, '\n'))
    }
}

plot_pep_intensities <- function(DT, output_dir, output_filename) {
  n_samples <- length(unique(DT$Sample))
  g <- ggplot(DT, aes(x=Sample, y=log10(Intensity))) +
    geom_violin(trim=FALSE, fill="steelblue")+
    theme_classic()+
    labs(fill = "",x="",y='Log10 Peptide Intensity')+
    theme(axis.text.x = element_text( angle=90))+
    geom_boxplot(width=0.1)+
    geom_hline( color='red', linetype='dashed',  aes(yintercept=quantile(log10(DT$Intensity), 0.50)))

    if (n_samples>50){
      ggsave(plot = g,filename = paste0(output_dir, output_filename),width = n_samples/10,height =6)
      cat(paste0('   -> ', output_dir, output_filename, '\n'))
    }else{
      ggsave(plot = g,filename = paste0(output_dir, output_filename),width = 8,height = 6)
      cat(paste0('   -> ', output_dir, output_filename, '\n'))
    }
}

## correlation
plot_correlation_heatmap <- function(DT.corrs, output_dir, output_filename) {
  n_samples <- length(unique(DT.corrs[,SampleA]))
  max_limit <- max(DT.corrs$Spearman)
  min_limit <- min(DT.corrs$Spearman)
  mid_limit <- as.numeric(format(((max_limit + min_limit) / 2), digits=3))
  g <- ggplot(DT.corrs, aes(x=SampleA, y=SampleB, fill=Spearman, label=Spearman)) +
    geom_tile() +
    geom_text(color='gray10') +
    theme_classic() +
    scale_fill_gradient2(low = "skyblue", high = "tomato1", mid = "white",
                         midpoint = mid_limit, limit = c(min_limit,max_limit),
                         space = "Lab", breaks=c(min_limit, mid_limit, max_limit),
                         name="Spearman\nCorrelation\n") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
  cat(paste0('   -> ', output_dir, output_filename, '\n'))
  ggsave(g,
        filename=paste0(output_dir, output_filename),
        height=1.2*n_samples,
        width=2*n_samples,
        units='cm',limitsize = FALSE)
  cat(paste0('   -> ', output_dir, output_filename, '\n'))

}

get_spearman <- function(DT.original) {
  DT <- DT.original
  #### Pairwise correlations between sample columns
  dt.samples <- DT[,-c(1:2)]     # Ignore info columns (subset to only intensity values)
  dt.corrs <- cor(log2(as.matrix(na.omit(dt.samples))+1), method='spearman')

  # Format correlations as 3 digits
  dt.corrs <- data.table(melt(dt.corrs, measure.vars=dt.corrs[,rn], value.name='Spearman'))
  dt.corrs <- dt.corrs[! is.na(Spearman)]
  setnames(dt.corrs, c('Var1', 'Var2'), c('SampleA','SampleB'))
  dt.corrs[, Spearman := as.numeric(format(Spearman, digits=3))]

  return(dt.corrs[])
}

get_pearson_matrix <- function(DT.original) {
    DT <- DT.original
    #### Pairwise correlations between sample columns
    dt.samples <- DT[,-c(1:2)]     # Ignore info columns (subset to only intensity values)
    dt.corrs <- cor(as.matrix(na.omit(dt.samples)), method='pearson')

    dt.corrs <- as.data.table(dt.corrs, keep.rownames=T)
    #dt.corrs <- melt(dt.corrs, measure.vars=dt.corrs[,rn], value.name='Spearman')
    #dt.corrs <- dt.corrs[! is.na(Spearman)]
    #setnames(dt.corrs, c('rn', 'variable'), c('SampleA','SampleB'))

    # Format correlations as 3 digits
    dt.corrs[, Spearman := as.numeric(format(Spearman, digits=3))]

    # Adjust levels such that both axes plot samples in the same order
    dt.corrs.levels <- sort(as.character(unique(dat.long$Sample)))
    dt.corrs[, SampleA := factor(SampleA, levels=dt.corrs.levels)]
    dt.corrs[, SampleB := factor(SampleB, levels=dt.corrs.levels)]
    return(dt.corrs[])
}

## PCA
get_PCs <- function(DT) {
  out <- list()
  ##cluster data(na=0)
  cluster_data <- DT %>%
    select_if(is.numeric) %>%
    # Replace NA values with 0
    mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data %>%
    mutate(across(everything(), ~ log2(. + 1)))

  ##PCA and plot
  pca_data=t(log2_cluster_data)
  pca=prcomp(pca_data, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
  out$summary <- as.data.table(t(summary(pca)$importance), keep.rownames=T)
  setnames(out$summary, c('component','stdv','percent','cumulative'))
  out$summary$percent=round(out$summary$percent*100, digits = 2)
  pca_df = as.data.frame(pca$x)[,1:5]
  pca_df$Sample=rownames(pca_df)
  pca_df$Condition=gsub('_[0-9]+$','',rownames(pca_df))
  out$components <- pca_df
  return(out)
}

plot_PCs <- function(PCA, output_dir, output_filename) {
  p <- ggplot(PCA$components, aes(x = PC1, y = PC2, color = Condition)) +
        geom_point(size=4) +
        xlab(paste0("PC1","(",PCA$summary$percent[1],"%)")) +
        ylab(paste0("PC2","(",PCA$summary$percent[2],"%)")) +
        theme_classic()
    ggsave(p,filename=paste0(output_dir, output_filename), height = 4,width = 5)
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
}

## HC cluster
plot_hierarchical_cluster <- function(DT, output_dir) {
  cluster_data <- DT %>%
    select_if(is.numeric) %>%
    # Replace NA values with 0
    mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data %>%
    mutate(across(everything(), ~ log2(. + 1)))

  dist_mat <- dist(t(log2_cluster_data)) #
  hc_cluster <- hclust(dist_mat,method = "complete")
  g <- ggdendrogram(hc_cluster, rotate=TRUE) + labs(title='Hierarchical clustering')
  cat(paste0('   -> ', output_dir, '\n'))
  ggsave(g, filename=paste0(output_dir, 'hc_cluster_log2.pdf'),width = 8,height = 6)
  cat(paste0('   -> ', output_dir, 'hc_cluster_log2.pdf', '\n'))
}

## umap
get_umap <- function(DT, neighbors) {
  ##cluster data(na=0)
  cluster_data <- DT %>%
    select_if(is.numeric)%>%
    # Replace NA values with 0
    mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data %>%
    mutate(across(everything(), ~ log2(. + 1)))

    set.seed(100)
    DT.umap <- umap(t(log2_cluster_data), n_neighbors=neighbors)
    DT.out <- as.data.table(DT.umap$layout, keep.rownames=TRUE)
    setnames(DT.out, c('Sample', 'UMAP1', 'UMAP2'))
    DT.out$Condition=gsub('_[0-9]+$','',DT.out$Sample)
    return(DT.out[])
}

plot_umap <- function(DT, output_dir, output_filename) {
    g <- ggplot(DT, aes(x=UMAP1, y=UMAP2, color=Condition)) +
    geom_point(size=4) +
    theme_classic()
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(output_dir, output_filename), height = 4,width = 5)
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
}

## ttest
do_t_test = function(DT,col, design_matrix,DE_dir,EA_dir) {
  DT <- data.frame(DT)
  name=names(DT)[sapply(DT, function(x) !all(is.numeric(x)))]
  #comparison
  conditions <- unique(design_matrix$condition)
  for (treatment in conditions) {
    controls=unique(design_matrix$control[design_matrix$condition==treatment])
    for (control in controls) {
      print(paste0(treatment, ' vs ', control, ' Differential Analysis by ttest '))
      treatment_samples=grep(treatment,colnames(DT),value = T)
      control_samples=grep(control,colnames(DT),value = T)
      n_treatment <- length(treatment_samples)
      n_control <- length(control_samples)
      #ttest table
      # Initial data transformation
      DT_ttest <- DT[,c(col, treatment_samples, control_samples)]

      # Convert NA to 0 and add missing value calculations
      DT_ttest <- DT_ttest %>%
        # Convert NA to 0
        mutate(across(where(is.numeric), ~ coalesce(., 0))) %>%
        # Calculate missing values
        mutate(
          missing_value = rowSums(select(., -col) == 0),  # Total missing values per row, excluding PTM column
          missing_value_c = rowSums(select(., all_of(control_samples)) == 0),  # Missing values in control samples per row
          missing_value_t = rowSums(select(., all_of(treatment_samples)) == 0)  # Missing values in treatment samples per row
        ) %>%
        # Filter missing values
        filter(
          !(missing_value_t > (n_treatment / 2) & missing_value_t < n_treatment),
          !(missing_value_c > (n_control / 2) & missing_value_c < n_control),
          missing_value != (n_treatment + n_control)
        ) %>%
        # Select relevant columns
        select(-missing_value, -missing_value_c, -missing_value_t)%>%
        as.data.frame()
      # Restore row names from PTM column and remove the PTM column
      rownames(DT_ttest) <- DT_ttest[,col]
      DT_ttest[,col]=NULL

      # Perform t-test  on treatment and control columns
      t_test <- apply(DT_ttest, 1, function(x){
        a =factor(c(rep('treatment',n_treatment),
                    rep("control",n_control)),
                  levels = c('treatment',"control"))
        fvalue=var.test(x~a)
        if (!is.na(fvalue$p.value)){
          if (fvalue$p.value > 0.05){
            result <- t.test(x~a, var.equal = T)
          }else{
            result <- t.test(x~a, var.equal = F)
          }
        }
        treatment_estimate <- as.numeric(unlist(result$estimate[1]))
        control_estimate <- as.numeric(unlist(result$estimate[2]))
        return(data.table('P.Value'=result$p.value,
                          'treatment_estimate'=treatment_estimate,
                          'control_estimate'=control_estimate)
        )
      })
      t_test <- rbindlist(t_test)
      result_ttest <- cbind(DT_ttest, t_test)
      # Merge DT and result_ttest based on PTM column
      result_ttest <- merge(DT[, name, drop=FALSE], result_ttest, by.x=col, by.y='row.names')
      # Calculate logFC and adjust P-values using dplyr
      result_ttest <- result_ttest %>%
        mutate(logFC = log2((treatment_estimate+1) / (control_estimate + 1)),
               adj.P.Val = p.adjust(P.Value, method='BH'))
      ezwrite(result_ttest[order(result_ttest$adj.P.Val),], DE_dir, paste0(treatment, '_vs_', control, '.tsv'))
      plot_volcano(DT.original =result_ttest ,
                   lfc_threshold = opt$lfc_threshold,
                   fdr_threshold = opt$fdr_threshold,
                   out_dir = DE_dir,label_col = 'Genes',
                   output_filename =paste0(treatment, '_vs_', control, '_ttest'),
                   labelgene = opt$labelgene)
      enrich_pathway(DT.original = result_ttest,
                     treatment = treatment,
                     control = control,
                     outdir = EA_dir,
                     lfc_threshold = opt$lfc_threshold,
                     fdr_threshold = opt$fdr_threshold,
                     enrich_pvalue = opt$enrich_pvalue)
    }
  }
}

## plot_volcano
plot_volcano <- function(DT.original, out_dir, output_filename, label_col ,lfc_threshold, fdr_threshold, labelgene) {
  options(ggrepel.max.overlaps = Inf)
  DT <- DT.original

  # Set initial group to 'Others' and update based on thresholds
  DT <- DT %>%
    mutate(Group = 'Others',
           Group = if_else(logFC >= lfc_threshold, 'UP', Group),
           Group = if_else(logFC <= -lfc_threshold, 'DOWN', Group),
           Group = if_else(adj.P.Val >= fdr_threshold, 'Others', Group),
           labeltext = '')

  # If labelgene is provided, update labeltext accordingly
  if (!is.null(labelgene)) {
    label_gene=fread(labelgene)
    DT <- DT %>%
      mutate(labeltext = if_else(!!sym(label_col) %in% label_gene$gene, !!sym(label_col), labeltext))
  } else{
    # Select top 5 genes for UP and DOWN groups
    top5_gene <- DT %>%
      filter(Group == 'UP') %>%
      arrange(desc(logFC)) %>%
      slice_head(n = 5) %>%
      bind_rows(
        DT %>%
          filter(Group == 'DOWN') %>%
          arrange(logFC) %>%
          slice_head(n = 5)
      )

    # Label top 5 genes using the specified label column
    DT <- DT %>%
      mutate(labeltext = if_else(!!sym(label_col) %in% top5_gene[[label_col]], !!sym(label_col), ''))
  }




  g <- ggplot(DT, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Group)) +
    scale_color_manual(breaks = c("DOWN", "Others", "UP"),
                       values = c("#67a9cf", "#969696", "#ef8a62")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    geom_label_repel(
      data = subset(DT),
      aes(label = labeltext),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    geom_vline(xintercept = lfc_threshold, linetype = "dashed") +
    geom_vline(xintercept = -lfc_threshold, linetype = "dashed") +
    theme_classic()

  cat(paste0(out_dir, output_filename,'_vocanol.pdf'))
  ggsave(g, filename = paste0(out_dir, output_filename,'_vocanol.pdf'), width = 8, height = 8)
  cat(paste0('   -> ', out_dir, output_filename,'_vocanol.pdf', '\n'))
}

plot_volcano_pep <- function(DT.original, treatment, control, log_base, lfc_threshold, fdr_threshold, out_dir, labelgene) {
  options(ggrepel.max.overlaps=Inf)
  DT <- DT.original
  DT[, 'Group' := 'Others']
  DT[logFC >= lfc_threshold, 'Group' := 'UP']
  DT[logFC <= -lfc_threshold, 'Group' := 'DOWN']
  DT[adj.P.Val >= fdr_threshold, 'Group' := 'Others']
  DT[, labeltext := '']

  up_gene_5 <- DT[Group == 'UP', ]
  up_gene_5 <- up_gene_5[order(up_gene_5$logFC,decreasing = T)[1:5],]
  down_gene_5 <- DT[Group == 'DOWN', ]
  down_gene_5 <- down_gene_5[order(down_gene_5$logFC,decreasing = F)[1:5],]
  top5_gene <- rbind(up_gene_5,down_gene_5)
  DT[Peptide_Sequence %in% top5_gene$Peptide_Sequence, labeltext := Peptide_Sequence]

  if (!is.null(labelgene)) {
    DT[Peptide_Sequence %in% labelgene, labeltext := Peptide_Sequence]
  }
  g <- ggplot(DT, aes(x=logFC, y=-log10(adj.P.Val))) +
    geom_point(aes(color = Group)) +
    scale_color_manual(breaks = c("DOWN", "Others", "UP"),
                       values=c("#67a9cf", "#969696","#ef8a62"))+
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_label_repel(
      data = subset(DT),
      aes(label = labeltext),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))+
    geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed")+
    geom_vline(xintercept=lfc_threshold, linetype="dashed")+
    geom_vline(xintercept=-lfc_threshold, linetype="dashed")+
    theme_classic()

  output_filename <- paste0(treatment, '_vs_', control, '.pdf')
  cat(paste0('   -> ', out_dir, output_filename, '\n'))
  ggsave(g, filename=paste0(out_dir, output_filename),width = 8,height = 8)
  cat(paste0('   -> ', output_dir, output_filename, '\n'))
}

plot_volcano_from_raw <- function(DT.original, treatment, control, log_base, lfc_threshold, fdr_threshold, out_dir) {
  DT <- DT.original
  DT[, log_foldchange := log(ratio, base=log_base)]
  log_lfc_threshold <- log(lfc_threshold, base=log_base)
  log_fdr_threshold <- -1*log(fdr_threshold, base=log_base)
  DT[, labeltext := '']
  g <- ggplot(DT, aes(x=log_foldchange, y=-1*log10(q))) +
    geom_point() +
    theme_few() +
    geom_label_repel(aes(label=labeltext)) +
    geom_hline(yintercept=log_fdr_threshold, linetype='dashed', alpha=0.5) +
    geom_vline(xintercept=log_lfc_threshold, linetype='dashed', alpha=0.5) +
    geom_vline(xintercept=(-1*log_lfc_threshold), linetype='dashed', alpha=0.5) +
    labs(x=paste0('Log[', opt$log_base, '](Intensity) fold change'),
         y='-Log[10](q)',
         title=paste0(treatment, ' vs ', control)
    )
  output_filename <- paste0(treatment, '_vs_', control, '.png')
  cat(paste0('   -> ', out_dir, output_filename, '\n'))
  ggsave(g, filename=paste0(out_dir, output_filename), height=16, width=16, units='cm')
  cat(paste0('   -> ', output_dir, output_filename, '\n'))
}

###pathway analysis
enrich_pathway = function(DT.original, treatment, control, outdir, lfc_threshold, fdr_threshold, enrich_pvalue){
  ##create dir
  dir=paste0(outdir,'/', treatment, '_vs_', control)
  if (!dir.exists(dir)){
    dir.create(dir,recursive = T)
  }

  DT <- DT.original %>%
    mutate(Group = 'Others') %>%
    mutate(Group = ifelse(logFC >= lfc_threshold, 'UP', Group)) %>%
    mutate(Group = ifelse(logFC <= -lfc_threshold, 'DOWN', Group)) %>%
    mutate(Group = ifelse(adj.P.Val >= fdr_threshold, 'Others', Group))%>%
    data.frame()

  ##ENTREZID names########
  entrizid = data.frame(bitr(DT$Genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop=T))
  DT <- merge(DT,entrizid,by.x='Genes',by.y='SYMBOL') %>%
    na.omit() %>%
    arrange(desc(logFC))
  all_gene_vector=DT$logFC
  names(all_gene_vector)=DT$ENTREZID

  ## up and down regulated genes
  up_genes=DT[which(DT$logFC>=lfc_threshold&DT$adj.P.Val<=fdr_threshold),]
  down_genes=DT[which(DT$logFC<=(-lfc_threshold)&DT$adj.P.Val<=fdr_threshold),]
  print('Processing up-gene enrichment analysis')
  if (nrow(up_genes)>0){
    enrichAll(gene_id=up_genes$ENTREZID, all_gene_vector=all_gene_vector, out_prefix = 'up',outdir=dir,width = 12,height = 8,enrich_pvalue=enrich_pvalue)
  }
  print('Processing down-gene enrichment analysis')
  if (nrow(down_genes)>0){
    enrichAll(gene_id=down_genes$ENTREZID, all_gene_vector=all_gene_vector, out_prefix = 'down',outdir=dir,width = 12,height = 8,enrich_pvalue=enrich_pvalue)
  }
  print('Processing Gene Set Enrichment Analysis')
  if (nrow(down_genes)+nrow(up_genes) >0){
    gseGO <- gseGO(geneList=all_gene_vector,
                   ont ="ALL",
                   keyType = "ENTREZID",
                   pvalueCutoff = 0.05,
                   OrgDb = 'org.Hs.eg.db')
    gseGO=setReadable(gseGO, OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
    gseGO_res=gseGO@result
    if (nrow(gseGO_res) > 0) {
      p=dotplot(gseGO, showCategory=10, split=".sign") + facet_grid(.~.sign)
      ggsave(paste0(dir,'/','gseGO_dotplot.pdf'),p,width = 10,height = 10)
      cat(paste0('   -> ', dir,'/','gseGO_dotplot.pdf', '\n'))
      x2 <- pairwise_termsim(gseGO)
      p=emapplot(x2)
      ggsave(paste0(dir,'/', 'gseGO_emap.pdf'),p,width = 10,height = 10)
      cat(paste0('   -> ', dir,'/', 'gseGO_emap.pdf', '\n'))
      gseKEGG <- gseKEGG(geneList     = all_gene_vector,
                         organism     = 'hsa',
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "none",
                         keyType       = "kegg")
      gseKEGG=setReadable(gseKEGG, OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
      gseKEGG_res=gseKEGG@result
      gseKEGG_res$ONTOLOGY='KEGG'
      p=dotplot(gseKEGG, showCategory=10, split=".sign") + facet_grid(.~.sign)
      ggsave(paste0(dir,'/', 'gseKEGG_dotplot.pdf'),p,width = 10,height = 10)
      cat(paste0('   -> ', dir,'/', 'gseKEGG_dotplot.pdf', '\n'))
      x2 <- pairwise_termsim(gseKEGG)
      p=emapplot(x2)
      ggsave(paste0(dir,'/', 'gseKEGG_emap.pdf'),p,width = 10,height = 10)
      cat(paste0('   -> ', dir,'/', 'gseKEGG_emap.pdf', '\n'))
      gse_res=rbind(gseKEGG_res,gseGO_res)
      write.csv(gse_res,paste0(dir,'/',  'gse_res.csv'),row.names = F)
    }else{
      print("No enriched pathways found in the Gene Set Enrichment Analysis analysis")}
    }
}

enrichAll = function(gene_id, all_gene_vector, outdir='.', out_prefix = NULL, width=12, height=8, enrich_pvalue=1){

  all_gene_vector <- na.omit(all_gene_vector)
  if (!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }
  enrich_res = data.frame()
  ### GO
  print('Processing GO')
  GO <- enrichGO(gene          = gene_id,
                 universe      = names(all_gene_vector),
                 OrgDb         = 'org.Hs.eg.db',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = enrich_pvalue,
                 qvalueCutoff  = enrich_pvalue,
                 readable      = TRUE)
  GO_res=GO@result
  if (nrow(GO) > 0) {
    p=dotplot(GO, showCategory=10,split='ONTOLOGY')+
      facet_grid(ONTOLOGY~.,scales = "free",space = "free")

    g=barplot(GO, showCategory=10,split='ONTOLOGY')+
      facet_grid(ONTOLOGY~.,scales = "free",space = "free")
    if (nrow(GO) < 5) {
      ggsave(plot = p,
           filename = file.path(outdir,paste0(out_prefix,'_GO_res_dotplot.pdf')),
           width=8,height=nrow(GO)*1.3)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_res_dotplot.pdf', '\n'))
      ggsave(plot = g,
           filename = file.path(outdir,paste0(out_prefix,'_GO_res_barplot.pdf')),
           width=8,height=nrow(GO)*1.3)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_res_barplot.pdf', '\n'))
      }
    else{
      ggsave(plot = p,
             filename = file.path(outdir,paste0(out_prefix,'_GO_res_dotplot.pdf')),
             width=8,height=8)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_res_dotplot.pdf', '\n'))
      ggsave(plot = g,
             filename = file.path(outdir,paste0(out_prefix,'_GO_res_barplot.pdf')),
             width=8,height=8)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_res_barplot.pdf', '\n'))
      }
  } else {
    print("No enriched pathways found in the GO analysis")}
  # GO CC
  print('Processing GO CC')
  ego_cc <- enrichGO(gene          = gene_id,
                     universe      = names(all_gene_vector),
                     OrgDb         = org.Hs.eg.db,
                     ont           = "CC",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = enrich_pvalue,
                     qvalueCutoff  = enrich_pvalue,
                     readable      = TRUE)
    if (nrow(ego_cc)>0){
      barplot(ego_cc, showCategory=20,label_format = 100)
      ggsave(file.path(outdir,paste0(out_prefix,'_GO_CC_barplot.pdf')),width = width,height=height)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_CC_barplot.pdf', '\n'))
      dotplot(ego_cc, showCategory=20,label_format = 100)
      ggsave(file.path(outdir,paste0(out_prefix,'_GO_CC_dotplot.pdf')),width = width,height=height)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_CC_dotplot.pdf', '\n'))
    }else {
      print('NO GO CC')
    }
  # GO BP
  print('Processing GO BP')
  ego_bp <- enrichGO(gene          = gene_id,
                     universe      = names(all_gene_vector),
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = enrich_pvalue,
                     qvalueCutoff  = enrich_pvalue,
                     readable      = TRUE)
    if (nrow(ego_bp)>0){
      barplot(ego_bp, showCategory=20,label_format = 100)
      ggsave(file.path(outdir,paste0(out_prefix,'_GO_BP_barplot.pdf')),width=width,height=height)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_BP_barplot.pdf', '\n'))
      dotplot(ego_bp, showCategory=20,label_format = 100)
      ggsave(file.path(outdir,paste0(out_prefix,'_GO_BP_dotplot.pdf')),width=width,height=height)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_BP_dotplot.pdf', '\n'))
    }else {
      print('NO GO BP')
      }
  # GO MF
  print('Processing GO MF')
  ego_mf <- enrichGO(gene          = gene_id,
                     universe      = names(all_gene_vector),
                     OrgDb         = org.Hs.eg.db,
                     ont           = "MF",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = enrich_pvalue,
                     qvalueCutoff  = enrich_pvalue,
                     readable      = TRUE)
    if (nrow(ego_mf)>0){
      barplot(ego_mf, showCategory=20,label_format = 100)
      ggsave(file.path(outdir,paste0(out_prefix,'_GO_MF_barplot.pdf')),width=width,height=height)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_MF_barplot.pdf', '\n'))
      dotplot(ego_mf, showCategory=20,label_format = 100)
      ggsave(file.path(outdir,paste0(out_prefix,'_GO_MF_dotplot.pdf')),width=width,height=height)
      cat(paste0('   -> ', outdir, out_prefix,'_GO_MF_dotplot.pdf', '\n'))
    } else {
      print('NO GO MF')
      }

  # KEGG
  print('Processing KEGG')
  KEGG <- enrichKEGG(gene         = gene_id,
                     organism     = 'hsa',
                     universe      = names(all_gene_vector),
                     pvalueCutoff = enrich_pvalue,
                     qvalueCutoff  = enrich_pvalue)
  if (!is.null(KEGG)) {
    KEGG=setReadable(KEGG, OrgDb='org.Hs.eg.db', keyType = "ENTREZID")
    KEGG_res=KEGG@result
    if (nrow(KEGG) > 0) {
      dotplot(KEGG, showCategory=20)
      ggsave(file.path(outdir,paste0(out_prefix,'_KEGG_dotplot.pdf')),width=width,height=height)
      cat(paste0('   -> ', outdir,out_prefix,'_KEGG_dotplot.pdf', '\n'))
      barplot(KEGG, showCategory=20)
      ggsave(file.path(outdir,paste0(out_prefix,'_KEGG_barplot.pdf')),width=width,height=height)
      cat(paste0('   -> ', outdir,out_prefix,'_KEGG_barplot.pdf', '\n'))
    } else {
      print("No enriched pathways found in the KEGG object.")}
    }else {
      print("No enriched pathways found in the KEGG object.")
      KEGG_res=data.frame()
      }

  print('Save enrichment analysis results')
  if (nrow(KEGG_res) > 0 &nrow(GO_res) > 0) {
    GO_res$category=''
    GO_res$subcategory=''
    KEGG_res$ONTOLOGY="KEGG"
    enrich_res=rbind(GO_res,KEGG_res)
    ezwrite(enrich_res, outdir, paste0(out_prefix,'_enrich_res.tsv'))
  }
  if((nrow(KEGG_res) == 0)){
    ezwrite(GO_res, outdir, paste0(out_prefix,'_enrich_res.tsv'))
  }
  if((nrow(GO_res) == 0)){
    ezwrite(KEGG_res, outdir, paste0(out_prefix,'_enrich_res.tsv'))
  }

}

##limma
##limma DE
do_limma=function(Log2_DT, design_matrix,DE_dir,EA_dir,lfc_threshold,fdr_threshold,enrich_pvalue){
  #data
  name=colnames(Log2_DT)[1:2]
  rownames(Log2_DT)=Log2_DT$Protein_Group
  Log2_DT=data.frame(Log2_DT)
  #comparation
  conditions <- unique(design_matrix$condition)
  for (treatment in conditions) {
    controls=unique(design_matrix$control[design_matrix$condition==treatment])
    for (control in controls) {
      print(paste0(treatment, ' vs ', control, ' Differential Analysis by Limma '))
      treatment_samples=grep(treatment,colnames(Log2_DT),value = T)
      control_samples=grep(control,colnames(Log2_DT),value = T)
      DT_limma <- Log2_DT[,c(treatment_samples, control_samples)]
      n_treatment <- length(treatment_samples)
      n_control <- length(control_samples)
      # Convert NA to 0
      DT_limma[is.na(DT_limma)] <- 0
      if (n_treatment>3&n_control>3) {
        # Drop rows (protein groups) with > 50% missingness in samples
        DT_limma$missing_value= apply(DT_limma, 1, function(x) sum(x==0))
        DT_limma$missing_value_c= apply(DT_limma[,control_samples], 1, function(x) sum(x==0))
        DT_limma$missing_value_t= apply(DT_limma[,treatment_samples], 1, function(x) sum(x==0))
        DT_limma <- DT_limma[!(DT_limma$missing_value_t >(n_treatment/2)&DT_limma$missing_value_t < n_treatment),]
        DT_limma <- DT_limma[!(DT_limma$missing_value_c >(n_control/2)&DT_limma$missing_value_c <n_control),]
        DT_limma <- DT_limma[DT_limma$missing_value != (n_treatment+n_control),]
        DT_limma[,grep('missing_value',colnames(DT_limma))]=NULL

      } else{
        # Drop rows (protein groups) with missingness in samples
        DT_limma$missing_value= apply(DT_limma, 1, function(x) sum(x==0))
        DT_limma$missing_value_c= apply(DT_limma[,control_samples], 1, function(x) sum(x==0))
        DT_limma$missing_value_t= apply(DT_limma[,treatment_samples], 1, function(x) sum(x==0))
        DT_limma <- DT_limma[DT_limma$missing_value_t %in% c(0,n_treatment),]
        DT_limma <- DT_limma[DT_limma$missing_value_c %in% c(0,n_control),]
        DT_limma <- DT_limma[DT_limma$missing_value != (n_treatment+n_control),]
        DT_limma[,grep('missing_value',colnames(DT_limma))]=NULL

      }



      #design

      group_list <- factor(c(rep('treatment',n_treatment),
                             rep("control",n_control)),
                           levels = c('treatment',"control"))
      limma_design <- model.matrix(~0+group_list)
      colnames(limma_design) <- levels(group_list)
      rownames(limma_design) <- colnames(DT_limma)
      cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = limma_design)

      #limma
      fit <- lmFit(DT_limma, limma_design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2, trend=TRUE)

      result_limma <- topTable(fit2, coef=1,n=Inf)
      result_limma=merge(Log2_DT[,c(name,treatment_samples, control_samples)],result_limma,by.x='Protein_Group',by.y=0)
      ezwrite(result_limma[order(result_limma$adj.P.Val),], DE_dir, paste0(treatment, '_vs_', control, '.tsv'))
      plot_volcano(DT.original = result_limma,
                   out_dir = DE_dir,
                   output_filename = paste0(treatment, '_vs_', control, '_limma'),
                   label_col = 'Genes',
                   lfc_threshold =opt$lfc_threshold,
                   fdr_threshold = opt$fdr_threshold,
                   labelgene =  opt$labelgene)
      enrich_pathway(DT.original = result_limma,
                     treatment = treatment,
                     control = control,
                     outdir = EA_dir,
                     lfc_threshold = lfc_threshold,
                     fdr_threshold = fdr_threshold,
                     enrich_pvalue = enrich_pvalue)

  }
  }
}



##heatmap
plot_heatmap <- function(DT) {
  DT_heatmap <- DT %>%
    select_if(is.numeric)%>%
    # Replace NA values with 0
    mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    filter(rowSums(.) > 0)

  sample_anno <- data.frame(condition = as.factor(sort(names(DT)[sapply(DT, function(x) all(is.numeric(x)))])))
  row.names(sample_anno) <- sample_anno$condition
  sample_anno$condition=gsub("_[0-9]*$","",sample_anno$condition)
  DT_heatmap[is.na(DT_heatmap)]=0
  pheatmap(log2(DT_heatmap+1),
           scale='column',
           show_rownames = F,
           annotation_col = sample_anno,
           cluster_cols = F,
           treeheight_row=0,
           filename=paste0(opt$outdir,'/','heatmap_all.pdf'))
}

plot_heatmap_subset <- function(DT_heatmap,gene) {
  sample_anno <- data.frame(condition = as.factor(sort(colnames(DT_heatmap)[3:ncol(DT_heatmap)])))
  row.names(sample_anno) <- sample_anno$condition
  sample_anno$condition=gsub("_[0-9]*$","",sample_anno$condition)
  DT_heatmap[is.na(DT_heatmap)]=0

  subset_DT=data.frame(DT_heatmap[DT_heatmap$Genes %in% gene$Gene ,])
  subset_DT=subset_DT[!duplicated(subset_DT$Genes),]
  rownames(subset_DT)=subset_DT$Genes

  pheatmap(log2(subset_DT[,3:ncol(subset_DT)]+1),
           scale='column',
           show_rownames = T,
           annotation_col = sample_anno,
           cluster_cols = F,
           treeheight_row=0,
           border_color =NA,
           filename=paste0(opt$outdir,'/','heatmap_label_gene.pdf'))

}

#PPI##########
##APMS ttest
do_t_test_APMS <- function(DT, treatment_samples, control_samples) {

  DT_ttest <- DT
  n_treatment <- length(treatment_samples)
  n_control <- length(control_samples)

  # Retain first three columns plus all treatment and control columns
  DT_ttest <- DT_ttest[,c(colnames(DT_ttest)[1:2], treatment_samples, control_samples), with=F]

  # Convert NA to 0
  DT_ttest[is.na(DT_ttest)] <- 0

  # Drop rows (protein groups) with > 50% missingness in samples
  DT_ttest[,'missing_value':= apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_samples,treatment_samples) ]
  DT_ttest <- DT_ttest[missing_value <= (n_treatment+n_control)/2]
  DT_ttest[, 'missing_value' := NULL]

  # Perform t-test  on treatment and control columns
  t_test <- apply(DT_ttest[,-c(1:2)], 1, function(x){
    a =factor(c(rep('treatment',n_treatment),
                rep("control",n_control)),
              levels = c('treatment',"control"))
    fvalue=var.test(x~a)
    if (!is.na(fvalue$p.value)){
      if (fvalue$p.value > 0.05){
        result <- t.test(x~a, var.equal = T)
      }else{
        result <- t.test(x~a, var.equal = F)
      }
    }
    treatment_estimate <- as.numeric(unlist(result$estimate[1]))
    control_estimate <- as.numeric(unlist(result$estimate[2]))
    return(data.table('P.Value'=result$p.value,
                      'treatment_estimate'=treatment_estimate,
                      'control_estimate'=control_estimate)
    )
  })


  t_test <- rbindlist(t_test)
  t_test <- cbind(DT_ttest[,c(1:2)], t_test)   # add back protein group / gene info cols
  t_test[, logFC := log2(treatment_estimate / (control_estimate+1))]
  t_test[, adj.P.Val := p.adjust(P.Value, method='BH')]
  return(t_test[])
}

get_ppi <- function(protein,DT,lfc_threshold, fdr_threshold) {
    ##canidate genes
    print('canidate genes')
    DT[, 'Group' := 'Others']
    DT[logFC >= lfc_threshold, 'Group' := 'UP']
    DT[logFC <= -lfc_threshold, 'Group' := 'DOWN']
    DT[adj.P.Val >= fdr_threshold, 'Group' := 'Others']
    up_genes=DT[which(DT$logFC>=lfc_threshold&DT$adj.P.Val<=fdr_threshold),]
    ##string db
    print('string db')
    string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=200, network_type="functional", input_directory="")
    protein_map=string_db$mp(protein)
    up_genes <- up_genes %>% string_db$map("Genes",removeUnmappedRows = TRUE)
    print("PPI")
    interactions=string_db$get_interactions(c(protein_map,up_genes$STRING_id))
    pro_interactions_from=interactions[interactions$from==protein_map,2:3]
    colnames(pro_interactions_from)=c('STRING_id','ppi_score')
    pro_interactions_to=interactions[interactions$to==protein_map,c(1,3)]
    colnames(pro_interactions_to)=c('STRING_id','ppi_score')
    pro_interactions=rbind(pro_interactions_from,pro_interactions_to)
    pro_interactions=pro_interactions[!duplicated(pro_interactions, fromLast=TRUE),]
    up_genes=merge(up_genes,pro_interactions,by='STRING_id')
    return(up_genes)
}

plot_PPI_rank <- function(t_test,PPI_score,lfc_threshold, fdr_threshold, output_dir, output_filename) {
    print("Plot PPI rank")
    up_genes=t_test[which(t_test$logFC>=lfc_threshold&t_test$adj.P.Val<=fdr_threshold),]
    up_genes=up_genes[order(up_genes$logFC,decreasing = T),]
    up_genes$rank=1:nrow(up_genes)
    up_genes_plot=merge(up_genes,PPI_score,by="Genes",all=T)
    g <- ggplot(up_genes_plot) +
    geom_point(aes(x=logFC.x, y=rank),color="#999999") +
    geom_point(aes(x=logFC.y, y=rank),color='#E69F00')+
    theme_classic()+
    labs(x="logFC",y='Rank')+
    geom_label_repel(
        data = up_genes_plot,
        aes(x=logFC.y, y=rank,label = Genes),
        size = 5,
        max.overlaps=50,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"))+
      theme(axis.title.x = element_text(size=20, face="bold"),
            axis.title.y = element_text(size=20, face="bold"),
            axis.text.x = element_text(size=14, angle=0),
            axis.text.y = element_text(size=14, angle=0),
            legend.title=element_text(size=20),
            legend.text = element_text(size=14))

    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(output_dir, output_filename), height=10, width=12)
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
}

plot_PPI_ven <- function(protein,t_test,lfc_threshold, fdr_threshold, output_dir, output_filename) {
  print('canidate genes')
  t_test[, 'Group' := 'Others']
  t_test[logFC >= lfc_threshold, 'Group' := 'UP']
  t_test[logFC <= -lfc_threshold, 'Group' := 'DOWN']
  t_test[adj.P.Val >= fdr_threshold, 'Group' := 'Others']
  up_genes=t_test[which(t_test$logFC>=lfc_threshold&t_test$adj.P.Val<=fdr_threshold),]
  ##string db
  print('string db')
  string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=200, network_type="functional", input_directory="")
  protein_map=string_db$mp(protein)
  up_genes <- up_genes %>% string_db$map("Genes",removeUnmappedRows = TRUE)
  neighbors=string_db$get_neighbors(protein_map)
  print("PPI Ven")
  ven=list(
    DE=unique(up_genes$STRING_id),
    String_db=unique(neighbors))
  vd=euler(ven)
  pdf(file=paste0(output_dir, output_filename))
  par(mar = c(2, 5, 5, 2))
  print(plot(vd,
       factor_names = TRUE, labels=list(font=2, cex=1.2),
       counts = TRUE,
       key=TRUE,
       fills = list(fill = c("#fbb4ae", "#b3cde3", "#ccebc5")),
       cex=1,
       edges = FALSE,
       quantities = list(type = c("counts"))))
  dev.off()

}

#mhcflurry########
get_mhcflurry_input=function(hla_typing,dat.long,sample,MHC_dir){
  print(paste0("Generate the MHCflurry input file of ",sample))
  dat.long$Peptide_Sequence=gsub('_|_.[0-9]','',dat.long$Peptide_Sequence)
  sample_allele=hla_typing[hla_typing$sample_name==sample,]
  sample_pep=dat.long[grep(sample,dat.long$Sample),]
  sample_pep$lengh=nchar(sample_pep$Peptide_Sequence)
  sample_pep=sample_pep[sample_pep$lengh<16,]
  sample_allele_peptide=data.frame()

    for (allele in sample_allele$allele) {
        tmp <- data.frame('allele' = allele,
                        'peptide' = unique(sample_pep$Peptide_Sequence))
        sample_allele_peptide = rbind(sample_allele_peptide,tmp)
    }

    cat(paste0('   -> ', MHC_dir,sample, '_allele_peptide.csv', '\n'))
    write.csv(sample_allele_peptide,file = paste0(MHC_dir,sample, '_allele_peptide.csv'),row.names = F)
}

runing_mhcflurry <- function(sample) {
  print(paste0("Runing MHCflurry for ",sample))
  system('mhcflurry-downloads fetch models_class1_presentation')
  command <- paste0("mhcflurry-predict ",MHC_dir,sample, '_allele_peptide.csv'," --out " ,MHC_dir,sample,"_mhc_predictions.csv")
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_predictions.csv", '\n'))
  system(command)

}

plot_mhc_affinity <- function(sample) {
  print(paste0("Plot for the Peptide/MHC I binding affinity prediction of ",sample))
  sample_prediction=fread(paste0(MHC_dir,sample,"_mhc_predictions.csv"))
  sample_prediction=sample_prediction[sample_prediction$mhcflurry_affinity<200
                                      & sample_prediction$mhcflurry_affinity_percentile<2,]
  sample_prediction=sample_prediction[order(sample_prediction$mhcflurry_affinity),]
  sample_prediction$rank=1:nrow(sample_prediction)

  g <- ggplot(sample_prediction, aes(x=-log10(mhcflurry_affinity), y=-log2(mhcflurry_affinity_percentile+1))) +
    geom_point(aes(color = allele))  +
    theme_classic()+
    geom_label_repel(
      data = subset(sample_prediction[1:5,]),
      aes(label = peptide),
      size = 3,nudge_x=0.2,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))+
    xlab('-Log10(MHC Affinity)') +
    ylab('-Log2(MHC Affinity Percentile)')
  ggsave(g, filename=paste0(MHC_dir,sample,"_mhc_affinity_allele.pdf"),width = 6,height = 5)
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_affinity_allele.pdf", '\n'))


  p <- ggplot(sample_prediction, aes(y=-log2(mhcflurry_affinity), x=rank)) +
    geom_point(aes(color = allele))  +
    theme_classic()+
    geom_label_repel(
      data = subset(sample_prediction[1:5,]),
      aes(label = peptide),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))+
    xlab('Rank') +
    ylab('-Log2(MHC Affinity)')
  ggsave(p, filename=paste0(MHC_dir,sample,"_mhc_affinity_rank.pdf"),width = 5,height = 5)
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_affinity_rank.pdf", '\n'))


  label_text=data.frame()
  for (hla in unique(sample_prediction$allele)) {
    tmp=sample_prediction[sample_prediction$allele == hla,]
    tmp=tmp[1:3,]
    label_text=rbind(label_text,tmp)
  }
  f <- ggplot(sample_prediction, aes(x=-log10(mhcflurry_affinity), y=-log2(mhcflurry_affinity_percentile+1))) +
    geom_point(aes(color = allele))  +
    facet_wrap(~ allele, ncol=2,scales="free_y")+
    theme_classic()+
    geom_label_repel(
      data = subset(label_text),
      aes(label = peptide),
      size = 2)+
    xlab('-Log10(MHC Affinity)') +
    ylab('-Log2(MHC Affinity Percentile)')+
    theme(legend.position="top")

  ggsave(f, filename=paste0(MHC_dir,sample,"_mhc_affinity_allele_separated.pdf"),width = 6,height = 8)
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_affinity_allele_separated.pdf", '\n'))

  f <- ggplot(sample_prediction, aes(x=allele, y=log2(mhcflurry_affinity))) +
    geom_dotplot(binaxis='y', stackdir='center',binwidth = 1/50,color='#E69F00',fill="#E69F00")+
    theme_classic()+
    theme(axis.text.x = element_text( angle=90))+
    labs(fill = "",x="HLA alle",y='Log2(mhcflurry_affinity)')
  ggsave(f, filename=paste0(MHC_dir,sample,"_mhc_affinity_allele_dotplot.pdf"),width = 5,height = 5)
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_affinity_allele_dotplot.pdf", '\n'))


  number=sample_prediction[, .N, by=allele]
  p=ggplot(number, aes(x=allele, y=N)) +
    geom_bar(stat="identity", fill="steelblue")+
    theme_classic()+
    labs(fill = "",x="HLA alle",y='Number of Peptide')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(p, filename=paste0(MHC_dir,sample,"_mhc_affinity_peptide_count.pdf"),width = 5,height = 5)
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_affinity_peptide_count.pdf", '\n'))


}

#pSILAC functions#################################
#format data
psilac_standardize_format = function(DT) {
  if ('PG.ProteinGroups' %in% colnames(DT)) {
    setnames(DT, 'PG.ProteinGroups', 'Protein_Group')
  }
  if ('PG.Genes' %in% colnames(DT)) {
    setnames(DT, 'PG.Genes', 'Genes')
  }
  if ('PG.CV' %in% colnames(DT)) {
    setnames(DT, 'PG.CV', 'CV')
  }
  if ('PEP.IsProteinGroupSpecific' %in% colnames(DT)) {
    setnames(DT, 'PEP.IsProteinGroupSpecific', 'ProteinGroupSpecific')
  }
  if ('PEP.StrippedSequence' %in% colnames(DT)) {
    setnames(DT, 'PEP.StrippedSequence', 'Peptide_Sequence')
  }
  if ('EG.PrecursorId' %in% colnames(DT)) {
    setnames(DT, 'EG.PrecursorId', 'Precursor')
  }
  DT=as.data.frame(DT)
  DT=DT[,grep('Protein_Group|Genes|Channel|Peptide_Sequence|Precursor',colnames(DT))]
  # Remove trailing file extensions
  extensions = '.mzML$|.mzml$|.RAW$|.raw$|.dia$|.DIA$|_Intensity|raw.PG.MS2|DIA_|raw.PEP.MS2|raw.EG.TargetReference| \\([^()]*\\)'
  extension_samplenames =  colnames(DT)[colnames(DT) %like% extensions]
  trimmed_samplenames = gsub(extensions, '', extension_samplenames)
  setnames(DT, extension_samplenames, trimmed_samplenames)
  # trim leading [N]
  colnames_out = gsub(pattern="\\[.*\\] ", replacement='', x=colnames(DT))
  setnames(DT, colnames_out)
  DT=data.table(DT)
  return(DT[])
}

##plot protein group or peptide counts of all the sample
plot_silac_pg_counts = function(DT.long, output_dir,intensity_cutoff,name,height,width) {
  pg_count <- DT.long %>%
    filter(Intensity >intensity_cutoff)%>%
    group_by(Sample) %>%  # Group by the Sample column
    summarise(counts = n_distinct(Protein_Group)) %>%  # Count unique Protein_Group
    mutate(Channel = gsub('.*Channel', 'Channel', Sample)) %>%   # Apply gsub using mutate
    mutate(Condition = gsub('[0-9]*\\.Channel', 'Channel', Sample))

  ezwrite(pg_count, output_dir, paste0(name,'pg_count_over_',intensity_cutoff,'.tsv'))

  p=ggplot(pg_count, aes(x=Sample, y=counts,fill=Channel)) +
    geom_bar(stat="identity")+
    theme_classic()+
    labs(fill = "",x="",y='Protein Group Number')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir,name, 'pg_count_over_',intensity_cutoff,'.pdf'),height=height,width = width)
  cat(paste0('   -> ', output_dir,name, 'pg_count_over_',intensity_cutoff,'.pdf', '\n'))

  pg_summary_data <- pg_count %>%
    group_by(Condition) %>%
    summarize(mean = mean(counts), sd = sd(counts)) %>%
    mutate(Channel = gsub('.*Channel', 'Channel', Condition)) %>%
    mutate(Sample = gsub('_Channel.*', '', Condition)) %>%
    arrange(Channel)  # Ascending order

  p=ggplot(pg_summary_data, aes(x=as.factor(Sample), y=mean,fill=Channel)) +
    geom_bar(stat="identity",position=position_dodge())+
    theme_classic()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    labs(fill = "",x="",y='Number of Protein Group')+
    scale_x_discrete(guide = guide_axis(angle = 90))+
    facet_grid(.~ Channel)
  ggsave(plot = p,filename = paste0(output_dir,'/',name,'pg_count_sd_over_',intensity_cutoff,'.pdf'),height=height,width = width)
  cat(paste0('   -> ', output_dir, '/',name,'pg_count_sd_over_',intensity_cutoff,'.pdf', '\n'))
  return(pg_count[])
}




plot_silac_pep_counts = function(DT.long, output_dir,intensity_cutoff,name,height,width) {
  pep_count <- DT.long %>%
    filter(Intensity >intensity_cutoff)%>%
    group_by(Sample) %>%  # Group by the Sample column
    summarise(counts = n_distinct(Peptide_Sequence)) %>%  # Count unique Peptide_Sequence
    mutate(Channel = gsub('.*Channel', 'Channel', Sample))%>%   # Apply gsub using mutate
    mutate(Condition = gsub('[0-9]*\\.Channel', 'Channel', Sample))

  ezwrite(pep_count, output_dir, paste0(name,'pep_count_over_',intensity_cutoff,'.tsv'))

  p=ggplot(pep_count, aes(x=Sample, y=counts,fill=Channel)) +
    geom_bar(stat="identity")+
    theme_classic()+
    labs(fill = "",x="",y='Peptide Number')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir, name,'pep_count_over_',intensity_cutoff,'.pdf'),height=height,width = width)
  cat(paste0('   -> ', output_dir, name,'pep_count_over_',intensity_cutoff,'.pdf', '\n'))

  pep_summary_data <- pep_count %>%
    group_by(Condition) %>%
    summarize(mean = mean(counts), sd = sd(counts)) %>%
    mutate(Channel = gsub('.*Channel', 'Channel', Condition)) %>%
    mutate(Sample = gsub('_Channel.*', '', Condition)) %>%
    arrange(Channel)  # Ascending order

  p=ggplot(pep_summary_data, aes(x=as.factor(Sample), y=mean,fill=Channel)) +
    geom_bar(stat="identity",position=position_dodge())+
    theme_classic()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    labs(fill = "",x="",y='Number of Peptide')+
    scale_x_discrete(guide = guide_axis(angle = 90))+
    facet_grid(.~ Channel)
  ggsave(plot = p,filename = paste0(output_dir,'/',name,'pep_count_sd_over_',intensity_cutoff,'.pdf'),height=5,width = 4)
  cat(paste0('   -> ', output_dir, '/',name,'pep_count_sd_over_',intensity_cutoff,'.pdf', '\n'))
  return(pep_count[])
}

##plot protein group or peptide intensity of all the sample
plot_silac_pg_intensity = function(DT.long, output_dir,height,width) {
  DT.long$Chanel=gsub('.*Channel','Channel',DT.long$Sample)
  DT.long$Condition=gsub('.Channel.*','',DT.long$Sample)
  DT.long$Sample=gsub('_[0-9]+$','',DT.long$Condition)
  p=ggplot(DT.long, aes(y=Chanel, x=log10(Intensity),fill=Chanel)) +
    geom_density_ridges()+
    facet_wrap(. ~ Sample,ncol = 3)+
    theme_classic()+
    theme(axis.text.y = element_blank())
  ggsave(plot = p,filename = paste0(output_dir, 'silac_pg_intensity','_sample.pdf'),width = width,height =height)
  cat(paste0('   -> ', output_dir, 'silac_pg_intensity','_sample.pdf', '\n'))
  p=ggplot(DT.long, aes(y=Condition, x=log10(Intensity),fill=Condition)) +
    geom_density_ridges()+
    facet_wrap(. ~ Chanel)+
    theme_classic()
  ggsave(plot = p,filename = paste0(output_dir, 'silac_pg_intensity','_chanel.pdf'),width = width,height =height)
  cat(paste0('   -> ', output_dir, 'silac_pg_intensity','_chanel.pdf', '\n'))

  p=ggplot(DT.long, aes(x=Condition, y=log10(Intensity))) +
    geom_boxplot(outlier.shape = NA, fill="#67a9cf") +
    theme_classic() +
    labs(fill = "",x="",y='Log10 Protein Intensity') +
    theme(axis.text.x = element_text( angle=90)) +
    geom_boxplot(width=0.1) +
    facet_wrap(. ~ Chanel,nrow=3)+
    geom_hline(color='#ef8a62', linetype='dashed',  aes(yintercept=quantile(log10(DT.long$Intensity), 0.50)))
  ggsave(plot = p,filename = paste0(output_dir, 'silac_pg_intensity','_separated.pdf'),width = height,height =width)
  cat(paste0('   -> ', output_dir, 'silac_pg_intensity','_separated.pdf', '\n'))


}

plot_silac_pep_intensity = function(DT.long, output_dir,name,height,width) {
  DT.long$Chanel=gsub('.*Channel','Channel',DT.long$Sample)
  DT.long$Condition=gsub('.Channel.*','',DT.long$Sample)
  DT.long$Sample=gsub('_[0-9]+$','',DT.long$Condition)
  p=ggplot(DT.long, aes(y=Chanel, x=log10(Intensity),fill=Chanel)) +
    geom_density_ridges()+
    facet_wrap(. ~ Sample,ncol = 3)+
    theme_classic()+
    theme(axis.text.y = element_blank())
  ggsave(plot = p,filename = paste0(output_dir,name, 'pep_intensity_sample.pdf'),width = width,height =height)
  cat(paste0('   -> ', output_dir,name, 'pep_intensity_sample.pdf', '\n'))
  p=ggplot(DT.long, aes(y=Sample, x=log10(Intensity),fill=Sample)) +
    geom_density_ridges()+
    facet_wrap(. ~ Chanel)+
    theme_classic()
  ggsave(plot = p,filename = paste0(output_dir,name, 'pep_intensity_chanel.pdf'),width = width,height =height)
  cat(paste0('   -> ', output_dir, name,'pep_intensity_chanel.pdf', '\n'))
  p=ggplot(DT.long, aes(y=Condition, x=log10(Intensity),fill=Sample)) +
    geom_density_ridges()+
    facet_wrap(. ~ Chanel)+
    theme_classic()+ theme(legend.position="none")
  ggsave(plot = p,filename = paste0(output_dir,name, 'pep_intensity_separated.pdf'),width = width,height =height)
  cat(paste0('   -> ', output_dir, name,'pep_intensity_separated.pdf', '\n'))
}

##Cluster analysis

##plot pca
psilac_pca <- function(DT, output_dir, output_filename) {
  #got pca
  PCA <- get_PCs(DT)
  PCA$components$Channel=gsub('.*Channel','Channel', PCA$components$Sample)
  PCA$components$Condition = gsub('_[0-9]*\\.Channel.*', '', PCA$components$Sample)
  #write pca results
  ezwrite(PCA$components, output_dir, paste0(output_filename,'.tsv'))
  ezwrite(PCA$summary, output_dir, paste0(output_filename,'_summary.tsv'))
  #plot pca results
  if (length(unique(PCA$components$Sample))>8) {
    ncol=2
  }else{
    ncol=1
  }
  p <- ggplot(PCA$components, aes(x = PC1, y = PC2, color = Condition,shape=Channel)) +
    geom_point(size=4) +
    xlab(paste0("PC1","(",PCA$summary$percent[1],"%)")) +
    ylab(paste0("PC2","(",PCA$summary$percent[2],"%)")) +
    theme_classic()+
    guides(color = guide_legend(ncol =  ncol))
  ggsave(p,filename=paste0(output_dir, output_filename,'.pdf'), height = 8,width = 10)
}

psilac_umap <- function(DT, output_dir, output_filename) {
  #got umap
  umap <- get_umap(DT,neighbors = 15)
  umap$Channel=gsub('.*Channel','Channel', umap$Sample)
  umap$Condition = gsub('_[0-9]*\\.Channel.*', '', umap$Sample)

  #write umap results
  ezwrite(umap, output_dir, paste0(output_filename,'.tsv'))
  #plot umap results
  if (length(unique(umap$Sample))>8) {
    ncol=2
  }else{
    ncol=1
  }
  g <- ggplot(umap, aes(x=UMAP1, y=UMAP2, color = Condition,shape=Channel)) +
    geom_point(size=4) +
    theme_classic() +
    guides(color = guide_legend(ncol =  ncol))

  cat(paste0('   -> ', output_dir, output_filename,'.pdf', '\n'))
  ggsave(g, filename=paste0(output_dir, output_filename,'.pdf'), height = 6,width = 8)
}



##Caculate the remaining light peptide by L/(L+H)
#Rate of loss of the light isotope (kLoss) by modeling the relative isotope abundance (RIA)
L_remain_channel=function(DT,output_dir,output_filename){
  DT=data.frame(DT)
  #na as 0
  DT_0=DT
  DT_0[is.na(DT_0)]=0
  #use the Channel to calculate the RIA
  DT_RIA=DT_0[grep('Channel1',colnames(DT_0),value = T)]/
    (DT_0[grep('Channel1',colnames(DT_0),value = T)]+DT_0[grep('Channel3',colnames(DT_0),value = T)])
  DT_RIA=cbind(DT_0[,-grep('Channel|Ratio',colnames(DT_0))],DT_RIA)
  colnames(DT_RIA)=gsub('Channel1','RIA',colnames(DT_RIA))
  #write csv for RIA ratio
  ezwrite(DT_RIA, output_dir, paste0(output_filename,'.tsv'))
  return(DT_RIA)
}

H_remain_channel=function(DT,output_dir,output_filename){
  DT=data.frame(DT)
  #na as 0
  DT_0=DT
  DT_0[is.na(DT_0)]=0
  #use the Channel to calculate the RIA
  DT_RIA=DT_0[grep('Channel3',colnames(DT_0),value = T)]/
    (DT_0[grep('Channel3',colnames(DT_0),value = T)]+DT_0[grep('Channel1',colnames(DT_0),value = T)])
  DT_RIA=cbind(DT_0[,-grep('Channel|Ratio',colnames(DT_0))],DT_RIA)
  colnames(DT_RIA)=gsub('Channel3','RIA',colnames(DT_RIA))
  #write csv for RIA ratio
  ezwrite(DT_RIA, output_dir, paste0(output_filename,'.tsv'))
  ####Protein groups which average %old did not show a continuous decay over time were excluded from further analysis.
  DT_fillter$t1mean=rowMeans(DT_fillter[,grep("D22",colnames(DT_fillter))])
  DT_fillter$t3mean=rowMeans(DT_fillter[,grep("D24",colnames(DT_fillter))])
  DT_fillter$t7mean=rowMeans(DT_fillter[,grep("D28",colnames(DT_fillter))])
  DT_fillter=DT_fillter[which(DT_fillter$t0mean>DT_fillter$t1mean&
                                DT_fillter$t1mean>DT_fillter$t3mean&
                                DT_fillter$t3mean>DT_fillter$t7mean),]
  DT_fillter[grep('mean',colnames(DT_fillter))]=NULL
  return(DT_fillter)
  return(DT_RIA)
}


#phospho functions#################################
##format data
phospho_standardize_format = function(DT) {
  DT=as.data.frame(DT)
  #as number
  DT[,grep('.raw.',colnames(DT))]=as.data.frame(apply(DT[,grep('.raw.',colnames(DT))],2,as.numeric))
  #change the colnames
  DT=data.table(DT)
  if ('PG.ProteinGroups' %in% colnames(DT)) {
    setnames(DT, 'PG.ProteinGroups', 'Protein_Group')
  }
  if ('PG.Genes' %in% colnames(DT)) {
    setnames(DT, 'PG.Genes', 'Genes')
  }
  if ('EG.ProteinPTMLocations' %in% colnames(DT)) {
    setnames(DT, 'EG.ProteinPTMLocations', 'PTM_Location')
  }
  if ('EG.ModifiedSequence' %in% colnames(DT)) {
    setnames(DT, 'EG.ModifiedSequence', 'Modified_Sequence')
  }
  if ('EG.PrecursorId' %in% colnames(DT)) {
    setnames(DT, 'EG.PrecursorId', 'Precursor')
  }
  if ('PEP.StrippedSequence' %in% colnames(DT)) {
    setnames(DT, 'PEP.StrippedSequence', 'Stripped_Sequence')
  }
  DT=as.data.frame(DT)
  DT=DT[,grep('Protein_Group|Genes|PTM_Location|Precursor|Modified_Sequence|Stripped_Sequence|raw',colnames(DT))]
  # Remove trailing file extensions
  extensions = '.mzML$|.mzml$|.RAW$|.raw$|.dia$|.DIA$|_Intensity|.raw.*| \\([^()]*\\)'
  extension_samplenames =  colnames(DT)[colnames(DT) %like% extensions]
  trimmed_samplenames = gsub(extensions, '', extension_samplenames)
  setnames(DT, extension_samplenames, trimmed_samplenames)
  # trim leading [N]
  colnames_out = gsub(pattern="\\[.*\\] ", replacement='', x=colnames(DT))
  setnames(DT, colnames_out)
  return(DT[])
}

plot_phospho_enrichment = function(DT_total,DT_phospho, output_dir,height,width) {
  counts_1=DT_total[, .N, by=Sample]
  counts_1$class='Total'
  counts_2=DT_phospho[, .N, by=Sample]
  counts_2$class='Phospho'
  count=rbind(counts_1,counts_2)
  p=ggplot(count, aes(x=Sample, y=N,fill=class)) +
    geom_bar(stat="identity",position = position_dodge(width = 0))+
    theme_classic()+
    labs(fill = "",x="",y='Number of Peptide')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir,'/','phospho_enrich_peptide_count.pdf'),height=height,width = width)
  percentage=data.frame(sample=counts_1$Sample,
                        percent=counts_2$N/counts_1$N*100)
  ezwrite(percentage,output_dir,'phospho_enrich_percentage.tsv')
  p=ggplot(percentage, aes(x=sample, y=percent)) +
    geom_bar(stat="identity", fill="#67a9cf")+
    theme_classic()+
    labs(fill = "",x="",y='Phospho Percentage %')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir,'/','phospho_enrich_percentage','.pdf'),height=height,width = width)
}

##plot phospho site counts
plot_phospho_site_counts = function(DT_phospho_long, output_dir,height,width) {
  counts=DT_phospho_long[, .N, by=Sample]
  counts$Condition=as.factor(gsub('_[0-9]+$','',counts$Sample))
  ezwrite(counts, output_dir, paste0('phospho_site_count','.tsv'))
  p=ggplot(counts, aes(x=Sample, y=N,fill=Condition)) +
    geom_bar(stat="identity")+
    theme_classic()+
    labs(fill = "",x="",y='Number of Phospho-sites')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir, 'phospho_site_count','.pdf'),height=height,width = width)
  cat(paste0('   -> ', output_dir, 'phospho_site_count','.pdf', '\n'))
  #sd
  summary_data <- counts %>%
    group_by(Condition) %>%
    summarize(mean = mean(N), sd = sd(N)) %>%
    arrange(Condition)
  ezwrite(summary_data, output_dir, paste0('phospho_site_count_condition','.tsv'))
  p=ggplot(summary_data, aes(x=as.factor(Condition), y=mean)) +
    geom_bar(stat="identity",fill="#67a9cf", position=position_dodge())+
    theme_classic()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    labs(fill = "",x="",y='Number of Phospho-sites')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir,'/','phospho_sites_number_sd','.pdf'),height=5,width = 4)
  cat(paste0('   -> ', output_dir, '/','phospho_sites_number_sd','.pdf', '\n'))
}

##plot phospho site distribution
plot_phospho_site_intensity = function(DT_phospho_long, output_dir) {
  n_samples <- length(unique(DT_phospho_long$Sample))
  DT_phospho_long$Condition=gsub('_[0-9]+$','',DT_phospho_long$Sample)
  p=ggplot(DT_phospho_long, aes(y=Sample, x=log10(Intensity),fill=Condition)) +
    geom_density_ridges()+
    theme_classic()+
    theme(axis.text.y = element_blank())
  g <- ggplot(DT_phospho_long, aes(x=Sample, y=log10(Intensity),fill=Condition)) +
    geom_boxplot(outlier.shape = NA) +
    theme_classic() +
    labs(fill = "",x="",y='Log10 Peptide Intensity') +
    theme(axis.text.x = element_text( angle=90)) +
    geom_boxplot(width=0.1) +
    geom_hline(color='#ef8a62', linetype='dashed',  aes(yintercept=quantile(log10(DT_phospho_long$Intensity), 0.50)))
  if (n_samples>50){
    ggsave(plot = g,filename = paste0(output_dir, 'Phospho_site_intensity_boxplot.pdf'),width = n_samples/10,height =6)
    cat(paste0('   -> ', output_dir, 'Phospho_site_intensity_boxplot.pdf', '\n'))
    ggsave(plot = p,filename = paste0(output_dir, 'Phospho_site_intensity.pdf'),width = 6,height =n_samples/5)
    cat(paste0('   -> ', output_dir, 'Phospho_site_intensity.pdf', '\n'))
  }else{
    ggsave(plot = g,filename = paste0(output_dir, 'Phospho_site_intensity_boxplot.pdf'),width = 6,height = 5)
    cat(paste0('   -> ', output_dir, 'Phospho_site_intensity_boxplot.pdf', '\n'))
    ggsave(plot = p,filename = paste0(output_dir, 'Phospho_site_intensity.pdf'),width = 4,height =6)
    cat(paste0('   -> ', output_dir, 'Phospho_site_intensity.pdf', '\n'))
  }
}

##phospho site DE analysis
#ttest
phospho_ttest = function(DT, design_matrix,DE_dir,EA_dir,PTM_EA_dir) {
  name=names(DT)[sapply(DT, function(x) !all(is.numeric(x)))]
  rownames(DT)=DT$PTM
  #comparation
  conditions <- unique(design_matrix$condition)
  for (treatment in conditions) {
    controls=unique(design_matrix$control[design_matrix$condition==treatment])
    for (control in controls) {
      print(paste0(treatment, ' vs ', control, ' Differential Analysis by ttest '))
      treatment_samples=grep(treatment,colnames(DT),value = T)
      control_samples=grep(control,colnames(DT),value = T)
      n_treatment <- length(treatment_samples)
      n_control <- length(control_samples)
      #ttest table
      # Initial data transformation
      DT_ttest <- DT[,c("PTM", treatment_samples, control_samples)]

      # Convert NA to 0 and add missing value calculations
      DT_ttest <- DT_ttest %>%
        # Convert NA to 0
        mutate(across(where(is.numeric), ~ coalesce(., 0))) %>%
        # Calculate missing values
        mutate(
          missing_value = rowSums(select(., -PTM) == 0),  # Total missing values per row, excluding PTM column
          missing_value_c = rowSums(select(., all_of(control_samples)) == 0),  # Missing values in control samples per row
          missing_value_t = rowSums(select(., all_of(treatment_samples)) == 0)  # Missing values in treatment samples per row
        ) %>%
        # Filter missing values
        filter(
          !(missing_value_t > (n_treatment / 2) & missing_value_t < n_treatment),
          !(missing_value_c > (n_control / 2) & missing_value_c < n_control),
          missing_value != (n_treatment + n_control)
        ) %>%
        # Select relevant columns
        select(-missing_value, -missing_value_c, -missing_value_t)%>%
        as.data.frame()
      # Restore row names from PTM column and remove the PTM column
      rownames(DT_ttest) <- DT_ttest$PTM
      DT_ttest <- DT_ttest %>% select(-PTM)
      # Perform t-test  on treatment and control columns
      t_test <- apply(DT_ttest, 1, function(x){
        a =factor(c(rep('treatment',n_treatment),
                    rep("control",n_control)),
                  levels = c('treatment',"control"))
        fvalue=var.test(x~a)
        if (!is.na(fvalue$p.value)){
          if (fvalue$p.value > 0.05){
            result <- t.test(x~a, var.equal = T)
          }else{
            result <- t.test(x~a, var.equal = F)
          }
        }
        treatment_estimate <- as.numeric(unlist(result$estimate[1]))
        control_estimate <- as.numeric(unlist(result$estimate[2]))
        return(data.table('P.Value'=result$p.value,
                          'treatment_estimate'=treatment_estimate,
                          'control_estimate'=control_estimate)
        )
      })
      t_test <- rbindlist(t_test)
      result_ttest <- cbind(DT_ttest, t_test)
      # Merge DT and result_ttest based on PTM column
      result_ttest <- merge(DT[, name, drop=FALSE], result_ttest, by.x='PTM', by.y='row.names')
      # Calculate logFC and adjust P-values using dplyr
      result_ttest <- result_ttest %>%
        mutate(logFC = log2((treatment_estimate+1) / (control_estimate + 1)),
               adj.P.Val = p.adjust(P.Value, method='BH'))
      ezwrite(result_ttest[order(result_ttest$adj.P.Val),], DE_dir, paste0(treatment, '_vs_', control, '.tsv'))
      plot_volcano(DT.original =result_ttest ,
                   lfc_threshold = opt$lfc_threshold,
                   fdr_threshold = opt$fdr_threshold,
                   out_dir = DE_dir,label_col = 'PTM',
                   output_filename =paste0(treatment, '_vs_', control, '_ttest'),
                   labelgene = opt$labelgene)
      enrich_pathway(DT.original = result_ttest,
                     treatment = treatment,
                     control = control,
                     outdir = EA_dir,
                     lfc_threshold = opt$lfc_threshold,
                     fdr_threshold = opt$fdr_threshold,
                     enrich_pvalue = opt$enrich_pvalue)
      PTM_EA(DT.original = result_ttest,
             treatment = treatment,
             control = control,
             outdir = PTM_EA_dir,
             lfc_threshold = opt$lfc_threshold,
             fdr_threshold = opt$fdr_threshold,
             enrich_pvalue = opt$enrich_pvalue)
    }
  }
}


#limma
phospho_limma = function(Log2_DT, design_matrix,DE_dir,EA_dir,PTM_EA_dir) {
  name=names(Log2_DT)[sapply(Log2_DT, function(x) !all(is.numeric(x)))]
  rownames(Log2_DT)=Log2_DT$PTM
  #comparation
  conditions <- unique(design_matrix$condition)
  for (treatment in conditions) {
    controls=unique(design_matrix$control[design_matrix$condition==treatment])
    for ( control in controls) {
      print(paste0(treatment, ' vs ', control, ' Differential Analysis by Limma '))
      treatment_samples=grep(treatment,colnames(Log2_DT),value = T)
      control_samples=grep(control,colnames(Log2_DT),value = T)
      n_treatment <- length(treatment_samples)
      n_control <- length(control_samples)

      #DT_limma
      DT_limma <- Log2_DT[,c(treatment_samples, control_samples)]
      DT_limma <- DT_limma %>%
        # Convert NA to 0
        mutate(across(where(is.numeric), ~ coalesce(., 0))) %>%
        #missing value
        mutate(
          missing_value = rowSums(. == 0),  # Total missing values per row
          missing_value_c = rowSums(select(., all_of(control_samples)) == 0),  # Missing values in control samples per row
          missing_value_t = rowSums(select(., all_of(treatment_samples)) == 0)  # Missing values in treatment samples per row
        ) %>%
        #fillter missing value more than 50%, but keep all missing
        filter(
          !(missing_value_t > (n_treatment / 2) & missing_value_t < n_treatment),
          !(missing_value_c > (n_control / 2) & missing_value_c < n_control),
          missing_value != (n_treatment + n_control)
        ) %>%
        select(-missing_value, -missing_value_c, -missing_value_t)

      #design
      group_list <- factor(c(rep('treatment',n_treatment),
                             rep("control",n_control)),
                           levels = c('treatment',"control"))
      limma_design <- model.matrix(~0+group_list)
      colnames(limma_design) <- levels(group_list)
      rownames(limma_design) <- colnames(DT_limma)
      cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = limma_design)
      #limma
      fit <- lmFit(DT_limma, limma_design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2, trend=TRUE)

      result_limma <- topTable(fit2, coef=1,n=Inf)
      result_limma=merge(Log2_DT[,c(name,treatment_samples, control_samples)],result_limma,by.x='PTM',by.y=0)
      ezwrite(result_limma[order(result_limma$adj.P.Val),], DE_dir, paste0(treatment, '_vs_', control, '_limma.tsv'))
      plot_volcano(DT.original =result_limma ,
                   lfc_threshold = opt$lfc_threshold,
                   fdr_threshold = opt$fdr_threshold,
                   out_dir = DE_dir,label_col = 'PTM',
                   output_filename =paste0(treatment, '_vs_', control, '_limma'),
                   labelgene = opt$labelgene)
      enrich_pathway(DT.original = result_limma,
                     treatment = treatment,
                     control = control,
                     outdir = EA_dir,
                     lfc_threshold = opt$lfc_threshold,
                     fdr_threshold = opt$fdr_threshold,
                     enrich_pvalue = opt$enrich_pvalue)
      PTM_EA(DT.original = result_limma,
             treatment = treatment,
             control = control,
             outdir = PTM_EA_dir,
             lfc_threshold = opt$lfc_threshold,
             fdr_threshold = opt$fdr_threshold,
             enrich_pvalue = opt$enrich_pvalue)
    }
  }
}


##PTM Enrichment Analysis
PTM_EA = function(DT.original, treatment, control, outdir, lfc_threshold, fdr_threshold, enrich_pvalue){
  PTM_db=read.gmt('src/ptm.sig.db.all.uniprot.human.v2.0.0.gmt')
  PTM_db$gene=gsub('-p.*','',PTM_db$gene)
  ##create dir
  dir=paste0(outdir,'/', treatment, '_vs_', control)
  if (!dir.exists(dir)){
    dir.create(dir,recursive = T)
  }
  DT <- DT.original %>%
    mutate(Group = 'Others') %>%
    mutate(Group = ifelse(logFC >= lfc_threshold, 'UP', Group)) %>%
    mutate(Group = ifelse(logFC <= -lfc_threshold, 'DOWN', Group)) %>%
    mutate(Group = ifelse(adj.P.Val >= fdr_threshold, 'Others', Group))%>%
    data.frame()
  DT$pho_site=paste0(DT$Protein_Group,';',gsub('\\(|\\)','',DT$PTM_Location))

  all_PTM=DT$logFC
  names(all_PTM)=DT$pho_site
  all_PTM <- sort(all_PTM, decreasing = TRUE)

  ## up and down regulated PTMs
  up_PTM=DT[which(DT$logFC>=lfc_threshold&DT$adj.P.Val<=fdr_threshold),]
  down_PTM=DT[which(DT$logFC<=(-lfc_threshold)&DT$adj.P.Val<=fdr_threshold),]
  if (nrow(up_PTM)>0){
    print('Processing up_PTM enrichment analysis')
    em <- enricher(up_PTM$pho_site,
                   TERM2GENE=PTM_db,
                   qvalueCutoff =enrich_pvalue,
                   universe = names(all_PTM) )
    em_res=em@result
    write.csv(em_res,paste0(dir,'/',  'up_enrich_res.csv'),row.names = F)
    em_res=em_res[em_res$p.adjust<=enrich_pvalue, ]
    if (nrow(em) > 0) {
      p=dotplot(em, showCategory=10)
      ggsave(paste0(dir,'/','up_enrich_dotplot.pdf'),p,width = 10,height = 10)
    }else{
      print("No enriched pathways found in up PTM Enrichment Analysis analysis")}
  }
  if (nrow(down_PTM)>0){
    print('Processing down_PTM enrichment analysis')
    em <- enricher(down_PTM$pho_site,
                   TERM2GENE=PTM_db,
                   qvalueCutoff =enrich_pvalue,
                   universe = names(all_PTM) )
    em_res=em@result
    write.csv(em_res,paste0(dir,'/',  'down_enrich_res.csv'),row.names = F)
    em_res=em_res[em_res$p.adjust<=enrich_pvalue, ]
    if (nrow(em_res) > 0) {
      p=dotplot(em, showCategory=10)
      ggsave(paste0(dir,'/','down_enrich_dotplot.pdf'),p,width = 10,height = 10)
    }else{
      print("No enriched pathways found in down PTM Enrichment Analysis analysis")}

  }
  print('Processing Gene Set Enrichment Analysis')
  if (nrow(down_PTM)+nrow(up_PTM) >0){
    gse <- GSEA(geneList = all_PTM,
                TERM2GENE = PTM_db,
                pvalueCutoff = enrich_pvalue)
    gse_res=gse@result
    if (nrow(gse_res) > 0) {
      p=dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
      ggsave(paste0(dir,'/','gse_dotplot.pdf'),p,width = 10,height = 10)
      cat(paste0('   -> ', dir,'/','gse_dotplot.pdf', '\n'))
      x2 <- pairwise_termsim(gse)
      p=emapplot(x2)
      ggsave(paste0(dir,'/', 'gse_emap.pdf'),p,width = 10,height = 10)
      cat(paste0('   -> ', dir,'/', 'gse_emap.pdf', '\n'))
      write.csv(gse_res,paste0(dir,'/',  'gse_res.csv'),row.names = F)
    }else{
      print("No enriched pathways found in the Gene Set Enrichment Analysis analysis")}
  }
}



##ksea
ksea=function(dir,outdir,substrates_cutoff,ksea_fdr){
  files=list.files(dir)
  files=grep('.*vs.*tsv',files,value = T)
  all_KSEA_Scores=list()
  for (i in files){
    filename=gsub('_i.*_vs|.tsv','',i)
    KS_outdir=paste0(outdir, filename,'/')
    if (!dir.exists(KS_outdir)){
      dir.create(KS_outdir,recursive = T)
    }
    DE_res=as.data.frame(fread(paste0(dir,'/',i)))
    DT=DE_res %>%
      select(Protein_Group,Genes,Modified_Sequence,PTM_Location,adj.P.Val,logFC)
    colnames(DT)=c("Protein","Gene","Peptide","Residue.Both","p","FC")
    DT$Peptide=gsub('_|\\[.*\\]','',DT$Peptide)
    DT$FC=2^DT$FC
    DT$Residue.Both=gsub('\\(|\\)','',DT$Residue.Both)
    DT=na.omit(DT)
    data("KSData")
    KSLinks= KSEA.KS_table(KSData,
                           DT,
                           NetworKIN = T,
                           NetworKIN.cutoff = 5)
    ezwrite(KSLinks,KS_outdir,'Kinase_Substrate_Links.tsv')
    KSEA_Scores=KSEA.Scores(KSData,
                            DT,
                            NetworKIN = T,
                            NetworKIN.cutoff = 5)
    ezwrite(KSEA_Scores,KS_outdir,'KSEA_Kinase_Scores.tsv')
    drawdata <- KSEA_Scores %>%
      filter(m >= substrates_cutoff) %>%  # Filter rows where m is less than 5,m.cutoff a numeric value between 0 and infinity indicating the min. # of substrates
      arrange(z.score) %>%  # Arrange by z.score in ascending order
      slice(c(1:10, (n() - 9):n())) %>%  # Select the smallest top 10 negative and largest top 10 positive z.scores
      mutate(Group = 'Others',
             Group = if_else(z.score > 0, 'UP', Group),
             Group = if_else(z.score < 0, 'DOWN', Group),
             Group = if_else(p.value >= ksea_fdr, 'Others', Group))


    p=ggplot(drawdata,aes(y=reorder(Kinase.Gene,z.score),x=z.score,fill=Group))+
      scale_fill_manual(breaks = c("DOWN", "Others", "UP"),
                        values = c("#67a9cf", "#969696", "#ef8a62"),
                        labels = c("Neg_Significant", "Not", "Pos_Significant"))+
      theme_classic()+
      geom_bar(stat='identity',position = 'dodge')+
      ylab('Kinase Gene')+
      xlab('Enrichment Score')+
      theme(legend.position="top")
    ggsave(filename = paste0(KS_outdir,'/enrich_plot_pvalue.pdf'),plot = p,width = 8, height = 6)
    p=ggplot(drawdata,aes(x=reorder(Kinase.Gene,Enrichment),y=Enrichment,fill=z.score))+
      theme_bw()+
      geom_bar(stat='identity',position = 'dodge')+
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5))+
      xlab('Kinase Gene')+
      ylab('Enrichment Score')+
      scale_fill_gradient2(high='#D14237',mid='white',low = '#3763B1')
    ggsave(filename = paste0(KS_outdir,'/enrich_plot.pdf'),plot = p,width = 8, height = 6)
    all_KSEA_Scores=c(all_KSEA_Scores,list(KSEA_Scores))
  }
  setwd(paste0(outdir))
  sample.labels=c(gsub('i.*_v','v',files))
  sample.labels=gsub('.tsv','',sample.labels)
  KSEA.Heatmap(all_KSEA_Scores,
               sample.labels=sample.labels,
               stats = 'p.value', #可选p.value或是FDR
               m.cutoff=5,
               p.cutoff=0.05,
               sample.cluster=F)
}

#Soma functions##########
##format data
soma_all_output=function(DT,output_dir){
  anno=getAnalyteInfo(DT)
  DT=data.frame(DT)
  DT_dat=data.frame(DT)%>%
    filter(grepl("Sample", SampleType, ignore.case = TRUE))
  rownames(DT_dat)=DT_dat$SampleId
  DT_dat = DT_dat[, grep('seq\\.', colnames(DT_dat))] %>%
    t() %>%
    as.data.frame()

  Buffer_mean <- DT %>%
    filter(SampleType == "Buffer") %>%
    select(matches("seq\\.", ignore.case = TRUE)) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    t() %>%
    as.data.frame()%>%
    {colnames(.) <- 'Buffer'; .}

  Calibrator_mean <- DT %>%
    filter(SampleType == "Calibrator") %>%
    select(matches("seq\\.", ignore.case = TRUE)) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    t() %>%
    as.data.frame()%>%
    {colnames(.) <- 'Calibrator'; .}

  DT_combined <- cbind(Buffer_mean, Calibrator_mean, DT_dat)
  DT_out=merge(anno[,grep('AptName|UniProt|EntrezGeneSymbol|TargetFullName|Organism|Type',colnames(anno))], DT_combined,by.x='AptName',by.y=0)
  DT_out=DT_out%>%
    rename(Protein_Group= UniProt)%>%
    rename(Genes= EntrezGeneSymbol)
  write.csv(DT_out,paste0(output_dir,'/','Soma_output.csv'),row.names = F)
  cat(paste0('   -> ', output_dir,'/','Soma_output.csv', '\n'))
  return(DT_out)
}

soma_sample_out=function(DT){
  anno=getAnalyteInfo(DT)%>%
    filter(Organism == "Human") %>%
    filter(Type == "Protein")
  DT=as.data.frame(DT)
  DT_dat=DT%>%
    filter(grepl("Sample", SampleType, ignore.case = TRUE))
  rownames(DT_dat)=DT_dat$SampleId
  DT_dat=DT_dat%>%
    select(matches("seq\\.", ignore.case = TRUE))%>%
    t()
  DT_dat=merge(anno[,grep('AptName|UniProt|EntrezGeneSymbol|TargetFullName',colnames(anno))], DT_dat,by.x='AptName',by.y=0,all.x=T)
  DT_dat=DT_dat %>%
    filter(UniProt != "")%>%
    filter(EntrezGeneSymbol != "") %>%
    rename(Protein_Group= UniProt)%>%
    rename(Genes= EntrezGeneSymbol)
  return(DT_dat)
}

Buffer_filter=function(DT,output_dir){
  DT=as.data.frame(DT)
  DT_filter <- DT %>%
    mutate(across(
      .cols = -c(Protein_Group, Genes, Buffer,Calibrator),  # Exclude PG_group, genes, and Buffer
      .fns = ~ ifelse(. < Buffer, NA, .)  # Apply the condition
    ))
  write.csv(DT_filter,paste0(output_dir,'/','Soma_fillter_output.csv'),row.names = F)
  cat(paste0('   -> ', output_dir,'/','Soma_fillter_output.csv', '\n'))
  return(DT_filter)
}

##soma pg count with plot sd
soma_plot_counts = function(counts, output_dir) {
  counts <- counts %>%
    filter(!is.na(SampleGroup)) %>%
    mutate(Sample = factor(Sample, levels = c(
      Sample[SampleGroup == "Buffer"],
      Sample[SampleGroup == "Calibrator"],
      Sample[!SampleGroup %in% c("Buffer", "Calibrator")][order(SampleGroup[!SampleGroup %in% c("Buffer", "Calibrator")])]
    )))
  p=ggplot(counts, aes(x=as.factor(Sample), y=N,fill=SampleGroup)) +
    geom_bar(stat="identity", position=position_dodge())+
    theme_classic()+
    labs(fill = "",x="",y='Number of Protein Group')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir,'/','protein_group_counts_over_buffer.pdf'),height=4,width = 15)
  cat(paste0('   -> ', output_dir, '/','protein_group_counts_over_buffer.pdf', '\n'))
  #sd
  summary_data <- counts %>%
    group_by(SampleGroup) %>%
    summarize(mean = mean(N), sd = sd(N)) %>%
    arrange(SampleGroup)
  p=ggplot(summary_data, aes(x=as.factor(SampleGroup), y=mean,fill=SampleGroup)) +
    geom_bar(stat="identity", position=position_dodge())+
    theme_classic()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    labs(fill = "",x="",y='Number of Protein Group')
  ggsave(plot = p,filename = paste0(output_dir,'/','protein_group_counts_over_buffer_group_sd','.pdf'),height=5,width = 6)
  cat(paste0('   -> ', output_dir, '/','protein_group_counts_over_buffer_group_sd','.pdf', '\n'))
  if ("FLAG" %in% counts$RowCheck) {
    p=ggplot(counts, aes(x = as.factor(Sample), y = N, fill = RowCheck)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      theme_classic() +
      labs(fill = "", x = "", y = "Number of Protein Group") +
      scale_x_discrete(guide = guide_axis(angle = 90))
    ggsave(plot = p,filename = paste0(output_dir,'/','protein_group_counts_over_buffer.pdf'),height=4,width = 15)
    cat(paste0('   -> ', output_dir, '/','protein_group_counts_over_buffer.pdf', '\n'))
  }else{
    counts <- counts %>%
      merge(condition_file, by.x = 'Sample', by.y = 'SampleId') %>%
      #filter(SampleType == "Sample") %>%   # remove control samples
      filter(!is.na(SampleGroup))
    p=ggplot(counts, aes(x=as.factor(Sample), y=N,fill=SampleGroup)) +
      geom_bar(stat="identity", position=position_dodge())+
      theme_classic()+
      labs(fill = "",x="",y='Number of Protein Group')+
      scale_x_discrete(guide = guide_axis(angle = 90))
    ggsave(plot = p,filename = paste0(output_dir,'/','protein_group_counts_over_buffer.pdf'),height=4,width = 15)
    cat(paste0('   -> ', output_dir, '/','protein_group_counts_over_buffer.pdf', '\n'))
    #sd
    summary_data <- counts %>%
      group_by(SampleGroup) %>%
      summarize(mean = mean(N), sd = sd(N)) %>%
      arrange(SampleGroup)%>%
      na.omit()
    p=ggplot(summary_data, aes(x=as.factor(SampleGroup), y=mean,fill=SampleGroup)) +
      geom_bar(stat="identity", position=position_dodge())+
      theme_classic()+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                    position=position_dodge(.9))+
      labs(fill = "",x="",y='Number of Protein Group')
    ggsave(plot = p,filename = paste0(output_dir,'/','protein_group_counts_over_buffer_group_sd','.pdf'),height=5,width = 4)
    cat(paste0('   -> ', output_dir, '/','protein_group_counts_over_buffer_group_sd','.pdf', '\n'))

  }
}

# Calculate the number of NAs in each protein
miss_value_plot=function(DT,outdir){
  na_per_protein <- rowSums(is.na(DT))
  na_per_protein_df <- data.frame(protein = DT$AptName, NA_Count = na_per_protein)
  ezwrite(na_per_protein_df,outdir,'na_per_protein.tsv')
  # Calculate the number of NAs in each sampleumn
  DT =DT %>%
    select_if(is.numeric)
  na_per_sample=sapply(DT, function(x) sum(is.na(x)))
  na_per_sample_df <- data.frame(sampleumn = names(na_per_sample), NA_Count = na_per_sample)
  ezwrite(na_per_sample_df,outdir,'na_per_sample.tsv')

  # Plot NA distribution in sampleumns
  p=ggplot(na_per_sample_df, aes(x = NA_Count)) +
    geom_histogram(color="darkblue", fill="lightblue")+
    theme_classic()+
    geom_vline(aes(xintercept=mean(NA_Count)),
               color="blue", linetype="dashed", size=1)+
    ggtitle("Distribution of NAs in samples")
  ggsave(paste0(outdir,'NA_in_sampels.pdf'),plot = p,width =6, height = 5 )
  cat(paste0('   -> ', outdir, 'NA_in_sampels.pdf', '\n'))
  # Plot NA distribution in proteins
  p=ggplot(na_per_protein_df, aes(x = NA_Count)) +
    geom_histogram(color="darkblue", fill="lightblue")+
    theme_classic()+
    geom_vline(aes(xintercept=mean(NA_Count)),
               color="darkblue", linetype="dashed", size=1)+
    ggtitle("Distribution of NAs in proteins")
  ggsave(paste0(outdir,'NA_in_proteins.pdf'),plot = p,width =6, height = 5 )
  cat(paste0('   -> ', outdir, 'NA_in_proteins.pdf', '\n'))
}

##soma pca plot
<<<<<<< HEAD
soma_plot_PCs <- function(adat,DT,condition_file, output_dir) {
  #pca all
  cluster_data <- adat %>%
    select(matches("seq\\.", ignore.case = TRUE)) %>%
    mutate(across(where(is.numeric), ~ log2(.)))
  pca=prcomp(cluster_data, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
  summary <- as.data.table(t(summary(pca)$importance), keep.rownames=T)
  setnames(summary, c('component','stdv','percent','cumulative'))
  summary$percent=round(summary$percent*100, digits = 2)
  pca_df = as.data.frame(pca$x)[,1:5]
  pca_df=merge(pca_df,adat,by=0)
  p=ggplot(pca_df, aes(x = PC1, y = PC2, color = SampleType)) +
    geom_point(size=4) +
    xlab(paste0("PC1","(",summary$percent[1],"%)")) +
    ylab(paste0("PC2","(",summary$percent[2],"%)")) +
    theme_classic()
  ggsave(p,filename=paste0(output_dir, 'pca_all.pdf'), height = 4,width = 5)
  #pca sample
=======
soma_get_PCs=function(DT,condition_file,cluster_dir){
  out <- list()
  ##cluster data(na=0)
>>>>>>> 1e7a12a (functions)
  cluster_data <- DT %>%
    select_if(is.numeric) %>%
    mutate(across(everything(), ~ log2(. )))
  pca_data=t(cluster_data)
  pca=prcomp(pca_data, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
  summary <- as.data.table(t(summary(pca)$importance), keep.rownames=T)
  setnames(summary, c('component','stdv','percent','cumulative'))
  summary$percent=round(summary$percent*100, digits = 2)
  pca_df = as.data.frame(pca$x)[,1:5]
  condition_file=condition_file%>%
    filter(!is.na(SampleGroup))
  pca_df <- merge(pca_df, condition_file[,grep('SampleId|SampleGroup|RowCheck',colnames(condition_file))],
                  by.x = 0, by.y = 'SampleId',all=T)
  ezwrite(pca_df, output_dir, 'PCA.tsv')
  ezwrite(summary, output_dir, 'PCA_summary.tsv')
  p=ggplot(pca_df, aes(x = PC1, y = PC2, color =SampleGroup )) +
    geom_point(size=4) +
    xlab(paste0("PC1","(",summary$percent[1],"%)")) +
    ylab(paste0("PC2","(",summary$percent[2],"%)")) +
    theme_classic()
  ggsave(p,filename=paste0(output_dir, 'pca.pdf'), height = 4,width = 5)
  cat(paste0('   -> ', output_dir, '/','pca.pdf', '\n'))
  if ("FLAG" %in% pca_df$RowCheck) {
    p=ggplot(pca_df, aes(x = PC1, y = PC2, color =RowCheck )) +
      geom_point(size=4) +
      xlab(paste0("PC1","(",summary$percent[1],"%)")) +
      ylab(paste0("PC2","(",summary$percent[2],"%)")) +
      theme_classic()
    ggsave(p,filename=paste0(output_dir, 'pca_flag.pdf'), height = 4,width = 5)
    cat(paste0('   -> ', output_dir, '/','pca_flag.pdf', '\n'))
  }
}

#soma umap
soma_get_umap <- function(DT, condition_file) {
  ##cluster data(na=0)
  cluster_data <- DT %>%
    select_if(is.numeric)%>%
    # Replace NA values with 0
    mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))  %>%
    # Filter rows where the sum of numeric values is greater than 0
    filter(rowSums(.) > 0)

  # Apply log2 transformation
  log2_cluster_data <- cluster_data %>%
    mutate(across(everything(), ~ log2(. + 1)))

  set.seed(100)
  DT.umap <- umap(t(log2_cluster_data))
  DT.out <- as.data.table(DT.umap$layout, keep.rownames=TRUE)
  setnames(DT.out, c('Sample', 'UMAP1', 'UMAP2'))
  if (!all(is.na(condition_file$SampleGroup))) {
    condition_file <- condition_file %>%
      filter(!is.na(SampleGroup))%>%
      select(SampleId,SampleGroup)
  }
  DT.out=merge(DT.out,condition_file,by.x='Sample',by.y='SampleId')
  return(DT.out[])
}

soma_plot_umap <- function(DT, output_dir, output_filename) {
  if (all(is.na(DT$SampleGroup))){
    g <- ggplot(DT, aes(x=UMAP1, y=UMAP2, color=SampleType)) +
      geom_point(size=4) +
      theme_classic()

  }else{
  g <- ggplot(DT, aes(x=UMAP1, y=UMAP2, color=SampleGroup)) +
    geom_point(size=4) +
    theme_classic()
  }
  cat(paste0('   -> ', output_dir, output_filename, '\n'))
  ggsave(g, filename=paste0(output_dir, output_filename), height = 4,width = 5)
}

#soma ttest
soma_test=function(DT.original, out_dir, output_filename ,lfc_threshold, fdr_threshold, labelgene) {
  options(ggrepel.max.overlaps = Inf)
  DT <- DT.original
  # Set initial group to 'Others' and update based on thresholds
  DT <- DT %>%
    mutate(Group = 'Others',
           Group = if_else(logFC >= lfc_threshold, 'UP', Group),
           Group = if_else(logFC <= -lfc_threshold, 'DOWN', Group),
           Group = if_else(adj.P.Val >= fdr_threshold, 'Others', Group),
           labeltext = '')
  # Select top 5 genes for UP and DOWN groups
  top5_gene <- DT %>%
    filter(Group == 'UP') %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 5) %>%
    bind_rows(
      DT %>%
        filter(Group == 'DOWN') %>%
        arrange(logFC) %>%
        slice_head(n = 5)
    )
  # Label top 5 genes using the specified label column
  DT <- DT %>%
    mutate(labeltext = if_else(AptName %in% top5_gene$AptName, Genes, ""))

  g <- ggplot(DT, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Group)) +
    scale_color_manual(breaks = c("DOWN", "Others", "UP"),
                       values = c("#67a9cf", "#969696", "#ef8a62")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    geom_label_repel(
      data = subset(DT),
      aes(label = labeltext),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    geom_vline(xintercept = lfc_threshold, linetype = "dashed") +
    geom_vline(xintercept = -lfc_threshold, linetype = "dashed") +
    theme_classic()

  cat(paste0(out_dir, output_filename,'_vocanol.pdf'))
  ggsave(g, filename = paste0(out_dir, output_filename,'_vocanol.pdf'), width = 8, height = 8)
  cat(paste0('   -> ', out_dir, output_filename,'_vocanol.pdf', '\n'))
}

soma_plot_volcano=function(DT.original, out_dir, output_filename ,lfc_threshold, fdr_threshold, labelgene) {
  options(ggrepel.max.overlaps = Inf)
  DT <- DT.original
  # Set initial group to 'Others' and update based on thresholds
  DT <- DT %>%
    mutate(Group = 'Others',
           Group = if_else(logFC >= lfc_threshold, 'UP', Group),
           Group = if_else(logFC <= -lfc_threshold, 'DOWN', Group),
           Group = if_else(adj.P.Val >= fdr_threshold, 'Others', Group),
           labeltext = '')
  # Select top 5 genes for UP and DOWN groups
  top5_gene <- DT %>%
    filter(Group == 'UP') %>%
    arrange(desc(logFC)) %>%
    slice_head(n = 5) %>%
    bind_rows(
      DT %>%
        filter(Group == 'DOWN') %>%
        arrange(logFC) %>%
        slice_head(n = 5)
    )
  # Label top 5 genes using the specified label column
  DT <- DT %>%
    mutate(labeltext = if_else(AptName %in% top5_gene$AptName, Genes, ""))

  g <- ggplot(DT, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Group)) +
    scale_color_manual(breaks = c("DOWN", "Others", "UP"),
                       values = c("#67a9cf", "#969696", "#ef8a62")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom") +
    geom_label_repel(
      data = subset(DT),
      aes(label = labeltext),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")) +
    geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    geom_vline(xintercept = lfc_threshold, linetype = "dashed") +
    geom_vline(xintercept = -lfc_threshold, linetype = "dashed") +
    theme_classic()

  cat(paste0(out_dir, output_filename,'_vocanol.pdf'))
  ggsave(g, filename = paste0(out_dir, output_filename,'_vocanol.pdf'), width = 8, height = 8)
  cat(paste0('   -> ', out_dir, output_filename,'_vocanol.pdf', '\n'))
}

soma_plot_heatmap <- function(DT_heatmap,condition_file ){
  DT_heatmap=DT_heatmap %>%
    select_if(is.numeric)
  if (all(is.na(condition$SampleGroup))) {
    pheatmap(log10(DT_heatmap),
             scale='column',
             show_rownames = F,
             show_colnames = F,
             cluster_cols = T,
             treeheight_row=0,
             filename=paste0(opt$outdir,'/','heatmap_all.pdf'))
  }else{
    condition_file <- condition_file %>%
      filter(!is.na(SampleGroup))
    sample_anno=data.frame(condition = as.factor(condition_file$SampleGroup))
    row.names(sample_anno) <- condition_file$SampleId
    pheatmap(log10(DT_heatmap),
             scale='column',
             show_rownames = F,
             show_colnames = F,
             annotation_col = sample_anno,
             cluster_cols = T,
             treeheight_row=0,
             filename=paste0(opt$outdir,'/','heatmap_all.pdf'))
  }
}

soma_plot_heatmap_subset <- function(DT_heatmap,condition_file,gene_subset ){
  condition_file <- condition_file %>%
    filter(!is.na(SampleGroup))
  sample_anno=data.frame(condition = as.factor(condition_file$SampleGroup))
  row.names(sample_anno) <- condition_file$SampleId
  gene_subset=fread(gene_subset)
  # Check the existing column names
  existing_names <- names(gene_subset)
  # Create a vector of names to look for
  potential_names <- c("gene","genes", "Gene", "Genes", "GENE", "GENES")
  # Find matching columns
  columns_to_rename <- existing_names[existing_names %in% potential_names]

  # Rename matching columns to 'gene'
  if (length(columns_to_rename) > 0) {
    gene_subset <- gene_subset %>%
      rename_with(~ "gene", .cols = all_of(columns_to_rename)) %>%
      select(gene, everything())  # Ensure 'gene' is the first column
  }
  DT_heatmap <- DT_heatmap %>%
    filter(Genes %in% gene_subset$gene) %>%
    select_if(is.numeric)
  pheatmap(log10(DT_heatmap),
           scale='column',
           show_rownames = F,
           show_colnames = F,
           annotation_col = sample_anno,
           cluster_cols = T,
           cluster_rows = T,
           treeheight_row=0,
           filename=paste0(opt$outdir,'/','heatmap_subset.pdf'))
}

soma_box_plot=function(soma_adat,treatment,control,out_dir,gene,t_test,levels,color){
  if (!is.null(gene)){
    target_map <- t_test %>%
      filter(Genes %in% gene) %>%  # Filter by relevant genes
      group_by(Genes) %>%                          # Group by Genes
      slice_min(order_by = adj.P.Val, n = 1) %>%  # Select the row with the smallest adj_pvalue
      ungroup()                                  # Ungroup to perform further operations
  }else{
    target_map <- t_test %>%
      group_by(Genes) %>%                             # Group by Genes
      slice_min(order_by = adj.P.Val, n = 1) %>%      # Select the row with the smallest adj.P.Val
      ungroup() %>%                                   # Ungroup to perform further operations
      slice_min(order_by = adj.P.Val, n = 20)         # Select the top 20 rows based on the smallest adj.P.Val
  }

  plot_tbl <- as.data.frame(soma_adat)  %>%            # plot non-center/scale data
    filter(SampleType == "Sample")  %>%    # rm control samples
    filter(!is.na(SampleGroup))  %>%  # rm NAs if present
    filter(SampleGroup %in% c(treatment,control)) %>%
    select(SampleGroup, all_of(target_map$AptName))%>%
    mutate(across(where(is.numeric), log10))  %>%
    melt(id.vars = "SampleGroup", variable.name = "AptName", value.name = "RFU")%>%
    left_join(target_map, by = "AptName")
  library(ggpubr)
  if (!is.null(levels)){
    plot_tbl$SampleGroup <- factor(plot_tbl$SampleGroup,   levels= levels)
  }
  # Plot with or without custom colors
  plot <- plot_tbl |>
    ggplot(aes(x = SampleGroup, y = RFU, fill = SampleGroup)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_jitter(shape = 16, width = 0.1, alpha = 0.5) +
    facet_wrap(~ Genes, ncol = 5, scales = "free_y") +
    ggtitle("Boxplots of Top Analytes by t-test") +
    labs(y = "log10(RFU)") +
    theme(plot.title = element_text(size = 21, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.position = "top") +
    stat_compare_means(method = "t.test", label = "p.signif",
                       aes(label = ..p.signif..),
                       label.x.npc = "center", vjust = 1)

  # Apply custom colors if provided
  if (!is.null(color)) {
    plot <- plot + scale_fill_manual(values = color)
    ggsave(paste0(out_dir,out,"_soma_ttest_top_boxplot.pdf"),plot,width =10, height = 6 )
    cat(paste0('   -> ', out_dir,"soma_ttest_top_boxplot.pdf", '\n'))
  } else {
    ggsave(paste0(out_dir,treatment, '_vs_', control,"_soma_ttest_top_boxplot.pdf"),plot,width =10, height = 6 )
    cat(paste0('   -> ', out_dir,treatment, '_vs_', control,"_soma_ttest_top_boxplot.pdf", '\n'))
  }

}






#Olink functions##########
#plot_intensities
plot_olink_intensities =function(DT, output_dir) {
  n_samples <- length(unique(DT$SampleID))
  g <- ggplot(DT, aes(x=SampleID, y=(NPX))) +
    geom_boxplot(outlier.shape = NA,fill='#67a9cf') +
    theme_classic() +
    labs(fill = "",x="",y='NPX') +
    theme(axis.text.x = element_text( angle=90)) +
    geom_boxplot(width=0.1) +
    geom_hline(color='#ef8a62', linetype='dashed',  aes(yintercept=quantile((DT$NPX), 0.50, na.rm = TRUE)))

  ggsave(plot = g,filename = paste0(output_dir, 'all_intensities.pdf'),width = 10,height =6)
  cat(paste0('   -> ', output_dir, 'all_intensities.pdf', '\n'))
  if ("Warning" %in% DT$SampleQC) {
    p=ggplot(DT, aes(x=SampleID, y=(NPX), fill = SampleQC)) +
      geom_boxplot(outlier.shape = NA) +
      theme_classic() +
      labs(fill = "",x="",y='NPX') +
      theme(axis.text.x = element_text( angle=90)) +
      geom_boxplot(width=0.1) +
      geom_hline(color='#ef8a62', linetype='dashed',  aes(yintercept=quantile((DT$NPX), 0.50, na.rm = TRUE)))
    ggsave(plot = p,filename = paste0(output_dir,'/','all_intensities_flag.pdf'),height=6,width = 10)
    cat(paste0('   -> ', output_dir, '/','all_intensities_flag.pdf', '\n'))
  }
}
#NC to filtter
data_fillter_nc=function(npx,output_dir){
  NEGATIVE_CONTROL <- npx %>%
    filter(AssayType == "assay")%>%
    dplyr::select(SampleID,SampleType,WellID,PlateID, NPX,OlinkID) |>
    tidyr::pivot_wider(names_from = OlinkID, values_from = NPX)%>%
    filter(SampleType == "NEGATIVE_CONTROL") %>%
    select(matches("OID", ignore.case = TRUE)) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    t() %>%
    as.data.frame()%>%
    {colnames(.) <- 'NEGATIVE_CONTROL'; .}
  data_sample <- npx |>
    dplyr::filter(SampleType == "SAMPLE") |>
    dplyr::filter(AssayType == "assay") |>
    dplyr::select(SampleID, UniProt,Assay, OlinkID, NPX) |>
    tidyr::pivot_wider(names_from = SampleID, values_from = NPX) |>
    dplyr::rename(Protein_Group = UniProt, Genes = Assay)
  dat_all=merge(data_sample,NEGATIVE_CONTROL,by.x='OlinkID',by.y=0)
  DT_filter <- dat_all %>%
    mutate(across(
      .cols = -c(Protein_Group, Genes, OlinkID,NEGATIVE_CONTROL),  # Exclude PG_group, genes, and Buffer
      .fns = ~ ifelse(. < NEGATIVE_CONTROL, NA, .)  # Apply the condition
    ))
  DT_filter_long = melt(DT_filter) %>%
    filter(!is.na(value)) %>%
    data.table()
  pgcounts=DT_filter_long[, .N, by=variable]%>%
    filter(variable != "NEGATIVE_CONTROL")
  p=ggplot(pgcounts, aes(x=variable, y=N)) +
    geom_bar(stat="identity", fill="#67a9cf")+
    theme_classic()+
    labs(fill = "",x="",y='Number of Protein Groups')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(filename = paste0(output_dir, 'protein_group_counts_over_nc.pdf'),plot = p,width = 8,height = 6)
  cat(paste0('   -> ', output_dir, 'protein_group_counts_over_nc.pdf', '\n'))
}
#LOD
lod_fillter=function(npx,output_dir){
  if (length(unique(grep('NEG_CTRL',value = T,npx$SampleID)))>10) {
    npx_lod=olink_lod(npx, lod_method = "NCLOD")
  } else{
    npx_lod=olink_lod(npx, lod_file_path = '/Users/liz36/Documents/Brain_Sample/olink/Olink_lod.csv', lod_method = "FixedLOD")
  }
  npx_lod_fillter=npx_lod%>%
    mutate(value = ifelse(PCNormalizedNPX < LOD, NA, PCNormalizedNPX)) %>%
    filter(!is.na(value))
  npx_lod_fillter=data.table(npx_lod_fillter)
  pgcounts=npx_lod_fillter[, .N, by=SampleID][!grepl('CTRL', SampleID)]
  p=ggplot(pgcounts, aes(x=SampleID, y=N)) +
    geom_bar(stat="identity", fill="#67a9cf")+
    theme_classic()+
    labs(fill = "",x="",y='Number of Protein Groups')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(filename = paste0(output_dir, 'protein_group_counts_over_lod.pdf'),plot = p,width = 8,height = 6)
  cat(paste0('   -> ', output_dir, 'protein_group_counts_over_old.pdf', '\n'))
}
#plot_pca
olink_plot_PCs <- function(DT,DT_all,condition_file, output_dir) {
  #pca all
  cluster_data <- DT_all %>%
    select(where(is.numeric)) %>%
    filter(complete.cases(.))
  pca=prcomp(cluster_data, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
  summary <- as.data.table(t(summary(pca)$importance), keep.rownames=T)
  setnames(summary, c('component','stdv','percent','cumulative'))
  summary$percent=round(summary$percent*100, digits = 2)
  pca_df = as.data.frame(pca$x)[,1:5]
  rownames(pca_df)=DT_all$SampleID
  condition=DT_all %>%select(!where(is.numeric))
  pca_df=merge(pca_df,condition,by.x=0,by.y='SampleID',all=T)
  p=ggplot(pca_df, aes(x = PC1, y = PC2, color = SampleType)) +
    geom_point(size=4) +
    xlab(paste0("PC1","(",summary$percent[1],"%)")) +
    ylab(paste0("PC2","(",summary$percent[2],"%)")) +
    theme_classic()
  ggsave(p,filename=paste0(output_dir, 'pca_all.pdf'), height = 4,width = 5)
  #pca sample
  pca_data <- DT %>%
    select_if(is.numeric)
  pca=prcomp(pca_data, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
  summary <- as.data.table(t(summary(pca)$importance), keep.rownames=T)
  setnames(summary, c('component','stdv','percent','cumulative'))
  summary$percent=round(summary$percent*100, digits = 2)
  pca_df = as.data.frame(pca$x)[,1:5]
  rownames(pca_df)=DT$SampleID
  condition_file=condition_file%>%
    filter(!is.na(Sample_Group))
  pca_df <- merge(pca_df, condition_file, by.x = 0, by.y = 'Sample_ID',all=T)
  ezwrite(pca_df, output_dir, 'PCA.tsv')
  ezwrite(summary, output_dir, 'PCA_summary.tsv')
  p=ggplot(pca_df, aes(x = PC1, y = PC2, color =Sample_Group )) +
    geom_point(size=4) +
    xlab(paste0("PC1","(",summary$percent[1],"%)")) +
    ylab(paste0("PC2","(",summary$percent[2],"%)")) +
    theme_classic()
  ggsave(p,filename=paste0(output_dir, 'pca.pdf'), height = 4,width = 5)
  cat(paste0('   -> ', output_dir, '/','pca.pdf', '\n'))
  condition=DT %>%select(!where(is.numeric))
  pca_df=merge(pca_df,condition,by.x='Row.names',by.y='SampleID',all=T)
  if ("FLAG" %in% pca_df$SampleQC) {
    p=ggplot(pca_df, aes(x = PC1, y = PC2, color =SampleQC )) +
      geom_point(size=4) +
      xlab(paste0("PC1","(",summary$percent[1],"%)")) +
      ylab(paste0("PC2","(",summary$percent[2],"%)")) +
      theme_classic()
    ggsave(p,filename=paste0(output_dir, 'pca_flag.pdf'), height = 4,width = 5)
    cat(paste0('   -> ', output_dir, '/','pca_flag.pdf', '\n'))
  }
}
#
olink_ttest =function(DT, treatment_samples, control_samples ){
  n_treatment <- length(treatment_samples)
  n_control <- length(control_samples)
  name=names(DT)[sapply(DT, function(x) !all(is.numeric(x)))]
  DT_ttest <- DT[,colnames(DT) %in% c(treatment_samples, control_samples)]
  rownames(DT_ttest) <- DT$OlinkID
  t_test <- apply(DT_ttest, 1, function(x){
    a =factor(c(rep('treatment',n_treatment),
                rep("control",n_control)),
              levels = c('treatment',"control"))
    fvalue=var.test(x~a)
    if (!is.na(fvalue$p.value)){
      if (fvalue$p.value > 0.05){
        result <- t.test(x~a, var.equal = T)
      }else{
        result <- t.test(x~a, var.equal = F)
      }
    }
    treatment_estimate <- as.numeric(unlist(result$estimate[1]))
    control_estimate <- as.numeric(unlist(result$estimate[2]))
    return(data.table('P.Value'=result$p.value,
                      'treatment_estimate'=treatment_estimate,
                      'control_estimate'=control_estimate)
    )
  })
  t_test <- rbindlist(t_test)
  t_test$OlinkID=DT$OlinkID
  result_ttest <- merge(DT[, name], t_test, by='OlinkID')
  # Calculate logFC and adjust P-values using dplyr
  result_ttest <- result_ttest %>%
    mutate(FC = (treatment_estimate) / (control_estimate),
           adj.P.Val = p.adjust(P.Value, method='BH'))
  ezwrite(result_ttest[order(result_ttest$adj.P.Val),], DE_dir, paste0(treatment, '_vs_', control, '.tsv'))


}
