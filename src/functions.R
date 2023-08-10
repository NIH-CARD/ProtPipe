
#### FUNCTIONS #####################################################################################
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

standardize_format <- function(DT.original) {
    # Accepts an input protein group intensity data.table, whether spectronaut or DIA-NN format,
    # and restructures into one consistent style for downstream processing
    DT <- copy(DT.original)
    if("Protein.Ids" %in% colnames(DT)) {
      print("DIAnn input")
        DT[, 'Protein.Ids' := NULL]
        DT[, 'Protein.Names' := NULL]
        DT[, 'First.Protein.Description' := NULL]
        setnames(DT, 'Protein.Group', 'Protein_Group')
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
    colnames_out <- gsub(pattern="\\..*\\.PG\\.Quantity", replacement='', x=colnames_out)   # remove suffix
    return(colnames_out)
}

melt_intensity_table <- function(DT) {
    # Converts intensity data.table to long format
    # info_cols <- c('Protein_Group', 'Genes', 'First_Protein_Description')
    DT.long <- melt(DT, 
    measure.vars=colnames(DT[,-c(1:2)]),
    variable.name='Sample',
    value.name='Intensity')
    DT.long=data.table(DT.long)
    return(DT.long)
}


plot_pg_counts <- function(DT, output_dir, output_filename) {
  n_samples <- nrow(DT)
  if (n_samples > 20) {
    p=ggplot(DT, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="steelblue")+
      theme_classic()+
      labs(fill = "",x="Sample",y='Number of Protein Groups')+
      scale_x_discrete(guide = guide_axis(angle = 90))
  }
  if (n_samples < 20) {
    p=ggplot(DT, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="steelblue")+
      theme_classic()+
      labs(fill = "",x="Sample",y='Number of Protein Groups')+
      scale_x_discrete(guide = guide_axis(angle = 90))+ 
      geom_text(aes(label=N, y=N + (0.05*max(pgcounts$N))))
  }
  
  if (n_samples>50){
    ggsave(plot = p,filename = paste0(output_dir, output_filename),width = n_samples/10,height = 6)
  }else {
    ggsave(plot = p,filename = paste0(output_dir, output_filename),width = 8,height = 6)
    } 
  }

plot_pep_counts <- function(DT, output_dir, output_filename) {
  n_samples <- nrow(DT)
  if (n_samples > 20) {
    p=ggplot(DT, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="steelblue")+
      theme_classic()+
      labs(fill = "",x="Sample",y='Number of Peptides')+
      scale_x_discrete(guide = guide_axis(angle = 90))
  }
  if (n_samples < 20) {
    p=ggplot(DT, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="steelblue")+
      theme_classic()+
      labs(fill = "",x="Sample",y='Number of Peptides')+
      scale_x_discrete(guide = guide_axis(angle = 90))+ 
      geom_text(aes(label=N, y=N + (0.05*max(pep_counts$N))))
  }
  
  if (n_samples>50){
    ggsave(plot = p,filename = paste0(output_dir, output_filename),width = n_samples/10,height = 6)
  }else {
    ggsave(plot = p,filename = paste0(output_dir, output_filename),width = 8,height = 6)
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
}


plot_density <- function(DT.original, output_dir, output_filename) {
    # Currently UNUSED as the beeswarm function serves the purpose well
    DT <- copy(DT.original)
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
        geom_vline(xintercept=intensity_median, color='red') +
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
}


shift_normalize_intensity <- function(DT.original) {
    DT <- copy(DT.original)
    # Get global median of intensity values
    global_median <- median(DT[, Intensity])
    DT[, 'sample_median' := median(Intensity), by=Sample]
    DT[, 'global_median' := global_median]
    DT[, 'NormInt' := Intensity - (sample_median - global_median)]
    DT[, c('Intensity', 'sample_median', 'global_median') := NULL]
    setnames(DT, 'NormInt', 'Intensity')
    return(DT[])
}


scale_normalize_intensity <- function(DT.original) {
    DT <- copy(DT.original)
    # Get global median of intensity values
    global_median <- median(DT[, Intensity])
    DT[, 'sample_median' := median(Intensity), by=Sample]
    DT[, 'global_median' := global_median]
    DT[, 'NormInt' := Intensity * (global_median / sample_median)]
    DT[, c('Intensity', 'sample_median', 'global_median') := NULL]
    setnames(DT, 'NormInt', 'Intensity')
    return(DT[])
}


plot_pg_intensities <- function(DT, output_dir, output_filename, plot_title) {
    n_samples <- length(unique(DT$Sample))
    g <- ggplot(DT, aes(x=Sample, y=log10(Intensity))) + 
        geom_boxplot(outlier.shape = NA, fill="steelblue") +
        theme_classic() +
        labs(fill = "",x="",y='Log10 Protein Intensity') +
        theme(axis.text.x = element_text( angle=90)) +
        geom_boxplot(width=0.1) +
        geom_hline(color='red', linetype='dashed',  aes(yintercept=quantile(log10(DT$Intensity), 0.50)))
  
    if (n_samples>50){
        ggsave(plot = g,filename = paste0(output_dir, output_filename),width = n_samples/10,height =6)
    }else{
        ggsave(plot = g,filename = paste0(output_dir, output_filename),width = 8,height = 6)
    }
}


plot_pep_intensities <- function(DT, output_dir, output_filename, plot_title) {
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
    }else{
        ggsave(plot = g,filename = paste0(output_dir, output_filename),width = 8,height = 6)
    }
}



convert_log_to_raw <- function(DT.original, log_base) {
    DT <- copy(DT.original)
    samplenames <- colnames(DT[,-c(1:2)])
    DT[, (samplenames) := lapply(.SD, function(x) log_base^x), .SDcols=samplenames]
    return(DT[])
}



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
        units='cm')
  
} 


get_spearman <- function(DT.original) {
  DT <- copy(DT.original)
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
    DT <- copy(DT.original)
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


get_PCs <- function(DT) {
  out <- list()
  ##cluster data(na=0)
  cluster_data=DT[,-c(1:2)]
  cluster_data[is.na(cluster_data)]=0
  cluster_data=cluster_data[which(rowSums(cluster_data)>0),]
  log2_cluster_data=log2(cluster_data+1)
  
  ##PCA and plot
  pca_data=t(log2_cluster_data)
  pca=prcomp(pca_data, center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
  out$summary <- as.data.table(t(summary(pca)$importance), keep.rownames=T)
  setnames(out$summary, c('component','stdv','percent','cumulative'))
  out$summary$percent=round(out$summary$percent*100, digits = 2)
  pca_df = as.data.frame(pca$x)
  pca_df$Condition=gsub('_[1-9]*$','',rownames(pca_df))
  out$components <- pca_df
  return(out)
}

plot_PCs <- function(PCA, output_dir, output_filename) {
  p <- ggplot(PCA$components, aes(x = PC1, y = PC2, color = Condition)) +
        geom_point(size=4) +
        xlab(paste0("PC1","(",PCA$summary$percent[1],"%)")) +
        ylab(paste0("PC2","(",PCA$summary$percent[2],"%)")) +
        theme_classic()
    ggsave(p,filename=paste0(output_dir, output_filename), height = 7,width = 9)
}


plot_hierarchical_cluster <- function(DT, output_dir) {
    cluster_data <- DT[,-c(1:2)]
    cluster_data[is.na(cluster_data)]=0
    cluster_data=cluster_data[which(rowSums(cluster_data)>0),]
    log2_cluster_data=log2(cluster_data+1)
  
  dist_mat <- dist(t(log2_cluster_data)) #
  hc_cluster <- hclust(dist_mat,method = "complete")
  g <- ggdendrogram(hc_cluster, rotate=TRUE) + labs(title='Hierarchical clustering')
  cat(paste0('   -> ', output_dir, '\n'))
  ggsave(g, filename=paste0(output_dir, 'hc_cluster_log2.pdf'))
}


get_umap <- function(DT, neighbors) {
    cluster_data=DT[,-c(1:2)]
    cluster_data[is.na(cluster_data)]=0
    cluster_data=cluster_data[which(rowSums(cluster_data)>0),]
    log2_cluster_data=log2(cluster_data+1)
  
    set.seed(100)
    DT.umap <- umap(t(log2_cluster_data), n_neighbors=neighbors)
    DT.out <- as.data.table(DT.umap$layout, keep.rownames=TRUE)
    setnames(DT.out, c('Sample', 'UMAP1', 'UMAP2'))
    DT.out$condition=gsub('_[1-9]*$','',DT.out$Sample)
    return(DT.out[])
}


plot_umap <- function(DT, output_dir, output_filename) {
    g <- ggplot(DT, aes(x=UMAP1, y=UMAP2, color=condition)) +
    geom_point(size=4) +
    theme_classic() 
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(output_dir, output_filename), height = 7,width = 9)
}

# filter_pgs <- function(DT.all, treatment_sample_names, treatment, control)  {
#     control_sample_names <- colnames(DT.all)[colnames(DT.all) %like% control]
#     n_treatment <- length(treatment_sample_names)
#     n_controls <- length(control_sample_names)

#     info_cols <- colnames(DT.all)[1:2]
#     DT <- DT.all[, c(info_cols, treatment_sample_names, control_sample_names), with=F]

#     # Filter: protein groups with at least half of control and treatment samples having non-zero value
#     DT[, 'N_missing_treatment' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(treatment_sample_names)]
#     DT[, 'N_missing_control' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_sample_names)]
#     DT <- DT[N_missing_treatment < (n_treatment/2)]
#     DT <- DT[N_missing_control < (n_controls/2)]
#     DT[, 'N_missing_treatment' := NULL]
#     DT[, 'N_missing_control' := NULL]
#     return(DT[])
# }


do_t_test <- function(DT, treatment_samples, control_samples) {
  
  dat <- copy(DT)
  n_treatment <- length(treatment_samples)
  n_control <- length(control_samples)
  
  # Retain first three columns plus all treatment and control columns
  dat <- dat[,c(colnames(DT)[1:2], treatment_samples, control_samples), with=F]
  
  # Convert NA to 0
  dat[is.na(dat)] <- 0
  
  # Drop rows (protein groups) with > 50% missingness in samples
  dat[,'missing_value':= apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_samples,treatment_samples) ]
  dat <- dat[missing_value <= (n_treatment+n_control)/2]
  dat[, 'missing_value' := NULL]
  
  
  # Drop rows (protein groups) with 0 variance in treatment OR control group
  dat[, 'control_var' := apply(.SD, 1, var), .SDcols=c(control_samples)]
  dat[, 'treatment_var' := apply(.SD, 1, var), .SDcols=c(treatment_samples)]
  dat <- dat[control_var != 0]
  dat <- dat[treatment_var != 0]
  dat[, c('control_var','treatment_var') := NULL]
  
  # Perform t-test  on treatment and control columns
  t_test <- apply(dat[,-c(1:2)], 1, function(x){
    a =factor(c(rep('treatment',n_treatment),
                rep("control",n_control)),
              levels = c('treatment',"control"))
    fvalue=var.test(x~a)
    if (!is.na(fvalue$p.value)){ 
      if (fvalue$p.value > 0.05){
        result <- t.test(x~a, var.equal = T)
      } else {
        result <- t.test(x~a, var.equal = F)
      }
    }
    treatment_estimate <- as.numeric(unlist(result$estimate[1]))
    control_estimate <- as.numeric(unlist(result$estimate[2]))
    return(data.table('P_value'=result$p.value,
                      'treatment_estimate'=treatment_estimate,
                      'control_estimate'=control_estimate)
    )
  })
  
  
  t_test <- rbindlist(t_test)
  t_test <- cbind(dat[,c(1:2)], t_test)   # add back protein group / gene info cols
  t_test[, log2FC := log2(treatment_estimate / (control_estimate+1))]
  t_test[, p.adj := p.adjust(P_value, method='BH')]
  return(t_test[])
}

# run_contrast <- function(DT.all, treatment_sample_names, treatment, control)  {
#     control_sample_names <- colnames(DT.all)[colnames(DT.all) %like% control]
#     n_treatment <- length(treatment_sample_names)
#     n_controls <- length(control_sample_names)

#     info_cols <- colnames(DT.all)[1:3]
#     DT <- DT.all[, c(info_cols, treatment_sample_names, control_sample_names), with=F]

#     # Filter: protein groups with at least half of control and treatment samples having non-zero value
#     DT[, 'N_missing_treatment' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(treatment_sample_names)]
#     DT[, 'N_missing_control' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_sample_names)]
#     DT <- DT[N_missing_treatment < (n_treatment/2)]
#     DT <- DT[N_missing_control < (n_controls/2)]

#     DT[, 'control_var' := apply(.SD, 1, var), .SDcols=c(control_sample_names)]
#     DT[, 'treatment_var' := apply(.SD, 1, var), .SDcols=c(control_sample_names)]
#     DT[, 'overall_var' := apply(.SD, 1, var), .SDcols=c(control_sample_names, treatment_sample_names)]
#     DT <- DT[control_var != 0][treatment_var != 0][overall_var != 0]
#     DT.out <- copy(DT[, info_cols, with=F])
#     DT <- DT[, c(treatment_sample_names, control_sample_names), with=F]
#     # Run contrast
#     pvalue <- apply(DT, 1, function(x) {
#                                 a = factor(c(rep(treatment, n_treatment), rep(control, n_controls)),
#                                     levels=c(treatment, control))
#                                 fvalue = var.test(x~a)
#                                 if (!is.na(fvalue$p.value)){
#                                     if (fvalue$p.value > 0.05) {
#                                         t.test(x~a, var.equal = TRUE)
#                                     } else {
#                                         t.test(x~a, var.equal = FALSE)
#                                     }
#                                 }
#         }
#     )

#     DT.out[, p := as.numeric(unlist(lapply(pvalue,function(x) x$p.value)))]
#     DT.out[, q := p.adjust(p, method='BH')]       # BH method to convert P to FDR (q) value
#     DT.out[, ratio := as.numeric(unlist(lapply(pvalue,function(x) x$estimate[1]/(x$estimate[2])))) ]
#     return(DT.out[])
# }

# run_contrast_on_raw_count <- function(DT.all, treatment_sample_names, treatment, control)  {
#     control_sample_names <- colnames(DT.all)[colnames(DT.all) %like% control]
#     n_treatment <- length(treatment_sample_names)
#     n_controls <- length(control_sample_names)

#     info_cols <- colnames(DT.all)[1:3]
#     DT <- DT.all[, c(info_cols, treatment_sample_names, control_sample_names), with=F]
#     DT[, lapply(.SD, function(x) log_base^x), .SDcols=c(control_sample_names, treatment_sample_names)]

#     # Filter: protein groups with at least half of control and treatment samples having non-zero value
#     DT[, 'N_missing_treatment' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(treatment_sample_names)]
#     DT[, 'N_missing_control' := apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_sample_names)]
#     DT <- DT[N_missing_treatment < (n_treatment/2)]
#     DT <- DT[N_missing_control < (n_controls/2)]
#     DT.out <- copy(DT[, info_cols, with=F])
#     DT <- DT[, c(treatment_sample_names, control_sample_names), with=F]

#     # Run contrast
#     pvalue <- apply(DT, 1, function(x) {
#                                 a = factor(c(rep(treatment, n_treatment), rep(control, n_controls)),
#                                     levels=c(treatment, control))
#                                 fvalue = var.test(x~a)
#                                 if (!is.na(fvalue$p.value)){
#                                     if (fvalue$p.value > 0.05) {
#                                         t.test(x~a, var.equal = TRUE)
#                                     } else {
#                                         t.test(x~a, var.equal = FALSE)
#                                     }
#                                 }
#         }
#     )

#     DT.out[, p := as.numeric(unlist(lapply(pvalue,function(x) x$p.value)))]
#     DT.out[, q := p.adjust(p, method='BH', n=.N)]       # BH method to convert P to FDR (q) value
#     DT.out[, ratio := as.numeric(unlist(lapply(pvalue,function(x) x$estimate[1]/(x$estimate[2])))) ]
#     return(DT.out[])
# }

do_t_test_APMS <- function(DT, treatment_samples, control_samples) {
  
    dat <- copy(DT)
    n_treatment <- length(treatment_samples)
    n_control <- length(control_samples)

    # Retain first three columns plus all treatment and control columns
    dat <- dat[,c(colnames(DT)[1:2], treatment_samples, control_samples), with=F]

    # Convert NA to 0
    dat[is.na(dat)] <- 0

    # Drop rows (protein groups) with > 50% missingness in samples
    dat[,'missing_value':= apply(.SD, 1, function(x) sum(x==0)), .SDcols=c(control_samples,treatment_samples) ]
    dat <- dat[missing_value <= (n_treatment+n_control)/2]
    dat[, 'missing_value' := NULL]

    # Perform t-test  on treatment and control columns
    t_test <- apply(dat[,-c(1:2)], 1, function(x){
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
    return(data.table('P_value'=result$p.value,
                        'treatment_estimate'=treatment_estimate,
                        'control_estimate'=control_estimate)
    )
    })


    t_test <- rbindlist(t_test)
    t_test <- cbind(dat[,c(1:2)], t_test)   # add back protein group / gene info cols
    t_test[, log2FC := log2(treatment_estimate / (control_estimate+1))]
    t_test[, p.adj := p.adjust(P_value, method='BH')]
    return(t_test[])
}



plot_volcano <- function(DT.original, treatment, control, log_base, lfc_threshold, fdr_threshold, out_dir, labelgene) {
    options(ggrepel.max.overlaps=Inf)
    DT <- copy(DT.original)
    DT[, 'Group' := 'Others']
    DT[log2FC >= lfc_threshold, 'Group' := 'UP']
    DT[log2FC <= -lfc_threshold, 'Group' := 'DOWN']
    DT[p.adj >= fdr_threshold, 'Group' := 'Others']
    DT[, labeltext := '']

    up_gene_5 <- DT[Group == 'UP', ]
    up_gene_5 <- up_gene_5[order(up_gene_5$log2FC,decreasing = T)[1:5],]
    down_gene_5 <- DT[Group == 'DOWN', ]
    down_gene_5 <- down_gene_5[order(down_gene_5$log2FC,decreasing = F)[1:5],]
    top5_gene <- rbind(up_gene_5,down_gene_5)
    DT[Genes %in% top5_gene$Genes, labeltext := Genes]

    if (!is.null(labelgene)) {
        DT[Genes %in% labelgene, labeltext := Genes]
    }
    g <- ggplot(DT, aes(x=log2FC, y=-log10(p.adj))) +
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
}

plot_volcano_from_raw <- function(DT.original, treatment, control, log_base, lfc_threshold, fdr_threshold, out_dir) {
    DT <- copy(DT.original)
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
}



###pathway analysis
enrich_pathway = function(DT.original, treatment, control, outdir, lfc_threshold, fdr_threshold, enrich_pvalue){
    ##create dir
    dir=paste0(outdir,'/', treatment, '_vs_', control)
    if (!dir.exists(dir)){
        dir.create(dir,recursive = T)
    }

    DT <- copy(DT.original)
    DT[, 'Group' := 'Others']
    DT[log2FC >= lfc_threshold, 'Group' := 'UP']
    DT[log2FC <= -lfc_threshold, 'Group' := 'DOWN']
    DT[p.adj >= fdr_threshold, 'Group' := 'Others']

    ##ENTREZID names
    entrizid = data.frame(bitr(DT$Genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop=T))

    DT <- merge(DT,entrizid,by.x='Genes',by.y='SYMBOL')
    DT <- na.omit(DT)
    DT <- DT[order(DT$log2FC,decreasing = T),]
    all_gene_vector=DT$log2FC
    names(all_gene_vector)=DT$ENTREZID
    ## up and down regulated genes 
    up_genes=DT[which(DT$log2FC>=lfc_threshold&DT$p.adj<=fdr_threshold),]
    down_genes=DT[which(DT$log2FC<=(-lfc_threshold)&DT$p.adj<=fdr_threshold),]
    ## ------- Enrichment -------
    print('Processing up-gene enrichment analysis')
    if (nrow(up_genes)>0){
        enrichAll(gene_id=up_genes$ENTREZID, all_gene_vector=all_gene_vector, MSigDb_gmt_dir = MSigDb_gmt_dir, out_prefix = 'up',outdir=dir,width = 12,height = 8,enrich_pvalue=enrich_pvalue)
    }
    print('Processing down-gene enrichment analysis')
    if (nrow(down_genes)>0){
        enrichAll(gene_id=down_genes$ENTREZID, all_gene_vector=all_gene_vector, MSigDb_gmt_dir = MSigDb_gmt_dir, out_prefix = 'down',outdir=dir,width = 12,height = 8,enrich_pvalue=enrich_pvalue)
    }
}


enrichAll = function(gene_id, all_gene_vector, outdir='.', MSigDb_gmt_dir = NULL, out_prefix = NULL, width=12, height=8, enrich_pvalue=1){

    all_gene_vector <- na.omit(all_gene_vector)

    if (!dir.exists(outdir)){
        dir.create(outdir,recursive = T)
    }

    enrich_res = data.frame()
    ### GO 
    print('Processing GO')
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
    if (is.null(ego_cc)) {
    print('NO GO CC')
    } else {head(ego_cc)
    tmp_df=ego_cc@result
    tmp_df$group='go_cc'
    enrich_res=rbind(enrich_res,tmp_df)

    if (nrow(ego_cc)>0){
        barplot(ego_cc, showCategory=20,label_format = 100)
        ggsave(file.path(outdir,paste0(out_prefix,'.go_cc.pdf')),width = width,height=height)
    }

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
    if (is.null(ego_bp)) {
    print('NO GO BP')
    }else {
    head(ego_bp)
    tmp_df=ego_bp@result
    tmp_df$group='go_bp'
    enrich_res=rbind(enrich_res,tmp_df)
    if (nrow(ego_bp)>0){
        barplot(ego_bp, showCategory=20,label_format = 100)
        ggsave(file.path(outdir,paste0(out_prefix,'.go_bp.pdf')),width=width,height=height)
    }
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
    if (is.null(ego_mf)) {
    print('NO GO MF')
    }else {
    head(ego_mf)
    tmp_df=ego_mf@result
    tmp_df$group='go_mf'
    enrich_res=rbind(enrich_res,tmp_df)

    if (nrow(ego_mf)>0){
        barplot(ego_mf, showCategory=20,label_format = 100)
        ggsave(file.path(outdir,paste0(out_prefix,'.go_mf.pdf')),width=width,height=height)
    }
    }

    ### KEGG
    #print('Processing KEGG')
    #kk <- enrichKEGG(gene         = gene_id,
    #                 organism     = 'hsa',
    #                universe      = names(all_gene_vector),
    #                pvalueCutoff = enrich_pvalue)
    #head(kk)
    #tmp_df=kk@result
    #tmp_df$group='kegg'
    #enrich_res=rbind(enrich_res,tmp_df)

    #if (nrow(kk)>0){
    #  barplot(kk, showCategory=20)
    #  ggsave(file.path(outdir,paste0(out_prefix,'.kegg.pdf')),width=width,height=height)
    # }
    print('Save enrichment analysis results')
    write.table(enrich_res,file.path(outdir,paste0(out_prefix,'.enrich_res.txt')),quote = F,sep = '\t')
}

get_ppi <- function(protein,DT,lfc_threshold, fdr_threshold) {
    ##canidate genes
    print('canidate genes')
    DT[, 'Group' := 'Others']
    DT[log2FC >= lfc_threshold, 'Group' := 'UP']
    DT[log2FC <= -lfc_threshold, 'Group' := 'DOWN']
    DT[p.adj >= fdr_threshold, 'Group' := 'Others']
    up_genes=DT[which(DT$log2FC>=lfc_threshold&DT$p.adj<=fdr_threshold),]
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
    up_genes=t_test[which(t_test$log2FC>=lfc_threshold&t_test$p.adj<=fdr_threshold),]
    up_genes=up_genes[order(up_genes$log2FC,decreasing = T),]
    up_genes$rank=1:nrow(up_genes)
    up_genes_plot=merge(up_genes,PPI_score,by="Genes",all=T)
    g <- ggplot(up_genes_plot) +
    geom_point(aes(x=log2FC.x, y=rank),color="#999999") +
    geom_point(aes(x=log2FC.y, y=rank),color='#E69F00')+
    theme_classic()+
    labs(x="Log2FC",y='Rank')+
    geom_label_repel(
        data = up_genes_plot,
        aes(x=log2FC.y, y=rank,label = Genes),
        size = 5,
        max.overlaps=50,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"))

    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(output_dir, output_filename), height=10, width=12)
}

plot_PPI_ven <- function(protein,t_test,lfc_threshold, fdr_threshold, output_dir, output_filename) {
  print('canidate genes')
  t_test[, 'Group' := 'Others']
  t_test[log2FC >= lfc_threshold, 'Group' := 'UP']
  t_test[log2FC <= -lfc_threshold, 'Group' := 'DOWN']
  t_test[p.adj >= fdr_threshold, 'Group' := 'Others']
  up_genes=t_test[which(t_test$log2FC>=lfc_threshold&t_test$p.adj<=fdr_threshold),]
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


get_mhcflurry_input=function(hla_typing,dat.long,sample,MHC_dir){
    print(paste0("Generate the MHCflurry input file of ",sample))
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
  command <- paste0('source activate mhcflurry-env', '\n',"mhcflurry-predict ",MHC_dir,sample, '_allele_peptide.csv'," --out " ,MHC_dir,sample,"_mhc_predictions.csv")
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
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))+
    xlab('-log10(MHC Affinity)') +
    ylab('-log2(MHC Affinity Percentile))')
  ggsave(g, filename=paste0(MHC_dir,sample,"_mhc_affinity_allele.pdf"),width = 8,height = 8)
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
    ylab('-log2(MHC Affinity)')
  ggsave(p, filename=paste0(MHC_dir,sample,"_mhc_affinity_rank.pdf"),width = 8,height = 8)
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
    xlab('-log10(MHC Affinity)') +
    ylab('-log2(MHC Affinity Percentile))')
  
  ggsave(f, filename=paste0(MHC_dir,sample,"_mhc_affinity_allele_separated.pdf"),width = 8,height = 8)
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_affinity_allele_separated.pdf", '\n'))
  
  f <- ggplot(sample_prediction, aes(x=allele, y=log2(mhcflurry_affinity))) +
    geom_dotplot(binaxis='y',, stackdir='center',binwidth = 1/50,color='#E69F00',fill="#E69F00")+
    theme_classic()
  ggsave(f, filename=paste0(MHC_dir,sample,"_mhc_affinity_allele_dotplot.pdf"),width = 8,height = 8)
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_affinity_allele_dotplot.pdf", '\n'))
  
  
  number=sample_prediction[, .N, by=allele]
  p=ggplot(number, aes(x=allele, y=N)) +
    geom_bar(stat="identity", fill="steelblue")+
    theme_classic()+
    labs(fill = "",x="",y='Number of Peptide')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(p, filename=paste0(MHC_dir,sample,"_mhc_affinity_peptide_count.pdf"),width = 8,height = 8)
  cat(paste0('   -> ', MHC_dir,sample,"_mhc_affinity_peptide_count.pdf", '\n'))
  
  
}

