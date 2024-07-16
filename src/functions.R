
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

convert_log_to_raw <- function(DT.original, log_base) {
  DT <- copy(DT.original)
  samplenames <- colnames(DT[,-c(1:2)])
  DT[, (samplenames) := lapply(.SD, function(x) log_base^x), .SDcols=samplenames]
  return(DT[])
}

#total proteomics############
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

melt_intensity_table <- function(DT) {
    # Converts intensity data.table to long format
    # info_cols <- c('Protein_Group', 'Genes', 'First_Protein_Description')
  DT.long <- melt(DT, 
                  measure.vars=colnames(DT)[!(colnames(DT) %in% c("Protein_Group", "Genes", "PTM",
                                                                  "PTM_Location", "Precursor", "Modified_Sequence"))],
                  variable.name='Sample',
                  value.name='Intensity')
  DT.long=data.table(DT.long)
  return(DT.long)
}

plot_pg_counts <- function(DT.long, output_dir, output_filename) {
  pgcounts <- DT.long[, .N, by=Sample]
  # Order samples by ascending counts
  ezwrite(pgcounts, QC_dir, 'protein_group_nonzero_counts.tsv')
  plot_pg_counts(pgcounts, QC_dir, 'protein_group_nonzero_counts.pdf')
  n_samples <- nrow(DT.long)
  if (n_samples > 20) {
    p=ggplot(DT, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="#67a9cf")+
      theme_classic()+
      labs(fill = "",x="",y='Number of Protein Groups')+
      scale_x_discrete(guide = guide_axis(angle = 90))
  }
  if (n_samples <= 20) {
    p=ggplot(DT, aes(x=Sample, y=N)) +
      geom_bar(stat="identity", fill="#67a9cf")+
      theme_classic()+
      labs(fill = "",x="",y='Number of Protein Groups')+
      scale_x_discrete(guide = guide_axis(angle = 90))+ 
      geom_text(aes(label=N, y=N + (0.05*max(DT$N))))
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

####correlation
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

###PCA
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
}

###HC cluster
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

####umap
get_umap <- function(DT, neighbors) {
    cluster_data=DT[,-c(1:2)]
    cluster_data[is.na(cluster_data)]=0
    cluster_data=cluster_data[which(rowSums(cluster_data)>0),]
    log2_cluster_data=log2(cluster_data+1)
  
    set.seed(100)
    DT.umap <- umap(t(log2_cluster_data), n_neighbors=neighbors)
    DT.out <- as.data.table(DT.umap$layout, keep.rownames=TRUE)
    setnames(DT.out, c('Sample', 'UMAP1', 'UMAP2'))
    DT.out$condition=gsub('_[0-9]+$','',DT.out$Sample)
    return(DT.out[])
}

plot_umap <- function(DT, output_dir, output_filename) {
    g <- ggplot(DT, aes(x=UMAP1, y=UMAP2, color=condition)) +
    geom_point(size=4) +
    theme_classic() 
    cat(paste0('   -> ', output_dir, output_filename, '\n'))
    ggsave(g, filename=paste0(output_dir, output_filename), height = 4,width = 5)
}

###ttest
do_t_test <- function(DT, treatment_samples, control_samples) {
  
  DT_ttest <- copy(DT)
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
  
  
  # Drop rows (protein groups) with 0 variance in treatment OR control group
  DT_ttest[, 'control_var' := apply(.SD, 1, var), .SDcols=c(control_samples)]
  DT_ttest[, 'treatment_var' := apply(.SD, 1, var), .SDcols=c(treatment_samples)]
  DT_ttest <- DT_ttest[control_var != 0]
  DT_ttest <- DT_ttest[treatment_var != 0]
  DT_ttest[, c('control_var','treatment_var') := NULL]
  
  # Perform t-test  on treatment and control columns
  t_test <- apply(DT_ttest[,-c(1:2)], 1, function(x){
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
  t_test <- cbind(DT_ttest[,c(1:2)], t_test)   # add back protein group / gene info cols
  t_test[, log2FC := log2(treatment_estimate / (control_estimate+1))]
  t_test[, p.adj := p.adjust(P_value, method='BH')]
  return(t_test[])
}
###ttest plot_volcano
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

plot_volcano_pep <- function(DT.original, treatment, control, log_base, lfc_threshold, fdr_threshold, out_dir, labelgene) {
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
  DT[Peptide_Sequence %in% top5_gene$Peptide_Sequence, labeltext := Peptide_Sequence]
  
  if (!is.null(labelgene)) {
    DT[Peptide_Sequence %in% labelgene, labeltext := Peptide_Sequence]
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

###ttest pathway analysis
enrich_pathway = function(DT.original, treatment, control, outdir, lfc_threshold, fdr_threshold, enrich_pvalue){
  ##create dir
  dir=paste0(outdir,'/', treatment, '_vs_', control)
  if (!dir.exists(dir)){
    dir.create(dir,recursive = T)
  }
  
  DT <- data.table(DT.original)
  DT[, 'Group' := 'Others']
  DT[log2FC >= lfc_threshold, 'Group' := 'UP']
  DT[log2FC <= -lfc_threshold, 'Group' := 'DOWN']
  DT[p.adj >= fdr_threshold, 'Group' := 'Others']
  
  ##ENTREZID names########
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
      barplot(ego_cc, showCategory=20,label_format = 100)+
        scale_fill_gradient(low = "#ef8a62", high = "#67a9cf", na.value = NA)
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
      barplot(ego_bp, showCategory=20,label_format = 100)+
        scale_fill_gradient(low = "#ef8a62", high = "#67a9cf", na.value = NA)
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
      barplot(ego_mf, showCategory=20,label_format = 100)+
        scale_fill_gradient(low = "#ef8a62", high = "#67a9cf", na.value = NA)
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
    control=unique(design_matrix$control[design_matrix$condition==treatment])
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
    limma_plot_volcano(DT_limma_res = result_limma ,
                       lfc_threshold = lfc_threshold,
                       fdr_threshold = fdr_threshold,
                       out_dir = DE_dir,
                       output_filename =paste0(treatment, ' vs ', control),
                       labelgene=opt$labelgene)
    enrich_pathway_limma(DT.original = result_limma,
                         treatment = treatment,
                         control = control,
                         outdir = EA_dir,
                         lfc_threshold = lfc_threshold,
                         fdr_threshold = fdr_threshold,
                         enrich_pvalue = enrich_pvalue)
  }
}
###limma plot_volcano
limma_plot_volcano <- function(DT_limma_res, lfc_threshold, fdr_threshold, out_dir, output_filename,labelgene) {
  options(ggrepel.max.overlaps=Inf)
  DT_limma_res=data.table(DT_limma_res)
  DT_limma_res[, 'Group' := 'Others']
  DT_limma_res[logFC >= lfc_threshold, 'Group' := 'UP']
  DT_limma_res[logFC <= -lfc_threshold, 'Group' := 'DOWN']
  DT_limma_res[adj.P.Val >= fdr_threshold, 'Group' := 'Others']
  DT_limma_res[, labeltext := '']
  
  up_gene_5 <- DT_limma_res[Group == 'UP', ]
  up_gene_5 <- up_gene_5[order(up_gene_5$logFC,decreasing = T)[1:5],]
  down_gene_5 <- DT_limma_res[Group == 'DOWN', ]
  down_gene_5 <- down_gene_5[order(down_gene_5$logFC,decreasing = F)[1:5],]
  top5_gene <- rbind(up_gene_5,down_gene_5)
  top5_gene$labeltext=top5_gene$Genes
  DT_limma_res[Genes %in% top5_gene$Genes, labeltext := Genes]
  if (!is.null(labelgene)) {
    DT_limma_res[Genes %in% labelgene, labeltext := Genes]
  }
  g <- ggplot(DT_limma_res, aes(x=logFC, y=-log10(adj.P.Val))) +
    geom_point(aes(color = Group)) +
    scale_color_manual(breaks = c("DOWN", "Others", "UP"), 
                       values=c("#67a9cf", "#969696","#ef8a62"))+
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_label_repel(max.overlaps = Inf,
                     data = subset(top5_gene),
                     aes(label = labeltext),
                     size = 5,
                     box.padding = unit(0.35, "lines"),
                     point.padding = unit(0.3, "lines"))+
    geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed")+ 
    geom_vline(xintercept=lfc_threshold, linetype="dashed")+ 
    geom_vline(xintercept=-lfc_threshold, linetype="dashed")+
    theme_classic()
  cat(paste0('   -> ', out_dir, output_filename, '\n'))
  ggsave(g, filename=paste0(out_dir, '/',output_filename,'_vocanol.pdf'),width = 8,height = 8)
}
##limma pathway analysis
enrich_pathway_limma = function(DT.original, treatment, control, outdir, lfc_threshold, fdr_threshold, enrich_pvalue){
  ##create dir
  dir=paste0(outdir,'/', treatment, '_vs_', control)
  if (!dir.exists(dir)){
    dir.create(dir,recursive = T)
  }
  
  DT <- data.table(DT.original)
  DT[, 'Group' := 'Others']
  DT[logFC >= lfc_threshold, 'Group' := 'UP']
  DT[logFC <= -lfc_threshold, 'Group' := 'DOWN']
  DT[adj.P.Val >= fdr_threshold, 'Group' := 'Others']
  
  ##ENTREZID names
  entrizid = data.frame(bitr(DT$Genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop=T))
  
  DT <- merge(DT,entrizid,by.x='Genes',by.y='SYMBOL')
  DT <- na.omit(DT)
  DT <- DT[order(DT$logFC,decreasing = T),]
  all_gene_vector=DT$logFC
  names(all_gene_vector)=DT$ENTREZID
  ## up and down regulated genes 
  up_genes=DT[which(DT$logFC>=lfc_threshold&DT$adj.P.Val<=fdr_threshold),]
  down_genes=DT[which(DT$logFC<=(-lfc_threshold)&DT$adj.P.Val<=fdr_threshold),]
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

##heatmap
plot_heatmap <- function(DT_heatmap) {
  sample_anno <- data.frame(condition = as.factor(sort(colnames(DT_heatmap)[3:ncol(DT_heatmap)])))
  row.names(sample_anno) <- sample_anno$condition
  sample_anno$condition=gsub("_[0-9]*$","",sample_anno$condition)
  DT_heatmap[is.na(DT_heatmap)]=0
  pheatmap(log2(DT_heatmap[,3:ncol(DT_heatmap)]+1),
           scale='row',
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
           scale='row',
           show_rownames = T,
           annotation_col = sample_anno,
           cluster_cols = F,
           treeheight_row=0,
           border_color =NA,
           filename=paste0(opt$outdir,'/','heatmap_label_gene.pdf'))
  
}

#####PPI##########
##APMS ttest
do_t_test_APMS <- function(DT, treatment_samples, control_samples) {
  
  DT_ttest <- copy(DT)
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
    return(data.table('P_value'=result$p.value,
                      'treatment_estimate'=treatment_estimate,
                      'control_estimate'=control_estimate)
    )
  })
  
  
  t_test <- rbindlist(t_test)
  t_test <- cbind(DT_ttest[,c(1:2)], t_test)   # add back protein group / gene info cols
  t_test[, log2FC := log2(treatment_estimate / (control_estimate+1))]
  t_test[, p.adj := p.adjust(P_value, method='BH')]
  return(t_test[])
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
        point.padding = unit(0.3, "lines"))+
      theme(axis.title.x = element_text(size=20, face="bold"),
            axis.title.y = element_text(size=20, face="bold"),
            axis.text.x = element_text(size=14, angle=0),
            axis.text.y = element_text(size=14, angle=0),
            legend.title=element_text(size=20),
            legend.text = element_text(size=14))

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

####mhcflurry########
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

####pSILAC functions#################################
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
  if ('PG.Genes' %in% colnames(DT)) {
    setnames(DT, 'PEP.IsGeneSpecific', 'GeneSpecific')
  }
  if ('EG.PrecursorId' %in% colnames(DT)) {
    setnames(DT, 'EG.PrecursorId', 'Precursor')
  }
  DT=as.data.frame(DT)
  DT=DT[,grep('Protein_Group|Genes|Channel',colnames(DT))]
  # Remove trailing file extensions
  extensions = '.mzML$|.mzml$|.RAW$|.raw$|.dia$|.DIA$|_Intensity|raw.PG.MS2|DIA_|raw.PEP.MS2|raw.EG.TargetReference| \\([^()]*\\)'
  extension_samplenames =  colnames(DT)[colnames(DT) %like% extensions]
  trimmed_samplenames = gsub(extensions, '', extension_samplenames)
  setnames(DT, extension_samplenames, trimmed_samplenames)
  # trim leading [N] 
  colnames_out = gsub(pattern="\\[.*\\] ", replacement='', x=colnames(DT))
  setnames(DT, colnames_out)
  #as number
  DT[,grep('Channel',colnames(DT))]=as.data.frame(apply(DT[,grep('Channe',colnames(DT))],2,as.numeric))
  DT=data.table(DT)
  return(DT[])
}

# pg count
plot_silac_pg_counts = function(DT.long, output_dir,height,width) {
  counts=DT.long[, .N, by=Sample]
  counts$Chanel=gsub('.*Channel','Channel',counts$Sample)
  ezwrite(counts, output_dir, paste0('silac_chanel_pgcount','.tsv'))
  p=ggplot(counts, aes(x=Sample, y=N,fill=Chanel)) +
    geom_bar(stat="identity")+
    theme_classic()+
    labs(fill = "",x="",y='Protein Group number')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir, 'silac_chanel_pgcount','.pdf'),height=height,width = width)
}

#pep count
plot_silac_pep_counts = function(DT.long, output_dir,height,width) {
  counts=DT.long[, .N, by=Sample]
  counts$Chanel=gsub('.*Channel','Channel',counts$Sample)
  ezwrite(counts, output_dir, paste0('silac_chanel_pepcount','.tsv'))
  p=ggplot(counts, aes(x=Sample, y=N,fill=Chanel)) +
    geom_bar(stat="identity")+
    theme_classic()+
    labs(fill = "",x="",y='Peptide number')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir, 'silac_chanel_pepcount','.pdf'),height=height,width = width)
}

##plot protein group or peptide intensity of all the sample
plot_silac_pg_intensity = function(DT.long, output_dir,height,width) {
  DT.long$Chanel=gsub('.*Channel','Channel',DT.long$Sample)
  DT.long$Condition=gsub('.Channel.*','',DT.long$Sample)
  DT.long$Sample=gsub('_[0-9]$','',DT.long$Condition)
  p=ggplot(DT.long, aes(y=Chanel, x=log10(Intensity),fill=Chanel)) + 
    geom_density_ridges()+ 
    facet_wrap(. ~ Sample,ncol = 3)+
    theme_classic()+
    theme(axis.text.y = element_blank())
  ggsave(plot = p,filename = paste0(output_dir, 'silac_pg_intensity','_sample.pdf'),width = width,height =height)
  p=ggplot(DT.long, aes(y=Condition, x=log10(Intensity),fill=Condition)) + 
    geom_density_ridges()+ 
    facet_wrap(. ~ Chanel)+
    theme_classic()
  ggsave(plot = p,filename = paste0(output_dir, 'silac_pg_intensity','_chanel.pdf'),width = width,height =height)
}

plot_silac_pep_intensity = function(DT.long, output_dir,height,width) {
  DT.long$Chanel=gsub('.*Channel','Channel',DT.long$Sample)
  DT.long$Condition=gsub('.Channel.*','',DT.long$Sample)
  DT.long$Sample=gsub('_[0-9]$','',DT.long$Condition)
  p=ggplot(DT.long, aes(y=Chanel, x=log10(Intensity),fill=Chanel)) + 
    geom_density_ridges()+ 
    facet_wrap(. ~ Sample,ncol = 3)+
    theme_classic()+
    theme(axis.text.y = element_blank())
  ggsave(plot = p,filename = paste0(output_dir, 'silac_pep_intensity','_sample.pdf'),width = width,height =height)
  p=ggplot(DT.long, aes(y=Condition, x=log10(Intensity),fill=Condition)) + 
    geom_density_ridges()+ 
    facet_wrap(. ~ Chanel)+
    theme_classic()
  ggsave(plot = p,filename = paste0(output_dir, 'silac_pep_intensity','_chanel.pdf'),width = width,height =height)
}



#####phospho functions#################################
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
  DT=as.data.frame(DT)
  DT=DT[,grep('Protein_Group|Genes|PTM_Location|Precursor|Modified_Sequence|raw',colnames(DT))]
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

##plot phospho site counts
plot_phospho_site_counts = function(DT_phospho_long, output_dir,height,width) {
  counts=DT_phospho_long[, .N, by=Sample]
  counts$Condition=as.factor(gsub('_[0-9]*','',counts$Sample))
  ezwrite(counts, output_dir, paste0('phospho_site_count','.tsv'))
  p=ggplot(counts, aes(x=Sample, y=N,fill=Condition)) +
    geom_bar(stat="identity")+
    theme_classic()+
    labs(fill = "",x="",y='Number of Phospho-sites')+
    scale_x_discrete(guide = guide_axis(angle = 90))
  ggsave(plot = p,filename = paste0(output_dir, 'phospho_site_count','.pdf'),height=height,width = width)
  #sd
  summary_data <- counts %>%
    group_by(Condition) %>%
    summarize(mean = mean(N), sd = sd(N)) %>%
    arrange(Condition)
  p=ggplot(summary_data, aes(x=as.factor(Condition), y=mean)) +
    geom_bar(stat="identity",fill="#67a9cf", position=position_dodge())+
    theme_classic()+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9))+
    labs(fill = "",x="",y='Number of Phospho-sites')
  p
  ggsave(plot = p,filename = paste0(output_dir,'/','phospho_sites_number_sd','.pdf'),height=5,width = 4)
}

##plot phospho site distribution
plot_phospho_site_intensity = function(DT_phospho_long, output_dir) {
  n_samples <- length(unique(DT_phospho_long$Sample))
  DT_phospho_long$Condition=gsub('_[0-9]*$','',DT_phospho_long$Sample)
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
    ggsave(plot = p,filename = paste0(output_dir, 'Phospho_site_intensity.pdf'),width = 6,height =n_samples/5)
  }else{
    ggsave(plot = g,filename = paste0(output_dir, 'Phospho_site_intensity_boxplot.pdf'),width = 6,height = 5)
    ggsave(plot = p,filename = paste0(output_dir, 'Phospho_site_intensity.pdf'),width = 4,height =6)
  }
}

##phospho site DE analysis

#ttest
phospho_ttest = function(DT_phospho, treatment_samples, control_samples) {
    DT_ttest <- copy(DT)
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
    
    
    # Drop rows (protein groups) with 0 variance in treatment OR control group
    DT_ttest[, 'control_var' := apply(.SD, 1, var), .SDcols=c(control_samples)]
    DT_ttest[, 'treatment_var' := apply(.SD, 1, var), .SDcols=c(treatment_samples)]
    DT_ttest <- DT_ttest[control_var != 0]
    DT_ttest <- DT_ttest[treatment_var != 0]
    DT_ttest[, c('control_var','treatment_var') := NULL]
    
    # Perform t-test  on treatment and control columns
    t_test <- apply(DT_ttest[,-c(1:2)], 1, function(x){
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
    t_test <- cbind(DT_ttest[,c(1:2)], t_test)   # add back protein group / gene info cols
    t_test[, log2FC := log2(treatment_estimate / (control_estimate+1))]
    t_test[, p.adj := p.adjust(P_value, method='BH')]
    return(t_test[])
  }
###ttest plot_volcano
plot_phospho_ttest_volcano <- function(DT.original, treatment, control, log_base, lfc_threshold, fdr_threshold, out_dir, labelgene) {
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

###ttest pathway analysis
phospho_ttest_enrich_pathway = function(DT.original, treatment, control, outdir, lfc_threshold, fdr_threshold, enrich_pvalue){
  ##create dir
  dir=paste0(outdir,'/', treatment, '_vs_', control)
  if (!dir.exists(dir)){
    dir.create(dir,recursive = T)
  }
  
  DT <- data.table(DT.original)
  DT[, 'Group' := 'Others']
  DT[log2FC >= lfc_threshold, 'Group' := 'UP']
  DT[log2FC <= -lfc_threshold, 'Group' := 'DOWN']
  DT[p.adj >= fdr_threshold, 'Group' := 'Others']
  
  ##ENTREZID names########
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
      barplot(ego_cc, showCategory=20,label_format = 100)+
        scale_fill_gradient(low = "#ef8a62", high = "#67a9cf", na.value = NA)
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
      barplot(ego_bp, showCategory=20,label_format = 100)+
        scale_fill_gradient(low = "#ef8a62", high = "#67a9cf", na.value = NA)
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
      barplot(ego_mf, showCategory=20,label_format = 100)+
        scale_fill_gradient(low = "#ef8a62", high = "#67a9cf", na.value = NA)
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

#limma
phospho_limma = function(Log2_DT, design_matrix,DE_dir,lfc_threshold,fdr_threshold,enrich_pvalue) {
  name=names(Log2_DT)[sapply(Log2_DT, function(x) !all(is.numeric(x)))]
  rownames(Log2_DT)=Log2_DT$Precursor
  #comparation
  conditions <- unique(design_matrix$condition)
  for (treatment in conditions) {
    control=unique(design_matrix$control[design_matrix$condition==treatment])
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
    design <- model.matrix(~0+group_list)
    colnames(design) <- levels(group_list)
    rownames(design) <- colnames(DT_limma)
    cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = design)
    #limma
    fit <- lmFit(DT_limma, design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, trend=TRUE)
    
    result_limma <- topTable(fit2, coef=1,n=Inf)
    result_limma=merge(Log2_DT[,c(name,treatment_samples, control_samples)],result_limma,by.x='Precursor',by.y=0)
    ezwrite(result_limma[order(result_limma$adj.P.Val),], DE_dir, paste0(treatment, '_vs_', control, '.tsv'))
  }
  limma_plot_volcano(DT_limma_res = result_limma ,
                     lfc_threshold = lfc_threshold,
                     fdr_threshold = fdr_threshold,
                     out_dir = DE_dir,
                     output_filename =paste0(treatment, ' vs ', control) )
  enrich_pathway(DT.original = result_limma,
                 treatment = treatment,
                 control = control,
                 outdir = DE_dir,
                 lfc_threshold = lfc_threshold,
                 fdr_threshold = fdr_threshold,
                 enrich_pvalue = enrich_pvalue)
}
###limma plot_volcano

###limma pathway analysis
