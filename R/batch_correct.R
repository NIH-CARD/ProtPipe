#batch_correct
library("pvca")
library("Biobase")
library("limma")
library("data.table")
library("ggplot2")
library('dplyr')
library('ggpubr')



#
#
#
# # data_pro=fread(data)
# # data_pro=standardize_format(data_pro)
# # colnames(data_pro)=trim_colnames(data_pro)
# #
# # condition=read.csv('slam.muscle.samples.csv')
# # condition$Sample.Identification=gsub('S','',condition$Sample.Identification)
# # condition = condition[, -c(1)]
#
# format_data <- function(assay, phenotype, sample_col_name = "Sample.Identification"){
#   assay[] <- lapply(df, function(x) {
#     if (is.numeric(x)) {
#       x[is.nan(x)] <- 0
#     }
#     return(x)
#   })
#
#   # sort rows and cols so that they match up
#   assay <- assay[ ,sort(names(assay))]
#   phenotype <- phenotype[sort(rownames(phenotype)), ]
#
#   #remove samples that don't appear in both tables
#   matching_names <- intersect(rownames(phenotype), colnames(assay))
#   phenotype_matched <- phenotype[rownames(phenotype) %in% matching_names, ]
#   assay_matched <- assay[, colnames(assay) %in% matching_names, drop=FALSE]
#
#   return(list(assay = assay_matched, phenotype=phenotype_matched))
# }
#
# # combines assay and phenotype data into an expressionSet for PVCA
# combined_data <- function(assay, phenotype){
#   phenotype_annotated <- new("AnnotatedDataFrame", data = phenotype)
#   assay_matrix <- as.matrix(assay)
#   ExpressionSet(assay_matrix, phenotype_annotated)
# }
# ## get PVCA ####################################################################
# getPVCA <- function(assay, phenotype, batchcols, biocols, pct_threshold=0.6){
#   c <- combined_data(assay, phenotype)
#   pvcaObj <- pvcaBatchAssess (c, c(batchcols, biocols), pct_threshold)
#
#   #remove interaction terms
#   ind <- which(sapply(pvcaObj$label, function(x) x %in% c(batchcols, biocols, "resid")))
#
#   pvca_df <<- data.frame(
#     Category = pvcaObj$label[ind],
#     Value = c(pvcaObj$dat[ind])
#   )
#   return(pvca_df)
# }
#
# ## plot PVCA ####################################################################
# plotPVCA <- function(df){
#   g <- ggbarplot(df, x = "Category", y = "Value",
#                  fill = "Category", ylab = "Weighted Avegerage Proportion Variance")+
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#
#   return(g)
# }
#
# ## PCA ####################################################################
# plot_batch_pca <- function(assay, phenotype, category){
#
#   #add dummy cols because get_PCs gets rid of first 2 cols
#   columns_to_add <- assay[, c(1, 2)]
#   before_pca <- cbind(columns_to_add, assay)
#   pca <- get_PCs(before_pca)
#
#   #makes continuous legend for age (ex: 12M, 36M)
#   if (str_detect(phenotype[[category]][1], "^\\d+M$")){
#     d_batches <- as.numeric(gsub('M','',(phenotype[[category]])))
#   }else{
#     d_batches = sort(factor(phenotype[[category]]))
#   }
#
#   p <- ggplot(pca$components, aes(x = PC1, y = PC2, color = d_batches)) +
#     geom_point(size=4) +
#     xlab(paste0("PC1","(",pca$summary$percent[1],"%)")) +
#     ylab(paste0("PC2","(",pca$summary$percent[2],"%)")) +
#     labs(color = category) +
#     theme_classic()
#
#   return(p)
# }
#
# ## Batch correct ####################################################################
# correct_batches <- function(assay, phenotype,
#                             batchcols = c("Acquisition.Batch", "Digestion.Batch", "MS.Batch"),
#                             biocols = c("Biological.Condition1", "Biological.Condition2", "Biological.Condition3")){
#
#   formatted <- format_data(assay, phenotype)
#
#   #remove empty cols
#   emptycols = c()
#   for (col in c(batchcols, biocols)){
#     if (all(is.na(formatted$phenotype[[col]]))){
#       emptycols <- c(emptycols, col)
#     }
#   }
#   batchcols <<- batchcols[!(batchcols %in% emptycols)]
#   biocols <<- biocols[!(biocols %in% emptycols)]
#
#   #determine how to batch correct based on PVCA values
#   high_var_cols <- character()
#
#   assayyy <<- formatted$assay
#   phenotype <<- formatted$phenotype
#
#   pvca_df <- getPVCA(formatted$assay, formatted$phenotype, batchcols, biocols)
#   pvca_df <- pvca_df[order(pvca_df$Value), ]
#
#   bio_variance <- sum(pvca_df$Value[pvca_df$Category %in% biocols])
#   for (col in batchcols){
#     v <- pvca_df[pvca_df$Category == col, "Value"]
#     print(paste(col, " variance = ", v))
#     if (v > bio_variance){
#       high_var_cols <- c(high_var_cols, col)
#       print(paste("high var: ", high_var_cols))
#     }
#   }
#   if(length(high_var_cols) == 0){
#     print("no significant batch effect observed and no correction applied")
#     return(NULL)
#   }else{
#     print(paste("Applying batch correction for the following condition(s): ", high_var_cols))
#   }
#
#   combined_row <<- apply(formatted$phenotype[high_var_cols], 1, function(row) paste(row, collapse = ""))
#
#   # create pca and pvca plots
#   pvca_before <- plotPVCA(pvca_df)
#   pca_before <- list()  # Initialize an empty list to store the plots
#   for (col in c(batchcols, biocols)) {
#     pca_before[[col]] <- plot_batch_pca(formatted$assay, formatted$phenotype, col)
#   }
#
#   formula <- as.formula(paste("~", paste(biocols, collapse = " + ")))
#   experimental_matrix <- model.matrix(formula, data = formatted$phenotype)
#   #experimental_matrix <<- model.matrix((biocols), formatted$phenotype)
#   log2_dat <- data.frame(formatted$assay) %>%
#     mutate(across(, ~ as.numeric(.))) %>%
#     mutate(across(, ~ replace(., is.na(.), 0))) %>%
#     mutate(across(, ~ log2(. + 1)))
#
#   log2_dat_batch=removeBatchEffect(log2_dat,
#                                    batch = combined_row,
#                                    design = experimental_matrix
#   )
#
#   #test this part
#   dat_batch=as.data.frame(2^(log2_dat_batch-1))
#   print(paste("neg values: ", sum(apply(dat_batch, 2, function(col) sum(col < 0)))))
#
#   ##WHY is there X in front of all col names now????
#   colnames(dat_batch) <- sub("^X", "", colnames(dat_batch))
#
#   pvca_df_after <- getPVCA(dat_batch, formatted$phenotype, batchcols, biocols)
#   #order the categories in the same way as pvca before batch correct
#   pvca_df_after <- pvca_df_after[order(factor(pvca_df_after$Category, levels = pvca_df$Category)), ]
#   pvca_after <- plotPVCA(pvca_df_after)
#   pca_after <- list()
#   for (col in c(batchcols, biocols)) {
#     pca_after[[col]] <- plot_batch_pca(dat_batch, formatted$phenotype, col)
#   }
#
#   return(list(data = dat_batch, pvca_before = pvca_before, pca_before = pca_before,
#               pvca_after = pvca_after, pca_after = pca_after, batchcols = batchcols,
#               biocols = biocols))
# }
#
# ## Save Batch correct info ############################################################
# save_batch <- function(batch_stuff, outdir){
#   ggsave(filename = paste0(outdir, "pcva_before.pdf"), plot = batch_stuff$pvca_before, width = 8, height = 6, dpi = 300)
#   ggsave(filename = paste0(outdir, "pvca_after.pdf"), plot = batch_stuff$pvca_after, width = 8, height = 6, dpi = 300)
#   columns <- c(batch_stuff$batchcols, batch_stuff$biocols)
#   for (col in columns){
#     ggsave(batch_stuff$pca_before[[col]],filename=paste0(outdir, col, '_pca_before.pdf'), height = 6,width = 5)
#     ggsave(batch_stuff$pca_after[[col]],filename=paste0(outdir, col, '_pca_after.pdf'), height = 6,width = 5)
#   }
# }
