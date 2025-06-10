


server <- function(input, output, session) {

  #read file uploads
  intensity <- fileUploadServer("intensity")
  sample_condition <- fileUploadServer("sample_condition")
  gene_labels_file <- fileUploadServer("gene_labels")
  heatmap_labels <- fileUploadServer("heatmap_labels")

  #### Reactive functions ############################################################################################

  #if sample conditions are provided, the data is reformatted so that the intensity columns
  #match the condition rows

  intensity_file <- reactive({
    req(input$use_example || !is.null(intensity()))
    if(input$use_example){
      return(data.table::fread("www/iPSC.csv", data.table=FALSE))
    }else{
      return(data.table::fread(intensity()$datapath, data.table=FALSE))
    }
  })
  condition_file <- reactive({
    req(sample_condition())
    return(data.table::fread(sample_condition()$datapath, data.table=FALSE))
  })
  # Dynamically generate dropdowns for column range selection
  output$column_range_ui <- renderUI({
    req(intensity_file())
    df <- intensity_file()
    choices <- names(df)

    # get default intensity columns
    intensity_cols <- detect_intensity_cols(df)
    intensity_cols <<- intensity_cols
    first <- intensity_cols[[1]]
    last <- intensity_cols[[length(intensity_cols)]]

    tagList(
      selectInput("lower_col", "Intensity columns start at:", choices = choices, selected = choices[first]),
      selectInput("upper_col", "Intensity columns end at:", choices = choices, selected = choices[last])
    )
  })

  # Validate and report selection
  output$range_result <- renderPrint({
    df <- intensity_file()

    # Ensure both selections are made
    req(input$lower_col, input$upper_col)

    cols <- names(df)
    lower_idx <- match(input$lower_col, cols)
    upper_idx <- match(input$upper_col, cols)

    if (is.na(lower_idx) || is.na(upper_idx)) {
      return("❌ Column not found in file.")
    }

    # Ensure lower <= upper
    if (lower_idx > upper_idx) {
      return("❌ First column must come before or be the same as the last column.")
    }

    selected <- df[, lower_idx:upper_idx, drop = FALSE]

    if (!all(sapply(selected, is.numeric))) {
      return("❌ All selected columns must be numeric.")
    }

    paste("✅ Selected columns:", input$lower_col, "to", input$upper_col,
          "| Count:", ncol(selected))
  })

  prot_data <- reactive({
    req(intensity_file())
    df <- intensity_file()

    # Ensure both selections are made
    req(input$lower_col, input$upper_col)

    cols <- names(df)
    lower_idx <- match(input$lower_col, cols)
    upper_idx <- match(input$upper_col, cols)

    if (!is.null(sample_condition())) {
      PD <- ProtPipe::create_protdata(dat = intensity_file(), intensity_cols = c(lower_idx:upper_idx), condition = condition_file())
    }else{
      PD <- ProtPipe::create_protdata(dat = intensity_file(), intensity_cols = c(lower_idx:upper_idx))
    }
    PD0 <<- PD
    #1 outlier removal
    if(input$remove_outliers == TRUE){
      PD <- ProtPipe::remove_outliers(PD, sds = input$outlier_sds)
    }

    #2 normalization
    if(input$normalize == TRUE){
      print(paste("Normalizing using", input$normalize_method))
      tryCatch({
        if (input$normalize_method == "mean") {
          PD <- ProtPipe::mean_normalize(PD)
        } else if (input$normalize_method == "median") {
          PD <- ProtPipe::median_normalize(PD)
        }
        PD1 <<- PD
      }, error = function(e) {
        print("Normalization failed")
        print(e)
      })
    }

    #3 imputation
    if(input$impute == TRUE){
      if(input$imputation_method == "zero"){
        PD <- ProtPipe::impute(PD, 0)
      }else if(input$imputation_method == "minimum"){
        PD <- ProtPipe::impute_min(PD, 1)
      }else if(input$imputation_method == "left-shifted distribution"){
        PD <- ProtPipe::impute_minimal(PD)
      }
    }

    #4 batch correction
    if(!is.null(input$batch_correct_column) && input$batch_correct == TRUE){
      PD <- ProtPipe::batch_correct(PD, input$batch_correct_column)
    }
    PD2 <<- PD
    return(PD)

  })



  ### Color Pallete ############################################################################################


  ### Batch Correction ############################################################################################

  output$batch_correct_column <- renderUI({
    req(intensity_file())
    req(sample_condition())

    choices <- names(prot_data()@condition)

    selectInput("batch_correct_column", "select condition for correction:", choices = choices)
  })



  #### QC ############################################################################################

  #select condition
  output$quality_control_condition <- renderUI({
    req(intensity_file())
    #req(sample_condition())

    choices <- names(prot_data()@condition)

    selectInput("qc_condition", "select condition to group by:", choices = choices)
  })

  # CV plot
  output$cv_graph <- renderPlot({
    req(intensity_file())
    #req(sample_condition())
    req(input$qc_condition)
    ProtPipe::plot_CVs(prot_data(), condition = input$qc_condition)
  })

  output$download_intensity <- downloadHandler(
    filename = function(){
      paste("intensities.pdf")
    },
    content = function(file){
      p <- ProtPipe::plot_pg_intensities(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  # intensity graph
  output$intensity_graph <- renderPlot({
    req(intensity_file())
    ProtPipe::plot_pg_intensities(prot_data())
  })

  output$download_intensity <- downloadHandler(
    filename = function(){
      paste("intensities.pdf")
    },
    content = function(file){
      p <- ProtPipe::plot_pg_intensities(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  # protein group counts
  output$pgroup_graph <- renderPlot({
    req(intensity_file())
    ProtPipe::plot_pg_counts(prot_data())
  })

  output$download_pg <- downloadHandler(
    filename = function(){
      paste("protein_groups_nonzero_counts.pdf")
    },
    content = function(file){
      p <- ProtPipe::plot_pg_counts(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_pg_tsv <- downloadHandler(
    filename = function(){
      paste("protein_group_nonzero_counts.tsv")
    },
    content = function(file){
      pgcounts <- normalized_file()$long[, .N, by=Sample]
      simplewrite(pgcounts, file)
    }
  )

  #correlation heatmap
  output$correlation_graph <- renderPlot({
    req(intensity_file())
    ProtPipe::plot_correlation_heatmap(prot_data())
  })

  output$download_cor <- downloadHandler(
    filename = function(){
      paste("sample_correlation.pdf")
    },
    content = function(file){
      ProtPipe::plot_correlation_heatmap(prot_data())
    }
  )

  output$download_cor_tsv <- downloadHandler(
    filename = function(){
      paste("sample_correlation.tsv")
    },
    content = function(file){
      dat.correlations <- get_spearman(normalized_file()$dat)
      simplewrite(dat.correlations, file)
    }
  )

  #### Clustering ############################################################################################

  #select condition
  output$clustering_condition <- renderUI({
    req(intensity_file())
    #req(sample_condition())

    choices <- names(prot_data()@condition)

    selectInput("cluster_condition", "select condition to group by:", choices = choices)
  })

  #hierarchical clustering
  output$hcluster <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded
    ProtPipe::plot_hierarchical_cluster(prot_data())
  })

  output$download_hcluster <- downloadHandler(
    filename = function(){
      paste("hierarchical_clustering.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      p <-  ProtPipe::plot_hierarchical_cluster(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  #PCA
  output$pca <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded
    #req(sample_condition())
    ProtPipe::plot_PCs(prot_data(), condition = input$cluster_condition)
  })

  output$download_pca <- downloadHandler(
    filename = function(){
      paste("PCA.pdf")
    },
    content = function(file){
      p <- ProtPipe::plot_PCs(prot_data(), condition = input$cluster_condition)
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_pca_tsv <- downloadHandler(
    filename = function(){
      paste('PCA.tsv')
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      dat <- clean_dat()$dat
      pca <- get_PCs(dat)
      simplewrite(pca$components, file)
    }
  )
  output$download_pca_sum <- downloadHandler(
    filename = function(){
      paste('PCA_summary.tsv')
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      dat <- clean_dat()$dat
      pca <- get_PCs(dat)
      simplewrite(pca$summary, file)
    }
  )

  #UMAP
  output$umap <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded'
    #req(sample_condition())
    ProtPipe::plot_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster_condition)
  })

  output$download_umap <- downloadHandler(
    filename = function(){
      paste("umap.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      p <- ProtPipe::plot_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster$condition)
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_umap_tsv <- downloadHandler(
    filename = function(){
      paste("umap.tsv")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      dat <- clean_dat()$dat
      if ((ncol(dat)-3)>input$neighbors) {
        tryTo('INFO: running UMAP',{
          umap <- get_umap(dat, input$neighbors)
          simplewrite(umap, file)
        }, 'ERROR: failed!')
      }
    }
  )

  #### Heatmap ############################################################################################


  #select condition
  output$protein_label <- renderUI({
    req(intensity_file())
    choices <- names(prot_data()@prot_meta)
    selectInput("protein_label", "select the column used to label proteins:", choices = choices)
  })

  #heatmap subset
  prot_labels <- reactive({
    req(intensity_file())  # Ensure file is uploaded
    req(input$protein_label)  # Ensure file is uploaded
    if(is.null(heatmap_labels())){
      return(NULL)
    }
    dat <- data.table::fread(heatmap_labels()$datapath, data.table=FALSE)
    if('Gene' %in% names(dat)){
      print("The csv must contain a column called Gene containing the labels")
    }
    return(dat$Gene)
  })

  #complete heatmap
  output$h_map <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded
    p <- ProtPipe::plot_proteomics_heatmap(prot_data(), protmeta_col = input$protein_label, genes = prot_labels())
    print(p)
    # grid::grid.newpage()
    # grid::grid.draw(p$gtable)
  })

  output$download_hmap <- downloadHandler(
    filename = function(){
      paste("heatmap.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      p <- ProtPipe::plot_proteomics_heatmap(prot_data(), protmeta_col = input$protein_label, genes = prot_labels())
      ggsave(file, plot=p, device = "pdf")
    }
  )
  #### Differential Intensity ############################################################################################

  #select condition
  output$de_condition <- renderUI({
    req(intensity_file())
    choices <- names(prot_data()@condition)
    selectInput("de_condition", "select the column used to compare groups:", choices = choices)
  })

  #select groups
  output$de_groups <- renderUI({
    req(intensity_file())
    req(input$de_condition)
    groups <- unique(prot_data()@condition[[input$de_condition]])

    tagList(
      selectInput("control_condition", "select the control groups:", choices = groups),
      selectInput("treatment_condition", "select the treatment groups:", choices = groups)
    )
  })

  #select column to label the proteins
  output$label_col <- renderUI({
    req(intensity_file())
    choices <- names(prot_data()@prot_meta)
    selectInput("label_col", "select the column used to label proteins:", choices = choices)
  })

  #custom labels for volcano
  gene_labels <- reactive({
    req(intensity_file())  # Ensure file is uploaded
    if(is.null(gene_labels_file())){
      return(NULL)
    }
    dat <- data.table::fread(gene_labels_file()$datapath, data.table=FALSE)

    if (!'Gene' %in% names(dat)) {
      warning("The uploaded gene label file must contain a column called 'Gene'.")
      return(NULL)
    }

    return(dat$Gene)
  })

  dea <- reactive({
    df <- ProtPipe::log2_transform(prot_data())
    condition <- input$de_condition
    control_group <- input$control_condition
    treatment_group <- input$treatment_condition

    return(ProtPipe::do_limma_by_condition(df,condition = condition, control_group = control_group, treatment_group = treatment_group))
  })

  #volcano plot
  output$volcano <- renderPlot({
    req(intensity_file())
    ProtPipe::plot_volcano(dea(), label_col = input$label_col, labelgene = gene_labels(), fdr_threshold = input$pvalue, lfc_threshold = input$logfc)
  })

  output$download_volcano <- downloadHandler(
    filename = function(){
      paste("volcano.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      p <- ProtPipe::plot_volcano(dea(), label_col = input$label_col, labelgene = gene_labels(), fdr_threshold = input$pvalue, lfc_threshold = input$logfc)
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_DE_tsv <- downloadHandler(
    filename = function(){
      paste("differential_expression_results.tsv")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      write.table(dea(), file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  #select condition
  output$gene_col <- renderUI({
    req(intensity_file())
    choices <- names(prot_data()@prot_meta)
    selectInput("gene_col", "select the column containing official gene symbols (e.g., TP53):", choices = choices)
  })

  #go and kegg db for different organisms
  selected_org <- reactive({
    req(input$organism)
    organism_map[[input$organism]]
  })

  # Create a reactiveVal to store pathway enrichment results
  enrichment_result <- reactiveVal(NULL)

  observeEvent(input$run_enrichment, {
    if (isTRUE(input$run_enrichment)) {
      # Disable the checkbox/button (if checkboxInput used as button, or use actionButton)
      shinyjs::disable("run_enrichment")  # requires shinyjs package and call to use it in UI

      # Optionally show a notification
      showNotification("Running enrichment analysis, please wait...", duration = NULL, id = "enrich_msg")

      # Run the long function (blocking)
      result <- ProtPipe::enrich_pathways(dea(), lfc_threshold=input$logfc, fdr_threshold=input$pvalue, enrich_pvalue=input$enrich_pval, go_org = selected_org()$OrgDb, kegg_org = selected_org()$kegg, gene_col = input$gene_col)
      enrichment_result(result)

      # Remove notification
      removeNotification("enrich_msg")

      # Re-enable button so user can rerun if needed
      shinyjs::enable("run_enrichment")
    }
  })

  #pathway enrichment plots
  output$go_up_enrich <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$go_up_dotplot
  })

  output$kegg_up_enrich <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$kegg_up_dotplot
  })
  output$go_down_enrich <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$go_down_dotplot
  })

  output$kegg_down_enrich <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$kegg_down_dotplot
  })
  output$go_gsea <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$gse_go_dotplot
  })

  output$kegg_gsea <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$gse_kegg_dotplot
  })


  output$download_DI <- downloadHandler(
    filename = function(){
      "Differential Intensity.zip"
    },
    content = function(file){

      req(intensity_file())  # Ensure file is uploaded
      req(input$design)  # Ensure file is uploaded

      #create a temporary directory to save the volcano plots
      di_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(di_directory)


      #need to add enrich pvalue input
      d <- process_design(input$design$datapath)
      des <- d$design
      conditions <- d$conditions
      control <- d$control
      des <- exclude_design(des, exclusions())
      dat <- clean_dat()$dat
      print(class(dat))
      differential_analysis(conditions, control, input$de_method, des, dat, di_directory, input$log_base,
                            input$foldchange, input$fdr, labelgene(), input$enrich_pvalue)

      zip::zip(
        zipfile = file,
        files = dir(di_directory),
        root = di_directory
      )
    }
  )

  output$download_enrichment <- downloadHandler(
    filename = function(){
      "Enrichement_Analysis.zip"
    },
    content = function(file){

      req(intensity_file())
      req(enrichment_result())

      #create a temporary directory to save the volcano plots
      output_dir <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(output_dir)

      # Save data frames
      for (name in names(enrichment_result()$results)) {
        df <- enrichment_result()$results[[name]]
        if (!is.null(df)) {
          write.table(
            df,
            file = file.path(output_dir, paste0(name, ".tsv")),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE
          )
        }
      }

      # # Save plots
      # for (name in names(enrichment_result()$plots)) {
      #   p <- enrichment_result()$plots[[name]]
      #   if (!is.null(p)) {
      #     pdf(file.path(output_dir, paste0(name, ".pdf")), width = 1, height = 1)
      #     print(p)
      #     dev.off()
      #   }
      # }

      zip::zip(
        zipfile = file,
        files = dir(output_dir),
        root = output_dir
      )
    }
  )


}





