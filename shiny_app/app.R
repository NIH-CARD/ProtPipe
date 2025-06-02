devtools::load_all()
library(shiny)
library(bslib)

# increase file size limit
options(shiny.maxRequestSize = 30*1024^2)

package_list = c('ggplot2', 'data.table', 'corrplot', 'umap', 'magick', 'ggdendro', 'ecodist','ggbeeswarm',
                 'ggrepel', 'ggthemes', 'foreach','reshape2','org.Hs.eg.db','clusterProfiler','pheatmap')
if(all((lapply(package_list, require, character.only=TRUE)))) {
  cat("INFO: All packages successfully loaded\n")
} else {
  cat("ERROR: One or more packages not available. Are you running this within the container?\n")
}

ui <- page_sidebar(
  title = "ProtPipe",
  sidebar = sidebar(
    fluidPage(
      card(
        card_header("Select plots to view"),
        selectInput("select", label = h3("View"),
                    choices = list("Input parameters" = 0, "Quality Control" = 1, "Clustering" = 2,
                                   "Differential Intensity" = 3, "Heatmap" = 4),
                    selected = 0),
        hr(),
        fluidRow(column(3, verbatimTextOutput("value")))
      )
    )
  ),

  ### Parameter input screen ############################################################################################
  conditionalPanel(condition = "input.select == 0",
                   fluidPage(
                     fluidRow(
                       column(width = 6,
                              card(
                                card_header("Upload protein intensity estimates csv file"),
                                fileInput("intensity", label = NULL),
                                  uiOutput("column_range_ui"),
                                  verbatimTextOutput("range_result")

                              )), column(width = 6,
                                         card(
                                           card_header("Upload sample condition csv file"),
                                           p("make sure row names match the column names of the intensity file exactly"),
                                           fileInput("sample_condition", label = NULL)
                                         ),
                              )),
                     card(
                       h3("Parameters"),
                       selectInput("normalize", label = "Normalize", choices = c("none", "shift", "scale"), selected = "none"),
                       textInput("exclusions", label = "Enter a semicolon-separated string of files to exclude from analysis",
                                 value = ""),
                       numericInput("minintensity", label = "Enter Minimum LINEAR (not log) intensity.", value = 0),
                       numericInput("sds", label = "Filter out samples with protein group counts > N standard deviations
                       from the mean", value = 3)
                     )
                   )
  ),


  ### Quality control screen ############################################################################################
  conditionalPanel(condition = "input.select == 1",
                   h2("Quality Control Information"),
                   fluidPage(
                     uiOutput("quality_control_condition"),
                     card(card_header("Coefficients of Variation (requires condition file)"),
                          plotOutput("cv_graph"),
                          downloadButton("download_cv", "Download Plot as PDF")
                     ),
                     card(card_header("Intensities"),
                          plotOutput("intensity_graph"),
                          downloadButton("download_intensity", "Download Plot as PDF")
                     ),
                     card(card_header("Non-zero Protein Group Counts"),
                          plotOutput("pgroup_graph"),
                          downloadButton("download_pg", "Download Plot as PDF"),
                          downloadButton("download_pg_tsv", "Download data as tsv")
                     ),
                     card(card_header("Correlation Heatmap"),
                          plotOutput("correlation_graph"),
                          downloadButton("download_cor", "Download Plot as PDF"),
                          downloadButton("download_cor_tsv", "Download data as tsv")
                     )
                   )
  ),

  ### Clustering screen ############################################################################################
  conditionalPanel(condition = "input.select == 2",
                   h2("Clustering Information"),
                   fluidPage(
                     uiOutput("clustering_condition"),
                     card(card_header("heirarchial clustering"),
                          plotOutput("hcluster"),
                          downloadButton("download_hcluster", "Download Plot as PDF")
                     ),card(uiOutput("cluster_groups")
                     ),card(card_header("PCA (requires condition file)"),
                            plotOutput("pca"),
                            downloadButton("download_pca", "Download Plot as PDF"),
                            downloadButton("download_pca_tsv", "Download PCA as tsv"),
                            downloadButton("download_pca_sum", "Download PCA summary as tsv")
                     ),card(card_header("UMAP (requires condition file)"),
                            fluidRow(
                              column(width = 3, div(style = "display: flex; align-items: center; height: 100%;",
                                                    sliderInput("neighbors", label = h3("Select number of neighbors for UMAP"),
                                                                min = 0, max = 100, value = 15))
                              ),column(width = 9,
                                       plotOutput("umap")
                              )
                            ),
                            downloadButton("download_umap", "Download Plot as PDF"),
                            downloadButton("download_umap_tsv", "Download PCA as tsv")
                     )
                   )
  ),
  ### Differential Intensity ############################################################################################
  conditionalPanel(condition = "input.select == 3",
                   h2("Differential Expression"),
                   card(card_header("Options"),
                     fluidPage(
                       fluidRow(
                         column(width = 4,
                                uiOutput("de_condition"),
                                uiOutput("de_groups")),
                         column(width = 4,
                                numericInput("logfc", label = "Enter log2 fold-change cutoff", value = 1),
                                numericInput("pvalue", label = "Enter pvalue cutoff", value = 0.01)),
                         column(width = 4,
                                uiOutput("label_col"),
                                p("optional: upload csv file with genes to label in volcano plot.
                                   Make sure column name is Genes. This must contain values present in the column
                                  selected above"),
                                fileInput("gene_labels", label = NULL))
                       )
                      )
                     ),
                     card(card_header("Volcano Plot"),
                          plotOutput("volcano"),
                          downloadButton("download_volcano", "Download Plot as PDF"),
                          downloadButton("download_DE_tsv", "Download differential expression results as tsv")
                     ),
                     card(card_header("Pathway Enrichment Options"),
                          fluidRow(
                            column(width=3,
                                  numericInput("enrich_pval", label = "Enter enrichment pvalue cutoff", value = 1)
                                   ),
                            column(width=3,
                                   selectInput(
                                     inputId = "organism",
                                     label = "Select Organism:",
                                     choices = c(
                                       "Human",
                                       "Mouse",
                                       "Nematode(C. elegans)",
                                       "Fruit fly"
                                     ),
                                     selected = "Human"  # default selection
                                   )
                            ),
                            column(width=3,
                                   uiOutput("gene_col")
                            ),
                            column(width=3,
                                   checkboxInput("run_enrichment", "Perform GO and KEGG enrichment (this may take a few minutes)", value = FALSE)
                            )
                          )
                          ),
                     card(card_header("Pathway Enrichment"),
                          fluidRow(
                          column(width=6 ,card(card_header("Upregulated GO Pathways"),plotOutput("go_up_enrich"))),
                          column(width=6 ,card(card_header("Upregulated KEGG Pathways"),plotOutput("kegg_up_enrich")))
                          ),
                          fluidRow(
                            column(width=6 ,card(card_header("Downregulated GO Pathways"),plotOutput("go_down_enrich"))),
                            column(width=6 ,card(card_header("Downregulated KEGG Pathways"),plotOutput("kegg_down_enrich")))
                          ),
                          fluidRow(
                          column(width=6 ,card(card_header("GO Gene Set Enrichment"),plotOutput("go_gsea"))),
                          column(width=6 ,card(card_header("Pathway Enrichment"),plotOutput("kegg_gsea")))
                          ),
                          downloadButton("download_enrichment", "Download pathway enrichment results")
                     )

  ),

  ### Heatmap ############################################################################################
  conditionalPanel(condition = "input.select == 4",
                   h2("Heatmap"),
                   fluidPage(
                     uiOutput("protein_label"),
                     p("optional: upload csv file with genes to include in heatmap subset.
                                   Make sure column name is Genes"),
                     fileInput("heatmap_labels", label = NULL),
                     card(card_header("Heatmap"),
                          plotOutput("h_map"),
                          downloadButton("download_hmap", "Download Plot as PDF")
                     )
                   )
                 )
  )


server <- function(input, output, session) {


  #### Reactive functions ############################################################################################

  #if sample conditions are provided, the data is reformatted so that the intensity columns
  #match the condition rows

  intensity_file <- reactive({
    req(input$intensity)
    return(data.table::fread(input$intensity$datapath, data.table=FALSE))
  })
  condition_file <- reactive({
    req(input$sample_condition)
    return(data.table::fread(input$sample_condition$datapath, data.table=FALSE))
  })
  # Dynamically generate dropdowns for column range selection
  output$column_range_ui <- renderUI({
    req(input$intensity)
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
    req(input$intensity)
    df <- intensity_file()

    # Ensure both selections are made
    req(input$lower_col, input$upper_col)

    cols <- names(df)
    lower_idx <- match(input$lower_col, cols)
    upper_idx <- match(input$upper_col, cols)

    if (!is.null(input$sample_condition)) {
       return(ProtPipe::create_protdata(dat = intensity_file(), intensity_cols = c(lower_idx:upper_idx), condition = condition_file()))
    }else{
      return(ProtPipe::create_protdata(dat = intensity_file(), intensity_cols = c(lower_idx:upper_idx)))
    }
  })

  exclusions <- reactive({
    if(input$exclusions == ""){
      exclusions <- NULL
    }else{
      exclusions <- strsplit(input$exclusions, split=';')[[1]]
    }
    exclusions
  })

  foldchange <- reactive({
    input$foldchange})


  geneHeatmap <- reactive({
    if(!is.null(input$heatmap_labels)){
      file <- fread(input$heatmap_labels$datapath)
      #hmap_labels <- as.vector(file$PG.Genes)
    }else{
      file <- NULL
    }
    return(file)
  })



  ### Color Pallete ############################################################################################


  ### Batch Correction ############################################################################################





  #### QC ############################################################################################

  #select condition
  output$quality_control_condition <- renderUI({
    req(input$intensity)
    req(input$sample_condition)

    choices <- names(prot_data()@condition)

    selectInput("qc_condition", "select condition to group by:", choices = choices)
  })

  # CV plot
  output$cv_graph <- renderPlot({
    req(input$intensity)
    req(input$sample_condition)
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
    req(input$intensity)
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
    req(input$intensity)
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
    req(input$intensity)
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
    req(input$intensity)
    req(input$sample_condition)

    choices <- names(prot_data()@condition)

    selectInput("cluster_condition", "select condition to group by:", choices = choices)
  })

  #hierarchical clustering
  output$hcluster <- renderPlot({
    req(input$intensity)  # Ensure file is uploaded
    ProtPipe::plot_hierarchical_cluster(prot_data())
  })

  output$download_hcluster <- downloadHandler(
    filename = function(){
      paste("hierarchical_clustering.pdf")
    },
    content = function(file){
      req(input$intensity)  # Ensure file is uploaded
      p <-  ProtPipe::plot_hierarchical_cluster(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  #PCA
  output$pca <- renderPlot({
    req(input$intensity)  # Ensure file is uploaded
    req(input$sample_condition)
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
      req(input$intensity)  # Ensure file is uploaded
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
      req(input$intensity)  # Ensure file is uploaded
      dat <- clean_dat()$dat
      pca <- get_PCs(dat)
      simplewrite(pca$summary, file)
    }
  )

  #UMAP
  output$umap <- renderPlot({
    req(input$intensity)  # Ensure file is uploaded'
    req(input$sample_condition)
    ProtPipe::plot_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster_condition)
  })

  output$download_umap <- downloadHandler(
    filename = function(){
      paste("umap.pdf")
    },
    content = function(file){
      req(input$intensity)  # Ensure file is uploaded
      p <- ProtPipe::plot_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster$condition)
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_umap_tsv <- downloadHandler(
    filename = function(){
      paste("umap.tsv")
    },
    content = function(file){
      req(input$intensity)  # Ensure file is uploaded
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
    req(input$intensity)
    choices <- names(prot_data()@prot_meta)
    selectInput("protein_label", "select the column used to label proteins:", choices = choices)
  })

  #heatmap subset
  prot_labels <- reactive({
    req(input$intensity)  # Ensure file is uploaded
    req(input$protein_label)  # Ensure file is uploaded
    if(is.null(input$heatmap_labels)){
      return(NULL)
    }
    dat <- data.table::fread(input$heatmap_labels$datapath, data.table=FALSE)
    if('Gene' %in% names(dat)){
      print("The csv must contain a column called Gene containing the labels")
    }
    return(dat$Gene)
  })

  #complete heatmap
  output$h_map <- renderPlot({
    req(input$intensity)  # Ensure file is uploaded
    p <- ProtPipe::plot_proteomics_heatmap(prot_data(), protmeta_col = input$protein_label, genes = prot_labels())
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
  })

  output$download_hmap <- downloadHandler(
    filename = function(){
      paste("heatmap.pdf")
    },
    content = function(file){
      req(input$intensity)  # Ensure file is uploaded
      p <- ProtPipe::plot_proteomics_heatmap(prot_data(), protmeta_col = input$protein_label, genes = prot_labels())
      ggsave(file, plot=p, device = "pdf")
    }
  )
  #### Differential Intensity ############################################################################################

  #select condition
  output$de_condition <- renderUI({
    req(input$intensity)
    choices <- names(prot_data()@condition)
    selectInput("de_condition", "select the column used to compare groups:", choices = choices)
  })

  #select groups
  output$de_groups <- renderUI({
    req(input$intensity)
    req(input$de_condition)
    groups <- unique(prot_data()@condition[[input$de_condition]])

    tagList(
      selectInput("control_condition", "select the control groups:", choices = groups),
      selectInput("treatment_condition", "select the treatment groups:", choices = groups)
    )
  })

  #select column to label the proteins
  output$label_col <- renderUI({
    req(input$intensity)
    choices <- names(prot_data()@prot_meta)
    selectInput("label_col", "select the column used to label proteins:", choices = choices)
  })

  #custom labels for volcano
  gene_labels <- reactive({
    req(input$intensity)  # Ensure file is uploaded
    if(is.null(input$gene_labels)){
      return(NULL)
    }
    dat <- data.table::fread(input$gene_labels$datapath, data.table=FALSE)

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
    req(input$intensity)
    ProtPipe::plot_volcano(dea(), label_col = input$label_col, labelgene = gene_labels(), fdr_threshold = input$pvalue, lfc_threshold = input$logfc)
  })

  output$download_volcano <- downloadHandler(
    filename = function(){
      paste("volcano.pdf")
    },
    content = function(file){
      req(input$intensity)  # Ensure file is uploaded
      p <- ProtPipe::plot_volcano(dea(), label_col = input$label_col, labelgene = gene_labels(), fdr_threshold = input$pvalue, lfc_threshold = input$logfc)
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_DE_tsv <- downloadHandler(
    filename = function(){
      paste("differential_expression_results.tsv")
    },
    content = function(file){
      req(input$intensity)  # Ensure file is uploaded
      write.table(dea(), file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  #select condition
  output$gene_col <- renderUI({
    req(input$intensity)
    choices <- names(prot_data()@prot_meta)
    selectInput("gene_col", "select the column containing official gene symbols (e.g., TP53):", choices = choices)
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
      result <- ProtPipe::enrich_pathways(dea())
      enrichment_result(result)

      # Remove notification
      removeNotification("enrich_msg")

      # Re-enable button so user can rerun if needed
      shinyjs::enable("run_enrichment")
    }
  })

  #pathway enrichment plots
  output$go_up_enrich <- renderPlot({
    req(input$intensity)
    req(enrichment_result())
    enrichment_result()$plots$go_up_dotplot
  })

  output$kegg_up_enrich <- renderPlot({
    req(input$intensity)
    req(enrichment_result())
    enrichment_result()$plots$kegg_up_dotplot
  })
  output$go_down_enrich <- renderPlot({
    req(input$intensity)
    req(enrichment_result())
    enrichment_result()$plots$go_down_dotplot
  })

  output$kegg_down_enrich <- renderPlot({
    req(input$intensity)
    req(enrichment_result())
    enrichment_result()$plots$kegg_down_dotplot
  })
  output$go_gsea <- renderPlot({
    req(input$intensity)
    req(enrichment_result())
    enrichment_result()$plots$gse_go_dotplot
  })

  output$kegg_gsea <- renderPlot({
    req(input$intensity)
    req(enrichment_result())
    enrichment_result()$plots$gse_kegg_dotplot
  })


  output$download_DI <- downloadHandler(
    filename = function(){
      "Differential Intensity.zip"
    },
    content = function(file){

      req(input$intensity)  # Ensure file is uploaded
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

  output$download_EN <- downloadHandler(
    filename = function(){
      "Enrichement Analysis.zip"
    },
    content = function(file){

      req(input$intensity)  # Ensure file is uploaded
      req(input$design)  # Ensure file is uploaded

      #create a temporary directory to save the volcano plots
      ea_directory <- file.path(tempdir(), as.integer(Sys.time()))
      dir.create(ea_directory)


      #need to add enrich pvalue input
      d <- process_design(input$design$datapath)
      des <- d$design
      conditions <- d$conditions
      control <- d$control
      des <- exclude_design(des, exclusions())
      dat <- clean_dat()$dat
      enrichment_analysis(conditions, control, input$de_method, des, dat, ea_directory, input$log_base,
                          input$foldchange, input$fdr, labelgene(), input$enrich_pvalue)

      zip::zip(
        zipfile = file,
        files = dir(ea_directory),
        root = ea_directory
      )
    }
  )

}




shinyApp(ui, server)
