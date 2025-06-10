source("helpers.R")


ui <- page_sidebar(
  title = tagList(h1("ProtPipe", style = "margin-bottom: 0;"),
  tags$div("Shiny app made by Jacob Epstein", style = "font-size: 0.9em; color: #666; margin-top: 0.2em;")),
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
                                fileUploadUI("intensity", label = NULL),
                                checkboxInput("use_example", "Or use our iPSC to neuron differentiation example dataset", value = FALSE),
                                uiOutput("column_range_ui"),
                                verbatimTextOutput("range_result")

                              )), column(width = 6,
                                         card(
                                           card_header("Upload sample condition csv file"),
                                           p("make sure row names match the column names of the intensity file exactly"),
                                           fileUploadUI("sample_condition", label = NULL)
                                         ),
                              )),
                     card(
                       h3("Parameters"),
                       fluidRow(
                         column(width = 6,
                         card(h4("Normalization"),
                              checkboxInput("normalize", label = "normalize", value = FALSE),
                              selectInput("normalize_method", label = "normalize_method", choices = c("mean", "median"), selected = "median"))),
                         column(width = 6,
                         card(h4("Outlier Removal"),
                              checkboxInput("remove_outliers", label = "remove outliers", value = FALSE),
                              numericInput("outlier_sds", label = "Remove samples with protein groups outside n standard deviations from the mean", value = 3)))
                       ),
                       fluidRow(
                         column(width = 6,
                         card(h4("Batch Correction"),
                              checkboxInput("batch_correct", label = "batch correct", value = FALSE),
                              uiOutput("batch_correct_column"),)),
                         column(width = 6,
                         card(h4("Imputation"),
                              checkboxInput("impute", label = "impute", value = FALSE),
                              selectInput("imputation_method", label = "imputation method", choices = c("zero", "minimun", "left-shifted distribution"), selected = "zero")))
                       ),
                       downloadButton("download_data", "Download processed data")
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
                                   fileUploadUI("gene_labels", label = NULL))
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
                                 numericInput("enrich_pval", label = "Enter enrichment pvalue cutoff", value = 0.05)
                          ),
                          column(width=3,
                                 selectInput("organism", "Select Organism:", choices = names(organism_map), selected = "Human")
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
                     fileUploadUI("heatmap_labels", label = NULL),
                     card(card_header("Heatmap"),
                          plotOutput("h_map"),
                          downloadButton("download_hmap", "Download Plot as PDF")
                     )
                   )
  )
)

