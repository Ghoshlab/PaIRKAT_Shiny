library(shiny)
library(shinyBS)

library(dplyr)
library(igraph)
library(DT)
library(ggnetwork)
library(colourpicker)
library(ggplot2)
library(RColorBrewer)
library(latex2exp)
library(ggrepel)
library(vroom)
library(KEGGREST)

source("../helpers.R")

ui <- function(request) {
  fluidPage(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
    ),
    navbarPage("PaIRKAT",
               tabPanel("About",
                        h1("About PaIRKAT"),
                        p("PaIRKAT is model framework for assessing statistical relationships between networks of 
                          metabolites (pathways) and clinical outcome. PaIRKAT queries the KEGG database to determine interactions between metabolites
                          from which network connectivity is constructed."),
                        h3("How to Use"),
                        tags$ol(
                          tags$li("Upload Data",
                                  p("Three datasets are required to use PaIRKAT"),
                                  tags$ol(
                                    tags$li("Clinical Data: Contains clinical outcomes of interest and any meaningful covariates to be adjusted for. Rows should be subjects and columns should be variables."),
                                    tags$li("Metabolite Measurements: Contains measurments of metabolites for all subjects. Rows should be subjects and columns should be metabolite names. One column should have subject IDs matching subject IDs in clinical data."),
                                    tags$li("Pathway Data: Contains linkage data and pathway information. Rows are metabolites and columns are variables. Should contain a column with KEGG IDs.")
                                  ),
                          ),
                          tags$li("Gather Pathways"),
                          p("The Gather Pathways tab will guide you though defining data linkage logic and querying the KEGG database to collect pathway information and form networks of metabolites."),
                          tags$li("Run PaIRKAT"),
                          p("Define outcome of interest and clinical covaraites to control for. Output from this analysis is a list of significant pathways associated with the outcome of interest."),
                          tags$li("Explore Results"),
                          p("Visual tools to explore the results of PaIRKAT analysis. There are two tools provided with PaIRKAT."),
                          tags$ol(
                            tags$li("Network Graph: This tool allows you to view significant pathways and their connectivity both within and between pathways. Many options are available to color and size nodes
                                    to better emphasize features of the network."),
                            tags$li("Plot Builder: This flexible tool allows you to make plots using any of the information entered into the application or derived from its functions.")
                          ),
                        ),
                        h3("Data Cleaning Guide"),
                        p("PaIRKAT requires clean data to function properly. Uploaded files should be in the .csv file format and should abide by the following criteria:"),
                        tags$ul(
                          tags$li("Rows of clinical data and metabolite measurements with missing data should be removed or imputed."),
                          tags$li("Metabolite measurements should be normalized."),
                          tags$li("Metabolite measurements and clinical data should both have a subject ID column with the same name for merging."),
                          tags$li("Remove any special characters."),
                          tags$li("Remove duplicate rows or column names."),
                          tags$li("Variables names should be in the first row of the spreadsheet with no leading empty rows prior to the start of the data."),
                          tags$li("Categorical variables should be converted to one-hot-encodings (aka dummy variables) with 1 or 0 representing inclusion or non-inclusion in a group respectively."),
                        ),
                        h3("Methods"),
                        p("PaIRKAT is a tool for improving testing power on high dimensional data by including graph topography in the kernel machine regression setting. Studies on high dimensional data can struggle 
                          to include the complex relationships between variables. The semi-parametric kernel machine regression model is a powerful tool for capturing these types of relationships. 
                          They provide a framework for testing for relationships between outcomes of interest and high dimensional data such as metabolomic, genomic, or proteomic pathways. 
                          We propose PaIRKAT, a method for including known biological connections between high dimensional variables into the kernel machine by representing them as edges of 
                          'graphs' or 'networks.' It is common for nodes (e.g. metabolites) to be disconnected from all others within the graph, which leads to meaningful decreases in testing power 
                          whether or not the graph information is included. We include a graph regularization or 'smoothing' approach for managing this issue."),
                        h3("Citing PaIRKAT"),
                        "PaIRKAT: A pathway integrated regression-based kernel association test with applications to metabolomics and COPD phenotypes",
                        br(),
                        "Charlie M. Carpenter, Weiming Zhang, Lucas Gillenwater, Cameron Severn, Tusharkanti Ghosh, Russel Bowler, Katerina Kechris, Debashis Ghosh",
                        br(),
                        "bioRxiv 2021.04.23.440821; doi: https://doi.org/10.1101/2021.04.23.440821",
                        
                        h3("Acknowledgements"),
                        p("The project described was supported by Award Number U01 HL089897 and Award Number
                        U01 HL089856 from the National Heart, Lung, and Blood Institute, and U01 CA235488 from
                        the National Cancer Institute. The content is solely the responsibility of the authors and does not
                        necessarily represent the official views of the National Heart, Lung, and Blood Institute, National
                          Cancer Institute, or the National Institutes of Health"),
                        h3("License"),
                        p("PaIRKAT is released under the GNU General Public License version 3 (GPLv3)"),
               ),
               tabPanel("Data Input",
                        fillPage(
                          sidebarPanel(h4("Clinical Data"),
                                       h6("Clinical variables associated with trial's subjects.
                   Subject IDs should be the first column."),
                                       fileInput("clin", h6("File Select"), accept = c(".csv")),
                                       
                                       h4("Metabolome Data"),
                                       h6("Metabolome data of the trial. Subject IDs should be the
                   first column. Metabolites should be normalized and imputed."),
                                       fileInput("metab", h6("File Select"), accept = c(".csv")),
                                       
                                       h4("Pathway Data"),
                                       h6("Pathway data linking metabolite names to KEGG IDs"),
                                       fileInput("pathway", h6("File Select"), accept = c(".csv")),
                                       h5("Optional"),
                                       h4("PaIRKAT Session"),
                                       fileInput("allData",
                                                 h6("Load past PaIRKAT session"),
                                                 accept = ".RData"),
                                       downloadButton("downloadData", "Save current PaIRKAT session"), 
                          ),
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Clinical Data", DT::dataTableOutput("clinTab")),
                              tabPanel("Metabolome", DT::dataTableOutput("metabTab")),
                              tabPanel("Pathways", DT::dataTableOutput("pathTab"))
                            )
                          )
                        )
               ),
               tabPanel("Pathways",
                        sidebarPanel(h4("Gather Pathways"),
                                     helpText('Column in pathway data that contains KEGG IDs.'),
                                     selectInput("pathID", "KEGG IDs", choices = NULL),
                                     helpText('Organism from which data were collected'),
                                     selectizeInput("organism", "Organism", choices = NULL),
                                     helpText('Column in pathway data that contains column names of metabolite data.'),
                                     selectInput("pathCol", "Metabolite Names", choices = NULL),
                                     helpText('Column in metabolitte measurements and clinical data with subject IDs (should be the same column name in both).'),
                                     selectInput("SID", "Subject IDs", choices = NULL),
                                     helpText('Select a minimum pathway size. Pathways with fewer than the minimum size will be filtered from results'),
                                     numericInput("minSize", "Minimum Pathway Size",
                                                  value = 10, min = 0),
                                     helpText('This step will take some time.'),
                                     actionButton("pathButton", "Gather Pathways")
                        ),
                        
                        mainPanel(selectInput("plotPath", "Pathway to plot",
                                              choices = NULL),
                                  plotOutput("plotPathway", height = "80vh"))
               ),
               tabPanel("PaIRKAT",
                        sidebarPanel(
                          h4("PaIRKAT Parameters"),
                          selectInput("Y", "Outcome", choices = NULL),
                          bsTooltip(id = "Y", title = "Select a variable with which to associate metabolites", placement = "bottom", trigger = "hover",
                                    options = NULL),
                          selectInput("X", "Clinical Covariates", choices = NULL, multiple = T),
                          bsTooltip(id = "X", title = "Select all variable to control for", placement = "bottom", trigger = "hover",
                                    options = NULL),
                          radioButtons("outType", "Outcome Type",
                                       c("Continuous" = "C",
                                         "Dichotomous" = "D")),
                          bsTooltip(id = "outType", title = "Is the outcome of interest continuous or dichotomous?", placement = "bottom", trigger = "hover",
                                    options = NULL),
                          bsCollapsePanel("Advanced Options",
                                          sliderInput("tau", "\u03C4", value = 1, min = 0, max = 10, step = 1)
                          ),
                          bsTooltip(id = "tau", title = "Penalty term", placement = "bottom", trigger = "hover",
                                    options = NULL),
                          numericInput("alpha", "Filter \u03B1 <",
                                       value = 0.05, min = 0, max = 1),
                          bsTooltip(id = "alpha", title = "Remove pathways with p-values above \u03B1", placement = "bottom", trigger = "hover",
                                    options = NULL),
                          actionButton("pKat", "Run PaIRKAT"),
                          downloadButton("downloadPaIRKAT", "Download Pathway Results"), 
                          downloadButton("downloadMetabolite", "Download Metabolite Results"), 
                        ),
                        
                        mainPanel(
                          h3("Model Information"),
                          "Outcome: ",textOutput("pKatY"),
                          "Clinical Covariates: ",textOutput("pKatX"),
                          h3("Results"),
                          tabsetPanel(
                            tabPanel("Pathway", DT::dataTableOutput("pKatTab")),
                            tabPanel("Metabolite", DT::dataTableOutput("pKatlmTab"))
                          ),
                          
                        ),
               ),
               tabPanel("Visualizations",
                        tabsetPanel(
                          tabPanel("Network Graph",
                                   fluidRow(
                                     column(4,
                                            fluidRow(
                                              column(12,
                                                     wellPanel(
                                                       #div(style = 'overflow-y: scroll; height: 70vh',
                                                       h4("Plot Control"),
                                                       bsCollapse(id = "collapseContainer", open = "Data", multiple = T,
                                                                  bsCollapsePanel("Data", 
                                                                                  DT::dataTableOutput("pKatpaths"),
                                                                                  style = "info"),
                                                                  bsCollapsePanel("Size", 
                                                                                  selectInput("nodeSize", "Node Size",
                                                                                              c("Constant" = "constant",
                                                                                                "Degree Centrality" = "degreeCentrality",
                                                                                                "Significance" = 'significance'),
                                                                                              selected = "significance",
                                                                                  ),
                                                                                  sliderInput("nodeSizeMax", "Constant/Max Size", value = 15, min = 1, max = 40, step = 1),
                                                                                  sliderInput("nodeSizeMin", "Min Size", value = 5, min = 1, max = 40, step = 1),
                                                                                  style = "info"),
                                                                  bsCollapsePanel("Color", 
                                                                                  selectInput("colorBy", "Color Nodes By",
                                                                                              c("Constant" = "constant",
                                                                                                "Pathway" = "pathway",
                                                                                                "Effect on Outcome" = "effectsize"),
                                                                                              selected = "pathway",
                                                                                  ),
                                                                                  colourInput("nodeColor", "Constant Color", "blue"),
                                                                                  selectInput("colorScaleGraph", "Color Scale",
                                                                                              c("Spectral","RdYlGn", "RdYlBu", "RdGy", "RdBu", 
                                                                                                "PuOr", "PRGn", "PiYG", 
                                                                                                "BrBG"
                                                                                              )),
                                                                                  checkboxInput("flipScale","Flip Scale"),
                                                                                  sliderInput("nodeAlpha", "Node Transparency", value = 1, min = 0, max = 1, step = 0.1),
                                                                                  style = "info"),
                                                                  bsCollapsePanel("Labels", 
                                                                                  checkboxInput('nodeLabels', label = "Node Labels", value = TRUE),
                                                                                  checkboxInput('graphLegend', label = "Legend", value = TRUE),
                                                                                  textInput("graphTitle", label = "Graph Title", value = "", width = NULL, placeholder = "Graph Title"),
                                                                                  
                                                                                  style = "info")
                                                                  
                                                       ),
                                                       #),
                                                       #actionButton("updateplot", "Redraw Plot"),
                                                       downloadButton("downloadNetworkPlot", "Export Graph Image"),
                                                       
                                                     )
                                              ),
                                            ),
                                     ),
                                     column(8, 
                                            plotOutput("pathResultsPlot", height = "90vh"),    
                                     ),
                                   )
                          ),
                          tabPanel("Plot Builder",
                                   fluidRow(
                                     column(4,
                                            fluidRow(
                                              column(12,
                                                     wellPanel(
                                                       h4("Plot Control"),
                                                       bsCollapse(id = "collapseContainer2", open = "Data", multiple = T,
                                                                  bsCollapsePanel("Data", 
                                                                                  selectInput("plotData", "Data Set",
                                                                                              c(
                                                                                                "Pathway Results" = "pathways",
                                                                                                "Metabolite Results" = "metabolites",
                                                                                                "Clinical/Metabolite Data" = "clinical"
                                                                                              )),
                                                                                  selectInput("plotY", "Y Variable",
                                                                                              choices = NULL),
                                                                                  selectInput("plotX", "X Variable",
                                                                                              choices = NULL),
                                                                                  selectInput("plotColor", "Color by",
                                                                                              choices = NULL),
                                                                                  selectInput("plotSize", "Size by",
                                                                                              choices = NULL),
                                                                                  style = "info"),
                                                                  bsCollapsePanel("Labels/Text", 
                                                                                  textInput("plotTitle", label = "Plot Title", value = "", width = NULL, placeholder = ""),
                                                                                  textInput("labelY", label = "Y-axis Label", value = "", width = NULL, placeholder = ""),
                                                                                  textInput("labelX", label = "X-axis Label", value = "", width = NULL, placeholder = ""),
                                                                                  textInput("labelC", label = "Color Label", value = "", width = NULL, placeholder = ""),
                                                                                  textInput("labelS", label = "Size Label", value = "", width = NULL, placeholder = ""),
                                                                                  numericInput(inputId = "plotFontSize",
                                                                                               label = "Font Size",
                                                                                               value = 18),
                                                                                  checkboxInput('graphLabels', label = "Pathway/Metabolite Labels", value = T),
                                                                                  bsCollapsePanel("Label IF", 
                                                                                                  selectInput("labelIfX", "Label IF X is",
                                                                                                              choices = c("All X",
                                                                                                                          "greater than",
                                                                                                                          "less than",
                                                                                                                          "|greater than|",
                                                                                                                          "|less than|")),
                                                                                                  numericInput(inputId = "xCutoff",
                                                                                                               label = "X cutoff",
                                                                                                               value = 0),
                                                                                                  selectInput("labelIfY", "Label IF Y is",
                                                                                                              choices = c("All Y",
                                                                                                                          "greater than",
                                                                                                                          "less than",
                                                                                                                          "|greater than|",
                                                                                                                          "|less than|")),
                                                                                                  numericInput(inputId = "yCutoff",
                                                                                                               label = "Y cutoff",
                                                                                                               value = 0),
                                                                                                  selectInput("cutoffLines", "Cutoff lines",
                                                                                                              choices = c("None",
                                                                                                                          "Dashed" = "dashed",
                                                                                                                          "Solid" = "solid"),
                                                                                                              selected = "dashed"),
                                                                                                  style = "info"),
                                                                                  numericInput(inputId = "plotLabelSize",
                                                                                               label = "Label Size",
                                                                                               value = 4),
                                                                                  style = "info"),
                                                                  bsCollapsePanel("Aesthetics", 
                                                                                  selectInput("plotTheme", "Plot Theme",
                                                                                              c("Gray" = "theme_gray",
                                                                                                "BW" = "theme_bw",
                                                                                                "Line Draw" = "theme_bw",
                                                                                                "Light" = "theme_light",
                                                                                                "Minimal" = "theme_minimal",
                                                                                                "Classic" = "theme_classic",
                                                                                                "Void" = "theme_void"),
                                                                                              selected = "theme_light"
                                                                                  ),
                                                                                  selectInput("colorScale", "Color Scale",
                                                                                              c("Auto","Blues", "BuGn", "BuPu", "GnBu", 
                                                                                                "Greens", "Greys", "Oranges", 
                                                                                                "OrRd", "PuBu", "PuBuGn", "PuRd", 
                                                                                                "Purples", "RdPu", "Reds", "YlGn", 
                                                                                                "YlGnBu", "YlOrBr", "YlOrRd"
                                                                                              )),
                                                                                  numericInput(inputId = "plotPointSize",
                                                                                               label = "Constant Point Size",
                                                                                               value = 1,
                                                                                               min = 0,
                                                                                               max = 10,
                                                                                               step = 0.1),
                                                                                  style = "info")
                                                       ),
                                                       downloadButton("downloadBuilderPlot", "Export Image"),
                                                     ),
                                                     
                                              )
                                            )
                                     ),
                                     column(8, 
                                            plotOutput("plotBuilderPlot", height = "90vh")
                                     ),
                                   )
                          )
                        ),     
               )     
    )
    ,style = "max-height: 100vh; overflow-y: auto;")
}