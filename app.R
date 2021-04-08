###############################
##
## Project: MetaboGuru
##
## Purpose: Shiny app for PaIRKAT
##
## Author: Charlie Carpenter & Cameron Severn
## Email: charles.carpenter@cuanschutz.edu
##
## Date Created: 2020-08-17
##
## ---------------------------
## Notes:
##
##
## ---------------------------

library(shiny)
library(tidyverse)
library(magrittr)
library(igraph)
library(matrixcalc)
library(MASS)
library(diffusr)
library(Matrix)
library(KEGGREST)
library(DT)
library(readxl)
library(igraph)
library(visNetwork)
library(colourpicker)
library(shinyBS)
library(ggplot2)
library(RColorBrewer)
library(latex2exp)
library(ggrepel)
library(grid)
library(gridExtra)
library(vroom)
library(BiocManager)
library(XVector)



source('helpers.R')

# Define UI for application that draws a histogram
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
                                    tags$li("Clinical Data: Contains clinical outcomes of interest and any meaningful covariates to be adjusted for. Rows should be subjects and columns should be variables"),
                                    tags$li("Metabolite Measurements: Contains measurments of metabolites for all subjects. Rows should match order of subjects in Clinical Data and columns should be metabolite names"),
                                    tags$li("Pathway Metadata: Contains linkage data and pathway information. Rows are metabolites and columns are variables. Should contain a column with KEGG IDs.")
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
                          tags$li("Rows or columns with missing data should be removed or imputed"),
                          tags$li("Remove any special characters"),
                          tags$li("Remove duplicate rows or column names"),
                          tags$li("Variables names should be in the first row of the spreadsheet with no leading empty rows prior to the start of the data"),
                          tags$li("Categorical variables should be converted to one-hot-encodings (aka dummy variables) with 1 or 0 representing inclusion or non-inclusion in a group respectively"),
                        ),
                        h3("Methods"),
                        p("PaIRKAT is a tool for improving testing power on high dimensional data by including graph topography in the kernel machine regression setting. Studies on high dimensional data can struggle 
                          to include the complex relationships between variables. The semi-parametric kernel machine regression model is a powerful tool for capturing these types of relationships. 
                          They provide a framework for testing for relationships between outcomes of interest and high dimensional data such as metabolomic, genomic, or proteomic pathways. 
                          We propose PaIRKAT, a method for including known biological connections between high dimensional variables into the kernel machine by representing them as edges of 
                          'graphs' or 'networks.' It is common for nodes (e.g. metabolites) to be disconnected from all others within the graph, which leads to meaningful decreases in testing power 
                          whether or not the graph information is included. We include a graph regularization or 'smoothing' approach for managing this issue."),
                        h3("Citing PaIRKAT"),
                        p("Paper Citation"),
                        h3("Grant Information"),
                        h3("License"),
        
                        
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
                                       
                                       h4("Metabolome Meta Data"),
                                       fileInput("pathway", h6("File Select"), accept = c(".csv")),
                                       
                                       h5("Optional"),
                                       h4("PaIRKAT Session"),
                                       fileInput("allData",
                                                 h6("Select .RData file created
                                                  from previous PaIRKAT session."),
                                                 accept = ".RData")
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
                                     
                                     
                                     helpText('Column in meta data that contains KEGG IDs.'),
                                     selectInput("pathID", "KEGG IDs", choices = NULL),
                                     
                                     helpText('Column in meta data that contains HMDB IDs.'),
                                     selectInput("hmdbID", "HMDB IDs (optional)", choices = NULL),
                                     
                                     helpText('Column in meta data that contains column names of metabolite data.'),
                                     selectInput("pathCol", "Metabolite Names", choices = NULL),
                                     
                                     helpText('Select a minimum pathway size. Pathways with fewer than the minimum size will be filtered from results'),
                                     numericInput("minSize", "Minimum Pathway Size",
                                                  value = 10, min = 0),
                                     
                                     helpText('This step will take some time.'),
                                     actionButton("pathButton", "Gather Pathways"),
                                     selectInput("plotPath", "Pathway to plot",
                                                 choices = NULL),
                                     
                                     #downloadButton("downloadPathway", "Download PaIRKAT_Pathways")
                        ),
                        
                        mainPanel(visNetworkOutput("pathPlot", height = "600px"))
               ),
               tabPanel("PaIRKAT",
                        fluidRow(
                          column(4,
                                 fluidRow(
                                   column(12,
                                          wellPanel(
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
                                            sliderInput("tau", "\u03C4", value = 1, min = 0, max = 10, step = 1),
                                            bsTooltip(id = "tau", title = "Penalty term", placement = "bottom", trigger = "hover",
                                                      options = NULL),
                                            numericInput("alpha", "Filter \u03B1 <",
                                                         value = 0.05, min = 0, max = 1),
                                            bsTooltip(id = "alpha", title = "Remove pathways with p-values above \u03B1", placement = "bottom", trigger = "hover",
                                                      options = NULL),
                                            
                                            actionButton("pKat", "Run PaIRKAT"),
                                            
                                            downloadButton("downloadPaIRKAT", "Download PaIRKAT Results"), 
                                          )
                                          
                                   ),
                                   
                                 )
                                 
                          ),
                          
                          
                          column(8,
                                 h3("Model Information"),
                                 "Outcome: ",textOutput("pKatY"),
                                 "Clinical Covariates: ",textOutput("pKatX"),
                                 h3("Results"),
                                 DT::dataTableOutput("pKatTab"),
                          ),
                          
                          
                        )
               ),
               tabPanel("Visualizations",
                        tabsetPanel(
                          tabPanel("Network Graph",
                                   fluidRow(
                                     column(4,
                                            
                                            fluidRow(
                                              column(12,
                                                     wellPanel(
                                                       div(style = 'overflow-y: scroll; height: 70vh',
                                                           h4("Plot Control"),
                                                           
                                                           bsCollapse(id = "collapseContainer", open = "Data", multiple = T,
                                                                      bsCollapsePanel("Data", 
                                                                                      DT::dataTableOutput("pKatpaths"),
                                                                                      style = "info"),
                                                                      bsCollapsePanel("Size", 
                                                                                      selectInput("node_size", "Node Size",
                                                                                                  c("None" = "none",
                                                                                                    "Number of Edges" = "num_edges",
                                                                                                    "Interpathway Connectivity" = "interpathway_connectivity",
                                                                                                    "Significance" = 'significance'
                                                                                                  )),
                                                                                      
                                                                                      style = "info"),
                                                                      bsCollapsePanel("Color", 
                                                                                      selectInput("colorBy", "Color Nodes By",
                                                                                                  c("Constant" = "constant",
                                                                                                    "Pathway" = "pathway",
                                                                                                    "Effect on Outcome" = "effectsize")
                                                                                      ),
                                                                                      colourInput("nodeColor", "Constant Color", "blue"),
                                                                                      
                                                                                      style = "info"),
                                                                      bsCollapsePanel("Labels", 
                                                                                      checkboxInput('nodeLabels', label = "Node Labels", value = TRUE),
                                                                                      #checkboxInput('graphLegend', label = "Legend", value = TRUE),
                                                                                      textInput("graphTitle", label = "Graph Title", value = "", width = NULL, placeholder = "Graph Title"),
                                                                                      
                                                                                      style = "info"),
                                                                      bsCollapsePanel("Physics", 
                                                                                      checkboxInput('graphPhysics', label = "Physics", value = TRUE),
                                                                                      sliderInput('gravitationalConstant', 
                                                                                                  label = "Gravitational Constant",
                                                                                                  min = -10000,
                                                                                                  max = 0,
                                                                                                  value = -2000),
                                                                                      sliderInput('centralGravity', 
                                                                                                  label = "Central Gravity",
                                                                                                  min = 0,
                                                                                                  max = 1,
                                                                                                  value = 0.3),
                                                                                      sliderInput('springLength', 
                                                                                                  label = "Spring Length",
                                                                                                  min = 1,
                                                                                                  max = 300,
                                                                                                  value = 95),
                                                                                      sliderInput('springConstant', 
                                                                                                  label = "Spring Constant",
                                                                                                  min = 0,
                                                                                                  max = 0.2,
                                                                                                  value = 0.04),
                                                                                      sliderInput('damping', 
                                                                                                  label = "Damping",
                                                                                                  min = 0,
                                                                                                  max = 1,
                                                                                                  value = 0.09),
                                                                                      sliderInput('avoidOverlap', 
                                                                                                  label = "Avoid Overlap",
                                                                                                  min = 0,
                                                                                                  max = 1,
                                                                                                  value = 0),
                                                                                      style = "info")
                                                           )
                                                       ),
                                                           actionButton("updateplot", "Redraw Plot"),
                                                           downloadButton("downloadNetwork", "Export Graph"),
                                                           downloadButton("downloadData", "Download PaIRKAT Session"), 
                                                       
                                                     )
                                              ),
                                              
                                            ),
                                            
                                            
                                     ),
                                     
                                     
                                     column(7, 
                                            visNetworkOutput("pathResultsPlot", height = "90vh"),
                                            
                                     ),
                                     column(1,plotOutput("pathLegend", height = "90vh")
                                     ),
                                     
                                     
                                   )
                          ),
                          tabPanel("Plot Builder",
                                   fluidRow(
                                     column(4,
                                            fluidRow(
                                              column(12,
                                                     wellPanel(
                                                       bsCollapse(id = "collapseContainer2", open = "Data", multiple = T,
                                                                  bsCollapsePanel("Data", 
                                                                                  selectInput("plotData", "Data Set",
                                                                                              c(
                                                                                                "Pathway Results" = "pathways",
                                                                                                "Metabolites" = "metabolites",
                                                                                                "Clinical Data" = "clinical"
                                                                                                
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
                                                                                                "Void" = "theme_void"
                                                                                              )),
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
                                                     )
                                              )
                                            )
                                     ),
                                     column(8, 
                                            plotOutput("plotBuilder", height = "90vh")
                                     ),
                                   )
                          )
                        ),
                        
               )
               
    )
  ,style = "max-height: 100vh; overflow-y: auto;")
}

# Define server logic 
server <- function(input, output, session) {
  ## Increasing max input file size
  options(shiny.maxRequestSize=30*1024^2)
  
  
  reac <- reactiveValues("clinDat" = F,
                         "metabDat" = F,
                         "pathDat" = F,
                         "networks" = F,
                         "pKatRslt" = F,
                         "comb_net" = F)
  
  
  # Data importing
  importedData <- reactive({
    req(input$allData)
    ext <- tools::file_ext(input$allData$datapath)
    validate(need(ext == "RData", "Please upload a valid RData file"))
    readRDS(input$allData$datapath)
  })
  
  observe({
    req(importedData())
    importNames <- names(importedData())
    
    
    if ("clinDat" %in% importNames){
      reac$clinDat <- T
    }
    if ("metabDat" %in% importNames){
      reac$metabDat <- T
    }
    if ("pathDat" %in% importNames){
      reac$pathDat <- T
    }
    if ("networks" %in% importNames){
      reac$networks <- T
    }
    if ("pKatRslt" %in% importNames){
      reac$pKatRslt <- T
    }
    if ("comb_net" %in% importNames){
      reac$comb_net <- T
    }
    
  })
  
  
  ####### Data Objects #########
  
  ## Data Input
  
  values <- reactiveValues(
    pullClinFrom = NULL,
    pullMetabFrom = NULL,
    pullPathFrom = NULL
  )
  
  observe({
    lapply(c('clin','allData'), function(x) {
      observe({
        input[[x]]
        values$pullClinFrom <- x
      })
    })
  })
  
  observe({
    lapply(c('metab','allData'), function(x) {
      observe({
        input[[x]]
        values$pullMetabFrom <- x
      })
    })
  })
  
  observe({
    lapply(c('pathway','allData'), function(x) {
      observe({
        input[[x]]
        values$pullPathFrom <- x
      })
    })
  })
  
  clinDat <- reactive({
    req(values$pullClinFrom)
    if (values$pullClinFrom == 'clin'){
      req(input$clin)
      ext <- tools::file_ext(input$clin$datapath)
      validate(need(ext == "csv", "Please upload a csv file"))
      vroom::vroom(input$clin$datapath,
                   .name_repair = 'minimal',
                   col_types = cols())
    }
    
    else if (values$pullClinFrom == 'allData' & reac$clinDat){
      importedData()$clinDat
      
    }
    
  })
  
  metabDat <- reactive({
    req(values$pullMetabFrom)
    if (values$pullMetabFrom == 'metab'){
      req(input$metab)
      ext <- tools::file_ext(input$metab$datapath)
      validate(need(ext == "csv", "Please upload a csv file"))
      vroom::vroom(input$metab$datapath,
                   .name_repair = 'minimal',
                   col_types = cols())
    }
    
    else if (values$pullMetabFrom == 'allData' & reac$metabDat){
      importedData()$metabDat
      
    }
    
  })
  
  pathDat <- reactive({
    req(values$pullPathFrom)
    if (values$pullPathFrom == 'pathway'){
      req(input$pathway)
      ext <- tools::file_ext(input$pathway$datapath)
      validate(need(ext == "csv", "Please upload a csv file"))
      vroom::vroom(input$pathway$datapath,
                   .name_repair = 'minimal',
                   col_types = cols())
    }
    
    else if (values$pullPathFrom == 'allData' & reac$pathDat){
      importedData()$pathDat
      
    }
    
  })
  
  
  ## Network Objects
  
  networks <- eventReactive({
    input$pathButton
    reac$networks
  },{
    if (input$pathButton > 0) {
      req(input$pathButton)
      kegg_format <- grepl("^[C,D]", unique(na.omit(pathDat()[[input$pathID]])))
      
      validate(
        need(sum(kegg_format)/length(kegg_format) > 0.8, "Invalid KEGG IDs, IDs should start with C or D"),
        need(sum(pathDat()[[input$pathCol]] %in% colnames(metabDat())) > 1, "Column containing metabolite names should match column names from metabolite data")
        
      )
      
      withProgress(message = 'Gathering Pathways', value = 0, {
        setProgress(0, detail = "Connecting to Database")
        .cr <- keggLink("compound", "reaction")
        hsapath <- unique(keggLink("pathway", "hsa"))
        
        .pID <- as.character(input$pathID)
        
        cr <- substr(.cr, 5, nchar(.cr))
        reactions <- names(cr)
        compId <- pathDat()[,.pID]
        compId <- unlist(strsplit(compId[!is.na(compId)], "[,]"))
        
        results <- data.frame(keggPath = hsapath,
                              stringsAsFactors=FALSE)
        np <- length(hsapath); comps <- list()
        
        for(i in 1:np){
          comps[[i]] <- keggGet(hsapath[i])[[1]]
          
          if(i < np){
            setProgress(i/np, message = 'Gathering Pathways',
                        detail = paste(i, "pathways gathered"))
          }
          if(i == np){
            setProgress((i-1)/np, message = 'Gathered Pathways',
                        detail = "Forming Networks")
          }
        }
        
        names(comps) <- hsapath
        comp <- sapply(comps, function(p) names(p$COMPOUND))
        compNames <- sapply(comps, function(p) p$NAME)
        results$inpathway <- sapply(comp, function(co) sum(compId %in% co))
        
        co <- sub(" - Homo sapiens (human)", "",
                  compNames[names(compNames) %in% results$keggPath],
                  fixed = T)
        
        results <- merge(results,
                         data.frame(keggPath = names(co),
                                    pathwayNames = co),
                         by = "keggPath")
        testPaths <- results[results$inpathway >= input$minSize, ]
        
        
        pdat <- list(testPaths = testPaths,
                     comps = comps[names(comps) %in% testPaths$keggPath],
                     # path.dat = pathDat()[!is.na(pathDat()$KEGG) & !duplicated(pathDat()$KEGG), ],
                     compoundReaction = cr)
        #### End of `pathList` function
        
        nn <- getNetworks(pathDat = pathDat(), metab = metabDat(),
                          database = input$pathwayDatabase, pdat = pdat,
                          pathCol = as.character(input$pathCol),
                          pathID = .pID)
        
      })
      nn
    }
    else if (reac$networks){
      importedData()$networks
      
    }
  })
  
  
  ## PairKAt Results
  
  pKatRslt <- eventReactive({
    input$pKat
    reac$pKatRslt
  },{
    
    
    if (input$pKat > 0){
      req(input$pKat)
      
      validate(
        need(is.numeric(clinDat()[[input$Y]]),"Outcome variable should be numeric"),
        need(input$X > 0, "Please select at least 1 clinical covariate. Unadjusted models coming soon."),
        need(!(input$Y %in% input$X), "The outcome should not be included in covariates")
      )
      
      .formula <- suppressWarnings(formula_fun(input$X))
      mm <- model.matrix(.formula, data = clinDat())
      
      npath <- nrow(networks()$pdat$testPaths)
      pKat.rslt <- data.frame(Pathway = character(npath),
                              `Pathway Size` = numeric(npath),
                              `Score Statistic` = numeric(npath),
                              pValue = numeric(npath))
      
      set.seed(input$seed)
      withProgress(message = 'Running PaIRKAT', value = 0, {
        for (i in 1:npath) {
          z <- PaIRKAT(G = networks()$networks[[i]],
                       out.type = input$outType,
                       Y = clinDat()[, input$Y], model = mm,
                       tau = input$tau, metab = metabDat())
          
          pKat.rslt[i,] <- c(networks()$pdat$testPaths$pathwayNames[i],
                             networks()$pdat$testPaths$inpathway[i],
                             z['Q'], z['pVal'], z['ka'], z['nu'])
          
          incProgress(1/npath,
                      message = "Running PaIRKAT",
                      detail = paste("Completed",
                                     networks()$pdat$testPaths$pathwayNames[i]))
          
        }
      })
      pKat.rslt$pValueFDR <- p.adjust(pKat.rslt$pValue, method = "BH")
      pKat.rslt$neg.log10.FDR.pValue <- -log10(pKat.rslt$pValueFDR)
      
      ## Linear model of single metabolites
      sig.path <- pKat.rslt$pValueFDR < input$alpha
      
      sig.net <- list(networks = networks()$networks[sig.path],
                      testPaths = networks()$pdat$testPaths[sig.path,])
      
      metab.lm <- metabMod(sig.net, Y = clinDat()[, input$Y],
                           clinDat = clinDat(), metab =  metabDat(),
                           .formula = .formula, out.type = input$outType)
      
      pKat.rslt$Pathway.Size <- as.numeric(pKat.rslt$Pathway.Size)
      
      pKat.rslt$Score.Statistic <- as.numeric(pKat.rslt$Score.Statistic)

      
      list(pKat.rslt = pKat.rslt[sig.path,], metab.lm = metab.lm, y = input$Y, X = input$X)
    }
    
    else if (reac$pKatRslt){
      importedData()$pKatRslt
    }
    
    
  })
  
  ## Combined Network
  
  comb_net <- eventReactive( {
    input$pKatpaths_rows_selected
    input$updateplot
    input$node_size
    input$nodeLabels
    input$colorBy
  },
  {
    req(pKatRslt())
    if (length(input$pKatpaths_rows_selected > 0)){
      
      nn <- networks()
      nn$networks <- rename_network_vertices(nn$networks)
      
      for (i in 1:length(nn$networks)){
        nn$networks[[i]] <- convert_network(nn$networks[i])
      }
      
      important_networks <- pKatRslt()$pKat.rslt$Pathway[input$pKatpaths_rows_selected]
      
      print(important_networks)
      
      comb_net <- combine_networks(nn$networks[important_networks])
      
      print(comb_net)
      
      value <- NULL
      
      if (input$node_size == "interpathway_connectivity"){
        for (i in 1:length(comb_net$nodes$group)){
          value <- c(value,length(strsplit(comb_net$nodes$group[i], ",")[[1]]))
        }
      } else if (input$node_size == "num_edges"){
        for (i in comb_net$nodes$label){
          value <- c(value,sum(comb_net$edges$from == i) +  sum(comb_net$edges$to == i))
        }
      } else if (input$node_size == "significance"){
        for (i in comb_net$nodes$label){
          value <- c(value,pKatRslt()$metab.lm$neg.log10.FDR.pVal[pKatRslt()$metab.lm$metab == i])
        }
      }
      else{
        value <- 1
      }
      
      
      
      
      colorVals <- NULL
      estimates <- NULL
      ps <- NULL
      
      abs_max <- abs(max(pKatRslt()$metab.lm$Estimate))
      color_scale <- brewer.pal(11,"RdYlGn")
      
      
      for (i in comb_net$nodes$label){
        
        estimate_i <- pKatRslt()$metab.lm$Estimate[pKatRslt()$metab.lm$metab == i]
        estimates <- c(estimates, estimate_i)
        
        p_i <- pKatRslt()$metab.lm$FDR.pVal[pKatRslt()$metab.lm$metab == i]
        ps <- c(ps, p_i)
        
      }
      
      comb_net$nodes$estimate <- estimates
      comb_net$nodes$p_FDR <- ps
      
      if (input$colorBy == "effectsize"){
        for (i in comb_net$nodes$label){
          print(i)
          estimate_i <- pKatRslt()$metab.lm$Estimate[pKatRslt()$metab.lm$metab == i]
          print(estimate_i)
          
          colorVals <- c(colorVals,color_scale[trunc((estimate_i/abs_max+1)*10.9999/2)])
          print(colorVals)
        }
        comb_net$nodes$color <- colorVals
      }
      
      comb_net$nodes$title <- paste("<p><strong>",comb_net$nodes$label,"</strong>",
                                    "<p> Pathway(s):",comb_net$nodes$group,
                                    "<p> \u00DF = ",signif(comb_net$nodes$estimate,4),
                                    "<p> FDR p = ",signif(comb_net$nodes$p_FDR,4)
      )
      
      if (input$nodeLabels == F) {
        comb_net$nodes$label <- ""
      }
      if (input$colorBy != "pathway"){
        comb_net$nodes$group <- NA
      }
      comb_net$nodes$value <- value
      comb_net$edges$id <- paste(comb_net$edges$from, comb_net$edges$to, sep = "_")
      comb_net
    }
  })
  
  ## Format network plot positions for exporting
  nodes_positions <- reactive({
    positions <- input$network_positions
    if(!is.null(positions)){
      nodes_positions <- do.call("rbind", lapply(positions, function(x){ data.frame(x = x$x, y = x$y)}))
      nodes_positions$id <- names(positions)
      nodes_positions
    } else {
      NULL
    }
  })
  
  ## Data store for session file saving/loading
  
  data_list <- reactive({
    data_list <- list()
    try({
      data_list$clinDat <- clinDat()
    })
    try({
      data_list$networks <- networks()
    })
    try({
      data_list$metabDat <- metabDat()
    })
    try({
      data_list$pathDat <- pathDat()
    })
    try({
      data_list$pKatRslt <- pKatRslt()
    })
    try({
      data_list$comb_net <- comb_net()
    })
    # try({
    #     data_list$nodes_positions <- nodes_positions()
    # })
    
    data_list
  })
  
  ############ Observers #############
  
  observeEvent(clinDat(),
               updateSelectInput(session, "X",
                                 choices = names(clinDat()) ))
  
  observeEvent(clinDat(),
               updateSelectInput(session, "Y",
                                 choices = names(clinDat()) ))
  
  observeEvent(pathDat(),
               updateSelectInput(session, "pathID",
                                 choices = names(pathDat()) ))
  
  observeEvent(pathDat(),
               updateSelectInput(session, "hmdbID",
                                 choices = names(pathDat()) ))
  
  observeEvent(pathDat(),
               updateSelectInput(session, "pathCol",
                                 choices = names(pathDat()) ))
  
  observeEvent(networks(),
               updateSelectInput(session, "plotPath",
                                 choices = names(networks()$networks),
                                 selected = names(networks()$networks)[1]))
  
  observeEvent(pKatRslt()$y,
               updateTextInput(session, "plotTitle",
                                 value = pKatRslt()$y))
  
  observeEvent(pKatRslt()$y,
               updateTextInput(session, "graphTitle",
                               value = pKatRslt()$y))
  
  observeEvent({
    clinDat()
    pKatRslt()
    input$plotData
  },{
    if (input$plotData == "clinical"){
      updateSelectInput(session, "plotX",
                        choices = names(clinDat()),
                        selected = names(clinDat())[5])
    }
    else if (input$plotData == "pathways"){
      updateSelectInput(session, "plotX",
                        choices = names(pKatRslt()$pKat.rslt),
                        selected = names(pKatRslt()$pKat.rslt)[2])
    }
    else if (input$plotData == "metabolites"){
      updateSelectInput(session, "plotX",
                        choices = names(pKatRslt()$metab.lm),
                        selected = names(pKatRslt()$metab.lm)[1])
    }
    
  })
  
  observeEvent({
    clinDat()
    pKatRslt()
    input$plotData
  },{
    if (input$plotData == "clinical"){
      updateSelectInput(session, "plotY",
                        choices = names(clinDat()),
                        selected = names(clinDat())[4])
    }
    else if (input$plotData == "pathways"){
      updateSelectInput(session, "plotY",
                        choices = names(pKatRslt()$pKat.rslt),
                        selected = names(pKatRslt()$pKat.rslt)[6])
    }
    else if (input$plotData == "metabolites"){
      updateSelectInput(session, "plotY",
                        choices = names(pKatRslt()$metab.lm),
                        selected = names(pKatRslt()$metab.lm)[7])
    }
    
  })
  
  observeEvent({
    clinDat()
    pKatRslt()
    input$plotData
  },{
    if (input$plotData == "clinical"){
      updateSelectInput(session, "plotColor",
                        choices = c("None",names(clinDat())),
                        selected = "None")
    }
    else if (input$plotData == "pathways"){
      updateSelectInput(session, "plotColor",
                        choices = c("None",names(pKatRslt()$pKat.rslt)),
                        selected = "None")
    }
    else if (input$plotData == "metabolites"){
      updateSelectInput(session, "plotColor",
                        choices = c("None",names(pKatRslt()$metab.lm)),
                        selected = "None")
    }
    
  })
  
  observeEvent({
    clinDat()
    pKatRslt()
    input$plotData
  },{
    if (input$plotData == "clinical"){
      updateSelectInput(session, "plotSize",
                        choices = c("None",names(clinDat())),
                        selected = "None")
    }
    else if (input$plotData == "pathways"){
      updateSelectInput(session, "plotSize",
                        choices = c("None",names(pKatRslt()$pKat.rslt)),
                        selected = "None")
    }
    else if (input$plotData == "metabolites"){
      updateSelectInput(session, "plotSize",
                        choices = c("None",names(pKatRslt()$metab.lm)),
                        selected = "None")
    }
    
  })
  
  
  observe({
    #req(comb_net())
    if (length(input$pKatpaths_rows_selected > 0)){
      print(comb_net()$nodes)
      try(
        visNetworkProxy("pathResultsPlot") %>%
          visSetData(nodes = comb_net()$nodes, edges = comb_net()$edges) %>%
          visLayout(randomSeed = 64209)
      )
      
      
      
    }
    
  })
  
  observe({
    visNetworkProxy("pathResultsPlot") %>%
      visPhysics(enabled = input$graphPhysics,
                 barnesHut = list(gravitationalConstant = input$gravitationalConstant,
                                  centralGravity = input$centralGravity,
                                  springLength = input$springLength,
                                  springConstant = input$springConstant,
                                  damping = input$damping,
                                  avoidOverlap = input$avoidOverlap))
  })
  
  observe({
    visNetworkProxy("pathResultsPlot") %>%
      visNodes(color = input$nodeColor)
  })
  
  observe({
    visNetworkProxy("pathResultsPlot") %>%
      visSetTitle(main = input$graphTitle)
  })
  
  observeEvent(input$updategraph,{
    visNetworkProxy("pathResultsPlot") %>% visSetData(nodes = comb_net()$nodes, edges = comb_net()$edges)
  })
  
  # get position info
  observeEvent(input$store_position, {
    visNetworkProxy("pathResultsPlot") %>% visGetPositions()
  })
  
  
  ########### Outputs ############
  
  output$clinTab <- DT::renderDataTable({
    req(clinDat())
    DT::datatable(clinDat(),options = list(
      pageLength = 100,
      scrollY = "70vh", 
      scrollX = T)) %>% formatSignif(columns = names(clinDat())[sapply(clinDat(), is.decimal)], digits = 4)
  })
  
  output$metabTab <- DT::renderDataTable({
    req(metabDat())
    DT::datatable(metabDat(),options = list(
      pageLength = 100,
      scrollY = "70vh", 
      scrollX = T)) #%>% formatSignif(columns = names(metabDat())[sapply(metabDat(), is.decimal)], digits = 4)
  })
  
  output$pathTab <- DT::renderDataTable({
    req(pathDat())
    DT::datatable(pathDat(),options = list(
      pageLength = 100,
      scrollY = "70vh", 
      scrollX = T)) #%>% formatSignif(columns = names(pathDat())[sapply(pathDat(), is.decimal)], digits = 4)
  })
  
  output$pathPlot <- renderVisNetwork({
    req(networks())
    req(input$plotPath)
    nn <- networks()
    
    nn$networks <- rename_network_vertices(nn$networks)
    g <- toVisNetworkData(nn$networks[[as.character(input$plotPath)]])
    visNetwork(g$nodes, g$edges, height = "70vh") %>%
      visPhysics(stabilization = FALSE) %>%
      visLayout(randomSeed = 64209) %>%
      visLegend()
  })
  
  output$pKatTab <- DT::renderDataTable({
    req(pKatRslt())
    DT::datatable(pKatRslt()$pKat.rslt %>% arrange(pValueFDR), options = list(
      pageLength = 100,
      caption = pKatRslt()$y,
      scrollY = "60vh")) %>% formatSignif(c("pValue","pValueFDR","Score.Statistic","neg.log10.FDR.pValue"), 2)
  })
  
  output$pKatY <- reactive({
    req(pKatRslt()$y)
    pKatRslt()$y
  })
  
  output$pKatX <- reactive({
    req(pKatRslt()$X)
    paste(pKatRslt()$X, sep = ",")
  })
  
  output$pKatpaths <- DT::renderDataTable({
    validate(
      need(req(pKatRslt()), "Run PaIRKAT before visualizing results")
    )
    
    DT::datatable(pKatRslt()$pKat.rslt[,c(1,2,5)]%>% arrange(pValueFDR), options = list(
      pageLength = 50,
      scrollY = "20vh")) %>% formatSignif("pValueFDR", 2)
  })
  
  output$pathResultsPlot <- renderVisNetwork({
    nodes <- data.frame("id"=c(0),"label"=c(0),"group"=c(0),"value"=c(0))
    edges <- data.frame("id"=c(0),"from"=c(0),"to"=c(0))
    
    visNetwork(nodes, edges) %>%
      visLayout(randomSeed = 64209) %>%
      visPhysics(enabled = TRUE,
                 stabilization = TRUE)
  })
  
  output$pathLegend <- renderPlot({
    
    req(comb_net())
    
    mid_rescaler <- function(mid = 0) {
      function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
        scales::rescale_mid(x, to, from, mid)
      }
    }
    
    if (input$colorBy == "effectsize"){
      df <- data.frame(x = c(0,1),y = c(0,1),color = c(-max(abs(pKatRslt()$metab.lm$Estimate)),max(abs(pKatRslt()$metab.lm$Estimate))))
      
      p <- ggplot(data = df, aes(x = x, y = y, color = color)) +
        geom_point() +
        #scale_color_gradient2(low = "#67000D", mid = "#FFFFFF", high = "#00441B") +
        scale_color_distiller(palette = "RdYlGn",
                              rescaler = mid_rescaler(),
                              direction = 1,
                              guide = guide_colorbar(barheight = 20, barwidth = 2)) +
        labs(color = "Estimated Effect")
      
      legend <- cowplot::get_legend(p)
      
      grid.draw(legend)
    }
  }
  
  )
  
  output$plotBuilder <- renderPlot(
    {
      req(clinDat())
      if (input$plotData == "clinical"){
        df <- clinDat()
      }
      else if (input$plotData == "pathways"){
        df <- pKatRslt()$pKat.rslt
      }
      
      else if (input$plotData == "metabolites"){
        df <- pKatRslt()$metab.lm
      }
      
      p <- ggplot(data = df, aes_string(y = input$plotY, x = input$plotX)) +
        ggtitle(label = input$plotTitle)
      
      if (input$plotSize == "None"){
        p <- p + geom_point(size = input$plotPointSize)
      } else {
        p <- p + geom_point()
      }
      
      if (input$plotColor != 'None'){
        
        
        if (class(df[[input$plotColor]]) == "numeric"){
          
          p <- p + aes_string(color=input$plotColor)
          if (input$colorScale != "Auto"){
            p <- p + scale_color_distiller(palette = input$colorScale, direction = 1)
          }
          
        }
        else if (class(df[[input$plotColor]]) == "character" | class(df[[input$plotColor]]) == "factor"){
          
          p <- p + aes_string(color=input$plotColor)
          #p <- p + scale_color_brewer(palette = input$colorScale)
        }
        
      }
      
      if (input$plotSize != 'None'){
        p <- p + aes_string(size=input$plotSize)
      }
      
      
      if (input$graphLabels){
        if (input$labelIfX != "All X"){
          validate(
            need(is.numeric(df[[input$plotX]]),"To use conditional labeling along the X axis, X must be a numeric variable")
          )
          if (input$labelIfX == "greater than"){
            
            labelDat <- df[df[input$plotX] > input$xCutoff,]
          } else if (input$labelIfX == "less than"){
            
            labelDat <- df[df[input$plotX] < input$xCutoff,]
          } else if (input$labelIfX == "|greater than|"){
            
            labelDat <- df[abs(df[input$plotX]) > input$xCutoff,]
          }else if (input$labelIfX == "|less than|"){
            
            labelDat <- df[abs(df[input$plotX]) < input$xCutoff,]
          }
        } else {
          labelDat <- df
        }
        
        if (input$labelIfY != "All Y"){
          validate(
            need(is.numeric(df[[input$plotY]]),"To use conditional labeling along the Y axis, Y must be a numeric variable")
          )
          if (input$labelIfY == "greater than"){
            
            labelDat <- labelDat[labelDat[input$plotY] > input$yCutoff,]
          } else if (input$labelIfY== "less than"){
            
            labelDat <- labelDat[labelDat[input$plotY] < input$yCutoff,]
          } else if (input$labelIfY == "|greater than|"){
            
            labelDat <- labelDat[abs(labelDat[input$plotY]) > input$yCutoff,]
          }else if (input$labelIfY == "|less than|"){
            
            labelDat <- labelDat[abs(labelDat[input$plotY]) < input$yCutoff,]
          }
        }
        
        if (input$plotData == "metabolites"){
          p <- p + geom_text_repel(data = labelDat, aes_string(label = "metab"), size = input$plotLabelSize)
        }
        else if (input$plotData == "pathways"){
          p <- p + geom_text_repel(data = labelDat, aes_string(label = "Pathway"), size = input$plotLabelSize)
        }
        
        if (input$cutoffLines != "None"){
          if (input$labelIfY == "greater than" | input$labelIfY == "less than"){
            p <- p + geom_hline(yintercept = input$yCutoff, linetype = input$cutoffLines)
          } else if (input$labelIfY == "|greater than|" | input$labelIfY == "|less than|"){
            p <- p + geom_hline(yintercept = input$yCutoff, linetype = input$cutoffLines)+
              geom_hline(yintercept = -as.numeric(input$yCutoff), linetype = input$cutoffLines)
          }
          
          if (input$labelIfX == "greater than" | input$labelIfX == "less than"){
            p <- p + geom_vline(xintercept = input$xCutoff, linetype = input$cutoffLines)
          } else if (input$labelIfX == "|greater than|" | input$labelIfX == "|less than|"){
            p <- p + geom_vline(xintercept = input$xCutoff, linetype = input$cutoffLines)+
              geom_vline(xintercept = -as.numeric(input$xCutoff), linetype = input$cutoffLines)
          }
        }
        
      }
      
      
      
      p <- p + get(input$plotTheme)()  +
        theme(text = element_text(size = input$plotFontSize))
      
      p
    }
  )
  
  
  output$downloadData <- downloadHandler(
    filename = "PaIRKAT_Data.RData",
    content = function(file) {
      
      saveRDS(data_list(), file = file)
    }
  )
  
  output$downloadPaIRKAT <- downloadHandler(
    filename = "PaIRKAT_Results.csv",
    content = function(file) {
      write.csv(pKatRslt()$pKat.rslt, file = file)
    }
  )
  
  output$downloadNetwork <- downloadHandler(
    filename = function() {
      paste('network-', Sys.Date(), '.html', sep='')
    },
    content = function(con) {
      nodes_positions <- nodes_positions()
      if(!is.null(nodes_positions)){
        nodes_save <- merge(comb_net()$nodes, nodes_positions, by = "id", all = T)
      } else  {
        nodes_save <- comb_net()$nodes
      }
      
      visNetwork(nodes = nodes_save, edges = comb_net()$edges, height = "800px") %>%
        visOptions(highlightNearest = TRUE) %>% visExport() %>%
        visPhysics(enabled = FALSE) %>% visEdges(smooth = TRUE) %>% visSave(con)
    }
  )
}

# Run the application
#enableBookmarking(store = "server")
shinyApp(ui = ui, server = server)
