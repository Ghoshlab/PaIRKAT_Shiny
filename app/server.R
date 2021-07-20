

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
  
  # Organism server-side selection
  updateSelectizeInput(session, "organism", choices = keggList("organism")[,3], server = T)
  
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
      metab_names <- pathDat()[[input$pathCol]]
      validate(
        need(sum(kegg_format)/length(kegg_format) > 0.8, "Invalid KEGG IDs, IDs should start with C or D"),
        need(sum(pathDat()[[input$pathCol]] %in% colnames(metabDat())) > 1, "Column containing metabolite names should match column names from metabolite data"),
        need((input$SID %in% names(metabDat())), "The subject ID varible name should be the same in both clinical and metabolite datasets."),
        need(length(metab_names[duplicated(metab_names)]) == 0, paste("Found duplicate metabolite names in Pathway Data. Please remove duplicate entries for: ",
                                                                      metab_names[duplicated(metab_names)]))
      )
      withProgress(message = 'Gathering Pathways', value = 0, {
        setProgress(0, detail = "Connecting to Database")
        species_lookup <- as.data.frame(keggList("organism")[,2:3])
        species_id <- species_lookup[species_lookup$species == input$organism,1]
        .cr <- keggLink("compound", "reaction")
        hsapath <- unique(keggLink("pathway", species_id))
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
        co <- sub(paste(" -", input$organism), "",
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
                          database = input$pathwayDatabase, 
                          pdat = pdat,
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
        need(is.numeric(clinDat()[[input$Y]]),"Outcome variable should be numeric."),
        # need(input$X > 0, "Please select at least 1 clinical covariate. Unadjusted models coming soon."),
        need(!(input$Y %in% input$X), "The outcome should not be included in covariates."),
        need(!(input$SID %in% input$X), "The subject ID should not be included in covariates."),
        need({
          if((input$outType == "D" & length(unique(clinDat()[[input$Y]])) != 2)){
            FALSE
          }
          else{
            TRUE
          }
        }, "Dichotomous outcome values should be 0 and 1")
      )
      ## Ordering both data sets by input subject id
      c.ord <- order(clinDat()[,input$SID])
      m.ord <- order(metabDat()[,input$SID])
      
      .formula <- formula_fun(input$Y, input$X)
      
      npath <- nrow(networks()$pdat$testPaths)
      pKat.rslt <- data.frame(Pathway = character(npath),
                              `Pathway Size` = numeric(npath),
                              `Score Statistic` = numeric(npath),
                              pValue = numeric(npath))
      set.seed(input$seed)
      withProgress(message = 'Running PaIRKAT', value = 0, {
        for (i in 1:npath){
          z <- PaIRKAT(formula.H0 = .formula, data = clinDat()[c.ord, ],
                       G = networks()$networks[[i]], metab = metabDat()[m.ord, ],
                       out.type = input$outType, tau = input$tau)
          
          pKat.rslt[i,] <- c(networks()$pdat$testPaths$pathwayNames[i],
                             networks()$pdat$testPaths$inpathway[i],
                             z$Q.adj, z$p.value)
          
          incProgress(1/npath,
                      message = "Running PaIRKAT",
                      detail = paste("Completed",
                                     networks()$pdat$testPaths$pathwayNames[i]))
        }
      })
      
      pKat.rslt$pValue[pKat.rslt$pValue == 0] <- 1e-9
      
      pKat.rslt$pValueFDR <- p.adjust(pKat.rslt$pValue, method = "BH")
      
      ## Fixing 0s for log transform
      #pKat.rslt$pValueFDR[pKat.rslt$pValueFDR == 0] <- 0.0001
      pKat.rslt$neg.log10.FDR.pValue <- -log10(pKat.rslt$pValueFDR)
      

      
      ## Linear model of single metabolites
      sig.path <- pKat.rslt$pValueFDR < input$alpha
      validate(
        need(sum(sig.path) > 0, "No significant pathways for chosen alpha.")
      )
      sig.net <- list(networks = networks()$networks[sig.path],
                      testPaths = networks()$pdat$testPaths[sig.path,])
      
      metab.lm <- metabMod(sig.net, formula.H0 = .formula,
                           data = clinDat()[c.ord, ], metab =  metabDat()[m.ord, ],
                           out.type = input$outType)
      pKat.rslt$Pathway.Size <- as.numeric(pKat.rslt$Pathway.Size)
      pKat.rslt$Score.Statistic <- as.numeric(pKat.rslt$Score.Statistic)
      metab.lm <- merge(metab.lm, pathDat(), by.x = "metab", by.y = networks()$pathCol)
      pKat.rslt <- pKat.rslt[sig.path,]
      pKat.rslt <- pKat.rslt %>% arrange(desc(neg.log10.FDR.pValue))
      list(pKat.rslt = pKat.rslt, metab.lm = metab.lm, y = input$Y, X = input$X)
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
    input$colorScaleGraph
  },
  {
    req(pKatRslt())
    if (length(input$pKatpaths_rows_selected > 0)){
      nn <- networks()
      nn$networks <- rename_network_vertices(nn$networks)
      
      for (i in 1:length(nn$networks)){
        V(nn$networks[[i]])$pathway <- names(nn$networks[i])
      }
      important_networks <- pKatRslt()$pKat.rslt$Pathway[input$pKatpaths_rows_selected]
      comb_net <- nn$networks[[important_networks[1]]]
      if (length(important_networks) > 1){
        netString<- NULL
        for (i in 1:length(important_networks)){
          netString <- c(netString,paste0("nn$networks[[important_networks[",i,"]]]"))
        }
        comb_net <- eval(parse(text = paste0("igraph::union(",
                                             paste(text = netString, 
                                                   collapse = ","),
                                             ")", 
                                             collapse = "")))
      }
      
      set.seed(1)
      comb_net <- ggnetwork(comb_net)
      comb_net <- merge(comb_net,count(comb_net,name))
      comb_net
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
  observeEvent(pathDat(),
               updateSelectInput(session, "SID",
                                 choices = names(clinDat()) ))
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
                        choices = c("None",unique(c(names(clinDat()),names(metabDat())))),
                        selected = names(clinDat())[5])
    }
    else if (input$plotData == "pathways"){
      updateSelectInput(session, "plotX",
                        choices = c("-Log10(p-value[FDR])" = "neg.log10.FDR.pValue",
                                    "Pathway Size" = "Pathway.Size",
                                    "Score Statistic" = "Score.Statistic",
                                    "p-value" = "pValue",
                                    "p-value[FDR]" = "pValueFDR"
                        ),
                        selected = "Score.Statistic")
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
                        choices = c("None",unique(c(names(clinDat()),names(metabDat())))),
                        selected = names(clinDat())[4])
    }
    else if (input$plotData == "pathways"){
      updateSelectInput(session, "plotY",
                        choices = c("-Log10(p-value[FDR])" = "neg.log10.FDR.pValue",
                                    "Pathway Size" = "Pathway.Size",
                                    "Score Statistic" = "Score.Statistic",
                                    "p-value" = "pValue",
                                    "p-value[FDR]" = "pValueFDR"
                        ),
                        selected = "neg.log10.FDR.pValue")
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
                        choices = c("None",unique(c(names(clinDat()),names(metabDat())))),
                        selected = "None")
    }
    else if (input$plotData == "pathways"){
      updateSelectInput(session, "plotColor",
                        choices = c("None" = "None",
                                    "-Log10(p-value[FDR])" = "neg.log10.FDR.pValue",
                                    "Pathway Size" = "Pathway.Size",
                                    "Score Statistic" = "Score.Statistic",
                                    "p-value" = "pValue",
                                    "p-value[FDR]" = "pValueFDR"
                        ),
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
                        choices = c("None",unique(c(names(clinDat()),names(metabDat())))),
                        selected = "None")
    }
    else if (input$plotData == "pathways"){
      updateSelectInput(session, "plotSize",
                        choices = c("None" = "None",
                                    "-Log10(p-value[FDR])" = "neg.log10.FDR.pValue",
                                    "Pathway Size" = "Pathway.Size",
                                    "Score Statistic" = "Score.Statistic",
                                    "p-value" = "pValue",
                                    "p-value[FDR]" = "pValueFDR"
                        ),
                        selected = "Pathway.Size")
    }
    else if (input$plotData == "metabolites"){
      updateSelectInput(session, "plotSize",
                        choices = c("None",names(pKatRslt()$metab.lm)),
                        selected = "None")
    }
  })
  
  ########### Outputs ############
  output$clinTab <- DT::renderDataTable({
    req(clinDat())
    DT::datatable(clinDat(),options = list(
      pageLength = 100,
      scrollY = "70vh", 
      scrollX = T)) #%>% formatSignif(columns = names(clinDat())[sapply(clinDat(), is.decimal)], digits = 4)
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
  
  output$pKatTab <- DT::renderDataTable({
    req(pKatRslt())
    formatted <- pKatRslt()$pKat.rslt %>% arrange(pValueFDR)
    names(formatted) <- c("Pathway", "Pathway Size", "Score Statistic","p-value","FDR p-value","-log10(FDR p-value)")
    DT::datatable(formatted, options = list(
      pageLength = 100,
      caption = pKatRslt()$y,
      scrollY = "60vh", 
      scrollX = T)) %>% formatSignif(c("Score Statistic","p-value","FDR p-value","-log10(FDR p-value)"), 2)
  })
  output$pKatlmTab <- DT::renderDataTable({
    req(pKatRslt())
    formatted <- pKatRslt()$metab.lm %>% dplyr::select(metab, Estimate, Std..Error,t.value,pVal,FDR.pVal,neg.log10.FDR.pVal) %>% arrange(FDR.pVal)
    names(formatted) <- c("Metabolite", "Effect Estimate", "SE","t-value","p-value","FDR p-value","-log10(FDR p-value)")
    DT::datatable(formatted, options = list(
      pageLength = 100,
      caption = pKatRslt()$y,
      scrollY = "60vh", 
      scrollX = T)) %>% formatSignif(c("Effect Estimate","SE","t-value","p-value","FDR p-value", "-log10(FDR p-value)"), 2)
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
  
  
  
  plotNetwork <- reactive({
    req(comb_net())
    df <- comb_net()
    df <- merge(df,pKatRslt()$metab.lm, by.x = "name", by.y = "metab")
    pathwayCols <- names(df)[grepl("^pathway",names(df))]
    mid_rescaler <- function(mid = 0) {
      function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
        scales::rescale_mid(x, to, from, mid)
      }
    }
    if(input$flipScale){
      direction <- -1
    }
    else{
      direction <- 1
    }
    
    p <- ggplot(df, aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = "black") +
      theme_blank() +
      scale_size_continuous(range = c(input$nodeSizeMin, input$nodeSizeMax)) +
      ggtitle(input$graphTitle)
    
    if (input$nodeSize == "significance"){
      if (input$colorBy == "effectsize"){
        p <- p + geom_nodes(aes_string(size = "neg.log10.FDR.pVal", 
                                       fill = "Estimate"),
                            alpha = input$nodeAlpha, 
                            shape = 21) +
          scale_fill_distiller(palette = input$colorScaleGraph, 
                               rescaler = mid_rescaler(),
                               direction = direction) +
          labs(fill = "Effect on Outcome", 
               size = "-log10(p-value[FDR])")
      }
      else if (input$colorBy == "pathway"){
        for (pathway in pathwayCols){
          p <- p + geom_nodes(data = df[!is.na(df[pathway]),], 
                              aes_string(size = "neg.log10.FDR.pVal", 
                                         fill = pathway),
                              alpha = input$nodeAlpha, 
                              shape = 21) +
            labs(size = "-log10(p-value[FDR])", fill = "Pathway")
        }
      }
      else{
        p <- p + geom_nodes(aes_string(size = "neg.log10.FDR.pVal"), 
                            fill = input$nodeColor,
                            alpha = input$nodeAlpha, 
                            shape = 21) +
          labs(size = "-log10(p-value[FDR])")
      }
    }
    
    else if (input$nodeSize == "degreeCentrality"){
      if (input$colorBy == "effectsize"){
        p <- p + geom_nodes(aes_string(size = "n", 
                                       fill = "Estimate"),
                            alpha = input$nodeAlpha, 
                            shape = 21) +
          scale_fill_distiller(palette = input$colorScaleGraph, 
                               rescaler = mid_rescaler(),
                               direction = direction) +
          labs(fill = "Effect on Outcome", size = "Connections")
      }
      else if (input$colorBy == "pathway"){
        for (pathway in pathwayCols){
          p <- p + geom_nodes(data = df[!is.na(df[pathway]),], 
                              aes_string(size = "n", 
                                         fill = pathway),
                              alpha = input$nodeAlpha) +
            labs(size = "Degree Centrality", fill = "Pathway")
        }
      }
      else{
        p <- p + geom_nodes(aes_string(size = "n"), 
                            fill = input$nodeColor,
                            alpha = input$nodeAlpha, 
                            shape = 21) +
          labs(size = "Connections")
      }
    }
    
    else{
      if (input$colorBy == "effectsize"){
        p <- p + geom_nodes(aes_string(fill = "Estimate"), 
                            size = input$nodeSizeMax,
                            alpha = input$nodeAlpha, 
                            shape = 21) +
          scale_fill_distiller(palette = input$colorScaleGraph,
                               rescaler = mid_rescaler(),
                               direction = direction) +
          labs(fill = "Effect on Outcome")
      }
      else if (input$colorBy == "pathway"){
        for (pathway in pathwayCols){
          p <- p + geom_nodes(data = df[!is.na(df[pathway]),], 
                              aes_string(fill = pathway),
                              size = input$nodeSizeMax,
                              alpha = input$nodeAlpha, 
                              shape = 21) +
            labs(fill = "Pathway")
        }
      }
      else{
        p <- p + geom_nodes(size = input$nodeSizeMax, 
                            fill = input$nodeColor, 
                            alpha = input$nodeAlpha, 
                            shape = 21)
      }
    }
    
    if (input$graphLegend == F){
      p <- p + theme(legend.position = "none")
    }
    
    if (input$nodeLabels){
      p <- p + geom_nodelabel_repel(aes(label = name),
                                    nudge_y = -0.05,
                                    label.size = NA,
                                    segment.size = 0,
                                    force = 2,
                                    color = "black",
                                    segment.colour = "transparent",
                                    fill = "transparent")
    }
    
    p 
    
  })
  
  output$pathResultsPlot <- renderPlot({
    req(plotNetwork())
    print(plotNetwork())
  })
  
  plotBuilder <- reactive(
    {
      req(clinDat())
      if (input$plotData == "clinical"){
        clin <- clinDat()
        metab <- metabDat()
        df <- merge(clin,metab, by = input$SID)
      }
      else if (input$plotData == "pathways"){
        df <- pKatRslt()$pKat.rslt
      }   
      else if (input$plotData == "metabolites"){
        df <- pKatRslt()$metab.lm
      }  
      
      yVar <- paste0("`",input$plotY,"`")
      xVar <- paste0("`",input$plotX,"`")
      
      
      p <- ggplot(data = df, aes_string(y = yVar, x = xVar)) +
        ggtitle(label = input$plotTitle) 
      if (input$plotSize == "None"){
        p <- p + geom_point(size = input$plotPointSize)
      } else {
        p <- p + geom_point()
      }     
      if (input$plotColor != 'None'){  
        colVar <- paste0("`",input$plotColor,"`")
        if (class(df[[input$plotColor]]) == "numeric"){
          
          p <- p + aes_string(color=colVar)
          if (input$colorScale != "Auto"){
            p <- p + scale_color_distiller(palette = input$colorScale, direction = 1)
          }      
        }
        else if (class(df[[input$plotColor]]) == "character" | class(df[[input$plotColor]]) == "factor"){
          p <- p + aes_string(color=colVar)
          #p <- p + scale_color_brewer(palette = input$colorScale)
        } 
      }    
      if (input$plotSize != 'None'){
        sizeVar <- paste0("`",input$plotSize,"`")
        p <- p + aes_string(size=sizeVar)
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
        theme(text = element_text(size = input$plotFontSize))+
        labs(x= input$labelX, y = input$labelY, color= input$labelC, size = input$labelS)
      p
    }
  )
  
  output$plotBuilderPlot <- renderPlot({
    req(plotBuilder())
    print(plotBuilder())
  })
  
  output$plotPathway <- renderPlot({
    req(networks())
    ggplot(data = networks()$networks[[input$plotPath]], aes(x = x, y = y, xend = xend, yend = yend)) +
      geom_edges(color = "black") +
      geom_nodes(size = 15, color = "gray") +
      geom_nodetext_repel(aes(label=label), nudge_y = -0.05, force = 2)+
      theme_blank()
  })
  
  output$downloadData <- downloadHandler(
    filename = "PaIRKAT_Data.RData",
    content = function(file) {
      saveRDS(data_list(), file = file)
    }
  )
  output$downloadPaIRKAT <- downloadHandler(
    filename = "PaIRKAT_Pathway_Results.csv",
    content = function(file) {
      write.csv(pKatRslt()$pKat.rslt, file = file)
    }
  ) 
  output$downloadMetabolite <- downloadHandler(
    filename = "PaIRKAT_Metabolite_Results.csv",
    content = function(file) {
      write.csv(pKatRslt()$metab.lm, file = file)
    }
  ) 
  output$downloadNetworkPlot <- downloadHandler(
    filename = "network.png",
    content = function(file) {
      ggsave(file,plotNetwork(), width = 16, height = 9, dpi = 300, limitsize = F)
    }
  )
  output$downloadBuilderPlot <- downloadHandler(
    filename = "plot.png",
    content = function(file) {
      ggsave(file,plotBuilder(), width = 16, height = 9, dpi = 300, limitsize = F)
    }
  )
}