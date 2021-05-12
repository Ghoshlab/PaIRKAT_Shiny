# Helpful functions
`%nin%` <- Negate(`%in%`)

########## Data Set Up #############

# pathVar = "KEGG",
getNetworks <- function(pathDat, metab, database = "KEGG", 
                        pdat, pathCol, pathID){
  
  networks <- lapply(pdat$testPaths$keggPath,
                     getNetwork, .comps = pdat$comps, 
                     .metab = metab, .pathDat = pathDat,
                     compoundReaction = pdat$compoundReaction,
                     .pathCol = pathCol, .pathID = pathID)
  
  ## Naming list
  names(networks) <- pdat$testPaths$pathwayNames[!is.na(pdat$testPaths$pathwayNames)]
  
  ## Removing networks without connections (getNetwork returns -1)
  keepNet <- !sapply(networks, is.double)
  networks <- networks[keepNet]
  pdat$testPaths <- pdat$testPaths[keepNet, ]
  pdat$comps <- pdat$comps[keepNet]
  
  return(list(networks = networks, pdat = pdat))
}

# pathVar,
pathList <- function(pathDat, database = "KEGG", min.size, organism = "hsa"){
  species_lookup <- as.data.frame(keggList("organism")[,2:3])
  species_name <- species_lookup[species_lookup$organism == organism,2]
  
  if(database == "KEGG"){
    .cr <- keggLink("compound", "reaction")
    hsapath <- unique(keggLink("pathway", species))
  } 
  
  cr <- substr(.cr, 5, nchar(.cr))
  reactions <- names(cr)
  
  compId <- as.character(pathDat[, pathID])
  compId <- unlist(strsplit(compId[!is.na(compId)], "[,]"))
  
  results <- data.frame(keggPath = hsapath, 
                        stringsAsFactors=FALSE)
  
  comps <- sapply(hsapath, function(p) keggGet(p)[[1]])
  comp <- sapply(comps, function(p) names(p$COMPOUND))
  compNames <- sapply(comps, function(p) p$NAME)
  results$inpathway <- sapply(comp, function(co) sum(compId %in% co))
  
  testPaths <- results[results$inpathway >= min.size, ]
  co <- sub(paste(" -", species_name), "",
            compNames[names(compNames) %in% testPaths$keggPath],
            fixed = T)
  
  testPaths <- merge(testPaths,
                     data.frame(keggPath = names(co), 
                                pathwayNames = co),
                     by = "keggPath")
  
  return(list(testPaths = testPaths, 
              comps = comps[names(comps) %in% testPaths$keggPath],
              pathDat = pathDat[!is.na(pathDat[, pathID]) & !duplicated(pathDat[, pathID]), ],
              compoundReaction = cr)
  )
}

## Calculates Laplacian of metabolite pathway ## pathVar,
getNetwork <- function(pathId, .comps, .metab, .pathDat, 
                       compoundReaction, .pathCol, .pathID){
  
  
  target_compound <- names(.comps[[pathId]]$COMPOUND)
  
  cvnames <- .pathDat[.pathDat[[.pathID]] %in% target_compound, ]
  varnames <- cvnames[[.pathCol]]
  
  path.v <- varnames[varnames %in% names(.metab)]
  path.v <- path.v[path.v %in% names(.metab)[names(.metab) %in% varnames ]]
  
  cnames <- cvnames[varnames %in% names(.metab), .pathID, drop = TRUE]
  
  ncomp <- length(cnames)
  reactions <- names(compoundReaction)
  
  A <- matrix(-1, length(cnames), length(cnames))
  diag(A) <- 0
  for(i in 1:ncomp){
    r1 <- reactions[which(compoundReaction==cnames[i])]
    j <- 1
    while (j < i){
      r2 <- reactions[which(compoundReaction==cnames[j])]
      common <- intersect(r1, r2)
      
      if(length(common) > 0) {
        A[i, j] <- A[j, i] <- 1
      } else A[i, j] <- A[j, i] <- 0
      
      j <- j + 1
    }
  }
  
  if(sum(A)==0) return(-1)
  G <- igraph::graph_from_adjacency_matrix(A, mode="undirected")
  V(G)$label <- path.v
  return(G)
}

########## Models #############

## Kernel test including network information through laplacian
PaIRKAT <- function(G, out.type, Y, model, tau = 1, metab){
  
  varnames <- V(G)$label
  ZZ <- scale(metab[, varnames[varnames %in% names(metab)]] )
  
  ## normalized Laplacian
  L <- graph.laplacian(G, normalized = T)
  rho <- median(dist(ZZ))
  Z <- ZZ %*% solve(diag(nrow(L)) + tau*L)
  K <- Gaussian_kernel(rho, Z)
  
  if(out.type == "C"){
    pp <- kernelScoreC(K, Y=as.matrix(Y), X=model)
  }
  
  if(out.type == "D"){
    pp <- kernelScoreD(K, Y=as.matrix(Y), X=model)
  }
  
  pp
}

## Calculates Laplacian of metabolite pathway
getLaplacian <- function(path_id, metabo2, tau){
  target_compound <- names(keggGet(path_id)[[1]]$COMPOUND)
  
  cvnames <- metabo2[metabo2$KEGG %in% target_compound, c("KEGG", "BIOCHEMICAL") ]
  varnames <- cvnames$BIOCHEMICAL
  cnames <- cvnames[which(varnames %in% names(analysis)), "KEGG"]
  ncomp <- length(cnames)
  
  A <- matrix(-1, length(cnames), length(cnames))
  diag(A) <- 0
  for(i in 1:ncomp){
    r1 <- reactions[which(compoundReaction==cnames[i])]
    j <- 1
    while (j < i){
      r2 <- reactions[which(compoundReaction==cnames[j])]
      common <- intersect(r1, r2)
      
      if(length(common) > 0) {
        A[i, j] <- A[j, i] <- 1
      } else A[i, j] <- A[j, i] <- 0
      
      j <- j + 1
    }
  }
  
  if(sum(A)==0) return(-1)
  L <- graph_from_adjacency_matrix(A, mode="undirected") %>%
    laplacian_matrix(normalized = T)
  # L <- normalize.laplacian(A)
  
  solve(diag(nrow(A)) + tau*L)
}

## Function for making formula from "..." in functions
formula_fun <- function(covs){
  cc <- character(0)
  if (length(covs) > 1){
    for(i in 2:length(covs)) cc <- paste(cc, paste0("`",covs[i],"`"), sep = "+")
    formula( paste0("~ ", paste0("`",covs[1],"`"), cc) )  ## pasting for final formula
  } else if (length(covs) == 1){
    formula(paste("~",paste0("`",covs,"`")))
  } else{
    formula(paste("~ 1"))
  }
  
}

Gaussian_kernel <- function(rho, Z){
  exp(-(1/rho)*as.matrix(dist(Z, method = "euclidean", upper = T)^2))
}

## Calculates scale param "ka" (kappa) and df "nu"
scaleChi <- function(P, K){ 
  ## Pieces of Itilde
  Itt <- (tr(P%*%K%*%P%*%K))/2
  Its <- tr(P%*%K%*%P)/2
  Iss <- tr(P)/2
  
  Itilde <- Itt - Its %*% solve(Iss) %*% t(Its)
  e <- tr(P%*%K)/2
  
  c(ka = Itilde/(2*e), ## Scale
    nu = (2*e^2)/Itilde) ## Degrees of Freedom
  
  ## Original paper (Liu, Lin, Gosh 2008)
  ## has Iss = tr(P^2)/2 but P is idempotent
}

kernelScoreC <- function(K, Y, X){
  ## Projection matrix
  P <- diag(nrow(X)) - X %*% (solve(t(X)%*%X) %*% t(X))
  R <- P%*%Y
  
  s2 <- sum(R^2)/(nrow(X) - ncol(X)) ## MSE = SS/rdf
  sc <- ( t(R)%*%K%*%R )/(2*s2)
  s_chi <- scaleChi(P, K)
  
  Q <- sc/s_chi['ka'] ## Scaled chisq stat
  pVal <- pchisq(Q, df = s_chi['nu'], lower.tail = FALSE)
  
  c(Q = Q, pVal = pVal, s_chi['ka'], s_chi['nu'])
}

kernelScoreD <- function(K, Y, X){
  ## Projection matrix
  m0 <- glm.fit(x=X, y=Y, family = binomial())$fitted.values
  D <- diag(m0*(1-m0))
  
  P <- t(X) %*% D
  P <- solve(t(X) %*% D %*% X) %*% P
  P <- D - (D %*% X %*% P)
  R <- Y-m0
  
  sc <- ( t(R)%*%K%*%R )/2
  
  s_chi <- scaleChi(P, K)
  Q <- sc/s_chi['ka'] ## Scaled chisq stat
  pVal <- pchisq(Q, df = s_chi['nu'], lower.tail = FALSE)
  
  out <- c(Q = Q, pVal = pVal, s_chi['ka'], s_chi['nu'])
  out
}

## modeling 1 metabolite at a time
metabMod <- function(sig.net, Y, clinDat, metab, .formula, out.type = "C"){
  V.labs <- lapply(sig.net$networks, function(G) V(G)$label)
  varnames <- unique(unlist(V.labs))
  ZZ <- metab[, varnames[varnames %in% names(metab)]]
  
  metab.lm <- data.frame(Estimate = numeric(),
                         `Std. Error` = numeric(),
                         `t value` = numeric(),
                         pVal = numeric())
  
  
  for(i in 1:ncol(ZZ)){
    #dd <- data.frame(Y, clinDat, ZZ[,i])
    dd <- tibble(Y, clinDat, ZZ[,i], .name_repair = "minimal")
    
    names(dd)[length(names(dd))] <- colnames(ZZ)[i]
    fp <- paste(colnames(dd)[1],
                paste(.formula, paste0("`",names(ZZ)[i],"`"), sep = "+")[2], sep = "~")
    
    .form <- as.formula(fp)
    if(out.type == "C"){
      mMod <- lm(.form, data = dd)
    } 
    if(out.type == "D") {
      mMod <- glm(.form, data = dd, family = binomial())
    }
    cf <- coef(summary(mMod))
    metab.lm[i, ] <- cf[nrow(cf), ]
    
  }
  metab.lm$metab <- names(ZZ)
  metab.lm$FDR.pVal <- p.adjust(metab.lm$pVal, method = "BH")
  metab.lm$neg.log10.FDR.pVal <- -log10(metab.lm$FDR.pVal)
  metab.lm$estimate_sign <- ifelse(metab.lm$Estimate > 0, "Positive","Negative")
  metab.lm
}


########## Plotting #############

count_shared_edges <- function(adjTrue1, adjTrue2, directed = T) {
  if(directed) se <- sum((adjTrue1 + adjTrue2)==2) / 2
  else se <- sum((adjTrue1 + adjTrue2)==2)
  se
}

## Trace of a matrix
tr <- function(x) sum(diag(x))

# Rename vertices to molecule names
rename_network_vertices <- function(networks){
  for (i in 1:length(networks)){
    g <- networks[[i]]
    V(g)$name <- V(g)$label
    networks[[i]] <- g
  }
  networks
}

# convert igraph to visNetwork object

convert_network <- function(network){
  n <- ggnetwork(network[[1]])
  n$pathway <- names(network[1])
  n
}

# count number of networks a given molecule is part of
count_nets <- function(network){
  length(strsplit(network$nodes$group[i], ",")[[1]])
}

# check if any nodes are common to at least any 2 groups and then combine group names to create a unique list of nodes

# is.decimal <- function(x){
#   if (is.numeric(x)){
#     mod(sum(x[1:3]),1) != 0
#   }
#   else{
#     FALSE
#   }
# }