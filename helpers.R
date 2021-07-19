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
  
  return(list(networks = networks, pdat = pdat, pathCol = pathCol))
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
PaIRKAT <- function(formula.H0, data, G, K=NULL, metab, out.type = "C", tau = 1){
  
  varnames <- V(G)$label
  ZZ <- scale(metab[, varnames[varnames %in% names(metab)]] )
  
  ## normalized Laplacian
  L <- graph.laplacian(G, normalized = T)
  rho <- median(dist(ZZ))
  Z <- ZZ %*% solve(diag(nrow(L)) + tau*L)
  if(is.null(K)) K <- Gaussian_kernel(rho, Z)
  
  if(out.type == "C"){
    pp <- SKAT.c(formula.H0, .data = data, K = K)
  }
  
  if(out.type == "D"){
    pp <- SKAT.b(formula.H0, .data = data, K = K)
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
formula_fun <- function(Y, covs){
  cc <- character(0)
  if (length(covs) > 1){
    ff <- paste("~", paste(paste0("`", covs,"`"), collapse = "+"))
  } else if (length(covs) == 1){
    ff <- paste("~", paste0("`", covs,"`"))
  } else{
    ff <- "~ 1"
  }
  paste(paste0("`",Y,"`"), ff)
}

Gaussian_kernel <- function(rho, Z){
  exp(-(1/rho)*as.matrix(dist(Z, method = "euclidean", upper = T)^2))
}

# Davies Test -------------------------------------------------------------

# The following code is taken from: https://github.com/jchen1981/SSKAT/blob/main/R/SSKAT.R
# Manuscript can be found at https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.21934

#Compute the tail probability of 1-DF chi-square mixtures
KAT.pval <- function(Q.all, lambda, acc=1e-9,lim=1e6){
  pval = rep(0, length(Q.all))
  i1 = which(is.finite(Q.all))
  for(i in i1){
    tmp <- davies(Q.all[i],lambda,acc=acc,lim=lim)
    pval[i] = tmp$Qq
    
    if(tmp$ifault>0) warning(paste("ifault =", tmp$ifault))
    # pval[i] = Sadd.pval(Q.all[i],lambda)
  }
  return(pval)
}

SKAT.c <- function(formula.H0, .data = NULL, K,
                   acc = 0.00001, lim = 10000, tol = 1e-10) {
  
  formula.H0 <- formula(formula.H0)
  print("continuous formula")
  print(formula.H0)
  m0 <- lm(formula.H0, data=.data)
  mX <- model.matrix(formula.H0, data=.data)
  
  res <- resid(m0); df <- nrow(mX)-ncol(mX)
  s2 <- sum(res^2)/df
  
  P0  <- diag(nrow(mX)) - mX %*% (solve(t(mX) %*% mX) %*% t(mX))
  PKP <- P0 %*% K %*% P0
  q <- as.numeric(res %*% K %*% res / s2)
  
  ee <- eigen(PKP - q * P0, symmetric = T) 
  lambda <- ee$values[abs(ee$values) >= tol]
  
  p.value <- KAT.pval(0, lambda=sort(lambda, decreasing=T), acc = acc, lim = lim)
  
  return(list(p.value=p.value, Q.adj=q))
}

SKAT.b <- function(formula.H0, .data = NULL, K,
                   acc = 0.00001, lim = 10000, tol = 1e-10) {
  
  formula.H0 <- formula(formula.H0)
  print("binary formula")
  print(formula.H0)
  X1 <- model.matrix(formula.H0, .data)
  lhs <- formula.H0[[2]]
  y <- eval(lhs, .data)
  
  y <- factor(y)
  
  
  if (nlevels(y) != 2) {
    stop('The phenotype is not binary!\n')
  } else {
    y <- as.numeric(y) - 1
  }
  
  glmfit <- glm(y ~ X1 - 1, family = binomial)
  
  betas <- glmfit$coef
  mu  <- glmfit$fitted.values
  eta <- glmfit$linear.predictors
  res.wk <- glmfit$residuals
  res <- y - mu
  
  w   <- mu * (1-mu)
  sqrtw <- sqrt(w)
  
  adj <- sum((sqrtw * res.wk)^2) 
  
  DX12 <- sqrtw * X1
  
  qrX <- qr(DX12)
  Q <- qr.Q(qrX)
  Q <- Q[, 1:qrX$rank, drop=FALSE]
  
  P0 <- diag(nrow(X1)) - Q %*% t(Q)
  
  DKD <- tcrossprod(sqrtw) * K
  tQK <- t(Q) %*% DKD
  QtQK <- Q %*% tQK 
  PKP <- DKD - QtQK - t(QtQK) + Q %*% (tQK %*% Q) %*% t(Q)
  q <- as.numeric(res %*% K %*% res) / adj
  ee <- eigen(PKP - q * P0, symmetric = T, only.values=T)  		
  lambda <- ee$values[abs(ee$values) >= tol]
  
  p.value <- KAT.pval(0, lambda=sort(lambda, decreasing = T), acc = acc, lim = lim) 
  
  return(list(p.value=p.value, Q.adj = q))
}

# Saddle pVal functions ---------------------------------------------------

saddle = function(x,lambda){
  d = max(lambda)
  lambda = lambda/d
  x = x/d
  k0 = function(zeta) -sum(log(1-2*zeta*lambda))/2
  kprime0 = function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
  kpprime0 = function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
  n = length(lambda)
  if (any(lambda < 0)) {
    lmin = max(1/(2 * lambda[lambda < 0])) * 0.99999
  } else if (x>sum(lambda)){
    lmin = -0.01
  } else {
    lmin = -length(lambda)/(2*x)
  }
  lmax = min(1/(2*lambda[lambda>0])) * 0.99999
  hatzeta = uniroot(function(zeta) kprime0(zeta) - x, lower = lmin, upper = lmax, tol = 1e-08)$root
  w = sign(hatzeta)*sqrt(2*(hatzeta * x-k0(hatzeta)))
  v = hatzeta*sqrt(kpprime0(hatzeta))
  if(abs(hatzeta)<1e-4){
    return(NA)
  } else{
    return( pnorm(w+log(v/w)/w, lower.tail=FALSE) )
  }
}

Liu.pval = function(Q, lambda){
  c1 = rep(0,4); for(i in 1:4){ c1[i] = sum(lambda^i) }
  muQ = c1[1];  sigmaQ = sqrt(2 *c1[2])
  s1 = c1[3]/c1[2]^(3/2);  s2 = c1[4]/c1[2]^2
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2));  d = s1 *a^3 - a^2;  l = a^2 - 2*d
  } else {
    l = 1/s2;  a = sqrt(l);  d = 0
  }
  muX = l+d;  sigmaX = sqrt(2)*a
  
  Q.Norm = (Q - muQ)/sigmaQ
  Q.Norm1 = Q.Norm*sigmaX + muX
  pchisq(Q.Norm1, df = l,ncp=d, lower.tail=FALSE)
}

Sadd.pval = function(Q.all,lambda){
  sad = rep(1,length(Q.all))
  if(sum(Q.all>0)>0){
    sad[Q.all>0] = sapply(Q.all[Q.all>0],saddle,lambda=lambda)
  }
  id = which(is.na(sad))
  if(length(id)>0){
    sad[id] = Liu.pval(Q.all[id], lambda)
  }
  return(sad)
}

# Univariate Test ---------------------------------------------------------

sqrt.inv <- function (V2) {
  eig.obj <- eigen(V2, symmetric = TRUE)
  vectors <- eig.obj$vectors
  values <- eig.obj$values
  ind <- values >= 1e-10
  values <- values[ind]
  vectors <- vectors[, ind]
  
  temp <- t(vectors) / (values)
  Vi2 <- vectors  %*% temp
  
  temp <- t(vectors) / sqrt(values)
  Vi <- vectors  %*% temp
  
  return(list(Vi = Vi, Vi2 = Vi2, rank = length(values)))
}

## modeling 1 metabolite at a time
metabMod <- function(sig.net, formula.H0, data, metab, out.type = "C"){
  V.labs <- lapply(sig.net$networks, function(G) V(G)$label)
  varnames <- unique(unlist(V.labs))
  ZZ <- metab[, varnames[varnames %in% names(metab)]]
  
  metab.lm <- data.frame(Estimate = numeric(),
                         `Std. Error` = numeric(),
                         `t value` = numeric(),
                         pVal = numeric())
  
  for(i in 1:ncol(ZZ)){
    dd <- tibble(data, ZZ[,i], .name_repair = "minimal")
    
    names(dd)[length(names(dd))] <- colnames(ZZ)[i]

    
    if ("1" %in% as.character(formula.H0)){
      fp <- paste(as.character(formula.H0)[2],paste0("`", names(ZZ)[i],"`"), sep = "~")

    }
    else {
      fp <- paste(formula.H0, paste0("`", names(ZZ)[i],"`"), sep = "+")
    }
    
    
    
    .form <- as.formula(fp)
    print("univariate formula")
    print(.form)
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

# Old Saterthwaitte Code --------------------------------------------------

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