########################################################################################
##### FLEXIBLE COPULA ESTIMATION
##### By Jonathan Siverskog
##### 2015-09-06
#####
##### CODED FOLLOWING ZIMMER (2012) AND HO ET AL. (2015)
#####
########################################################################################

install <- function(packages) {
  for(i in 1:length(packages)) {
    if(!is.element(packages[i], installed.packages()[,1])) install.packages(packages[i])
    library(packages[i], character.only = TRUE)
  }
  warnings()
}

packages <- c("np", "boot", "copula", "rugarch", "fGarch", "foreach", "doParallel")

install(packages)

source("https://raw.githubusercontent.com/siverskog/r_functions/master/diagnostics.R")

########################################################################################
##### SELECT BEST GARCH-SPECIFICATION ##################################################
########################################################################################

garch.spec <- function(data, ar.order = 0:4, garch.model = c("sGARCH", "eGARCH", "gjrGARCH"), error.dist = c("norm", "std", "ged"), q.lag = 10, signif.lvl = 0.05) {
  
  ##### SETUP LIST ALL SPECIFICATIONS #####
  
  X <- list()
  count <- 0
  
  for(n in 1:ncol(data)) {
    for(i in 1:length(garch.model)) {
      for(j in 1:length(ar.order)) {
        for(k in 1:length(error.dist)) {
          count <- count+1
          X[[count]] <- c(n, garch.model[i], ar.order[j], error.dist[k])
        }
      }
    }
  }
  
  X <- do.call(rbind, X)
  
  ##### FUNCTION TO RUN THROUGH FOREACH LOOP #####
  
  garch.test <- function(X, i) {
    
    spec <- ugarchspec(variance.model = list(model = X[i,2], garchOrder = c(1,1)),
                       mean.model = list(armaOrder= c(as.numeric(X[i,3]), 0)),
                       distribution.model = X[i,4])
    
    fit <- ugarchfit(spec, data[,as.numeric(X[i,1])], solver = "hybrid")
    resid <- fit@fit$residuals/fit@fit$sigma
    
    q <- Box.test(resid, type = "Ljung-Box", lag = q.lag, fitdf = as.numeric(X[i,3]))$p.value
    q2 <- Box.test(resid^2, type = "Ljung-Box", lag = q.lag, fitdf = as.numeric(X[i,3]))$p.value
    llh <- fit@fit$LLH
    
    return(c(q, q2, llh))
    
  }
  
  ##### FOREACH LOOP #####
  
  print("SETTING UP CLUSTER...")
  cluster <- makeCluster(detectCores()-1)
  registerDoParallel(cluster)
  
  print("ESTIMATING ALL MODELS. THIS MAY TAKE SOME TIME...")
  test <- foreach(i = 1:nrow(X), .packages = c("rugarch", "stats")) %dopar% {
    
    garch.test(X = X, i = i)
    
  }
  
  stopCluster(cluster)
  test <- do.call(rbind, test)
  
  ##### SELECT BEST MODEL #####
  
  best <- as.data.frame(matrix(NA, ncol = 3, nrow = ncol(data)))
  colnames(best) <- c("Model", "AR", "ErrorDist")
  
  for(i in 1:ncol(data)) {
    
    good <- X[,1]==i & test[,1]>=signif.lvl & test[,2]>=signif.lvl
    
    if(all(!good)) next
    
    best[i,] <- X[test[,3]==max(test[good,3]),-1]
    
  }
  
  return(best)
  
}

########################################################################################
##### APPLY AR-GARCH-FILTER TO DATA ####################################################
########################################################################################

garch.filter <- function(x, type, error.dist, garch, arma, package) {
  
  if(package=="rugarch") {
    
    spec <- ugarchspec(variance.model = list(model = type, garchOrder = garch),
                       mean.model = list(armaOrder= arma),
                       distribution.model = error.dist)
    
    if(is.vector(x)) {
      fit <- ugarchfit(spec, x)
      res <- fit@fit$residuals/fit@fit$sigma
    } else {
      res <- list()
      for(i in 1:ncol(x)) {
        fit <- ugarchfit(spec, x[,i])
        res[[i]] <- residuals(fit)/sigma(fit)
      }
    }
  } else {
    
    form <- as.formula(paste("~", "arma", "(", arma[1], ",", arma[2], ")", "+", "garch", "(", garch[1], ",", garch[2], ")", sep = ""))
    
    if(is.vector(x)) {
      fit <- garchFit(formula = form, data = x, trace = FALSE)
      res <- fit@residuals/fit@sigma.t
    } else {
      res <- list()
      for(i in 1:ncol(x)) {
        fit <- garchFit(formula = form, data = x[,i], trace = FALSE)
        res[[i]] <- fit@residuals/fit@sigma.t
      }
    }
  }
  
  if(!is.vector(x)) {
    
    res <- data.frame(do.call(cbind, res))
    colnames(res) <- colnames(x)
    rownames(res) <- rownames(x)
    
  }
  
  return(res)
  
}

########################################################################################
##### DESCRIPTIVE STATISTICS ###########################################################
########################################################################################

desc.stat <- function(df, dec = 2, dlog = TRUE, obsperyear = 260, only.stars = TRUE) {
  
  result <- as.data.frame(matrix(NA, nrow = ncol(df), ncol = 9))
  colnames(result) <- c("Obs", "Mean", "StdDev", "Skewness", "Kurtosis", "JB", "Q(10)", "$Q^2$(10)", "ARCH(10)")
  
  for(i in 1:ncol(df)) {
    
    x <- as.numeric(na.exclude(df[,i]))
    
    if(length(x)>50) {
      result[i,"Obs"] <- format(round(length(x), 0), nsmall = 0)
      result[i,"Skewness"] <- format(round(.skew(x), dec), nsmall = dec)
      result[i,"Kurtosis"] <- format(round(.kurt(x), dec), nsmall = dec)
      
      result[i,"JB"] <- sign(jb.test(x), digits = dec, only.stars = only.stars)
      result[i,"Q(10)"] <- sign(q.test(x, lag = 10), digits = dec, only.stars = only.stars)
      result[i,"$Q^2$(10)"] <- sign(q.test(x, lag = 10, sq = TRUE), digits = dec, only.stars = only.stars)
      result[i,"ARCH(10)"] <- sign(arch.test(x, lag = 10), digits = dec, only.stars = only.stars)
    } else {
      result[i,] <- rep(NA, ncol(result))
      result[i,"Obs"] <- 0
    }
    
    if(dlog) {
      result[i,"Mean"] <- format(round(mean(x)*obsperyear*100, dec), nsmall = dec)
      result[i,"StdDev"] <- format(round(sd(x)*sqrt(obsperyear)*100, dec), nsmall = dec)
    } else{
      result[i,"Mean"] <- format(round(mean(x), dec), nsmall = dec)
      result[i,"StdDev"] <- format(round(sd(x), dec), nsmall = dec)
    }
    
  }
  
  
  
  if(dlog) {rownames(df)[2:3] <- c("Mean (%)", "StdDev (%)")}
  
  rownames(result) <- colnames(df)
  
  return(result)
  
}

########################################################################################
##### COMPUTE BANDWIDTH FOR CDF ########################################################
########################################################################################

compute.bandwidth <- function(data, x, y, ckertype = "gaussian", nmulti = 30, bwtype = "adaptive_nn") {
  
  bw <- list()
  
  for(i in 1:length(x)) {
    print(paste("COMPUTING CDF BANDWIDTH FOR PAIR", i, sep = " "))
    form <- as.formula(paste("~",colnames(data)[x[i]], "+", colnames(data)[y[i]], sep = " "))
    bw[[i]] <- npudistbw(form, data = data, ckertype = ckertype, bwmethod = "cv.cdf", nmulti = nmulti, bwtype = bwtype)
  }
  
  return(bw)
  
}

########################################################################################
##### COPULA FUNCTION ##################################################################
########################################################################################

nonparametric.copula <- function(data, x, y, bw, bwtype = "adaptive_nn", grid = seq(-2,2,0.1), as.vector = TRUE) {
  
  prob <- list()
  
  for(i in 1:length(bw)) {
    
    ucdf.x <- npudistbw(as.formula(paste("~",colnames(data)[x[i]], sep = " ")), data = data, bws = bw[[i]]$bw[1], bandwidth.compute = FALSE, bwtype = bwtype)
    ucdf.y <- npudistbw(as.formula(paste("~",colnames(data)[y[i]], sep = " ")), data = data, bws = bw[[i]]$bw[2], bandwidth.compute = FALSE, bwtype = bwtype)
    
    z <- data.frame(grid)
    colnames(z) <- colnames(data)[x[i]]
    ux <- fitted(npudist(bws = ucdf.x,newdata = data.frame(z)))
    
    colnames(z) <- colnames(data)[y[i]]
    uy <- fitted(npudist(bws = ucdf.y,newdata = data.frame(z)))
    
    uxy <- cbind(ux, uy)
    copula <- npcopula(bws = bw[[i]], data = data, u = uxy)
    C <- diag(matrix(copula$copula,length(grid),length(grid)))
    prob[[i]] <- c((C/uy)[1:(ceiling(length(grid)/2))],((1-ux-uy+C)/(1-uy))[(ceiling(length(grid)/2)):length(grid)])
    
  }
  
  if(as.vector) {
    return(unlist(prob))
  } else {
    result <- do.call(cbind, prob)
    colnames(result) <- paste(colnames(data)[x], colnames(data)[y], sep = ".")
    return(result)
  }
}

########################################################################################
##### BOOTSTRAP ########################################################################
########################################################################################

boot.copula <- function(data, x, y, grid, rep = 5, block.size = 20, sim = "fixed", bw, bwtype = "adaptive_nn", ...) {
  
  env <- environment()
  #count <- 0
  #pb <- winProgressBar(title = "Initializing", min = 0, max = rep, width = 300)
  
  ### FUNCTION TO PASS THROUGH BOOTSTRAP ###
  
  func <- function(data, x, y, grid, bw, bwtype) {
    #curVal <- get("count", envir = env) + detectCores()-1
    #assign("count", curVal, envir = env)
    #setWinProgressBar(get("pb", envir = env), curVal, title = paste("Bootstrap:", round(curVal/rep*100, 0), "% done"))
    res <- garch.filter(data, ...)
    result <- nonparametric.copula(data = res, x = x, y = y, bw = bw, grid = grid)
    return(result)
  }
  
  ### BOOTSTRAP ###
  
  print("CREATING CLUSTER...")
  cluster <- makeCluster(detectCores()-1)
  clusterExport(cluster, varlist = list("env", "data", "func", "rep", "block.size", "sim", "x", "y", "grid", "bw", "bwtype"), envir = environment())
  clusterExport(cluster, varlist = list("garch.filter", "nonparametric.copula"))
  
  clusterCall(cluster, function() library(np))
  clusterCall(cluster, function() library(rugarch))
  clusterCall(cluster, function() library(fGarch))

  print("STARTING BOOTSTRAP...")

  bc <- tsboot(data,
               func,
               R = rep-1,
               l = block.size,
               sim = sim,
               n.sim = nrow(data),
               x = x,
               y = y,
               grid = grid,
               bw = bw,
               bwtype = bwtype,
               parallel = "snow",
               ncpus = detectCores()-1,
               cl = cluster)
  
  stopCluster(cluster)
  
  mu <- as.data.frame(matrix(apply(bc$t, 2, mean), ncol = length(bw), nrow = ceiling(length(grid)/2)*2, dimnames = list(NULL, paste(colnames(data)[x], colnames(data)[y], sep = "."))))
  sd <- as.data.frame(matrix(apply(bc$t, 2, sd), ncol = length(bw), nrow = ceiling(length(grid)/2)*2, dimnames = list(NULL, paste(colnames(data)[x], colnames(data)[y], sep = "."))))
  
  return(list(mu = mu, sd = sd, grid = grid))
  
}

########################################################################################
##### EMPIRICAL COPULA #################################################################
########################################################################################

empirical.copula <- function(data, x, y, grid) {
  
  prob <- list()
  
  for(i in 1:length(x)) {
    
    #res <- garch.filter(data)
    res <- data
    ucdf.x <- ecdf(res[,x[i]])
    ucdf.y <- ecdf(res[,y[i]])
    ux <- ucdf.x(grid)
    uy <- ucdf.y(grid)
    
    U <- pobs(res[,c(x[i], y[i])])
    u <- as.matrix(expand.grid(ux,uy))
    Cn <- C.n(u = u, U = U)
    C <- diag(matrix(Cn, length(grid), length(grid)))
    
    prob[[i]] <- c((C/uy)[1:(ceiling(length(grid)/2))],((1-ux-uy+C)/(1-uy))[(ceiling(length(grid)/2)):length(grid)])
    
  }
  
  prob <- do.call(cbind, prob)
  colnames(prob) <- paste(colnames(data)[x], colnames(data)[y], sep = ".")
  
  return(list(prob = prob, grid = grid))
  
}

########################################################################################
##### PLOT #############################################################################
########################################################################################

plot.copula <- function(c = NULL, bc, ec = NULL, mfrow, w = 200, h = 200, print = FALSE, file = "copula.png") {
  
  if(print) png(filename = file, width = w, height = h, units = "mm", res = 600)
  par(mfrow = mfrow, mar = c(2,3,2,1))
  
  if(!is.null(c)) bc$mu <- as.data.frame(c)
  
  neg.grid <- bc$grid[1:ceiling(length(bc$grid)/2)]
  pos.grid <- bc$grid[ceiling(length(bc$grid)/2):length(bc$grid)]
  neg <- 1:(nrow(bc$mu)/2)
  pos <- (nrow(bc$mu)/2 + 1):nrow(bc$mu)
  
  for(i in 1:ncol(bc$mu)) {
    
    neg.grid <- bc$grid[1:ceiling(length(bc$grid)/2)]
    pos.grid <- bc$grid[ceiling(length(bc$grid)/2):length(bc$grid)]
    
    neg <- 1:(nrow(bc$mu)/2)
    pos <- (nrow(bc$mu)/2 + 1):nrow(bc$mu)
    
    hi.neg <- bc$mu[neg,i]+1.96*bc$sd[neg,i]
    lo.neg <- bc$mu[neg,i]-1.96*bc$sd[neg,i]
    hi.pos <- bc$mu[pos,i]+1.96*bc$sd[pos,i]
    lo.pos <- bc$mu[pos,i]-1.96*bc$sd[pos,i]
    
    plot(NULL,
         xlim = c(bc$grid[1], bc$grid[length(bc$grid)]),
         ylim = c(min(lo.pos, lo.neg), max(hi.pos, hi.neg)),
         xlab = NA,
         ylab = NA,
         main = colnames(bc$mu)[i],
         las = 1)
    
    polygon(c(neg.grid, rev(neg.grid)),c(lo.neg, rev(hi.neg)),col="lightgray",border=FALSE)
    polygon(c(rev(pos.grid), pos.grid),c(rev(hi.pos), lo.pos),col="lightgray",border=FALSE)
    points(neg.grid, bc$mu[neg,i], type = "l", lwd = 2, col = "red")
    points(pos.grid, bc$mu[pos,i], type = "l", lwd = 2, col = "red")
    abline(v=0)
    
    if(!is.null(ec)) {
      points(neg.grid, ec$prob[neg,i], pch = "+")
      points(pos.grid, ec$prob[pos,i], pch = "+") 
    }
  }
  
  if(print) dev.off()
  
}
