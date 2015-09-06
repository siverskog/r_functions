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

packages <- c("np", "boot", "copula", "rugarch", "fGarch")

install(packages)

########################################################################################
##### APPLY AR-GARCH-FILTER TO DATA ####################################################
########################################################################################

garch.filter <- function(x, type = "sGARCH", garch = c(1,1), arma = c(1,0), package = "fGarch") {
  
  if(package=="rugarch") {
    
    spec <- ugarchspec(variance.model = list(model = type, garchOrder = garch), mean.model = list(armaOrder= arma))
    
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

copula <- function(data, x, y, bw, bwtype = "adaptive_nn", grid = seq(-2,2,0.1), as.vector = TRUE) {
  
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

boot.copula <- function(data, x, y, grid, rep = 5, block.size = 20, sim = "fixed", bw, bwtype = "adaptive_nn") {
  
  assign("boot.num", 0, envir = globalenv())
  
  ### FUNCTION TO PASS THROUGH BOOTSTRAP ###
  
  func <- function(data, x, y, grid, bw, bwtype) {
    assign("boot.num", boot.num+1, envir = globalenv())
    print(paste("STRAPPIN DEM BOOTS", boot.num, "OF", rep, sep = " "))
    res <- garch.filter(data)
    result <- copula(data = res, x = x, y = y, bw = bw, grid = grid)
    return(result)
  }
  
  ### BOOTSTRAP ###
  
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
               bwtype = bwtype)
  
  mu <- as.data.frame(matrix(apply(bc$t, 2, mean), ncol = length(bw), nrow = ceiling(length(grid)/2)*2, dimnames = list(NULL, paste(colnames(data)[x], colnames(data)[y], sep = "."))))
  sd <- as.data.frame(matrix(apply(bc$t, 2, sd), ncol = length(bw), nrow = ceiling(length(grid)/2)*2, dimnames = list(NULL, paste(colnames(data)[x], colnames(data)[y], sep = "."))))
  
  return(list(mu = mu, sd = sd, grid = grid))
  
}

########################################################################################
##### WRAP FUNCTION ####################################################################
########################################################################################

copula.function <- function(data, x, y, grid, bwtype = "adaptive_nn", ckertype = "gaussian", nmulti = 5, rep = 5, sim = "fixed", block.size = 20) {
  
  res <- garch.filter(data)
  bw <- compute.bandwidth(res, x, y, ckertype, nmulti, bwtype)
  bc <- boot.copula(data, x, y, grid, rep, block.size, sim, bw)
  return(bc)
  
}

########################################################################################
##### EMPIRICAL COPULA #################################################################
########################################################################################

empirical.copula <- function(data, x, y, grid) {
  
  prob <- list()
  
  for(i in 1:length(x)) {
    
    res <- garch.filter(data)
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
    points(neg.grid, bc$mu[neg,i], type = "l", lwd = 2)
    points(pos.grid, bc$mu[pos,i], type = "l", lwd = 2)
    abline(v=0)
    
    if(!is.null(ec)) {
      points(neg.grid, ec$prob[neg,i], pch = "+", type = "b", lty = 3)
      points(pos.grid, ec$prob[pos,i], pch = "+", type = "b", lty = 3) 
    }
  }
  
  if(print) dev.off()
  
}
