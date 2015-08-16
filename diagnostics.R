########################################################################################
##### DIAGNOSTICS TESTING ##############################################################
##### Jonathan Siverskog  ##############################################################
##### 2015-08-13          ##############################################################
########################################################################################

#source("https://raw.githubusercontent.com/siverskog/r_functions/master/diagnostics.R")

########################################################################################
##### CRITICAL VALUE FOR ADF-TEST ######################################################
########################################################################################

.adf.pval <- function(stat, size, type) {
  
  ### SIMULATED CRITICAL VALUES FOR DF-STATISTIC ###
  
  adf.array <- array(NA, c(6,3,3), dimnames = list("size" = c(25, 50, 100, 250, 500, Inf),
                                                   "sign" = c(0.01, 0.05, 0.10),
                                                   "type" = c("ct", "c", "nc")))
  
  adf.array[,,"ct"] <- as.matrix(data.frame("0.01" = c(-4.38, -4.15, -4.04, -3.99, -3.98, -3.96),
                             "0.05" = c(-3.60, -3.50, -3.45, -3.43, -3.42, -3.41),
                             "0.10" = c(-3.24, -3.18, -3.15, -3.13, -3.13, -3.12)))
  
  adf.array[,,"c"] <- as.matrix(data.frame("0.01" = c(-3.75, -3.58, -3.51, -3.46, -3.44, -3.43),
                            "0.05" = c(-3.00, -2.93, -2.89, -2.88, -2.87, -2.86),
                            "0.10" = c(-2.63, -2.60, -2.58, -2.57, -2.57, -2.57)))
  
  adf.array[,,"nc"] <- as.matrix(data.frame("0.01" = c(-2.66, -2.62, -2.60, -2.58, -2.58, -2.58),
                             "0.05" = c(-1.95, -1.95, -1.95, -1.95, -1.95, -1.95),
                             "0.10" = c(-1.60, -1.61, -1.61, -1.62, -1.62, -1.62)))
  
  ### INTERPOLATE CRITICAL VALUES FOR SAMPLE SIZE ###
  
    if(size>500) {
    tmp <- adf.array[6,,type]
  } else if(size<25) {
    tmp <- adf.array[1,,type]
  } else {
    
    func <- function(x) {
      approx(names(x[-6]), x[-6], xout = size)$y
    }
    
    tmp <- apply(adf.array[,,type], 2, func)
    
  }
  
  ### INTERPOLATE P-VALUE BASED ON DF-STATISTIC ###
  
  if(abs(stat)>=abs(tmp[1])) {
    p <- 0.01
  } else if(abs(stat)<=abs(tmp[3])) {
    p <- 0.5
  } else {
    p <- exp(approx(tmp, log(as.numeric(names(tmp))), xout = stat)$y)
  }
  
  return(p)
  
}

########################################################################################
##### LAG FUNCTION #####################################################################
########################################################################################

.lag <- function(x, lag) {
  
  if(lag>0) {
    return(c(rep(NA, lag), x[-((length(x)-lag+1):length(x))]))
  } else {
    return(x)
  }
  
}

########################################################################################
##### LAG FUNCTION FOR MATRICES ########################################################
########################################################################################

.lag.mat <- function(x, lags = rep(1, ncol(x)), contemp = FALSE) {
  
  x <- as.matrix(x)
  if(is.null(colnames(x))) {colnames(x) <- rep("", ncol(x))}
  if(all(lags==0) & !contemp) {return(NULL)}
  
  z <- matrix(NA, ncol = sum(lags+contemp), nrow = nrow(x), dimnames = list(rownames(x), 1:sum(lags+contemp)))
  count <- 0
  
  for(i in 1:length(lags)) {
    for(j in ((1-contemp):lags[i])) {
      
      count <- count+1
      z[,count] <- .lag(x[,i], j)
      
      if(j>0) {
        colnames(z)[count] <- paste(colnames(x)[i], ".", "l", j, sep = "")
      } else {
        colnames(z)[count] <- colnames(x)[i]
      }
      
    }
  }
  
  return(z)
  
}

########################################################################################
##### ADF FUNCTION #####################################################################
########################################################################################

.adf.function <- function(x, lag, type = "ct") {
  
  dx <- c(NA, base::diff(x)) # FIRST DIFFERENCE OF X
  x.l1 <- .lag(x,1) # LAGGED X
  dx.l <- .lag.mat(dx,lag) # LAGS OF FIRST DIFFERENCE
  tt <- 1:length(x) # TIME TREND
  c <- rep(1, length(x)) # CONSTANT
  
  if(type == "ct") {
    z <- na.exclude(cbind(dx, x.l1, dx.l, c, tt))
  } else if(type == "c") {
    z <- na.exclude(cbind(dx, x.l1, dx.l, c))
  } else {
    z <- na.exclude(cbind(dx, x.l1, dx.l))
  }
  
  return(lm(z[,1] ~ z[,-1] -1))
  
}

########################################################################################
##### SCHWARZ INFORMATION CRITERION ####################################################
########################################################################################

.sic <- function(lm.fit) {
  
  L <- logLik(lm.fit)[1]
  k <- length(lm.fit$coefficients)
  n <- length(lm.fit$fitted.values)
  
  return(((-2*L)/n)+((k*log(n))/n))
  
}

########################################################################################
##### AKAIKE INFORMATION CRITERION #####################################################
########################################################################################

.aic <- function(lm.fit) {
  
  L <- logLik(lm.fit)[1]
  k <- length(lm.fit$coefficients)
  n <- length(lm.fit$fitted.values)
  
  return(((-2*L)/n)+((2*k)/n))
  
}

########################################################################################
##### ADD STARS FOR SIGNIFICANCE #######################################################
########################################################################################

sign <- function(x, digits = 2, only.stars = FALSE) {
  
  if(class(x)!="diag.test") {warning("x must be of class 'diag.test'")}
  stat <- x$stat
  p <- x$p
  
  stat <- format(round(stat, digits), nsmall = digits) # FORMAT THE NUMBER OF DIGITS
  if(x$test=="adf") { stat <- paste(stat, "(", x$lag, ")", sep = "") }
  if(only.stars){ stat <- "" }
  
  if (p<=0.01) {
    result <- paste(stat, "***", sep = "")
  } else if (p<=0.05){
    result <- paste(stat, "**", sep = "")
  } else if (p<=0.1) {
    result <- paste(stat, "*", sep = "")
  } else {
    result <- paste(stat, sep = "")
  }
  
  #names(result) <- x$name
  
  return(result)
  
}

########################################################################################
##### ADF TEST #########################################################################
########################################################################################

#' Augmented Dickey-Fuller test
#' 
#' @param x A vector of a time series.
#' @param lag Scalar setting the lag-length. If \code{optim = TRUE} \code{lag} sets the maximum lag length.
#' @param type Sets constant and trend in the regression. Possible choices are \code{c("ct", "c", "nc")}.
#' @param optim Logical indicating whether the Schwarz Information Criterion should be used to select the optimum lag-length.
#' @return The Dickey-Fuller test statistic.
#' @examples
#' 

adf.test <- function(x, lag = 10, type = "ct", optim = TRUE) {
  
  x <- na.exclude(x)
  
  if(optim) { # LOOPS OVER ALL POSSBLE LAG-LENGTHS
      
      z <- NULL
      
      for(i in 0:lag) {
        tmp <- .adf.function(x, i, type)
        z <- c(z, .sic(tmp)) # SAVES THE SIC
        }
      
      lag <- which.min(z)-1 # SELECTS LAG BASED ON MINIMUM SIC
      lm.fit <- .adf.function(x, lag, type)
      
    } else {
    lm.fit <- .adf.function(x, lag, type)
    }
  
  df <- summary(lm.fit)$coef[1,1]/summary(lm.fit)$coef[1,2] # DICKEY-FULLER TEST STATISTIC
  t <- length(resid(lm.fit)) # SAMPLE SIZE

  p <- .adf.pval(df, t, type)
  
  test <- list(test = "adf", stat = df, p = p, type = type, lag = lag, name = paste("ADF", "(", type, ")", sep = ""))
  class(test) <- "diag.test"
  
  return(test)
  
}

########################################################################################
##### SKEWNESS #########################################################################
########################################################################################

.skew <- function(x, sample = TRUE) {
  
  if(sample) {
    return(mean((x-mean(x))^3)/sd(x)^3)}
  else {
    return( mean((x-mean(x))^3) / (sum((x-mean(x))^2)/length(x))^(3/2) )
  }
  
}

########################################################################################
##### KURTOSIS #########################################################################
########################################################################################

.kurt <- function(x, sample = TRUE) {
  
  if(sample) {
    return(mean((x-mean(x))^4)/sd(x)^4)
  } else {
    return( mean((x-mean(x))^4) / (sum((x-mean(x))^2)/length(x))^2 )
  }
  
  return(result)
  
}

########################################################################################
##### JARQUE-BERA TEST #################################################################
########################################################################################

jb.test <- function(x, sample = TRUE) {
  
  x <- na.exclude(x)
  s <- .skew(x, sample)
  k <- .kurt(x, sample)
  n <- length(x)
  
  jb <-  n * ( (s^2)/6 + ((k-3)^2)/24 )
  p <- pchisq(jb, 2, lower.tail = F)
  
  test <- list(test = "jb", stat = jb, p = p, name = "JB")
  class(test) <- "diag.test"
  
  return(test)
  
}

########################################################################################
##### LJUNG-BOX Q-TEST #################################################################
########################################################################################

q.test <- function(x, lag = 20, sq = FALSE) {
  
  x <- na.exclude(x)
  
  lm.fit <- lm(x ~ 1)
  x <- lm.fit$residuals
  if(sq){ x <- x^2 }
  
  n <- length(x)
  rho <- as.vector(acf(x, lag.max = lag, plot = FALSE)$acf)[-1]
  q <- n * (n + 2) * sum(rho^2 / (n-seq(1,lag)))
  p <- pchisq(q, lag, lower.tail = FALSE)
  
  if(sq) {
    test <- list(test = "q2", stat = q, p = p, name = paste("Q^2", "(", lag, ")", sep = "")) 
  } else {
    test <- list(test = "q", stat = q, p = p, name = paste("Q", "(", lag, ")", sep = ""))
  }
  
  
  class(test) <- "diag.test"
  
  return(test)
  
}

########################################################################################
##### ARCH-LM ##########################################################################
########################################################################################
# Modified from code by Bernard Pfaff

arch.test <- function (x, lag = 20, demean = FALSE) {
  x <- as.vector(x)
  if(demean) x <- scale(x, center = TRUE, scale = FALSE)
  lags <- lag + 1
  mat <- embed(x^2, lags)
  arch.lm <- summary(lm(mat[, 1] ~ mat[, -1]))
  stat <- arch.lm$r.squared * length(resid(arch.lm))
  df <- lags - 1
  p <- pchisq(stat, df = df, lower.tail = FALSE)
  
  test <- list(test = "arch-lm", stat = stat, p = p, name = paste("ARCH", "(", lag, ")", sep = ""))
  class(test) <- "diag.test"
  
  return(test)
}

