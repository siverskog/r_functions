########################################################################################
##### DIAGNOSTICS TESTING ##############################################################
##### Jonathan Siverskog  ##############################################################
##### 2015-08-13          ##############################################################
#####                     ##############################################################
#####                     ##############################################################

########################################################################################
##### CRITICAL VALUE FOR ADF-TEST ######################################################
########################################################################################

.adf.crit <- function(size, sign, type) {
  
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
  
  
  if(size>500) {
    return(adf.array[5,as.character(sign),type])
  } else if(size<25) {
    return(adf.array[1,as.character(sign),type])
  } else {
    return(approx(x = n, y = adf.array[-6,as.character(sign),type], xout = size)$y)
  }
  
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
##### OPTIMAL LAG FOR ADF ##############################################################
########################################################################################

.adf.optim <- function(x, max = 10, type = "ct") {
  
  z <- NULL
  
  for(i in 0:max) {
    
    tmp <- .adf.function(x, i, type)
    z <- cbind(z, .sic(tmp))
    
  }
  
  return(.adf.function(x, which.min(z)-1, type))
  
}

########################################################################################
##### ADD STARS FOR SIGNIFICANCE #######################################################
########################################################################################

.significance <- function(p, stat, only.stars = FALSE) {
  
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
#' @param print.lag Logical indicating whether lag-length should be returned
#' @param only.stars Logical indicating whether test statistic should not be returned
#' @param digits Scalar setting the number of digits to be returned.
#' @return The Dickey-Fuller test statistic.
#' @examples
#' 

adf.test <- function(x, lag = 10, type = "ct", optim = TRUE, print.lag = TRUE, only.stars = FALSE, digits = 2) {
  
  if(optim) {
    lm.fit <- .adf.optim(x, lag, type)
  } else {
    lm.fit <- .adf.function(x, lag, type)
  }
  
  df <- summary(lm.fit)$coef[1,1]/summary(lm.fit)$coef[1,2] # DICKEY-FULLER TEST STATISTIC
  t <- length(resid(lm.fit)) # SAMPLE SIZE
  sign <- c(0.01, 0.05, 0.1) # SIGNIFICANCE LEVELS
  p <- 0.9
  
  for(i in seq_along(sign)) {
    if(abs(df)>abs(.adf.crit(t, sign[i], type))) { # TEST AGAINST CRITICAL VALUES
      p <- sign[i]
      break()
    }
  }
  
  df <- format(round(df, digits), nsmall = digits) # FORMAT THE NUMBER OF DIGITS
  if(print.lag) { df <- paste(df, "(", length(x)-length(resid(lm.fit))-1, ")", sep = "") }
  df <- .significance(p, df, only.stars)
  
  return(df)
  
}

########################################################################################
##### ... ##############################################################################
########################################################################################