parameterValidation <- function(lambda, eta){
  # input validation
  if((lambda >= 1) | (lambda <= -1)){
    stop("error: -1 < lambda < 1")
  }
  
  if(eta <= 2){
    stop("error: 2 < eta < infty")
  }
}
#' @title
#' Returns the pdf of Hansen's (1994) 'skewed t' distribution
#'
#' @description
#' Returns the pdf at x of Hansen's (1994) 'skewed t' distribution
#' 
#' @usage 
#' dthanssen(x, lambda = 0, eta = 3, log.p = FALSE)
#'
#' @param x a vector of quantiles of the distribution
#' @param lambda the value of the skewness parameter 
#' @param eta the value of the degrees of freedom parameter 
#' @param log.p if 'TRUE' return the value in log form
#'
#' @returns a vactor of pdf at each element of x
#'
#' @note 
#' Based on tdis_pdf.m from Andrew Patton
#' 
#' @references 
#' Hansen, B. E. (1994). Autoregressive conditional density estimation. International Economic Review, 705-730.
#' 
#' @seealso [unique()], [stringi::stri_unique()] 
#'
#' @examples
#' x = c(0,-1, 1)
#' pdf = dthanssen(x, lambda = 0, eta = 3, log.p = FALSE)
#' pdf
#' 
#' @export
#' 
dthanssen <- function(x, lambda = 0, eta = 3, log.p = FALSE){
  
  # input validation
  parameterValidation(lambda, eta)
  
  # Calculating hanssen variables
  c = gamma((eta+1)/2)/(((pi*(eta-2))^0.5) * gamma(eta/2))
  a = 4*lambda*c*((eta-2)/(eta-1))
  b = (1+ 3*lambda^2 - a^2)^0.5
  
  # initialize the pdf function
  f = as.numeric(NA)
  
  # determine teh value of lambda for each x provided
  L = ifelse(test = (x < -a/b), yes = -lambda, no = +lambda)
  
  # calculate the value of the pdf for each value of x provided
  f = b*c*(1 + (1/(eta-2)) * ((b*x+a)/(1+L))^2 )^(-(eta+1)/2)
  
  if(log.p){f= log(f)}
  
  return(f)
}


pthanssen <- function(q, lambda = 0, eta = 3, lower.tail = TRUE, log.p = FALSE){
  
  # input validation
  parameterValidation(lambda, eta)
  
  # Calculating hanssen variables
  c = gamma((eta+1)/2)/(((pi*(eta-2))^0.5) * gamma(eta/2))
  a = 4*lambda*c*((eta-2)/(eta-1))
  b = (1+ 3*lambda^2 - a^2)^0.5

  # initialize the pdf function
  cdf = as.numeric(NA)
  
  # determine the value of lambda for each x provided
  L = ifelse(test = (q < -a/b), yes = -lambda, no = +lambda)
  
  # For a closed solution the CDF rely in using the classical t distribution.
  # Determine the correction for value lower part of the integral.
  HC = ifelse(test = (q < -a/b), yes = 0, no = 0.5)
  TC = ifelse(test = (q < -a/b), yes = 0, no = ((1 - lambda) / 2))
  
  # Determine the corresponding quantile in a t distribution
  y <- (b * q + a) / (1 + L) * sqrt(eta / (eta - 2))
  
  # determine the value of the CDF
  cdf <- HC + ( (1 + L) * (pt(y, df = eta) - TC) )

  if(!lower.tail){cdf <- 1-cdf}
  
  if(log.p){cdf <- log(cdf)}
  
  return(cdf)
  
}

pthanssen_v1 <- function(q, lambda = 0, eta = 3, lower.tail = TRUE, log.p = FALSE){
  
  base = 0.0001
  x = seq(from = -100, to = q, by = base)
  
  y = dthanssen(x, lambda = lambda, eta = eta)
  
  LT = sum(y*base  + c(y[1], base::diff(y)) * base /2)
  
  if(lower.tail){
    f = LT
  } else { 
    f = 1-LT
  }
  
  if(log.p){f = log(f)}
  
  return(f)    
}


qthanssen <- function(p, lambda = 0, eta = 3, lower.tail = TRUE){
  
  # input validation
  parameterValidation(lambda, eta)
  
  # ToDo: valiudade probability between [0, 1]
  
  if(!lower.tail){p <- 1-p}
  
  # Calculating hanssen variables
  c = gamma((eta+1)/2)/(((pi*(eta-2))^0.5) * gamma(eta/2))
  a = 4*lambda*c*((eta-2)/(eta-1))
  b = (1+ 3*lambda^2 - a^2)^0.5
  
  # determine the probability is below or above the mean (remember that the
  # distribution is assimetric)
  L = ifelse(test = (p < ((1-lambda)/2)), yes = -lambda, no = +lambda)
  
  # For a closed solution the CDF rely in using the classical t distribution.
  # Determine the correction for value lower part of the integral.
  HC = ifelse(test = (p < ((1-lambda)/2)), yes = 0, no = 0.5)
  TC = ifelse(test = (p < ((1-lambda)/2)), yes = 0, no = ((1 - lambda) / 2))
  
  # initialize the pdf function
  inv = as.numeric(NA)
  
  # Calculate the inverse distribution
  inv = (1+L)/b * sqrt((eta-2)/eta) *qt( HC + (p - TC)/(1+L), df = eta, lower.tail = TRUE)-a/b;
  
  return(inv)
}


qthanssen_v1 <- function(p, lambda = 0, eta = 3, lower.tail = TRUE){
  # dnorm(x, mean = 0, sd = 1, log = FALSE)
  # pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # rnorm(n, mean = 0, sd = 1
  
  base = 0.0001
  x = seq(from = -100, to = 100, by = base)
  
  y = dthanssen(x, lambda = lambda, eta = eta)
  
  LT = cumsum(y*base  + c(y[1], base::diff(y)) * base /2)
  
  if(lower.tail) { 
    id = max(which(LT <= p))
  }
  else {
    id = min(which(LT >= 1-p))
  }
  
  return(x[id])
}

# benchmark_results <- microbenchmark(
#   qtV1 = qthanssen_v1(runif(1)),
#   qtV2 = qthanssen_v2(runif(1)),
#   times = 100L # Number of iterations
# )


rthanssen <- function(n, lambda = 0, eta = 3){
  rnd <- runif(n = n, min = 0, max = 1)
  return(qthanssen(p = rnd, lambda = lambda, eta = eta, lower.tail = TRUE))
}

rthanssen_v1 <- function(n, lambda = 0, eta = 3){
  # dnorm(x, mean = 0, sd = 1, log = FALSE)
  # pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  # rnorm(n, mean = 0, sd = 1
  
  sampleCDF <- runif(n)
  
  base = 0.0001
  x = seq(from = -100, to = 100, by = base)
  
  y = dthanssen(x, lambda = lambda, eta = eta)
  
  LT = cumsum(y*base  + c(y[1], base::diff(y)) * base /2)
  
  return(x[sapply(sampleCDF, function(x){max(which(LT <= x))})])
}

# rthanssen(n=100, lambda = 0, eta = 3)
# rthanssen_v1(n=100, lambda = 0, eta = 3)
  
# benchmark_results <- microbenchmark(
#   rndV1 = rthanssen(n=100, lambda = 0, eta = 3),
#   rndV2 = rthanssen_v1(n=100, lambda = 0, eta = 3),
#   times = 100L # Number of iterations
# )
# benchmark_results

hist(rthanssen(n=1000000, lambda = 0.5, eta = 3), breaks = "FD", probability = TRUE)


nll.thanssen <- function(param, data){
  lambda <- param[1]
  eta <- param[2]
  ll <- sum(dthanssen(x = data, lambda, eta, log.p = TRUE))
  return(-ll)
}

mle.thanssen <- function(data, lambda = 0, eta = 2.1){
  
  param <- c("lambda" = lambda, "eta" = eta)
  # param <- c("lambda" = 1, "eta" = 1.1)
  
  # Restriction coefficients 
  ui = matrix(data = c(1,0,-1,0,0,1), nrow = 3, ncol = 2, byrow = TRUE)
  
  # Restriction matrix
  ci = matrix(data = c(-0.9999999999,-0.9999999999,2.0000000001), nrow = 3, ncol = 1, byrow = TRUE)
  
  # theta = matrix(data = param, nrow = 2, ncol = 1, byrow = TRUE)
  # ui %*% param - ci >= 0.
  
  mOpt <- constrOptim(theta = param,
                      f = nll.thanssen,
                      grad=NULL, ui=ui, ci=ci,
                      control = list(),
                      outer.iterations = 100,
                      outer.eps = 1e-05, 
                      data=data)
  return(mOpt)    
}

