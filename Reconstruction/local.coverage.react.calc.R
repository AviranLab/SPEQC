manual.mean.var <- function(minus.stops, minus.pass, plus.stops, 
                            plus.pass, est.flag, norm.flag, log.flag, zero.flag) {
    # compares to local.coverage.react.calc settings: 1, 0, 0, 0. (Numerator of ML beta)
    minus.coverage <- minus.stops + minus.pass
    plus.coverage <- plus.stops + plus.pass
    
    beta <- (plus.stops / plus.coverage) - (minus.stops / minus.coverage)
    gamma <- minus.stops / minus.coverage
    
    # for difference estimate calcuations
    if (est.flag == 1) { 
    beta.mean <- beta - (beta * gamma)
    beta.var <- ((beta + gamma - beta * gamma) * (1 - beta + gamma - beta * gamma) / plus.coverage) +
                    ((gamma * (1 - gamma)) / minus.coverage)
    }
    
    # for ratio estimate calculations
    if (est.flag == 2) {
    beta.mean <- 1 + beta * ((1 / gamma) - 1)
    beta.var <- (((beta + gamma - beta * gamma) * (1 - beta - gamma + beta * gamma)) / (gamma ^ 2 * plus.coverage)) + 
      (((beta + gamma - beta * gamma) ^ 2 * gamma * (1- gamma))) / (gamma ^ 4 * minus.coverage)
    }
    
    #adjust values for normalized reactivities
    if (norm.flag == 1 | norm.flag == 2 |norm.flag == 3 | norm.flag == 4 | norm.flag == 5 | norm.flag == 6){
      un.normalized <- local.coverage.react.calc(minus.stops, minus.pass, 
                                                 plus.stops, plus.pass, est.flag, 
                                                 0, log.flag, zero.flag) # get un-normalized react
      normalizer <- get.normalizer(un.normalized, norm.flag)  
      beta.mean <- beta.mean / normalizer
      beta.var <- beta.var / normalizer ^ 2
    }
    
  return(data.frame( mean = beta.mean, variance = beta.var))
    
}

local.coverage.react.calc <- function(minus.stops, minus.pass, plus.stops, plus.pass,
                                      estimate.flag, norm.flag, log.flag, zero.flag) {
  
  # estimate flags: 1 = difference, 2 = ratio
  # norm.flag: 1 = 2%-8% normalization, 2 = boxplot normalization
  # log.flag: 1 = log counts and pass individually (before summing), 2 = log rates, 3 = log before normalize
    # 4 = log after normalize, 5 = sum coverage before log
  # zero.flag: 0 = keep negative reactivities, 1 = set negative reactivities to zero
  
  minus.coverage <- minus.stops + minus.pass
  plus.coverage <- plus.stops + plus.pass
  # log flag = 1 counts become log pseudocounts
  if (log.flag == 1) {
    minus.stops <- log.maker(minus.stops)
    minus.coverage <- log.maker(minus.coverage)
    plus.stops <- log.maker(plus.stops)
    plus.coverage <- log.maker(plus.coverage)
  }
  
 if (log.flag == 5) {
   new.log.df <- log.coverage.creater(data.frame(minus.counts = minus.stops, plus.counts = plus.stops))
   
   minus.stops <- new.log.df$minus.counts
   minus.pass <- new.log.df$minus.pass
   plus.stops <- new.log.df$plus.counts
   plus.pass <- new.log.df$plus.pass
   
   minus.coverage <- minus.stops + minus.pass
   plus.coverage <- plus.stops + plus.pass
   
 }
 
  #  for SHAPE-Seq
  minus.stoprate <- minus.stops / minus.coverage
  plus.stoprate <- plus.stops / plus.coverage
  

  
  # single-end coverage normalization
  # removed weird Ding averaging, because it just multiples each 
  # factor by the length
  
  # averaging total counts by length does not change roc scores at all
  if (estimate.flag == 4 | estimate.flag == 5) {
  minus.stoprate <- minus.stops / (sum(minus.stops))
  plus.stoprate <- plus.stops / (sum(plus.stops))
  }
  
  if (log.flag == 2) {
    minus.stoprate <- log.maker(minus.stoprate)
    plus.stoprate <- log.maker(plus.stoprate)
  }
  
  if (estimate.flag == 1 | estimate.flag == 4) {react <- plus.stoprate - minus.stoprate}
  if (estimate.flag == 2 | estimate.flag == 5) {react <- plus.stoprate / minus.stoprate - 1}
  if (estimate.flag == 3) {react <- (plus.stoprate - minus.stoprate) / (1 - minus.stoprate)}
  
  #weird bug, sometimes output is a 1 column df. This fixes it
  if (length(dim(react)) > 0) {react <- as.vector(react[, 1])}
  
  if(zero.flag == 1) {react[react < 0] <- 0}
  if (log.flag == 3) {
    react <- log(react)
    react <- log.postprocess(react)
  }
  
  if (norm.flag == 1) {react <- two.eight.normalize(react)}
  if (norm.flag == 2) {react <- boxplot.normalize(react)}
  if (norm.flag == 3) {react <- five.ten.normalize(react)}
  if (norm.flag == 4) {react <- ten.ten.normalize(react)} 
  if (norm.flag == 5) {react <- straight.ten.normalize(react)}
  if (norm.flag == 6) {react <- thirty.ten.normalize(react)}
  
  if (log.flag == 4) {
    react <- log(react)
    react <- log.postprocess(react)
  }
  
  return(react)
  
}




log.postprocess <- function(react) {
  real.react <- react[!is.infinite(react) & !is.na(react)]
  b <- min(real.react)
  new.react <- react - b
  #new.react <- react
#  new.react[is.infinite(react) | is.na(react)] <- NaN
  new.react[is.infinite(react) | is.na(react)] <- 0
  #new.react <- new.react + 1
  
  return(new.react)
}


log.coverage.creater <- function(counts.df) {
  # just sum the stops of every stop after each desired residue, but log first
  counts.df$minus.counts <- log.maker(counts.df$minus.counts)
  counts.df$plus.counts <- log.maker(counts.df$plus.counts)
  
  counts.df$minus.pass <- 0
  counts.df$plus.pass <- 0
  r <- length(counts.df$minus.pass)
  for (i in 2:r) {
    counts.df$minus.pass[i] <- sum(counts.df$minus.counts[1:i - 1])
    counts.df$plus.pass[i] <- sum(counts.df$plus.counts[1: i - 1])
  }
  return(counts.df)
}



log.maker <- function(counts) {
    mod.counts <- counts
    mod.counts[mod.counts == 0] <- 1
    counts.output <- log(mod.counts)
    return(counts.output)
}

