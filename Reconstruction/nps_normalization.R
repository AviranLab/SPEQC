#zero'ing of reactivity moved to the calculator

straight.ten.normalize <- function(raw.estimates) {
  sorted <- raw.estimates[order(raw.estimates)]
  normalize.range <- c(round(length(sorted) * .9), round(length(sorted) * 1))
  normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  normalized.reactivity <- raw.estimates / normalizer
  # normalized.reactivity[normalized.reactivity < 0] <- 0
  return(normalized.reactivity)
} 

two.eight.normalize <- function(raw.estimates) {
  sorted <- raw.estimates[order(raw.estimates)]
  normalize.range <- c(round(length(sorted) * .9), round(length(sorted) * .98))
  normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  normalized.reactivity <- raw.estimates / normalizer
  # normalized.reactivity[normalized.reactivity < 0] <- 0
  return(normalized.reactivity)
} 

five.ten.normalize <- function(raw.estimates) {
  sorted <- raw.estimates[order(raw.estimates)]
  normalize.range <- c(round(length(sorted) * .85), round(length(sorted) * .95))
  normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  normalized.reactivity <- raw.estimates / normalizer
  # normalized.reactivity[normalized.reactivity < 0] <- 0
  return(normalized.reactivity)
}

ten.ten.normalize <- function(raw.estimates) {
  sorted <- raw.estimates[order(raw.estimates)]
  normalize.range <- c(round(length(sorted) * .80), round(length(sorted) * .90))
  normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  normalized.reactivity <- raw.estimates / normalizer
  # normalized.reactivity[normalized.reactivity < 0] <- 0
  return(normalized.reactivity)
}

thirty.ten.normalize <- function(raw.estimates) {
  sorted <- raw.estimates[order(raw.estimates)]
  normalize.range <- c(round(length(sorted) * .60), round(length(sorted) * .70))
  normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  normalized.reactivity <- raw.estimates / normalizer
  # normalized.reactivity[normalized.reactivity < 0] <- 0
  return(normalized.reactivity)
}

boxplot.normalize <- function(raw.estimates) {
  outlier.limit <- 0
  non.zero.only <- 0
  
  sorted.estimates <- raw.estimates[order(raw.estimates)]
  #sorted.estimates[sorted.estimates < 0] <- 0
  box.data <- boxplot(sorted.estimates, plot = FALSE)
  trimmed <- sorted.estimates[! sorted.estimates %in% box.data$out]
  trimmed <- trimmed[!is.na(trimmed) & !is.infinite(trimmed)]
  if (non.zero.only == 1) {
    box.data.pos <- boxplot(sorted.estimates[sorted.estimates > 0], plot = F)
    trimmed <- sorted.estimates[! sorted.estimates %in% box.data.pos$out]
    trimmed <- trimmed[!is.na(trimmed) & !is.infinite(trimmed)] 
  }
  #optional 10% max threshold for outliers
  outlier.threshold <- .12
  if (outlier.limit == 1) {
    if (length(box.data$out) > (length(raw.estimates) * outlier.threshold)) {
      trimmed <- sorted.estimates[1:round(length(sorted.estimates) * 
                                            (1 - outlier.threshold))]
    }
  }
  normalize.range <- c(round(length(trimmed) * .9), length(trimmed) * 1)
  #print(normalize.range)
  #print(trimmed[normalize.range[1]:normalize.range[2]])
  normalizer <- mean(trimmed[normalize.range[1]:normalize.range[2]], na.rm = TRUE)
  normalized.reactivity <- raw.estimates / normalizer
  # normalized.reactivity[normalized.reactivity < 0] <- 0 
  return(normalized.reactivity)
}


boxplot.winsorize <- function(raw.estimates) {
  sorted.estimates <- raw.estimates[order(raw.estimates)]
  #sorted.estimates[sorted.estimates < 0] <- 0
  box.data <- boxplot(sorted.estimates, plot = FALSE)
  trimmed <- sorted.estimates[! sorted.estimates %in% box.data$out]
  top.limit <- max(trimmed)
  bottom.limit <- min(trimmed)
  print(c(top.limit, bottom.limit))
    
  normalize.range <- c(round(length(trimmed) * .9), length(trimmed) * 1)
  normalizer <- mean(trimmed[normalize.range[1]:normalize.range[2]], na.rm = TRUE)
    

}


get.normalizer <- function(raw.estimates, flag) {
  if (flag == 1) {
    sorted <- raw.estimates[order(raw.estimates)]
    normalize.range <- c(round(length(sorted) * .9), round(length(sorted) * .98))
    normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  if (flag == 2) {
    sorted.estimates <- raw.estimates[order(raw.estimates)]
    #sorted.estimates[sorted.estimates < 0] <- 0
    box.data <- boxplot(sorted.estimates, plot = FALSE)
    trimmed <- sorted.estimates[! sorted.estimates %in% box.data$out]
    
    #optional 10% max threshold for outliers
    outlier.limit <- 0
    outlier.threshold <- .1
    if (outlier.limit == 1) {
      if (length(box.data$out) > (length(raw.estimates) * outlier.threshold)) {
        trimmed <- sorted.estimates[1:round(length(sorted.estimates) * .9)]
      }
    }
    normalize.range <- c(round(length(trimmed) * .9), length(trimmed) * 1)
    #print(normalize.range)
    #print(trimmed[normalize.range[1]:normalize.range[2]])
    normalizer <- mean(trimmed[normalize.range[1]:normalize.range[2]], na.rm = TRUE)
  }
  
  if (flag == 3) {
      sorted <- raw.estimates[order(raw.estimates)]
      normalize.range <- c(round(length(sorted) * .85), round(length(sorted) * .95))
      normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  
  if (flag == 4) {
      sorted <- raw.estimates[order(raw.estimates)]
      normalize.range <- c(round(length(sorted) * .80), round(length(sorted) * .90))
      normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  
  if (flag == 5) {
      sorted <- raw.estimates[order(raw.estimates)]
      normalize.range <- c(round(length(sorted) * .90), round(length(sorted) * 1))
      normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  
  if (flag == 6) {
      sorted <- raw.estimates[order(raw.estimates)]
      normalize.range <- c(round(length(sorted) * .60), round(length(sorted) * .7))
      normalizer <- mean(sorted[normalize.range[1]:normalize.range[2]])
  }
  
  return(normalizer)
}