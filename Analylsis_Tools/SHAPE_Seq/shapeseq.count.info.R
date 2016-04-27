count.info <- function(count.list) {
  # reads shape2.coveragelist and generates information for each RNA replicate 
  #total counts, coverage rate, ratio of coverage in plus and minus, hit rate
  
  output.df <- data.frame(minus.coverage.rate = numeric(), plus.coverage.rate = numeric(),
                          total.minus.coverage = numeric(), total.plus.coverage = numeric(),
                          rate.ratio = numeric(), total.ratio = numeric(),
                          hit.rate = numeric())
  
  for (i in 1:length(count.list)) {
    rep1.minus <- count.list[[i]]$minus.counts[count.list[[i]]$residue > 0]
    rep1.plus <- count.list[[i]]$plus.counts[count.list[[i]]$residue > 0]
    rep1.minus.pass <- count.list[[i]]$minus.pass[count.list[[i]]$residue > 0]
    rep1.plus.pass <- count.list[[i]]$plus.pass[count.list[[i]]$residue > 0]
    
    rep2.minus <- count.list[[i]]$minus.counts.1[count.list[[i]]$residue > 0]
    rep2.plus <- count.list[[i]]$plus.counts.1[count.list[[i]]$residue > 0]
    rep2.minus.pass <- count.list[[i]]$minus.pass.1[count.list[[i]]$residue > 0]
    rep2.plus.pass <- count.list[[i]]$plus.pass.1[count.list[[i]]$residue > 0]
    
    rep3.minus <- count.list[[i]]$minus.counts.2[count.list[[i]]$residue > 0]
    rep3.plus <- count.list[[i]]$plus.counts.2[count.list[[i]]$residue > 0]
    rep3.minus.pass <- count.list[[i]]$minus.pass.2[count.list[[i]]$residue > 0]
    rep3.plus.pass <- count.list[[i]]$plus.pass.2[count.list[[i]]$residue > 0]
    
    rep1.length <- length(rep1.minus)
    rep1.minus.coverage <- sum(rep1.minus)
    rep1.plus.coverage <- sum(rep1.plus)
    rep1.total.coverage <- sum(rep1.minus.coverage + rep1.plus.coverage)
    
    rep1.hit.rate <- (rep1.plus / (rep1.plus + rep1.plus.pass)) - 
      (rep1.minus / (rep1.minus + rep1.minus.pass))
    rep1.hit.rate[rep1.hit.rate < 0] <- 0
    rep1.hit.rate <- sum(rep1.hit.rate) / rep1.length
    
    rep2.length <- length(rep2.minus)
    rep2.minus.coverage <- sum(rep2.minus)
    rep2.plus.coverage <- sum(rep2.plus)
    rep2.total.coverage <- sum(rep2.minus.coverage + rep2.plus.coverage)
    rep2.hit.rate <- sum((rep2.plus / (rep2.plus + rep2.plus.pass)) - 
                           (rep2.minus / (rep2.minus + rep2.minus.pass))) / rep2.length
    
    rep2.hit.rate <- (rep2.plus / (rep2.plus + rep2.plus.pass)) - 
      (rep2.minus / (rep2.minus + rep2.minus.pass))
    rep2.hit.rate[rep2.hit.rate < 0] <- 0
    rep2.hit.rate <- sum(rep2.hit.rate) / rep2.length
    
    
    rep3.length <- length(rep3.minus)
    rep3.minus.coverage <- sum(rep3.minus)
    rep3.plus.coverage <- sum(rep3.plus)
    rep3.total.coverage <- sum(rep3.minus.coverage + rep3.plus.coverage)
    rep3.hit.rate <- sum((rep3.plus / (rep3.plus + rep3.plus.pass)) - 
                           (rep3.minus / (rep3.minus + rep3.minus.pass))) / rep3.length
    
    rep3.hit.rate <- (rep3.plus / (rep3.plus + rep3.plus.pass)) - 
      (rep3.minus / (rep3.minus + rep3.minus.pass))
    rep3.hit.rate[rep3.hit.rate < 0] <- 0
    rep3.hit.rate <- sum(rep3.hit.rate) / rep3.length
    
    output.df <- rbind(output.df, 
                       data.frame(minus.coverage.rate = c(rep1.minus.coverage / rep1.length,
                                                          rep2.minus.coverage / rep2.length,
                                                          rep3.minus.coverage / rep3.length),
                                  plus.coverage.rate = c(rep1.plus.coverage / rep1.length,
                                                         rep2.plus.coverage / rep2.length,
                                                         rep3.plus.coverage / rep3.length),
                                  total.minus.coverage = c(rep1.minus.coverage,
                                                           rep2.minus.coverage,
                                                           rep3.minus.coverage),
                                  total.plus.coverage = c(rep1.plus.coverage,
                                                          rep2.plus.coverage,
                                                          rep3.plus.coverage),
                                  rate.ratio = c(rep1.plus.coverage / rep1.minus.coverage / rep1.length,
                                                 rep2.plus.coverage / rep2.minus.coverage / rep2.length,
                                                 rep3.plus.coverage / rep3.plus.coverage / rep3.length), 
                                  total.ratio = c(rep1.plus.coverage / rep1.minus.coverage,
                                                  rep2.plus.coverage / rep2.minus.coverage,
                                                  rep3.plus.coverage / rep3.plus.coverage),
                                  hit.rate = c(rep1.hit.rate, rep2.hit.rate, rep3.hit.rate)))
    
  }
  
  return(output.df)
}