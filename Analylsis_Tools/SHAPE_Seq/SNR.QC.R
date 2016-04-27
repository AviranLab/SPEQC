load("~/Box Sync/Information_Content_in_Probes/Code/Shape-Seq_data_analysis/shape2.coveragelist.RData")

source('~/Box Sync/Information_Content_in_Probes/Code/Shape-Seq_data_analysis/shape.seq.tools.R')
source('~/Box Sync/Information_Content_in_Probes/Code/shape-seq QC/SNR_syncable//shapeseq.count.info.R')

rna.range <- reactivity.quicksplit(react.1101, 1)
react.1101 <- reactivitylooper.shape2list(shape2.coveragelist, 1, 1, 0, 1)
count.summary.stats <- count.info(shape2.coveragelist)
bootstrap.mean.var <- read.csv('~/Box Sync/Information_Content_in_Probes/Code/shape-seq QC/mean variance/diff.MR_28_zeroed.csv')


snr.generator <- function (react, rna.range, snr.type.flag) {
  rnas <- 8
  output <- numeric()
  total.mean <- numeric()
  total.sd <- numeric()
  
  for (i in 1:rnas) {
    sub.react.df <- replicate.subsetter(react, rna.range, i)
    rep1 <- sub.react.df$rep1.react
    rep2 <- sub.react.df$rep2.react
    rep3 <- sub.react.df$rep3.react
    
    # get mean or median SNR for 3 replicates
    if (snr.type.flag == 1) { 
      snr.3rep <- apply(data.frame(rep1, rep2, rep3), 1, mean) /
        apply(data.frame(rep1, rep2, rep3), 1, sd)
      
      output <- c(output, (mean(snr.3rep, na.rm = T)))
    }
    
    # get pairwise SNR mean or median (12, 13, 23)
    if (snr.type.flag == 2) { 
      snr.12 <- apply(data.frame(rep1, rep2), 1, mean) /
        apply(data.frame(rep1, rep2), 1, sd)
      snr.13 <- apply(data.frame(rep1, rep3), 1, mean) /
        apply(data.frame(rep1, rep3), 1, sd)
      snr.23 <- apply(data.frame(rep2, rep3), 1, mean) /
        apply(data.frame(rep2, rep3), 1, sd)
      
      output <- c(output, (c(mean(snr.12, na.rm = T), mean(snr.13, na.rm = T), 
                             mean(snr.23, na.rm = T))))
      
      
    }

    # get residue SNR for 3 replicates
    if (snr.type.flag == 3) {
      rep.mean <- apply(data.frame(rep1, rep2, rep3), 1, mean)
      rep.sd <- apply(data.frame(rep1, rep2, rep3), 1, sd)
      snr.3rep <- rep.mean / rep.sd
      
      output <- c(output, snr.3rep)
      total.mean <- c(total.mean, rep.mean)
      total.sd <- c(total.sd, rep.sd)
    }
  }
  if (snr.type.flag == 3) {
    plot(total.mean, total.sd, log = 'xy')
    lines(c(1e-5, 1e1), c(1e-5, 1e1), col = 'red')
    hist(total.mean / total.sd, breaks = 100)
    hist(log(total.mean) / log(total.sd), breaks = 100)
    return(data.frame(total.mean, total.sd))
  }
  return(output)
}

pairwise.correlation.test <- function(react, rna.range, this.cor) {
  rnas <- 8
  output <- numeric()

  
  for (i in 1:rnas) {
    sub.react.df <- replicate.subsetter(react, rna.range, i)
    rep1 <- sub.react.df$rep1.react
    rep2 <- sub.react.df$rep2.react
    rep3 <- sub.react.df$rep3.react
    
    cor.12 <- cor(rep1, rep2, method = this.cor)
    cor.13 <- cor(rep1, rep3, method = this.cor)
    cor.23 <- cor(rep2, rep3, method = this.cor)
    
    output <- c(output, cor.12, cor.13, cor.23)
  }
  return(output)
}

count.data.parser <- function(count.data, pairwise.flag, snr.flag) {
  #pairwise flag either processes all 3 replicates as one, or pairwise 12, 23, 13
  #snr.flag allows for snr of the data values to be generated, rather than the absolute value
  
  rnas <- 8
  summary.mean <- data.frame(minus.coverage.rate = numeric(), plus.coverage.rate = numeric(), 
                             total.minus.coverage = numeric(),total.plus.coverage = numeric(), 
                             rate.ratio = numeric(), total.ratio = numeric(), hit.rate = numeric())
  summary.sd <- data.frame(minus.coverage.rate = numeric(), plus.coverage.rate = numeric(), 
                             total.minus.coverage = numeric(),total.plus.coverage = numeric(), 
                             rate.ratio = numeric(), total.ratio = numeric(), hit.rate = numeric())
  
  for (i in 1:rnas) {
    begin <- (i - 1) * 3 + 1
    middle <- (i - 1) * 3 + 2
    finish <- (i - 1) * 3 + 3

    if (pairwise.flag == 0) {
      sub.mean <- apply(count.data[begin:finish, ], 2, mean)
      sub.sd <- apply(count.data[begin:finish, ], 2, sd)
      summary.mean[i, ] <- sub.mean
      summary.sd[i, ] <- sub.sd
    }
    
    if (pairwise.flag == 1) {
      sub.mean.12 <- apply(count.data[begin:middle, ], 2, mean)
      sub.mean.13 <- apply(count.data[c(begin, finish), ], 2, mean)
      sub.mean.23 <- apply(count.data[middle:finish, ], 2, mean)
      
      sub.sd.12 <- apply(count.data[begin:middle, ], 2, sd)
      sub.sd.13 <- apply(count.data[c(begin, finish), ], 2, sd)
      sub.sd.23 <- apply(count.data[middle:finish, ], 2, sd)
      
      summary.mean[begin, ] <- sub.mean.12
      summary.mean[middle, ] <- sub.mean.13
      summary.mean[finish, ] <- sub.mean.23
      
      summary.sd[begin, ] <- sub.sd.12
      summary.sd[middle, ] <- sub.sd.13
      summary.sd[finish, ] <- sub.sd.23
    }
  }
 if (snr.flag == 0) {return(list(summary.mean, summary.sd))}
 if (snr.flag == 1) {return(summary.mean / summary.sd)}
}




replicate.subsetter <- function(react, rna.range, i) {
    rep1 <- (i - 1) * 3 + 1
    rep2 <- (i - 1) * 3 + 2
    rep3 <- (i - 1) * 3 + 3
    
    rep1.range <- rna.range[rep1, ]
    rep2.range <- rna.range[rep2, ]
    rep3.range <- rna.range[rep3, ]
    
    rep1.react <- react$reactivity[rep1.range[1]:rep1.range[2]]
    rep2.react <- react$reactivity[rep2.range[1]:rep2.range[2]]
    rep3.react <- react$reactivity[rep3.range[1]:rep3.range[2]]
    
    return(data.frame(rep1.react, rep2.react, rep3.react))
}

snr.bootrep <- function(boot.means, rna.range) {
    # find the mean SNR per residue for each RNA
    output <- numeric()
    for (i in 1:24) {
        sub.mean <- boot.means$reactivity[rna.range[i, 1]:rna.range[i, 2]]
        sub.sd <- sqrt(boot.means$react.var[rna.range[i, 1]:rna.range[i, 2]])
        output <- c(output, mean(sub.mean / sub.sd, na.rm = T))
        #output <- c(output, mean(sub.mean, na.rm = T))
    }
    return(output)
}


