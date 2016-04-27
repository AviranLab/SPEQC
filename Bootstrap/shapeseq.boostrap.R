source('~/Box Sync/Information_Content_in_Probes/Code/Shape-Seq_data_analysis/shape.seq.tools.R')
# boot.1x.20.1000 <- plot.boot.mean.var(get.boot.mean(shapeseq.react.looper(shape2boot.1x.20, 1, 0, 0, 0)))

multi.boot.react.looper <- function(multi.boot.list, estimate.flag, norm.flag, log.flag, zero.flag) {
  output.list <- list()
  
  for (i in 1:length(multi.boot.list)) {
    sub.boot <- multi.boot.list[[i]]
    react.data <- mean.list.to.df(get.boot.mean(shapeseq.react.looper(sub.boot, estimate.flag,
                                                                 norm.flag, log.flag, zero.flag)))
    output.list[[i]] <- react.data
  }
  return(output.list)
}





mean.boot.counts <- function(boot.list) {
  output.list <- list()
  
  for (i in 1:length(boot.list)) { # 8 rnas
    sub1 <- boot.list[[i]]
    return.df <- data.frame(residue = boot.list[[i]][[1]][[1]]$residue, 
                            sequence = boot.list[[i]][[1]][[1]]$sequence, pairing = boot.list[[i]][[1]][[1]]$pairing)
    for (j in 1:length(sub1)) { # 3 replicates
      sub2 <- sub1[[j]]
      
      minus.counts <- sub2[[1]]$minus.counts
      plus.counts <- sub2[[1]]$plus.counts
      minus.pass <- sub2[[1]]$minus.pass
      plus.pass <- sub2[[1]]$plus.pass
      
      for (k in 2:length(sub2)) {
        # sum counts in each category
        sub3 <- sub2[[k]]
        
        minus.counts <- minus.counts + sub3$minus.counts
        plus.counts <- plus.counts + sub3$plus.counts
        minus.pass <- minus.pass + sub3$minus.pass
        plus.pass <- plus.pass + sub3$plus.pass 
      }
      # average counts
      minus.counts <- minus.counts / length(sub2)
      plus.counts <- plus.counts / length(sub2)
      minus.pass <- minus.pass / length(sub2)
      plus.pass <- plus.pass / length(sub2)
      return.df <- data.frame(return.df, minus.counts = minus.counts, plus.counts = plus.counts,
                              minus.pass = minus.pass, plus.pass = plus.pass)
    }
    output.list[[i]] <- return.df 
  }
  return(output.list)
}


get.coverage <- function(shape2list) {
  minus.stop.names <- c('minus.counts', 'minus.counts.1', 'minus.counts.2')
  plus.stop.names <- c('plus.counts', 'plus.counts.1', 'plus.counts.2')
  minus.pass.names <- c('minus.pass', 'minus.pass.1', 'minus.pass.2')
  plus.pass.names <- c('plus.pass', 'plus.pass.1', 'plus.pass.2')  

  output.df <- data.frame(minus = numeric(), plus = numeric())
  for (i in 1:length(shape2list)) {
    sub1 <- shape2list[[i]]
    for (j in 1:3) {
       minus.cov <- sum(sub1[[minus.stop.names[j]]] + sub1[[minus.pass.names[j]]])
       plus.cov <- sum(sub1[[plus.stop.names[j]]] + sub1[[plus.pass.names[j]]])
       
       output.df <- rbind(output.df, data.frame(minus = minus.cov, plus = plus.cov))
    }
  }
  return(output.df)
}




mean.list.to.df <- function(mean.list) {
  output.df <- data.frame(pairing = numeric(), mean = numeric(), var = numeric(), zeroes = numeric())
  for (i in 1:length(mean.list)) {
    sub.1 <- mean.list[i]
    for(j in 1:length(sub.1)) {
      sub2 <- sub.1[[j]]
      for (k in 1:length(sub2)) {
      output.df <- rbind(output.df, sub2[[k]])
      }
    }
  }
  output.df$react.sd <- sqrt(output.df$var)
  names(output.df) <- c('pairing', 'reactivity', 'react.var', 'zeroes', 'react.sd')
  return(output.df)
} 


quick.snr.plot <- function(stats.input, window.length, append.df) {
  over.five.snr <- numeric()
  snr.average <- numeric()
  for (i in 1:length(stats.input)) {

    for (j in 1:length(stats.input[[i]])) {
      snr <- stats.input[[i]][[j]]$mean / sqrt(stats.input[[i]][[j]]$var)
      snr[is.na(snr)] <- 0 # to visualize which have mostly zeroes
#       plot(snr)
#       lines(seq(0,400,1), rep(5,401))
      
      #a plot with bar graphs of paired or unpaired, and averaged snr scores
      residue <- seq(1, length(snr), 1)
      new.pairing <- stats.input[[i]][[j]]$pairing 
      new.zeroes <- stats.input[[i]][[j]]$zeroes

      new.pairing[new.pairing != 0] <- 1
      new.pairing[new.pairing == 0] <- 5

      output.df <- data.frame(residue, snr, mean.snr = 0)
     
      for (k in 1:length(residue)) {
        if (k - window.length > 0) {
          start <- k - window.length
        } else {
          start <- 1
        }
        if (k + window.length < length(residue)) {
          end <- k + window.length
        } else {
          end <- length(residue)
        }
        output.df$mean.snr[k] <- mean(snr[start:end])        
      }
      plot(output.df$mean.snr, main = c(i,j), type = 'l')
      x <- sum(output.df$mean.snr > 5) / length(output.df$mean.snr)
      over.five.snr <- c(over.five.snr, x)
      y <- mean(output.df$mean.snr)
      snr.average <- c(snr.average, y)
#       print(over.five.snr)
#       print(snr.average)
            
      #lines(seq(0,400,1), rep(5,401))
      #lines(new.pairing, col = 'red')
      #lines(10 - new.zeroes/2, col = 'blue')
    }
  }
  summary.df <- data.frame(append.df, over.five = over.five.snr, snr.average = snr.average)

  return(summary.df)
}

plot.boot.mean.var <- function(sum.stat.list) {
  zero.limit <- 20
  
  plot.df <- data.frame(react.mean = numeric(), react.var = numeric(), sd = numeric(), zeroes = numeric())
  for (i in 1:length(sum.stat.list)) {
    #print(sum.stat.list[[i]])
    for(j in 1:length(sum.stat.list[[i]])) {
      x <- sum.stat.list[[i]][[j]]
      print(x)
      plot.df <- rbind(plot.df, data.frame(react.mean = x$mean[x$zeroes <= zero.limit], 
                                            react.var = x$var[x$zeroes <= zero.limit],
                                            react.sd = sqrt(x$var[x$zeroes <= zero.limit]),
                                            zeroes = x$zeroes[x$zeroes <= zero.limit]))
      
    }
  }
  y <- plot.df
  #plot(y$react.mean, y$react.var , main = deparse(substitute(sum.stat.list)), ylim = c(0,.5))
  plot(y$react.mean, y$react.var , main = deparse(substitute(sum.stat.list)))
  
  #print((y[y$react.mean > 5 & !is.na(y$react.var), ]))
  #hist(y$react.var[y$react.mean > 5 & !is.na(y$react.var)], breaks = seq(0,5, .000001))
  return(plot.df)
}

get.boot.mean <- function(react.output) {
  output.list <- list()
  for (i in 1:length(react.output)) {
    sub1 <- react.output[[i]]
    output.list[[i]] <- list()
    for (j in 1:length(sub1)) {
      sub2 <- sub1[[j]]
      
      combine.react <- data.frame(residue = sub2[[1]]$residue, pairing = sub2[[1]]$pairing)
      for(k in 1:length(sub2)) {
          combine.react <- data.frame(combine.react, react = sub2[[k]]$reactivity)
      }
      x <- (combine.react[, 3:dim(combine.react)[2]])
      react.mean <- apply(x, 1, mean, na.rm = T)
      react.var <- apply(x, 1, var, na.rm = T)
      neg.reacts <- rep(0, dim(x)[1])
      for (l in 1:dim(x)[1]){
        #neg.reacts[l] <- sum(is.na(x[l, ]))
        neg.reacts[l] <- sum(x[l, ] <= 0 | is.na(x[l, ]))
      }
      output.list[[i]][[j]] <- data.frame(pairing = combine.react$pairing, mean = react.mean, 
                                          var = react.var, zeroes = neg.reacts)

    }
  }
  return(output.list)
}


shapeseq.react.looper <- function(shape2.list, estimate.flag, norm.flag,
                                  log.flag, zero.flag) {
  output.list <- list()
  for (i in 1:length(shape2.list)) {
    output.list[[i]] <- list()
    sub1 <- shape2.list[[i]]
    for (j in 1:length(sub1)) {
      output.list[[i]][[j]] <- list()
      sub2 <- sub1[[j]]
      boot.target <- sub2 ### was this a bug?! It originally said sub2[[j]]. Actually, this isn't targeted at anything. Seems okay
      for (k in 1:length(sub2)) {
       sub3 <- sub2[[k]] 
       react <- local.coverage.react.calc(sub3$minus.counts, sub3$minus.pass, 
                          sub3$plus.counts, sub3$plus.pass, estimate.flag, 
                          norm.flag, log.flag, zero.flag)
       
       
        mean.var <- manual.mean.var(sub3$minus.counts, sub3$minus.pass, 
                        sub3$plus.counts, sub3$plus.pass, estimate.flag, norm.flag,
                        log.flag, zero.flag)
       
       
       output.list[[i]][[j]][[k]] <- data.frame(residue = sub3$residue, 
                                                pairing = sub3$pairing, reactivity = react,
                                                mean.var)
       output.list[[i]][[j]][[k]] <- output.list[[i]][[j]][[k]][
         output.list[[i]][[j]][[k]]$residue > 0, ]
       
      }
    }
  }
  return(output.list)
}


shapeseq.bootstrap <- function(residue, sequence, pairing, minus.counts, plus.counts, scale.it) {
  #creates a list in a list in a list. 
  # 1st list: RNA type
  # 2nd list: replicate experiment
  # 3rd list: bootstrap replicates of each experiment. in data.frames
  stop.list <- list(minus.counts, plus.counts)
  output.df <- data.frame(residue, sequence, pairing)
  for (j in 1:2) {
    stop.counts <- stop.list[[j]]
    reads.sim <- sum(stop.counts) * scale.it
    boot.it <- sample(length(stop.counts), reads.sim, replace = T, prob = stop.counts / sum(stop.counts))
    boot.table <- table(boot.it)
    boot.residue <- as.numeric(names(boot.table))
    boot.vector <- rep(0, length(stop.counts))
    for (i in 1:length(boot.residue)) {
      boot.vector[boot.residue[i]] <- boot.table[i]
      count.pass <- 0
      for (k in 2:length(boot.vector)) {
      count.pass[k] <- sum(boot.vector[1:k - 1])
      }
    }
  output.df <- data.frame(output.df, boot.vector, count.pass)
  }
  names(output.df) <- c('residue', 'sequence', 'pairing',
                        'minus.counts', 'minus.pass', 'plus.counts', 'plus.pass')
  return(output.df)
}



split.shape2list <- function(shape2.list, boot.times, scale.value) {
  ptm <- proc.time()
  name.transfer <- names(shape2.list)
  #boot.times <- 20
  minus.stop.names <- c('minus.counts', 'minus.counts.1', 'minus.counts.2')
  plus.stop.names <- c('plus.counts', 'plus.counts.1', 'plus.counts.2')
  minus.pass.names <- c('minus.pass', 'minus.pass.1', 'minus.pass.2')
  plus.pass.names <- c('plus.pass', 'plus.pass.1', 'plus.pass.2')
  output.list <- list()
  
  for (i in 1:length(shape2.list)) {
    output.list[[i]] <- list()
    for (j in 1:3){ 
      output.list[[i]][[j]] <- list()
      split.this <- shape2.list[[i]]
      reassemble <- data.frame(residue = split.this$residue, sequence = split.this$sequence, 
                             pairing = split.this$pairing)
      reassemble <- data.frame(reassemble, minus.stops = split.this[minus.stop.names[j]], 
                               minus.pass = split.this[minus.pass.names[[j]]], 
                               plus.stops = split.this[plus.stop.names[[j]]],
                               plus.pass = split.this[plus.pass.names[[j]]])
      for (k in 1:boot.times) {
        booted.df <- shapeseq.bootstrap(reassemble$residue, reassemble$sequence, reassemble$pairing, 
                         reassemble$minus.counts, reassemble$plus.counts, scale.value)
        output.list[[i]][[j]][[k]] <- booted.df
        print(proc.time() - ptm)
      }
    }
  }
  names(output.list) <- name.transfer
  return(output.list)
}

quick.combine <- function(stat.output) {
  output.df <- data.frame(pairing = numeric(), mean = numeric(), var = numeric(), zeroes = numeric())
  for (i in 1:length(stat.output)) {
    sub1 <- stat.output[[i]]
    for (j in 1:length(sub1)) {
      sub2 <- unlist(sub1[[j]])
      output.df <- rbind(output.df, sub2)
    }
  }
  return(output.df)
}

quick.combine.single <- function(stat.output) {
  output.df <- data.frame(residue = numeric(), pairing = numeric(), reactivity = numeric())
  for (i in 1:length(stat.output)) {
    sub1 <- stat.output[[i]]
    for (j in 1:length(sub1)) {
      sub2 <- unlist(sub1[[j]])
      output.df <- rbind(output.df, sub2)
    }
  }
  return(output.df)
}

downsample.stats <- function(down.none, down.10, down.100, down.1000) {
  boxplot(list(down.none, down.10, down.100, down.1000) )
  
  b.none <- boxplot(down.none, plot = F)
  begin <- b.none$stats[[1]]
  end <- b.none$stats[[5]]
  
  print(b.none$stats[[1]])
  print(b.none$stats[[3]])
  print(b.none$stats[[5]])
  print(length(b.none$out))
  print(length(down.none[down.none >= begin & down.none <= end]))
        
  sample.list <- list(down.10, down.100, down.1000)
  
  for (i in 1:length(sample.list)) {
    sample.this <- sample.list[[i]]
    sample.bp <- boxplot(sample.this, plot = F)
    print(sample.bp$stats[[1]])
    print(sample.bp$stats[[3]])
    print(sample.bp$stats[[5]])
    print(length(sample.bp$out))
    print(length(sample.this[sample.this >= begin & sample.this <= end]))
  }
  
}
