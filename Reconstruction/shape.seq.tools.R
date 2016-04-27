library('PearsonDS')




reactivity.quicksplit <- function(reactivity.list, replicates) {
  
  #'replicates' should be 1 or 3, so that you're isolating a single experiment or grouping each by RNA assayed
  
  start.stop.points <- c(which(reactivity.list$residue == 1), length(reactivity.list$residue) + 1)
  #print(start.stop.points)
  output.list <- list()
  
  rna.mat <- matrix(nrow = (length(start.stop.points) - 1) / replicates, ncol = 2)
  
  for (i in 1:((length(start.stop.points) - 1) / replicates)) {
    
    rna.mat[i, 1] <- start.stop.points[(i - 1) * replicates + 1]
    rna.mat[i, 2] <- start.stop.points[i * replicates + 1] - 1
    green <- start.stop.points[(i - 1) * replicates + 1]
    red <- start.stop.points[i * replicates + 1] - 1
    
    output.list[[i]] <- reactivity.list[green:red, ]
#     output.list[[i]] <- reactivity.list[start.stop.points[(i - 1) * replicates + 1]:start.stop.points[i * replicates + 1], ]
#     print(start.stop.points[(i - 1) * replicates + 1])
#     print(start.stop.points[i * replicates + 1] - 1)
  }
  #return(output.list)
  return(rna.mat)
}
  

reactivitylooper.shape2list <- function(shape2.list, estimate.flag, norm.flag, log.flag, zero.flag) {
  #shape-seq data
  output.df <- data.frame(residue = numeric(), paired = numeric(), reactivity = numeric())
  minus.stop.names <- c('minus.counts', 'minus.counts.1', 'minus.counts.2')
  plus.stop.names <- c('plus.counts', 'plus.counts.1', 'plus.counts.2')
  minus.pass.names <- c('minus.pass', 'minus.pass.1', 'minus.pass.2')
  plus.pass.names <- c('plus.pass', 'plus.pass.1', 'plus.pass.2')
  count <- 0
    
  
 for (i in 1:length(shape2.list)) {
#   for (i in 1:1){
    for (j in 1:3) {
    
      # for shapeseq
      react <- reactivity.test(shape2.list[[i]][minus.stop.names[j]], 
                         shape2.list[[i]][minus.pass.names[j]],
                         shape2.list[[i]][plus.stop.names[j]],
                      shape2.list[[i]][plus.pass.names[j]], 
                      estimate.flag, norm.flag, log.flag, zero.flag)
      #print(acf(react))
           meanvar.equation <- manual.mean.var(shape2.list[[i]][minus.stop.names[j]],
             shape2.list[[i]][minus.pass.names[j]],
             shape2.list[[i]][plus.stop.names[j]],
             shape2.list[[i]][plus.pass.names[j]], 
             estimate.flag, norm.flag, log.flag, zero.flag)
      
      sub.df <- data.frame(residue = shape2.list[[i]]['residue'],
                           paired = shape2.list[[i]]['pairing'], 
                           sequence = shape2.list[[i]]['sequence'],
                           reactivity = react, calc.mean = meanvar.equation[, 1],
                           calc.var = meanvar.equation[, 2])
      #print(sub.df)
#       sub.df$pairing <- RNA.pairing.sorter(sub.df$pairing)
#       sub.df$pairing[sub.df$pairing == 2] <- 1
      #remove all residues that are not in actual RNA (including complete pass throughs)
      sub.df <- sub.df[sub.df$residue > 0, ]
      output.df <- rbind(output.df, sub.df)
      #print(mean(react), na.rm =  T)
      count <- count + 1
      #print(count)
    }
 }
  return(output.df)
}

reactivity.test <- function(minus.stop, minus.pass, plus.stop, plus.pass,
                            estimate.flag, norm.flag, log.flag, zero.flag) {
  #react.calculator <- function(order, sequence, minus, plus, pairing, log.flag, 
  #                             estimate.flag, norm.flag, subset.flag, zero.flag) {
    # Flags: 
    # log.flag: 0 = no log, 1 = log 
    # estimate.flag: 1 = Ding Formula, 2 = Talkish Formula, 3 = Background Balance Formula
    # norm.flag: 1 = 2% 8% norm, 2 = boxplot norm
    # subset.flag: 0 = no subsetting, 1 = yes subsetting
    # zero.flag: 0 = no zero'ing of negative reactivities after normalization, 1 = yes zero'ing
  
  #output.df <- react.calculator(0, 0, minus.stop, plus.stop, 0, log.flag = 1, 1, 2, 0, 0) #ding
  #output.df <- react.calculator(0, 0, minus.stop, plus.stop, 0, log.flag = 0, 2, 2, 0, 1) # talkish normalized
 # output.df <- manual.ML.calc(0, minus.stop, minus.pass, plus.stop, plus.pass, 0)[, 1] 
  
  #for some reason, gets presented as dataframe. Converting to vector
#   minus.stop <- minus.stop[1:dim(minus.stop)[1], ]
#   plus.stop <- plus.stop[1:dim(plus.stop)[1], ]
#   minus.pass <- minus.pass[1:dim(minus.pass)[1], ]
#   plus.pass <- plus.pass[1:dim(plus.pass)[1], ]
 
  #estimate.flag, norm.flag, log.flag, zero.flag
  output.df <- local.coverage.react.calc(minus.stop, minus.pass, plus.stop, plus.pass, 
                                         estimate.flag, norm.flag, log.flag, zero.flag)
#   print(output.df)
  return(output.df)
}

likelihood.ratio <- function(shape.data, gap, beta.gap) {
  paired <- hist(shape.data$ding.react[shape.data$pairing != 0 & shape.data$residue > 0],
                 freq = F, ylim = c(0, 8), xlim = c(0, 2), breaks = seq(0,10, gap), 
                 main = 'Paired reactivities')
  unpaired <- hist(shape.data$ding.react[shape.data$pairing == 0 & shape.data$residue > 0],
                   freq = F, ylim = c(0, 8), xlim = c(0, 2), breaks = seq(0,10, gap), 
                   main = 'unpaired reactivities')
  plot(seq(0, 10 - gap, gap), paired$density / unpaired$density,
        main = 'Likelihood ratio of reactivities', xlim = c(0, 2))
  
  paired.beta <- hist(shape.data$beta[shape.data$pairing != 0 & shape.data$residue > 0],
                 freq = F, ylim = c(0, 10), xlim = c(0, .2), breaks = seq(0,10, beta.gap), 
                 main = 'Paired betas')
  unpaired.beta <- hist(shape.data$beta[shape.data$pairing == 0 & shape.data$residue > 0],
                   freq = F, ylim = c(0, 10), xlim = c(0, .2), breaks = seq(0,10, beta.gap), 
                   main = 'unpaired betas')
  print(unpaired.beta$density[1:20])
  print(paired.beta$density[1:20])
  plot(seq(0, 10 - beta.gap, beta.gap), paired.beta$density / unpaired.beta$density, 
        main = 'Likelihood ratio of betas', xlim = c(0, .1))
}

coverage.creater <- function(counts.df) {
# just sum the stops of every stop after each desired residue
  counts.df$minus.pass <- 0
  counts.df$plus.pass <- 0
#   counts.df$minus.pass.1 <- 0
#   counts.df$plus.pass.1 <- 0
#   counts.df$minus.pass.2 <- 0
#   counts.df$plus.pass.2 <- 0
  r <- length(counts.df$minus.pass)
  for (i in 2:r) {
    counts.df$minus.pass[i] <- sum(counts.df$minus.counts[1:i - 1])
    counts.df$plus.pass[i] <- sum(counts.df$plus.counts[1: i - 1])
#     counts.df$minus.pass.1[i] <- sum(counts.df$minus.counts.1[1:i - 1])
#     counts.df$plus.pass.1[i] <- sum(counts.df$plus.counts.1[1: i - 1])
#     counts.df$minus.pass.2[i] <- sum(counts.df$minus.counts.2[1:i - 1])
#     counts.df$plus.pass.2[i] <- sum(counts.df$plus.counts.2[1: i - 1])
  }
#   minus.names <- grep('minus', names(counts.df))
#   plus.names <- grep('plus', names(counts.df))
#   
#   for (j in 1:length(minus.names)) {
#     
#   }
#   
#   counts.df$ding.react <- react.calculator(counts.df$residue, counts.df$sequence, 
#                                            counts.df$minus.counts, counts.df$plus.counts,
#                                            counts.df$pairing, 1, 1, 1, 0, 1)
#   counts.df$ding.react.1 <- react.calculator(counts.df$residue, counts.df$sequence, 
#                                            counts.df$minus.counts.1, counts.df$plus.counts.1,
#                                            counts.df$pairing, 1, 1, 1, 0, 1)
#   counts.df$ding.react.2 <- react.calculator(counts.df$residue, counts.df$sequence, 
#                                            counts.df$minus.counts.2, counts.df$plus.counts.2,
#                                            counts.df$pairing, 1, 1, 1, 0, 1)
#   counts.df <- data.frame(counts.df, manual.ML.calc(residues = counts.df$residue, counts.df$minus.counts,
#                                                     counts.df$minus.pass, counts.df$plus.counts, 
#                                                     counts.df$plus.pass, 0))
  return(counts.df)
}

single.rna.maker <- function(stacked.rna, flag) {
  keep.names <- names(stacked.rna[, 1:5])
  if (flag == 1) {
    split.rna <- stacked.rna[, c(1,2,3,4,5)]
  }
  if (flag == 2) {
    split.rna <- stacked.rna[, c(1,2,3,6,7)]
  }
  if (flag == 3) {
    split.rna <- stacked.rna[, c(1,2,3,8,9)]
  }
  names(split.rna) <- keep.names
  split.rna <- coverage.creater(split.rna)
  return(split.rna)
}

# what happen if you perturb the real data in simulations, introduce random distortions, see how robust it is, does it distort the data. 
# what if one equation stabilizes variance well, but disturbs the general pictures? Reduce sequencing coverage
# Measure preturbations at certain threshold: does it change by 15%? What fraction of nucleotides stay in window after perturbation?
# A possible way to look at it: process the 0's and non-zeros as different statistics

count.ct.combiner <- function(counts, ct.file) {
  residue <- c(-999, ct.file[, 1])
  sequence <- c('0', as.character(ct.file[, 2]))
  pairing <- c(0, ct.file[, 3])
#   minus.counts <- counts$minus.counts + counts$minus.counts.1 + counts$minus.counts.2
#   plus.counts <- counts$plus.counts + counts$plus.counts.1 + counts$plus.counts.2
  output.df <- data.frame(residue, sequence, pairing, counts)
  return(output.df)
}

RNA.pairing.sorter <- function(native) {
  # reads native.ct files and determines if each residue is unpaired or paired
  # current scoring system is 1 = unpaired, 2 = helix end, 3 = stacked
  index <- 0
  for (i in 1:length(native)) {
    flag.5prime <- 0 #flags indicate whether or not the 5' or 3' ends are paired. 
    flag.3prime <- 0
    if (i > 1) {
      if (native[i - 1] != 0) {
        flag.5prime <- 1
      }   
    }
    if (i < length(native) - 1) {
      if (native[i + 1] != 0) {
        flag.3prime <- 1
      }
    }
    if (native[i] == 0) { #if unpaired
      index[i] <- 0
    } 
    else if (flag.5prime == 1 & flag.3prime == 1) { # if stacked, paired on either side
      index[i] <- 2
    }
    else { #everything else...which at this point is helix ends
      index[i] <- 1
    }
  }
  return(index)
}



ct.reader <- function(ct.file) {
  output.df <- read.table(ct.file, skip = 1)[c(1, 2, 5)]
  names(output.df) <- c('residue', 'sequence', 'paired')
  return(output.df)
}

rdat.looper <- function(rdat.namelist) {
  output.list <- list()
  for (i in 1:length(rdat.namelist)) {
    print(rdat.namelist[i])
    rdat.subset <- rdat.reader(rdat.namelist[i])
    output.list[[i]] <- rdat.subset 
    names(output.list)[i] <- rdat.namelist[i]
  }
  return(output.list)
}

rdat.reader <- function(rdat.name) {
  file.names <- list.files(pattern = rdat.name)
  print(file.names)
  for (i in 1:length(file.names)) {
    lines <- readLines(file.names[i])
    read.location <- grep('READS', lines)
    plus.counts <- as.numeric(unlist(strsplit(lines[read.location[1]], split = '\t'))[-1:-3])
    minus.counts <- as.numeric(unlist(strsplit(lines[read.location[2]], split = '\t'))[-1:-3])
    print(length(plus.counts))
    print(length(minus.counts))
    if (i == 1) {output.df <- data.frame(plus.counts, minus.counts)}
    else {output.df <- data.frame(output.df, plus.counts, minus.counts)}
  }

  return(output.df)
}

get.namelist <- function() {
  output <- strtrim(list.files(pattern = '_0001'), 7)
  return(output)
}




unlist <- function (x) {
  ## do.call(c, ...) coerces factor to integer, which is undesired
  # taken from Rsamtools documentation
  x1 <- x[[1L]]
  if (is.factor(x1)) {
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}