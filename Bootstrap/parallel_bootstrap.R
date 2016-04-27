###################################################
#Code for 10000 bootstraps
###################################################

downsampling_factors <- c(1,10,100,1000)
names_df <- c("full", "down_10", "down_100", "down_1000")
bundle_size <- 100
num_bundles <- 200

for (i in downsampling_factors) {
    for (j in 1:num_bundles) {
        save(bootstrap(bundle_size, shape2.coveragelist, i), file=paste(names_df[which(downsampling_factors==i)], "_","j", sep=""))
    }
}


#----------------------------------------------------------------------------------------
#A function for bootstrapping data.
#----------------------------------------------------------------------------------------
bootstrap <- function(num_bootstrapSamples,data,downsample_factor=1, random_primer = FALSE) {
    
    #a function to tell foreach how to combine lists of kind data as they are generated in parallel
    combine <- function(...) {
        
        all_lists <- list(...)
        list.combined <- list()
        
        pass_cal <- function(counts) return(c(0,cumsum(counts[1:(length(counts)-1)])))
        
        if (random_primer) {
            list.combined$minus <- lapply(all_lists, '[[', 'minus')
            list.combined$plus <- lapply(all_lists, '[[', 'plus')
        } else {
            
            list.combined$minus.counts <- sapply(all_lists, '[[', 'minus.counts')
            list.combined$minus.pass <- sapply(all_lists, function(x) pass_cal(x$minus.counts))
            #list.combined$minus.pass <- sapply(all_lists, '[[', 'minus.pass')
            list.combined$plus.counts <- sapply(all_lists, '[[', 'plus.counts')
            list.combined$plus.pass <- sapply(all_lists, function(x) pass_cal(x$plus.counts))
            #list.combined$plus.pass <- sapply(all_lists, '[[', 'plus.pass')
        }
        return(list.combined)
    }
    
    transcript_names <- names(data)
    num_transcripts <- length(transcript_names)
    boot <- list()
    
    for (i in 1:num_transcripts) {
        
        num_biological_replicates <- length( if (random_primer) data[[i]] else ls(data[[i]], pat= "minus.counts") )
        
        for (j in 1:num_biological_replicates) {
            
            if (random_primer) {
                reads <- data[[i]][[j]]
                sampleSize <- list(minus="", plus="")
                freq <- list(minus =reads$minus/sampleSize$minus, plus= reads$plus/sampleSize$plus)
            } else {
                reads <- data[[i]]
                freq <- list(minus.counts="", plus.counts="")
                sampleSize <- list(minus.counts="", plus.counts="")
                
                if (j==1) {
                    sampleSize$minus.counts <- floor(sum(reads$minus.counts)/downsample_factor)
                    #sampleSize$minus.pass <- sum(reads$minus.pass)
                    sampleSize$plus.counts <- floor(sum(reads$plus.counts)/downsample_factor)
                    #sampleSize$plus.pass <- sum(reads$plus.pass)
                    
                    ###############################
                    #####probable error, sampleSize$minus.counts should be replaced with sum(reads$minus.counts)
                    ###############################
                    freq$minus.counts <- reads$minus.counts/sum(reads$minus.counts)
                    #freq$minus.pass <- reads$minus.counts/sampleSize$minus.pass
                    freq$plus.counts <- reads$plus.counts/sum(reads$plus.counts)
                    #freq$plus.pass <- reads$minus.counts/sampleSize$plus.pass
                } else {
                    eval( parse( text=paste("sampleSize$minus.counts <- floor(sum(reads$minus.counts.", j-1,")","/downsample_factor)", sep="")))
                    #eval( parse( text=paste("sampleSize$minus.pass <- sum(reads$minus.pass.", j-1,")", sep="")))
                    eval( parse( text=paste("sampleSize$plus.counts <- floor(sum(reads$plus.counts.", j-1,"/downsample_factor))", sep="")))
                    #eval( parse( text=paste("sampleSize$plus.pass <- sum(reads$plus.pass.", j-1,")", sep="")))
                    
                    eval( parse( text= paste("freq$minus.counts <- reads$minus.counts.", j-1,"/sum(reads$minus.counts)",sep="")))
                    eval( parse( text= paste("freq$plus.counts <- reads$plus.counts.", j-1, "/sum(reads$plus.counts)",sep="")))
                }
                
            }
            
            boot_current <- foreach(i=1:num_bootstrapSamples , .combine = combine, .multicombine=TRUE) %dopar% SampleGenerator(freq, sampleSize, random_primer)
            eq <- paste("boot","$",transcript_names[i],"$","rep", j, " <- ", "boot_current",sep="" )
            eval(parse(text= eq))
            
            rm(boot_current)
        }
    }
    
    return(boot)
}

#----------------------------------------------------------------------------------------
#A function to generate minus channel and plus channel fragment counts
#Requires frequency of k-fragments in minus and plus channel and target sample sizes for paired end reads.
#Requires stops and passes for minus and plus channel and sample sizes for single end reads.
#----------------------------------------------------------------------------------------
SampleGenerator <- function (freq,sampleSize, random_primer) {
    
    print(sampleSize)
    if (random_primer) {
        
        n <- nrow(freq$minus)
        m <- ncol(freq$minus)
        
        frag_count_minus <- matrix(, nrow= n, ncol= m)
        frag_count_plus <- matrix(, nrow= n, ncol= n)
        
        freq_vector_minus <- as.vector(freq$minus)
        freq_vector_plus <- as.vector(freq$plus)
        frag_map_sites_minus <- sample(1:(m*n),sampleSize$minus,rep=TRUE,prob=freq_vector_minus)
        frag_map_sites_plus <- sample(1:(m*n),sampleSize$plus,rep=TRUE,prob=freq_vector_plus)
        
        jk_frag_counts_minus <- tabulate(frag_map_sites_minus, nbins= m*n)
        jk_frag_counts_plus <- tabulate(frag_map_sites_plus, nbins= m*n)
        
        frag_count_minus <- matrix(jk_frag_counts_minus, nrow= n, ncol= m)
        frag_count_plus <- matrix(jk_frag_counts_plus, nrow= n, ncol= m)
        
        sample_gen <- list(minus= frag_count_minus, plus= frag_count_plus)
        
    } else {
        n <- length(freq$minus.counts)
        sample_gen <- list(minus.counts="", plus.counts="")
        sample_gen$minus.counts <- tabulate(sample(n, sampleSize$minus.counts, rep= TRUE, prob= freq$minus.counts), nbins=n)
        #sample_gen$minus.pass <- tabulate(sample(n, sampleSize$minus.pass, rep= TRUE, prob= freq$minus.pass), nbins=n)
        sample_gen$plus.counts <- tabulate(sample(n, sampleSize$plus.counts, rep= TRUE, prob= freq$plus.counts), nbins=n)
        #sample_gen$plus.pass <- tabulate(sample(n, sampleSize$plus.pass, rep= TRUE, prob= freq$plus.pass), nbins=n)
    }
    
    return(sample_gen)
}

#----------------------------------------------------------------------------
#for the variance vs mean for various flag combinations for 10000x bootstrap 09/15/15
jknife_subset_sizes <- c(10,100,1000,10000)
downsampling_cases <- c("full", "down_10", "down_100", "down_1000")

estimate.flags <- 1
norm.flags <- 1
log.flags <- 0
zero.flags <- 1

for (down_by in 1:length(downsampling_cases)) {
    file_names <- list.files("/Users/aviran/Desktop/new_bootstrap_91815", pattern=downsampling_cases[down_by])
    
    for (bundle in 1:length(file_names)) {
        
        load(paste("/Users/aviran/Desktop/new_bootstrap_91815/", file_names[bundle], sep=""))
        data <- x
        transcript_names <- names(data)
        num_transcripts <- length(transcript_names)
        num_biological_replicates <- length( data[[1]] )
        reactivities <- list()
        for (p in 1:num_transcripts) {
            for (q in 1:num_biological_replicates) {
                eval(parse(text=paste("reactivities$",transcript_names[p],"_",q, "<- cbind()", sep="")))
            }
        }
        
        for (i in 1:num_transcripts) {
            
            num_biological_replicates <- length( data[[i]] )
            #n <- nrow(reads)
            
            for (j in 1:num_biological_replicates) {
                
                reads <- data[[i]][[j]]
                n <- nrow(reads$plus.pass)
                
                for (bundle_size in 1:100) {
                    eval(parse(text=paste("reactivities$",transcript_names[i],"_",j, "<- cbind(reactivities$",transcript_names[i],"_",j, ", matrix(local.coverage.react.calc(reads$minus.counts[5:n, bundle_size],
                    reads$minus.pass[5:n, bundle_size],
                    reads$plus.counts[5:n, bundle_size],
                    reads$plus.pass[5:n, bundle_size],
                    estimate.flags,norm.flags,
                    log.flags,zero.flags), nrow=n-4, ncol=1))", sep="")))
                    
                }
                
            }
            
        }
        
    }
    # }
    
    #}
    #}
    #}
    
    save(reactivities, file=paste("/Users/aviran/Desktop/new_bootstrap_91815/reactivities_" ,downsampling_cases[down_by], ".Rdata", sep=""))
    
}

num_of_subsets <- 100
num_transcripts <- length(reactivities_full)
for (down_by in 1:length(downsampling_cases)) {
    load(paste("/Users/aviran/Desktop/new_bootstrap_91815/reactivities_", downsampling_cases[down_by], ".Rdata", sep=""))
    
    for ( subset_sizes in c(10,100,1000,10000) ) {
        mean_var <- c()
        var_var <- c()
        for (i in 1:num_transcripts) {
            n <- nrow(reactivities[[i]])
            variances <- matrix(, nrow= n)
            for(num_subsets in 1:num_of_subsets) {
                
                subset_columns <- floor(runif(subset_sizes, 1, 30000))
                subset <- reactivities[[i]][,subset_columns]
                
                variances <- cbind(variances, matrix(apply(subset,1,var, na.rm= TRUE), nrow=n, ncol=1))
            }
            mean_var <- c(mean_var, apply(variances, 1, mean, na.rm=TRUE))
            var_var <- c(var_var, apply(variances, 1, var, na.rm= TRUE))
        }
        setEPS()
        postscript(file=paste("/Users/aviran/Desktop/new_bootstrap_91815/", downsampling_cases[down_by], "_subset_size_",
        subset_sizes,"_linear_scale.eps", sep=""))
        plot(mean_var, var_var)
        dev.off()
        
        setEPS()
        postscript(file=paste("/Users/aviran/Desktop/new_bootstrap_91815/", downsampling_cases[down_by], "_subset_size_",
        subset_sizes,"_log_scale.eps", sep=""))
        plot(mean_var, var_var, log="xy")
        dev.off()
    }
}





####################################################
#Code for getting best fit to bootsrap mean vs variance data
####################################################
method <- "lm"
#method <- "loess"
cutoff <- -4
data<- data.frame(means = read.csv(file="diff.MR_28_zeroed.csv", header=TRUE, sep=",")$reactivity,
sdevs = read.csv(file="diff.MR_28_zeroed.csv", header=TRUE, sep=",")$react.sd)
data <- log10(data[with(data, order(means)), ])
data <- data[-which(is.infinite(as.numeric(data[,1])) | as.numeric(data[,1]) < cutoff), ]

eval(parse(text= paste("fit <- ", method, "(sdevs ~ means, data= data)", sep="")))
#fit  <- loess(sdevs ~ means, data= data)
plot(data)
lines(data$means, predict(fit), col = "blue")
fit

#use following to get predicted variance if using loess
#predict(fit, ...) #replace ellipsis with your mean value.

# abline(fit)
# 
# nls_fit <- nls(sdevs ~ sqrt(means + b * means^(2)), data, start = list(b = -1))
# lines(Data$x, predict(nls_fit), col = "red")
# pr <- predict(fit)


