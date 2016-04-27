####################################################
#Code for Coverage Quality Index- testing- #v2.0 three parameters for low, medium and high reactivities
####################################################
library(xlsx)
library(foreach)
library(parallel)
library(doParallel)
registerDoParallel(cores=detectCores(all.tests=TRUE))

#1. Define ranges for low, medium and high reactivities
# a. Range definitions from Kevin Weeks
#2. Calculate requirement factor for reactivities in all three ranges
#3. In every range, find 90 percentile as the predicted index
#4. If a single index is desired, use the index for medium range of reactivities.

low_range <- c(0.1, 0.3)
medium_range <- c(0.3, 0.7)
high_range <- 0.7
#what percentile to take for giving indices
index_criteria <- c(0.95)

data <- shape2.coveragelist
accuracy <- list(confidence=c(), percent_var=c(), good_predictions=c(), index_criteria=c())
#for ( min_reactivity in c(0.001, 0.005, 0.01) ) {
min_reactivity <- 0.1
for (criteria in index_criteria) {
    for (confidence in c(0.95)) {
        for (percent_var in c(0.4)) {
            #min_reactivity <- 0.001
            #       confidence <- 0.95
            z_value <- qnorm(1- ((1-confidence)/2))
            #       percent_var <- 0.3
            
            #1. calculate demand for all transcripts
            #2. generate a bootstrap sample which satisfies the demand
            #3. bootstrap the sample to test and see how many give the output we want
            
            #1. calculate demand for all transcripts
            transcript_names <- names(data)
            num_transcripts <- length(transcript_names)
            demand <- list()
            supply <- list()
            per_residue_info <- list()
            per_residue_req <- list()
            betas <- list()
            gammas <- list()
            formula_gammas <- "minus.counts[5:n]/(minus.counts[5:n] + minus.pass[5:n])"
            #formula_betas <- "((plus.counts[5:n]/(plus.counts[5:n] + plus.pass[5:n])) - g)/(1- g)"
            formula_betas <- "local.coverage.react.calc(minus.counts[5:n], minus.pass[5:n], plus.counts[5:n], plus.pass[5:n], 1, 0, 0, 1)"
            #formula_norm <- "get.normalizer(b, 1)"
            formula_req <- "ceiling((((b + g - (b * g)) * (1- (b + g - (b * g)))) + (g*(1-g)/size_factor))/((percent_var*b/z_value)^2 ))"
            #     req_low <- c()
            #     req_medium <- c()
            #     req_high <- c()
            index_summary <- list()
            norm_constants <- c()
            req <- c()
            
            for (i in 1:num_transcripts) {
                num_biological_replicates <- length(ls(data[[i]], pat= "minus.counts"))
                #print(num_biological_replicates)
                reads <- data[[i]]
                n <- nrow(reads)
                for (j in 1:num_biological_replicates) {
                    
                    if (j==1) {
                        minus.counts <- matrix(reads$minus.counts, nrow=n, ncol=1)
                        minus.pass <- matrix(reads$minus.pass, nrow=n, ncol=1)
                        plus.counts <- matrix(reads$plus.counts, nrow=n, ncol=1)
                        plus.pass <- matrix(reads$plus.pass, nrow=n, ncol=1)
                        
                    } else {
                        eval(parse(text= paste("minus.counts <- matrix(reads$minus.counts.", j-1,", nrow=n, ncol=1)", sep="")))
                        eval(parse(text= paste("minus.pass <- matrix(reads$minus.pass.", j-1,", nrow=n, ncol=1)", sep="")))
                        eval(parse(text= paste("plus.counts <- matrix(reads$plus.counts.", j-1,", nrow=n, ncol=1)", sep="")))
                        eval(parse(text= paste("plus.pass <- matrix(reads$plus.pass.", j-1,", nrow=n, ncol=1)", sep="")))
                    }
                    
                    eval(parse(text= paste("supply$", transcript_names[i], j, "_minus <- sum(minus.counts)", sep="")))
                    eval(parse(text= paste("supply$", transcript_names[i], j, "_plus <- sum(plus.counts)", sep="")))
                    eval(parse(text= paste("size_factor <- supply$", transcript_names[i], j,
                    "_minus/supply$", transcript_names[i], j, "_plus", sep="")))
                    eval(parse(text= paste("gammas$", transcript_names[i], j, " <- matrix(", formula_gammas,", nrow=n-4, ncol=1)", sep="")))
                    eval(parse(text= paste("g <- gammas$", transcript_names[i], j, sep="")))
                    #     print(j)
                    #     print(dim(plus.counts))
                    #     print(dim(g))
                    eval(parse(text= paste("betas$", transcript_names[i], j, " <- matrix(", formula_betas, ", nrow=n-4, ncol=1)", sep="")))
                    eval(parse(text= paste("b <- betas$", transcript_names[i], j, sep="")))
                    eval(parse(text= paste("per_residue_info$", transcript_names[i], j, " <- matrix(plus.counts[5:n] + plus.pass[5:n])", sep="")))
                    eval(parse(text= paste("per_residue_req$", transcript_names[i], j, " <- matrix(", formula_req, ", nrow=n-4, ncol=1)", sep="")))
                    #index where the value of beta is just greater than min_reactivity specified by the user
                    #index_min_reactivity <- which(b >= min_reactivity)
                    #index_min_reactivity <- which(b == min(b[which(b > min_reactivity)]))
                    #print(b[index_min_reactivity])
                    norm <- get.normalizer(b, 1)
                    norm_constants <- c(norm_constants, norm)
                    index_low_reactivity <- which(b/norm >= low_range[1] & b/norm < low_range[2])
                    index_medium_reactivity <- which(b/norm >= medium_range[1] & b/norm < medium_range[2])
                    index_high_reactivity <- which(b/norm >= high_range)
                    #min_reactivity <- c(min_reactivity, 0.001/norm)
                    eval(parse(text= paste( "rf <- per_residue_req$", transcript_names[i], j, "/per_residue_info$", transcript_names[i], j, sep="")))
                    print(index_medium_reactivity)
                    rf_low <- rf[index_low_reactivity]
                    rf_medium <- rf[index_medium_reactivity]
                    rf_high <- rf[index_high_reactivity]
                    
                    sorted_low <- rf_low[order(rf_low)]
                    sorted_medium <- rf_medium[order(rf_medium)]
                    sorted_high <- rf_high[order(rf_high)]
                    
                    #         range_low <- c(round(length(sorted_low) * .9), round(length(sorted_low) * .98))
                    #         range_medium <- c(round(length(sorted_medium) * .9), round(length(sorted_medium) * .98))
                    #         range_high <- c(round(length(sorted_high) * .9), round(length(sorted_high) * .98))
                    
                    req_factor_low <- quantile(sorted_low, criteria)
                    req_factor_medium <- quantile(sorted_medium, criteria)
                    req_factor_high <- quantile(sorted_high, criteria)
                    
                    req_factor <- max(req_factor_low, req_factor_medium, req_factor_high, na.rm=T)
                    req <- c(req, req_factor)
                    #req_factor <- mean(sorted[range[1]:range[2]])
                    #     eval(parse(text= paste("req_factor <- max(per_residue_req$", transcript_names[i], j,
                    #                            "[index_min_reactivity]/per_residue_info$", transcript_names[i], j,"[index_min_reactivity])", sep="")))
                    #     eval(parse(text= paste("req_factor <- per_residue_req$", transcript_names[i], j,
                    #                            "[index_min_reactivity]/per_residue_info$", transcript_names[i], j,"[index_min_reactivity]", sep="")))
                    #print(req_factor)
                    #         req_low <- c(req_low, req_factor_low)
                    #         req_medium <- c(req_medium, req_factor_medium)
                    #         req_high <- c(req_high, req_factor_high)
                    
                    eval(parse(text=paste("index_summary$",transcript_names[i], j, "<- c(req_low = req_factor_low, req_medium =  req_factor_medium,
                    req_high = req_factor_high)", sep="")))
                    #         eval(parse(text=paste("index_summary <- list(transcript=c(index_summary$transcript, paste(", transcript_names[i], j,
                    #                               ", sep=\"\" )), req_low = c(index_summary$req_low, req_factor_low),
                    #                               req_medium = c(index_summary$req_medium, req_factor_medium),
                    #                               req_high = c(index_summary$req_high, req_factor_high))", sep="")))
                    eval(parse(text= paste("demand$", transcript_names[i], j,
                    "_plus <- ceiling(req_factor*supply$", transcript_names[i], j, "_plus)", sep="")))
                    #     eval(parse(text= paste("print(demand$", transcript_names[i], j, "_plus, supply$",
                    #                            transcript_names[i], j, "_plus)", sep="")))
                    eval(parse(text= paste("demand$", transcript_names[i], j,
                    "_minus <- ceiling(size_factor* demand$", transcript_names[i], j, "_plus)", sep="")))
                    rm(g, b, size_factor)
                    
                    
                }
            }
            
            write.xlsx(index_summary, paste("/Users/aviran/Box Sync/N7G_Analysis/Nathan's thing/supply_demand_thing/summary_indices",percent_var,"_", confidence,"_", criteria,".xlsx" ,sep=""))
            
            #2. bootstrap with demand calculated above
            
            combine <- function(...) {
                
                all_lists <- list(...)
                list.combined <- list()
                
                pass_cal <- function(counts) return(c(0,cumsum(counts[1:(length(counts)-1)])))
                
                list.combined$minus.counts <- sapply(all_lists, '[[', 'minus.counts')
                list.combined$minus.pass <- sapply(all_lists, function(x) pass_cal(x$minus.counts))
                #list.combined$minus.pass <- sapply(all_lists, '[[', 'minus.pass')
                list.combined$plus.counts <- sapply(all_lists, '[[', 'plus.counts')
                list.combined$plus.pass <- sapply(all_lists, function(x) pass_cal(x$plus.counts))
                #list.combined$plus.pass <- sapply(all_lists, '[[', 'plus.pass')
                
                return(list.combined)
            }
            
            SampleGenerator <- function (freq,sampleSize) {
                
                n <- length(freq$minus.counts)
                sample_gen <- list(minus.counts="", plus.counts="")
                sample_gen$minus.counts <- tabulate(sample(n, sampleSize$minus.counts, rep= TRUE, prob= freq$minus.counts), nbins=n)
                sample_gen$plus.counts <- tabulate(sample(n, sampleSize$plus.counts, rep= TRUE, prob= freq$plus.counts), nbins=n)
                
                return(sample_gen)
            }
            
            transcript_names <- names(data)
            num_transcripts <- length(transcript_names)
            booty <- list()
            
            for (i in 1:num_transcripts) {
                
                num_biological_replicates <- length(ls(data[[i]], pat= "minus.counts") )
                
                for (j in 1:num_biological_replicates) {
                    
                    reads <- data[[i]]
                    freq <- list(minus.counts="", plus.counts="")
                    sampleSize <- list(minus.counts="", plus.counts="")
                    
                    if (j==1) {
                        freq$minus.counts <- reads$minus.counts/sum(reads$minus.counts)
                        freq$plus.counts <- reads$plus.counts/sum(reads$plus.counts)
                    } else {
                        eval( parse( text= paste("freq$minus.counts <- reads$minus.counts.", j-1,"/sum(reads$minus.counts)",sep="")))
                        eval( parse( text= paste("freq$plus.counts <- reads$plus.counts.", j-1, "/sum(reads$plus.counts)",sep="")))
                    }
                    
                    eval(parse(text=paste("sampleSize$minus.counts <- demand$", transcript_names[i], j,"_minus" , sep="")))
                    eval(parse(text=paste("sampleSize$plus.counts <- demand$", transcript_names[i], j,"_plus" , sep="")))
                    boot_current <- foreach(i=1:100 , .combine = combine, .multicombine=TRUE) %dopar% SampleGenerator(freq, sampleSize)
                    eq <- paste("booty","$",transcript_names[i],"$","rep", j, " <- ", "boot_current",sep="" )
                    eval(parse(text= eq))
                    
                    rm(boot_current)
                }
            }
            
            
            
            
            #3. bootstrap the sample to test and see how many give the output we want
            
            data1 <- booty
            transcript_names <- names(data1)
            num_transcripts <- length(transcript_names)
            mean_confidence<- list(transcript=c(), flag=c(), mean=c())
            demand_performance <- c()
            for (i in 1:num_transcripts) {
                
                num_biological_replicates <- length(data1[[i]])
                #print(num_biological_replicates)
                for (j in 1:num_biological_replicates) {
                    #print(c(i,j))
                    boot <- data1[[i]][[j]]
                    num_boot_replicates <- ncol(boot[[1]])
                    n <- nrow(boot[[1]])
                    
                    estimate.flags=1
                    norm.flags=1
                    log.flags=0
                    zero.flags=1
                    
                    
                    reactivities <- cbind()
                    
                    for (p in 1:num_boot_replicates) {
                        reactivities <- cbind(reactivities, matrix(local.coverage.react.calc(boot$minus.counts[5:n,p], boot$minus.pass[5:n,p], boot$plus.counts[5:n,p],boot$plus.pass[5:n,p],
                        estimate.flags,norm.flags,log.flags,zero.flags), nrow=n-4, ncol=1))
                    }
                    
                    reliable_counts <- c()
                    #real_counts <- 0
                    #     eval(parse(text=paste("b <- betas$",i,"]]",)))
                    for (p in 1:(n-4)) {
                        
                        if (betas[[(3*(i-1))+j]][p]/norm_constants[(3*(i-1))+j] >= min_reactivity) reliable_counts <- c(reliable_counts, length(which(abs(reactivities[p,] - (betas[[(3*(i-1))+j]][p]/norm_constants[(3*(i-1))+j]) ) <=  percent_var*betas[[(3*(i-1))+j]][p]/norm_constants[(3*(i-1))+j] )))
                        #reliable_counts <- c(reliable_counts, length(which(abs(reactivities[p,] - betas[[(3*(i-1))+j]][p]) <=  percent_var*betas[[(3*(i-1))+j]][p] & reactivities[p,] >= min_reactivity)))
                        #if (betas[[(3*(i-1))+j]][p] >= min_reactivity) real_counts <- real_counts+1
                    }
                    print(reliable_counts)
                    #     setEPS()
                    #     postscript(file=paste("/Users/aviran/Box Sync/N7G_Analysis/Nathan's thing/supply_demand_thing/",
                    #                           transcript_names[i],j,"_",estimate.flags,norm.flags,log.flags,zero.flags,".eps" ,sep=""))
                    #     barplot(reliable_counts, names.arg=c(1:(n-4)))
                    #     dev.off()
                    #print(mean_confidence$transcript)
                    mean_confidence <- list(transcript = c(mean_confidence$transcript, paste(transcript_names[i],j,"_",sep="")),
                    flag= c(mean_confidence$flag, paste("_",estimate.flags, norm.flags, log.flags,zero.flags, sep="")),
                    mean=c(mean_confidence$mean, mean(reliable_counts))) #, out_of = c(mean_confidence$out_of, real_counts) )
                    
                    
                    #transcript <- c(paste(transcript_names[i],j, sep=""), transcript)
                }
                
            }
            print(mean_confidence)
            demand_performance <- mean_confidence$mean
            write.xlsx(mean_confidence, paste("/Users/aviran/Box Sync/N7G_Analysis/Nathan's thing/supply_demand_thing/summary_demand",percent_var, confidence,criteria, ".xlsx" ,sep="_"))
            
            #4. Compare with the mean_confidence from supply
            transcript_names <- names(data)
            num_transcripts <- length(transcript_names)
            booty_supply <- list()
            
            for (i in 1:num_transcripts) {
                
                num_biological_replicates <- length(ls(data[[i]], pat= "minus.counts") )
                
                for (j in 1:num_biological_replicates) {
                    
                    reads <- data[[i]]
                    freq <- list(minus.counts="", plus.counts="")
                    sampleSize <- list(minus.counts="", plus.counts="")
                    
                    if (j==1) {
                        freq$minus.counts <- reads$minus.counts/sum(reads$minus.counts)
                        freq$plus.counts <- reads$plus.counts/sum(reads$plus.counts)
                        sampleSize$minus.counts <- sum(reads$minus.counts)
                        sampleSize$plus.counts <- sum(reads$plus.counts)
                        
                    } else {
                        eval( parse( text= paste("freq$minus.counts <- reads$minus.counts.", j-1,"/sum(reads$minus.counts)",sep="")))
                        eval( parse( text= paste("freq$plus.counts <- reads$plus.counts.", j-1, "/sum(reads$plus.counts)",sep="")))
                        eval( parse( text=paste("sampleSize$minus.counts <- sum(reads$minus.counts.", j-1,")", sep="")))
                        eval( parse( text=paste("sampleSize$plus.counts <- sum(reads$plus.counts.", j-1,")", sep="")))
                        
                    }
                    
                    #     eval(parse(text=paste("sampleSize$minus.counts <- demand$", transcript_names[i], j,"_minus" , sep="")))
                    #     eval(parse(text=paste("sampleSize$plus.counts <- demand$", transcript_names[i], j,"_plus" , sep="")))
                    boot_current <- foreach(i=1:100 , .combine = combine, .multicombine=TRUE) %dopar% SampleGenerator(freq, sampleSize)
                    eq <- paste("booty_supply","$",transcript_names[i],"$","rep", j, " <- ", "boot_current",sep="" )
                    eval(parse(text= eq))
                    
                    rm(boot_current)
                }
            }
            
            #5. bootstrap the sample from supply to test and see how many give the output we want
            
            data2 <- booty_supply
            transcript_names <- names(data2)
            num_transcripts <- length(transcript_names)
            mean_confidence<- list(transcript=c(), flag=c(), mean=c())
            
            for (i in 1:num_transcripts) {
                
                num_biological_replicates <- length(data2[[i]])
                #print(num_biological_replicates)
                for (j in 1:num_biological_replicates) {
                    #print(c(i,j))
                    boot <- data2[[i]][[j]]
                    num_boot_replicates <- ncol(boot[[1]])
                    n <- nrow(boot[[1]])
                    
                    estimate.flags=1
                    norm.flags=1
                    log.flags=0
                    zero.flags=1
                    
                    
                    reactivities <- cbind()
                    
                    for (p in 1:num_boot_replicates) {
                        reactivities <- cbind(reactivities, matrix(local.coverage.react.calc(boot$minus.counts[5:n,p], boot$minus.pass[5:n,p], boot$plus.counts[5:n,p],boot$plus.pass[5:n,p],
                        estimate.flags,norm.flags,log.flags,zero.flags), nrow=n-4, ncol=1))
                    }
                    
                    reliable_counts <- c()
                    #real_counts <- 0
                    #     eval(parse(text=paste("b <- betas$",i,"]]",)))
                    for (p in 1:(n-4)) {
                        if (betas[[(3*(i-1))+j]][p]/norm_constants[(3*(i-1))+j]  >= min_reactivity) reliable_counts <- c(reliable_counts, length(which(abs(reactivities[p,] - (betas[[(3*(i-1))+j]][p]/norm_constants[(3*(i-1))+j]) ) <=  percent_var*betas[[(3*(i-1))+j]][p]/norm_constants[(3*(i-1))+j] )))
                        #reliable_counts <- c(reliable_counts, length(which(abs(reactivities[p,] - betas[[(3*(i-1))+j]][p]) <=  percent_var*betas[[(3*(i-1))+j]][p] & reactivities[p,] >= min_reactivity)))
                        #if (betas[[(3*(i-1))+j]][p] >= min_reactivity) real_counts <- real_counts+1
                    }
                    
                    #     setEPS()
                    #     postscript(file=paste("/Users/aviran/Box Sync/N7G_Analysis/Nathan's thing/supply_demand_thing/",
                    #                           transcript_names[i],j,"_",estimate.flags,norm.flags,log.flags,zero.flags,".eps" ,sep=""))
                    #     barplot(reliable_counts, names.arg=c(1:(n-4)))
                    #     dev.off()
                    #print(mean_confidence$transcript)
                    mean_confidence <- list(transcript = c(mean_confidence$transcript, paste(transcript_names[i],j,"_",sep="")),
                    flag= c(mean_confidence$flag, paste("_",estimate.flags, norm.flags, log.flags,zero.flags, sep="")),
                    mean=c(mean_confidence$mean, mean(reliable_counts))) #, out_of= c(mean_confidence$out_of, real_counts))
                    
                    
                    #transcript <- c(paste(transcript_names[i],j, sep=""), transcript)
                }
                
            }
            supply_performance <- mean_confidence$mean
            write.xlsx(mean_confidence, paste("/Users/aviran/Box Sync/N7G_Analysis/Nathan's thing/supply_demand_thing/summary_supply",min_reactivity,percent_var, confidence,criteria,".xlsx" ,sep="_"))
            
            #find the accuracy of predictions
            good_predictions <- 0
            for (num in 1:24) {
                if (round(demand_performance[num])/100 >= confidence & req[num]<=1) good_predictions <- good_predictions + 1
                if (round(demand_performance[num])/100 >= confidence & (req[num] >1 & round(supply_performance[num])/100 < confidence) ) {
                    good_predictions <- good_predictions + 1
                }
            }
            #print(good_predictions/24)
            accuracy <- list(confidence= c(accuracy$confidence, confidence),
            percent_var= c(accuracy$percent_var, percent_var), good_predictions= c(accuracy$good_predictions, good_predictions),
            index_criteria= c(accuracy$index_criteria, criteria))
        }
    }
}
#}

write.xlsx(accuracy, "/Users/aviran/Box Sync/N7G_Analysis/Nathan's thing/supply_demand_thing/summary_accuracy.xlsx")
#------------------------------------------------------------------


