library(stringr)
library(mgcv)
library(GA)
library(ggplot2)

#-------------------INPUT loci of interest for reference ----------------------

loci_of_interest <- readRDS("loci_of_interest.RDS")


#-----------------FORMATTING function------------------

formating_ampCov <- function(amplicon_coverage, loci_of_interest) {
  
  amplicon_coverage <- amplicon_coverage[!grepl("neg", amplicon_coverage$SampleID, ignore.case = TRUE), ] #remove neg controls
  
  # Create "amplicon_length" column
  numbas <- str_extract(amplicon_coverage$Locus, "(?<=_v[1234]-).*?(?=-[^-]*$)")
  numbas <- data.frame(do.call("rbind", str_split(numbas, "-")))
  numbas[] <- lapply(numbas, as.numeric)
  amplicon_coverage$amplicon_length <- abs(numbas$X1 - numbas$X2)
  
  # Create "amplicon.group" column
  last_dash <- gsub(".*-(.*)$", "\\1", amplicon_coverage$Locus)
  numbas2 <- gsub("[^0-9]", "", last_dash)
  numbas2 <- gsub("12", "1", numbas2)
  amplicon_coverage$amplicon.group <- numbas2
  
  # Create columns for each loci of interest filled with 0's
  loci_names <- names(loci_of_interest)
  
  for (locus in loci_names) {
    amplicon_coverage[[locus]] <- 0
    amplicon_coverage[[locus]][amplicon_coverage$Locus %in% loci_of_interest[[locus]]] <- 1
  }
  
  return(amplicon_coverage)
}

#-------------------------------estCNV function--------------------------

estCNV <- function(data, sample.name="<sample>", plot.gam=F, verbose=F, k.gam=4,
                   loci.of.interest=c("HRP2","HRP3","MDR1","MDR2","PM1","PM2","PM3","PM4"), 
                   amp.to.remove=c(""), # amp.to.remove=c("Pf3D7_05_v3-960810-961048-1B", "Pf3D7_03_v3-85507-85756-1A")
                   threshold = 50){
  
  if(median(data[,"OutputPostprocessing"] < threshold)){ ##MODIFY AS DESIRED (2500 SEEMS GOOD AFTER NORMALIZATION OF SAMPLES ACROSS miseq/nextseq RUNS)
    if(verbose){
      print(paste("Insufficient reads for sample",sample.name))
    }
    return(NA)
  } else {
    data <- data[data[,"amplicon_length"] < 275,]
    data <- data[!(data[,"Locus"] %in% amp.to.remove),]
    data$amplicon.group <- as.factor(data$amplicon.group)
    if(length(unique(data$amplicon.group))>1){
      form <- as.formula(paste("OutputPostprocessing ~ s(amplicon_length, k=k.gam, by=(amplicon.group)) + amplicon.group +",
                               paste(loci.of.interest, collapse=" + ")) ) 
    } else {
      form <- as.formula(paste("OutputPostprocessing ~ s(amplicon_length, k=k.gam) +",
                               paste(loci.of.interest, collapse=" + ")) ) 
    }
    
    fit <- gam(form, data=data, family="poisson")
    
    data$locus.of.interest <- apply(data[,loci.of.interest], 1, function(x) sum(x)>0)
    lengths <- min(data$amplicon_length):max(data$amplicon_length)  
    pred.frame <- cbind(data.frame(amplicon_length=lengths),  
                        data.frame(matrix(data=0, nrow=length(lengths), ncol=length(loci.of.interest))))
    names(pred.frame) <- c("amplicon_length",loci.of.interest)
    pred.frame.1 <- cbind(pred.frame, data.frame(amplicon.group=unique(data$amplicon.group)[1]))  
    pred.1 <- predict(fit, newdata = pred.frame.1, type="response")
    
    if(length(unique(data$amplicon.group))>1){
      pred.frame.2 <- cbind(pred.frame, data.frame(amplicon.group=unique(data$amplicon.group)[2]))  
      pred.2 <- predict(fit, newdata = pred.frame.2, type="response")
    }
    
    amp.of.interest <- data[data$locus.of.interest==1,c("Locus","OutputPostprocessing","amplicon_length","amplicon.group")]
    # amp.of.interest <- data[,c("Locus","OutputPostprocessing","amplicon_length","amplicon.group")] #RETURNS ALL AMPLICONS
    amp.of.interest <- cbind(amp.of.interest, pred.frame[1:nrow(amp.of.interest), 2:ncol(pred.frame)])
    amp.of.interest$expected.reads <- predict(fit, newdata = amp.of.interest, type="response")
    fc.amplicon <- amp.of.interest$OutputPostprocessing / amp.of.interest$expected.reads
    names(fc.amplicon) <- amp.of.interest$Locus
    
    coef <- summary(fit)$p.coef[loci.of.interest]
    se <- summary(fit)$se[loci.of.interest]
    
    if(plot.gam){
      plot(data$amplicon_length, data$OutputPostprocessing+0.5, log = "y", xlab="Amplicon length (bp)", ylab="Reads", main=sample.name,
           col=as.numeric(data$amplicon.group), pch=1+(data$locus.of.interest*15), las=1)
      points(lengths, pred.1, type="l", lwd=2)
      if(length(unique(data$amplicon.group))>1){
        points(lengths, pred.2, type="l", lwd=2, col=2)
      }
    }
    return(list(fc.locus=exp(coef), coef=coef, se=se, fc.amplicon=fc.amplicon))
  }
}


### ---------------EXPECTED FOLD CHANGE FOR EACH CONTROL SAMPLE FROM THE DATASET --------------------
expected_foldchanges <- read.csv("expected_foldchanges.csv")

expected_foldchanges_filepaths <-  paste0("CNV_runs_sample_coverage/", expected_foldchanges$filename)

expected_foldchanges$filepaths <- expected_foldchanges_filepaths
  


### INIT LOOP
#------------------------------FORMATING------------------------------------

results_list <- list()

# Loop through each file
for (i in seq_along(expected_foldchanges_filepaths)) {
  # Read amplicon_coverage file
  filepath <- expected_foldchanges_filepaths[i]
  filename <- basename(filepath)
  
  sample_name_ <- expected_foldchanges$control_name[i]
  iteration_name <- paste0(filepath, "___", sample_name_)
  
  #this is the amplicon table that will be referenced by the GA. 1 = used, 0 = left out
  #amp_table <- data.frame(amplicons = unique(amplicon_coverage_formatted$Locus), used = 1)
  unique_amplicons <- unique(amplicon_coverage_formatted$Locus)
  
  
  # Check if the file is already processed
  if (iteration_name %in% names(results_list)) {
    cat(paste("\nSkipping control:", iteration_name, "as it's already processed\n"))
    next  # Skip to the next iteration
  }
  
  amplicon_coverage <- read.table(filepath, header = TRUE)
  
  # Format input for estCNV
  amplicon_coverage_formatted <- formating_ampCov(amplicon_coverage = amplicon_coverage, loci_of_interest = loci_of_interest)
  
  #---------------------------FOLD CHANGE ESTIMATION--------------------------------
  
  ###############C GENETIC ALGO ####################3
  
  # SUBSET AMPLICONS FUNCTION
  subsetAmplicons <- function(amplicon_indices, all_amplicons) {
    selected_amplicons <- all_amplicons[amplicon_indices == 1]
    return(selected_amplicons)
  }
  
  # FITNESS FUNCTION
  fitness_function <- function(amplicon_indices, sample_name = sample_name_) {
    selected_amplicons <- subsetAmplicons(amplicon_indices, unique_amplicons)
    
    # Ensure at least one amplicon from each locus of interest is selected
    for (locus in names(loci_of_interest)) {
      locus_amplicons <- loci_of_interest[[locus]]
      if (all(!(locus_amplicons %in% selected_amplicons))) {
        # If none of the locus amplicons are selected, select one randomly
        selected_amplicons <- c(selected_amplicons, sample(locus_amplicons, 1))
      }
    }
    
    # Subset sample data
    sample_subset <- amplicon_coverage_formatted[amplicon_coverage_formatted$SampleID == sample_name, ]
    
    # Exclude amplicons not in the selected set
    excluded_amplicons <- setdiff(unique_amplicons, selected_amplicons)
    sample_subset <- sample_subset[!(sample_subset$Locus %in% excluded_amplicons),]
    
    # Check if sample_subset is empty
    if (nrow(sample_subset) == 0) {
      print("No rows left after amplicon subset. Returning Inf RMSE.")
      return(Inf)
    }
    
    # Run estCNV
    result_CNV <- estCNV(data = sample_subset, sample.name = sample_name, plot.gam = F)
    
    # Extract and format the $fc.locus elements corresponding to the loci from expected_foldchanges_loci
    expected_foldchanges_loci <- expected_foldchanges[expected_foldchanges$control_name == sample_name, ][3:4]
    
    fc <- as.data.frame(result_CNV$fc.locus[expected_foldchanges_loci$locus])
    observed_foldchanges_loci <- data.frame(loci = rownames(fc), observed_foldchange = fc[,1], row.names=NULL)
    
    # Calculate RMSE
    rmse <- sqrt(mean((expected_foldchanges_loci$expected_foldchange - observed_foldchanges_loci$observed_foldchange)^2))
    #print(rmse)
    
    # Return the inverse of RMSE as fitness
    return(1/rmse)
  }
  
  # PLOT FUNCTION (to be used by GA)
  plot_fitness <- function(obj) {
    plot(obj)
  }
  
  # Define GA parameters
  pop_size <- 100
  generations <- 40
  mutation_prob <- 0.1
  crossover_prob <- 0.8
  elitism <- 10
  chrom_length <- length(unique_amplicons)
  
  # Run GA with real-time plotting
  ga_result <- ga(type = "binary", fitness = fitness_function, nBits = chrom_length,
                  popSize = pop_size, maxiter = generations, pmutation = mutation_prob,
                  pcrossover = crossover_prob, elitism = elitism, keepBest = TRUE,
                  run = 100, monitor = plot_fitness, seed = 420)
  
  # SOLUTION FROM GA
  used_amplicons <- as.numeric(ga_result@solution[1,]) #amplicons
  best_solution  <- as.data.frame(cbind(amplicons = unique_amplicons, used_amplicons = used_amplicons))
  best_solution$used_amplicons <- as.numeric(best_solution$used_amplicons)
  
  #print(best_solution)
  
  n_amplicons <- sum(best_solution$used_amplicons)
  best_fitness <- max(ga_result@fitness)**-1
  
  cat("\n")
  print(iteration_name)
  print(paste("Optimal # of amplicons =", n_amplicons))
  print(paste("Lowest RMSE =", round(best_fitness, 5)))
  cat("\n")
  
  # Check if at least one amplicon from each locus of interest was used
  used_amplicons <- best_solution$amplicons[best_solution$used_amplicons == 1]
  loci_used <- sapply(loci_of_interest, function(locus_amplicons) any(locus_amplicons %in% used_amplicons))
  
  if (sum(loci_used) == length(names(loci_used))){
    print("At least 1 amplicon of each loci of interest was used")
  }else{
    print(paste("No", names(loci_used[loci_used == FALSE]), "amplicons were used"))
  }
  
  results_list[[iteration_name]] <- best_solution

}


# Merge resulting dataframes
for (i in seq_along(results_list)) {
  col_name <- paste0("used_amplicons_", i)
  results_list[[i]] <- setNames(results_list[[i]], c("amplicon", col_name))
}

merged_GA_result <- results_list[[1]]
for (i in 2:length(results_list)) {
  merged_GA_result <- merge(merged_GA_result, results_list[[i]], by = "amplicon", all = TRUE)
}

colnames(merged_GA_result)[-1] <- basename(names(results_list))

#checkpoint
#write.csv(merged_GA_result, "merged_GA_result.csv") ################################

#---------------------------- ANALYZE RESULTS ---------------------------

# 1) pca

# Exclude the first column (amplicon names) for PCA
pca_data <- merged_GA_result[, -1]
pca_data  <- t(pca_data)

# Perform Principal Component Analysis
pca_result <- prcomp(pca_data, scale. = TRUE)

# Access the results
summary(pca_result)

pca_df <- as.data.frame(pca_result$x)

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 4) +
  labs(title = "PCA of merged_GA_result",
       x = "Principal Component 1",
       y = "Principal Component 2")



# 2) barplot of percentage of controls that used each amplicon

#percentage of controls that used each amplicon
amplicon_results <- data.frame(amplicon = merged_GA_result$amplicon)
amplicon_results$percentage_used <- rowSums(merged_GA_result[, -1])/ length(merged_GA_result[, -1])
amplicon_results$loci <- NA  # Initialize the column with NAs

for (f in 1:length(loci_of_interest)) {
  matching_amplicons <- amplicon_results$amplicon %in% loci_of_interest[[f]]
  amplicon_results$loci[matching_amplicons] <- names(loci_of_interest)[f]
}

amplicon_results <- amplicon_results[order(-amplicon_results$percentage_used), ]

# Create a sorted barplot using ggplot2
ggplot(amplicon_results, aes(x = reorder(amplicon, -percentage_used), y = percentage_used, fill = loci)) +
  geom_bar(stat = "identity") +
  labs(title = "Percentage of Amplicons Used by GA",
       x = "Amplicon",
       y = "Percentage Used") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5)) +
  # Add orange color to bars corresponding to amplicons in all_loci_amplicons
  scale_fill_manual(values = c("red", "blue", "green", "purple", "yellow", "pink", "orange", "black"))


#testing common amplicons between 2 runs:

# used_amplicons2 <- as.numeric(ga_result@solution[1,]) #amplicons
# best_solution2  <- as.data.frame(cbind(amplicons = unique_amplicons, used_amplicons = used_amplicons2))
# best_solution2$used_amplicons <- as.numeric(best_solution2$used_amplicons)
# 
# print(best_solution2)
# 
# n_amplicons <- sum(best_solution2$used_amplicons)
# best_fitness <- max(ga_result@fitness)**-1
# 
# print(paste("Optimal # of amplicons =", n_amplicons))
# print(paste("Lowest RMSE =", round(best_fitness, 5)))
# 


# best_solution$used_amplicons
# best_solution2$used_amplicons

# common<-best_solution[best_solution$used_amplicons + best_solution2$used_amplicons == 2,]
# sum(common$used_amplicons)



FOLDCHANGE <- data.frame(rep(NA, 8)) #number of rows = amount of loci of interest

for (name in unique(amplicon_coverage_formatted$SampleID)){
  tryCatch({
    result <- estCNV(amplicon_coverage_formatted[amplicon_coverage_formatted$SampleID == name, ], plot.gam = FALSE, sample.name = name)
    tmp_df <- do.call(cbind.data.frame, result[1]) 
    colnames(tmp_df) <- name
    FOLDCHANGE <- cbind(FOLDCHANGE, tmp_df)
    
  }, error = function(e) {
    #message("Error: ", e$message)
  })
}

FOLDCHANGE_final<-FOLDCHANGE[,-1]
FOLDCHANGE_final<-t(FOLDCHANGE_final)
FOLDCHANGE_final<-as.data.frame(FOLDCHANGE_final)
