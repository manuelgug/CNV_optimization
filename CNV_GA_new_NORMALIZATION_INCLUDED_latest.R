library(stringr)
library(mgcv)
library(GA)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(Rtsne)
library(corrplot)
library(factoextra)

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
  
  amplicon_coverage <- read.table(filepath, header = TRUE)
  
  # Format input for estCNV
  amplicon_coverage_formatted <- formating_ampCov(amplicon_coverage = amplicon_coverage, loci_of_interest = loci_of_interest)
  
  #this is the amplicon table that will be referenced by the GA. 1 = used, 0 = left out
  #amp_table <- data.frame(amplicons = unique(amplicon_coverage_formatted$Locus), used = 1)
  unique_amplicons <- unique(amplicon_coverage_formatted$Locus)
  
  
  # Check if the file is already processed
  if (iteration_name %in% names(results_list)) {
    cat(paste("\nSkipping control:", iteration_name, "as it's already processed\n"))
    next  # Skip to the next iteration
  }

  
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
    
    # Expected loci fc
    expected_foldchanges_loci <- expected_foldchanges[expected_foldchanges$control_name == sample_name, ][3:4]
    
    ###############################################################################
    #Subset single-copy control data for normalization
    controls <- amplicon_coverage_formatted[!grepl("(?i)Dd2|PM|HB3", amplicon_coverage_formatted$SampleID) & grepl("(?i)3D7", amplicon_coverage_formatted$SampleID), ]
    
    # Exclude amplicons not in the control set
    excluded_amplicons <- setdiff(unique_amplicons, selected_amplicons)
    controls <- controls[!(controls$Locus %in% excluded_amplicons),]
    
    #estCNV for controls (calculate normalization factor)
    unique_controls <- unique(controls$SampleID)
    fc_controls <- data.frame(nrow =1)
    
    # Loop estCNV through unique_controls
    for (control in unique_controls){
      
      controls_CNV <- estCNV(data = amplicon_coverage_formatted[amplicon_coverage_formatted$SampleID == control,], sample.name = control, plot.gam = FALSE)
      fc_controls_ <- as.data.frame(controls_CNV$fc.locus[expected_foldchanges_loci$locus])
      
      # Append new columns to fc_controls
      fc_controls <- cbind(fc_controls, fc_controls_)
    }
    
    fc_controls
    fc_controls <- fc_controls[,-1]
    
    #mean fc of single-copy controls
    fc_controls_norm_factor <- as.data.frame(rowMeans(fc_controls))
    colnames(fc_controls_norm_factor) <- "norm_factor"
    
    fc_controls_norm_factor[fc_controls_norm_factor < 0.001] <- 0
    ###############################################################################
    
    # Ensure at least one amplicon from each locus of interest is selected
    for (locus in names(loci_of_interest)) {
      locus_amplicons <- loci_of_interest[[locus]]
      if (all(!(locus_amplicons %in% selected_amplicons))) {
        # If none of the locus amplicons are selected, select one randomly
        selected_amplicons <- c(selected_amplicons, sample(locus_amplicons, 1))
      }
    }
    
    # Exclude amplicons not in the selected set
    excluded_amplicons <- setdiff(unique_amplicons, selected_amplicons)
    amplicon_coverage_formatted <- amplicon_coverage_formatted[!(amplicon_coverage_formatted$Locus %in% excluded_amplicons),]
    
    # Subset sample data
    sample_subset <- amplicon_coverage_formatted[amplicon_coverage_formatted$SampleID == sample_name, ]
    
    # Check if sample_subset is empty
    if (nrow(sample_subset) == 0) {
      print("No rows left after amplicon subset. Returning Inf RMSE.")
      return(Inf)
    }
    
    # Run estCNV on CNV sample
    result_CNV <- estCNV(data = sample_subset[sample_subset$SampleID == sample_name,], sample.name = sample_name, plot.gam = F)
    
    fc <- as.data.frame(result_CNV$fc.locus[expected_foldchanges_loci$locus])
    observed_foldchanges_loci <- data.frame(loci = rownames(fc), observed_foldchange = fc[,1], row.names=NULL)
    
   observed_foldchanges_loci <- as.data.frame(observed_foldchanges_loci)
   rownames(observed_foldchanges_loci) <- observed_foldchanges_loci[,1]
   observed_foldchanges_loci <- observed_foldchanges_loci[-1]
    

    
    #apply normalization factor
    observed_foldchanges_loci <-  as.data.frame(observed_foldchanges_loci$observed_foldchange/fc_controls_norm_factor$norm_factor)

    colnames(observed_foldchanges_loci) <- "observed_foldchange"
    rownames(observed_foldchanges_loci) <- rownames(fc_controls_norm_factor)
    
    #Change Inf values to 0. Infs appear as a result of 0/0 division during normalization
    #observed_foldchanges_loci[observed_foldchanges_loci$observed_foldchange == Inf] <- 0

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
  mutation_prob <- 0.2
  crossover_prob <- 0.8
  elitism <- 10
  chrom_length <- length(unique_amplicons)
  
  # Run GA with real-time plotting
  ga_result <- ga(type = "binary", fitness = fitness_function, nBits = chrom_length,
                  popSize = pop_size, maxiter = generations, pmutation = mutation_prob,
                  pcrossover = crossover_prob, elitism = elitism, keepBest = TRUE,
                  run = 100, monitor = plot_fitness, seed = 420, parallel = 10)
  
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
  print(paste("Lowest RMSE =", best_fitness))
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
write.csv(merged_GA_result, "merged_GA_result_NORMALIZED.csv") ################################



#---------------------------- ANALYZE RESULTS ---------------------------

##################inputs##################
merged_GA_result <- read.csv("merged_GA_result.csv", row.names = 1) 
loci_of_interest <- readRDS("loci_of_interest.RDS")
##########################################

## data formatting

#percentage of controls that used each amplicon
amplicon_results <- data.frame(amplicon = merged_GA_result$amplicon)
amplicon_results$percentage_used <- rowSums(merged_GA_result[, -1])/ length(merged_GA_result[, -1])

amplicon_results$loci <- NA  # Initialize the column with NAs

for (f in 1:length(loci_of_interest)) {
  matching_amplicons <- amplicon_results$amplicon %in% loci_of_interest[[f]]
  amplicon_results$loci[matching_amplicons] <- names(loci_of_interest)[f]
}

melted_merged_GA_result <- melt(merged_GA_result)
melted_merged_GA_result <- separate(melted_merged_GA_result, variable, into = c("run", "control"), sep = "___")

unique_variables <- length(unique(melted_merged_GA_result$control))


##### exploration

# percentage of controls that used each amplicon in their optimal solution, also check for the loci of interest

ggplot(amplicon_results, aes(x = reorder(amplicon, -percentage_used), y = percentage_used, fill = loci)) +
  geom_bar(stat = "identity") +
  labs(title = "Amplicons Used by GA for optimal Fold Change Calculation on Controls",
       x = "Amplicons",
       y = "Percentage of Controls") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5)) +
  # Add orange color to bars corresponding to amplicons in all_loci_amplicons
  scale_fill_manual(values = c("red", "green", "blue", "orange", "purple", "yellow", "black", "violet"))

#are the runs of the controls a factor that influences the use of amplicons by the GA?

amps_used <- melted_merged_GA_result %>% 
  group_by(run,control) %>%
  summarize(amps_used = sum(value))

ggplot(amps_used, aes(x = reorder(control, -amps_used), y = amps_used, fill = run )) +
  geom_bar(stat = "identity")+
  labs(title = "",
       x = "Control strain",
       y = "Number of Amplicons Used") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5))+
  guides(fill = guide_legend(ncol = 1))


#### DO AMPLICONS CO-OCURR IN MOST OPTIMAL SOLUTIONS? (aka are there "GOOD" and "BAD" amplicons?). IS THERE A GENERALIZABLE SET OF AMPLICONS FOR CNV CALCULATION?

# 1) kmeans
optimal_k <- 2:5
plot_list <- list()
clusters <- data.frame(matrix(nrow = dim(merged_GA_result)[1], ncol = 0)) 

# Perform k-means clustering and plot for each optimal k
for (k in optimal_k) {
  # Perform k-means clustering with the optimal k
  set.seed(69)
  kmeans_result <- kmeans(merged_GA_result[, -1], centers = k, nstart = 10000, iter.max = 10000)
  
  # Add cluster assignment as a new column to clusters
  clusters[[paste0("cluster_", k)]] <- kmeans_result$cluster
  
}

# 2) pca
color_data <- data.frame(
  loci = sort(na.omit(unique(amplicon_results$loci))),
  color = c("red", "green", "blue", "orange", "purple", "yellow", "black", "violet")
)

pca_data <- merged_GA_result[, -1]

pca_result <- prcomp(pca_data, scale. = TRUE)

#add metadata
pca_df <- as.data.frame(cbind(pca_result$x, amplicon_results))
pca_df <- merge(pca_df, color_data, by = "loci", all.x = TRUE)
clusters$amplicon <- merged_GA_result$amplicon
pca_df <- merge(pca_df, clusters, by =c("amplicon"))

variance_explained <- summary(pca_result)$importance["Proportion of Variance", ]

ggplot(pca_df, aes(PC1, PC2, color = ifelse(!is.na(loci), loci, NA), fill = percentage_used, shape = factor(cluster_3))) +
  geom_point(size = 6, alpha = ifelse(!is.na(pca_df$loci), 1, 0.4), stroke = 1.5) +
  labs(title = "PCA of co-occurrence of amplicons in GA solutions + Clustering",
       x = paste("PC1 (", round(variance_explained[1] * 100, 2), "%)", sep = ""),
       y = paste("PC2 (", round(variance_explained[2] * 100, 2), "%)", sep = "")) +
  scale_fill_gradient(low = "black", high = "cyan") +
  scale_shape_manual(values = c(21, 22, 23)) +  # Specify shapes 21 to 23
  scale_color_manual(values = setNames(color_data$color, color_data$loci), na.value = "white") +
  theme_minimal()


OPTIMIZED_SET_OF_AMPLICONS <- pca_df[pca_df$cluster_3 == 1,]$amplicon

print(paste0("OPTIMIZED_SET_OF_AMPLICONS = ", length(OPTIMIZED_SET_OF_AMPLICONS)))

#RESULTS:
# Amplicons used in more controls are also grouped together more often, so there's such thing as good/bad amplicons for CNV calculation
# Cluster 1 of k = 3 includes the most used amplicons, and also includes at least 1 amplicon of each locus of interest; k = 2 is too broad and k = 4 doesn't include all loci of interest in the high abundance cluster
# Cluster 1 of K = 3 it's a good candidate as an optimized generalizable amplicon set for CNV estimation


### BENCHMARKING OPTIMIZED AMPLICON SET: CROSS VALIDATION

# 1) format inputs for estCNV

amp_cov_inputs <- unique(expected_foldchanges_filepaths)

estcnv_inputs <- list()

for (i in seq_along(amp_cov_inputs)) {
  # Read amplicon_coverage file
  filepath <- amp_cov_inputs[i]
  filename <- basename(filepath)
  
  sample_name_ <- expected_foldchanges$control_name[i]
  #iteration_name <- paste0(filepath, "___", sample_name_)
  
  amplicon_coverage <- read.table(filepath, header = TRUE)
  
  # Format input for estCNV
  amplicon_coverage_formatted <- formating_ampCov(amplicon_coverage = amplicon_coverage, loci_of_interest = loci_of_interest)
  
  estcnv_inputs[[filepath]] <- amplicon_coverage_formatted
}


# 2) get optimized inputs for estCNV
subset_by_locus <- function(df) {
  df_subset <- subset(df, Locus %in% OPTIMIZED_SET_OF_AMPLICONS)
  return(df_subset)
}

estcnv_inputs_OPTIMIZED <- lapply(estcnv_inputs, subset_by_locus)



# 3) function tu run estCNV on input lists
run_estCNV <- function(INPUT_LIST) {
  RESULTS_LIST <- list()  # Initialize empty list to store results
  FOLDCHANGE <- data.frame(rep(NA, 8))  # Initialize FOLDCHANGE dataframe outside the loop
  
  # Loop through each element in INPUT_LIST
  for (j in seq_along(INPUT_LIST)) {
    amplicon_coverage_formatted <- INPUT_LIST[[j]]  # Get dataframe
    
    # Initialize FOLDCHANGE dataframe
    FOLDCHANGE <- data.frame(rep(NA, 8))
    
    # Loop through unique SampleIDs
    for (name in unique(amplicon_coverage_formatted$SampleID)) {
      tryCatch({
        result <- estCNV(amplicon_coverage_formatted[amplicon_coverage_formatted$SampleID == name, ], plot.gam = FALSE, sample.name = name)
        tmp_df <- do.call(cbind.data.frame, result[1]) 
        colnames(tmp_df) <- name
        FOLDCHANGE <- cbind(FOLDCHANGE, tmp_df)
      }, error = function(e) {
        # Handle errors if necessary
        #message("Error: ", e$message)
      })
    }
    
    # Process FOLDCHANGE_final
    FOLDCHANGE_final <- FOLDCHANGE[, -1]
    FOLDCHANGE_final <- t(FOLDCHANGE_final)
    FOLDCHANGE_final <- as.data.frame(FOLDCHANGE_final)
    
    
    ##########################  NORMALIZE WITH SINGLE COPY CONTROLS ###########################
    
    controls <- FOLDCHANGE_final[!grepl("(?i)Dd2|PM|HB3", rownames(FOLDCHANGE_final)) & grepl("(?i)3D7", rownames(FOLDCHANGE_final)), ]
    control_means <- colMeans(controls) #MEAN FOLD CHANGE OF CONTROLS
    
    #Round everything belos 0.001 to 0 both in control means and samples fold change to avoid erroneous normalization 
    control_means[control_means < 0.001] <- 0
    control_means[control_means == 1] <- 0 
    
    # create an empty data frame for normalized data
    FOLDCHANGE_final_NORMALIZED <- data.frame(matrix(0, nrow = nrow(FOLDCHANGE_final), ncol = ncol(FOLDCHANGE_final)))
    colnames(FOLDCHANGE_final_NORMALIZED) <- colnames(FOLDCHANGE_final)
    
    # normalize FOLDCHANGE_final by control_means
    for (i in 1:ncol(FOLDCHANGE_final)) {
      FOLDCHANGE_final_NORMALIZED[, i] <- FOLDCHANGE_final[, i] / control_means[i]
    }
    
    #Change Inf values to 0. Infs appear as a result of 0/0 division during normalization
    FOLDCHANGE_final_NORMALIZED[FOLDCHANGE_final_NORMALIZED == Inf] <- 0
    
    rownames(FOLDCHANGE_final_NORMALIZED) <- rownames(FOLDCHANGE_final)
    
    ########################################################################################### 
    
    # Append FOLDCHANGE_final to RESULTS_LIST list
    RESULTS_LIST[[j]] <- FOLDCHANGE_final_NORMALIZED
  }
  
  return(RESULTS_LIST)
}


# 4) run estCNV with all amplicons (estcnv_inputs)
estcsv_inputs_RESULTS <- run_estCNV(estcnv_inputs)

lapply(estcsv_inputs_RESULTS[[3]], median)


# 5) run estCNV with OPTIMIZED_SET_OF_AMPLICONS (estcnv_inputs_OPTIMIZED)
estcsv_inputs_RESULTS_OPTIMIZED <- run_estCNV(estcnv_inputs_OPTIMIZED)

lapply(estcsv_inputs_RESULTS_OPTIMIZED[[3]], median)



### CROSS VALIDATION OF SUBSET OF AMPLICONS

# 0) develop a way of benchmarking results of estCNV against expected results ( is it on fitness function already? i think so )
# 1) subset amplicons of cluster 1 of k = 3
# 2) use subset to calculate fold change in all controls with known genotype
# 3) evaluate agains using all amplicons. is it better?



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
