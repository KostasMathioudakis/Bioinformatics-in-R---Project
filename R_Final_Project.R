###
### Kostas Mathioudakis CSD3982
### Michail Toutoudakis CSD4054
###

### Libraries
library(gplots)

# File paths
# We assume that the files are in the same folder as the exercise wanted
param_path = "pars_final.txt"
obs_path = "ms_obs_final.out"
sim_path = "ms_sim_final.out"


# Function for reading datasets
read_datasets <- function(file_path) {
  lines <- readLines(file_path)
  datasets <- list()
  dataset <- NULL
  
  for (line in lines) {
    if (line != "") {
      vector <- as.numeric(strsplit(line, split = "")[[1]])
      if (is.null(dataset)) {
        dataset <- vector
      } else {
        dataset <- rbind(dataset, vector)
      }
    } else {
      datasets[[length(datasets) + 1]] <- dataset
      dataset <- NULL
    }
  }
  # If the last dataset was not added to the list, add it
  if (!is.null(dataset)) {
    datasets[[length(datasets) + 1]] <- dataset
  }
  
  return(datasets)
}

# Function to compute average pairwise differences
average_pairwise_diff <- function(sequences) {
  num_sequences <- nrow(sequences)
  total_diff <- 0
  num_pairs <- 0
  
  for (i in 1:(num_sequences - 1)) {
    for (j in (i + 1):num_sequences) {
      total_diff <- total_diff + sum(sequences[i,] != sequences[j,])
      num_pairs <- num_pairs + 1
    }
  }
  
  return(total_diff / num_pairs)
}

calculate_w <- function(sequences) {
  S <- ncol(sequences)
  n <- nrow(sequences)
  a_1 <- 0
  for (i in 1:(n-1)) {
    a_1 <- a_1 + (1/i)
  }
  w = S/a_1
  return (w)
}

tajimas_d <- function(sequences, k, w) {
  S <- ncol(sequences)
  n <- nrow(sequences)
  a1 <- 0
  a2 <- 0
  b1 <- 0
  b2 <- 0
  c1 <- 0
  c2 <- 0
  e1 <- 0
  e2 <- 0
  de <- 0 #denominator
  D <- 0
  
  
  for (i in 1:(n-1)){
    a1 <- a1 + (1/i)
  }
  for (i in 1:(n-1)){
    a2 <- a1 + (1/(i^2))
  }
  
  b1 <- (n+1) / (3*(n-1))
  b2 <- ( 2*(n^2+n+3) ) / ( 9*n*(n-1) ) 
  
  c1 <- b1 - (1/a1)
  c2 <- b2 - (n+2)/(a1*n) + a2/a1^2
 
  e1 <- c1/a1
  e2 <- c2/(a1^2+a2)
  
  de <- sqrt(e1*S + e2*S*(S-1))
  D <- (k-w)/de
  return (D)
}

standardize <- function(vec) {
  vec_mean <- mean(vec)
  vec_sd <- sd(vec)
  norm_vec <- (vec-vec_mean)/vec_sd
  return (norm_vec)
}

standardize_obs <- function(obs, sims_vec){
  sims_mean <- mean(sims_vec)
  sims_sd <- sd(sims_vec)
  norm_obs <- (obs - sims_mean) / sims_sd
  return (norm_obs)
}

euclidean_distance <- function(obs_vec, sims_vec) {
  distance <- sqrt(sum((obs_vec - sims_vec)^2))
  return(distance)
}

### Execution ###

# Reading data
param_data = read.table(param_path)
sim_datasets <- read_datasets("ms_sim_final.out")
obs_dataset <- read_datasets(obs_path)[[1]]

# Using function to calculate k
obs_k <- average_pairwise_diff(obs_dataset)
sims_k <- lapply(sim_datasets, average_pairwise_diff)

# Using function to calculate w
obs_w <- calculate_w(obs_dataset)
sims_w <- lapply(sim_datasets, calculate_w)

# Using function to calculate Tajima's D
obs_D <- tajimas_d(obs_dataset, obs_k, obs_w)
sims_D <- mapply(tajimas_d, sim_datasets, sims_k, sims_w)

# Make them numeric
sims_D <- as.numeric(sims_D)
sims_k <- as.numeric(sims_k)
sims_w <- as.numeric(sims_w)

# Standardize vectors
sims_k_std <- standardize(sims_k)
sims_w_std <- standardize(sims_w)
sims_D_std <- standardize(sims_D)
obs_k_std <- standardize_obs(obs_k, sims_k)
obs_w_std <- standardize_obs(obs_w, sims_w)
obs_D_std <- standardize_obs(obs_D, sims_D)

obs_vec <- c(obs_D_std, obs_w_std, obs_k_std)
sim_vec <- c(sims_D_std, sims_w_std, sims_k_std)

# Calculate Euclidean distances
distances <- sapply(1:length(sims_D_std), function(i) {
  sims_vec <- c(sims_D_std[i], sims_w_std[i], sims_k_std[i])
  euclidean_distance(obs_vec, sims_vec)
})

# Find the 500 smallest distances
min_distance_indices <- order(distances)[1:500]

min_distance_params <- param_data[min_distance_indices, ]

params_mean <- mean(min_distance_params)
params_median <- median(min_distance_params)

hist(min_distance_params, breaks=20,  main="Histogram of Min Distance Params", xlab="Params", col="blue", border="black")
plot(density(min_distance_params), main="Density plot of Min Distance Params", xlab="Params", ylab="Density")
print(params_mean)
print(params_median)
