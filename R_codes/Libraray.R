#--------------------------------------------------------
#----------------- General Functions --------------------
#--------------------------------------------------------
#
kendall_tau_matrix_to_pearson_correlation_matrix <- function(kendall_tau_matrix) {
  sin((pi/2)*kendall_tau_matrix)
}


positive_semi_definite_maker <- function(A) {
  # making positive definite
  minimum_eigen_value <- min(eigen(A)$values)
  d <- nrow(A)
  if(minimum_eigen_value < 0) {
    A <- A - diag(minimum_eigen_value, nrow = d, ncol = d)
  }
  A
}


# Estimate pearson correlation using the mapping of the empirical kendall's tau.
## The random_set_ratio arguments contorols the ratio of samples incorporated in estimation.
## It is used in validation step for parameter tuning.
weighted_rank_based_pearson_correlation_estimator <- function(datasets, random_set_ratio=1) {
  number_of_all_cases <- sum(unlist(lapply(datasets, nrow)))
  number_of_nodes <- ncol(datasets[[1]])
  empirical_kendall <- matrix(data = 0, nrow = number_of_nodes, ncol = number_of_nodes)
  
  for(p in 1:length(datasets)) {
    number_of_samples <- nrow(datasets[[p]])
    # Create random sub set of samples
    indices <- sample(1:number_of_samples, size = random_set_ratio*number_of_samples, replace = F)
    # Estimate empirical kendall's tau using weighted mean of empirical kandall's tau of datasets
    empirical_kendall <- empirical_kendall +
      (number_of_samples/number_of_all_cases)*cor(datasets[[p]][indices,], method = "kendall");
  }
  empirical_pearson = kendall_tau_matrix_to_pearson_correlation_matrix(empirical_kendall)
  # make it positive semi definite
  empirical_pearson <- positive_semi_definite_maker(empirical_pearson)
  empirical_pearson
}


#--------------------------------------------------------
#-------- Functions in Real Data Experiments ------------
#--------------------------------------------------------
#
instability_evaluator_of_solution_pathes <- function(solution_pathes, number_of_nodes){
  number_of_random_set <- length(solution_pathes)
  
  instability_results <- data.frame()
  set_index <- rep(1, number_of_random_set)
  while(T){
    max_index = NA;
    max_lambda_value = 0
    edge_frequency = rep(0, number_of_nodes^2)
    number_of_results <- 0;
    for(i in 1:number_of_random_set){
      edge_existence = rep(0, number_of_nodes^2)
      if(set_index[i] <= length(solution_pathes[[i]])) {
        if(set_index[i] > 1) {
          edge_existence[(solution_pathes[[i]][[set_index[i]-1]]$active_set+1)] <- 1
        }
        edge_frequency <- edge_frequency + edge_existence
        
        number_of_results = number_of_results+1
        
        if(max_lambda_value < solution_pathes[[i]][[set_index[i]]]$knots_lambdas){
          max_index = i
          max_lambda_value = solution_pathes[[i]][[set_index[i]]]$knots_lambdas
        }
      }
    }
    if(number_of_results < number_of_random_set){
      break;
    }
    
    if(max_lambda_value == 0){
      break;
    }
    
    Adj <- matrix(edge_frequency/number_of_random_set, nrow = number_of_nodes)
    edge_existence_probability <- Adj[upper.tri(Adj)]
    instability <- sum(2*(1-edge_existence_probability)*edge_existence_probability)/choose(number_of_nodes,2)
    
    instability_results <- rbind(instability_results,
                                 data.frame(lambda=max_lambda_value,
                                            instability=instability))
    
    set_index[max_index] = set_index[max_index]+1;
  }
  instability_results
}


library(stats)
fisher_test <- function(all_genes, hub_genes, functional_genes) {
  functional_index <- all_genes %in% functional_genes
  hub_index <- all_genes %in% hub_genes
  
  pDNA_test <-
    matrix(c(sum(functional_index&hub_index), sum(!functional_index&hub_index),
             sum(functional_index&!hub_index), sum(!functional_index&!hub_index)),
           nrow = 2,
           dimnames = list(inDR = c("T", "F"),
                           isHub = c("T", "F")))
  fisher.test(pDNA_test, alternative = "greater")$p.value
}



#--------------------------------------------------------
#------ Functions in Synthetic Data Experiments ---------
#--------------------------------------------------------
#
library("GGMselect")
library("igraph")
# create two complete precision matrix
generate_reference_models <- function(numberOfNodes, numberOfSamples,
                                      type="ScaleFree", densityOfGraph = 0.2,
                                      power = 1, numberOfChanges, mult=1)
{
  #######################################
  ########### Generate Model ############
  #######################################
  library("GGMselect")
  library("igraph")

  # zero matrix generation
  A <- matrix(rep(0,numberOfNodes*numberOfNodes), nrow=numberOfNodes);
  B <- matrix(rep(0,numberOfNodes*numberOfNodes), nrow=numberOfNodes);
  totalPossibleEdges <- choose(numberOfNodes, 2);

  # create common model
  randomWeights <- rnorm(totalPossibleEdges);
  randomWeights <- randomWeights + mult*sign(randomWeights);


  # define base of precision matrix
  A[lower.tri(A)] <- randomWeights;
  A <- t(A)
  A[lower.tri(A)] <- randomWeights;

  B <- A;


  changeMask <- matrix(data = 0, nrow = numberOfNodes, ncol = numberOfNodes)
  indexes <- 1:choose(numberOfNodes,2)
  changeMask[lower.tri(changeMask)] <- indexes
  changeMask <- t(changeMask)
  changeMask[lower.tri(changeMask)] <- indexes
  mask <- sample(indexes,numberOfChanges)
  changeMask[changeMask %in% mask] <- -1
  changeMask[changeMask != -1] <- 0
  # create structure
  if(type == "Full") {
    A <- A*matrix(1, nrow = numberOfNodes, ncol = numberOfNodes)
    B <- B*(matrix(1, nrow = numberOfNodes, ncol = numberOfNodes)+changeMask)
  }  else if(type == "Erdos") {
    library(igraph)
    g <- erdos.renyi.game(numberOfNodes, densityOfGraph, loops = F, directed = F)
    graphMatrix <- as_adjacency_matrix(g, sparse = F);
    A <- A*graphMatrix;
    B <- B*(graphMatrix+changeMask);
  }  else if(type == "ScaleFree") {
    g <- barabasi.game(numberOfNodes, power, directed = F);
    graphMatrix <- as_adjacency_matrix(g, sparse = F);
    A <- A*graphMatrix;
    B <- B*(graphMatrix+changeMask);
  } else if(type == "Simple") {
    for(i in 1:numberOfNodes)
      for(j in 1:numberOfNodes)
        A[i,j] = decay_value^(abs(i-j))
    B = A
    for(i in 1:numberOfChanges)
      B[i,numberOfNodes-numberOfChanges+i] = B[numberOfNodes-numberOfChanges+i, i] = decay_value
  } else {
    errorCondition(message = "The type value is not standard!")
  }

  if(type != "Simple") {
    # making positive definite
    minimumEigenValue <- min(c(eigen(A)$values, eigen(B)$values))

    if(minimumEigenValue < 1)
    {
      A <- A + diag(1-minimumEigenValue, nrow = numberOfNodes, ncol = numberOfNodes)
      B <- B + diag(1-minimumEigenValue, nrow = numberOfNodes, ncol = numberOfNodes)
    }
  }


  #########################################
  ########### Generate Samples ############
  #########################################
  covarianceMatrixA <- solve(A);
  XA <- rmvnorm(numberOfSamples, mean=rep(0,numberOfNodes), sigma=covarianceMatrixA);

  covarianceMatrixB <- solve(B);
  XB <- rmvnorm(numberOfSamples, mean=rep(0,numberOfNodes), sigma=covarianceMatrixB);
  list(samplesA = XA, precisionMatrixA = A, covarianceMatrixA = covarianceMatrixA,
       samplesB = XB, precisionMatrixB = B, covarianceMatrixB = covarianceMatrixB);
}

# Estimated differential structure evaluator
result_evaluator <- function(real_differences, estimated_differences){
  real_differences_support <- sign(real_differences)
  estimated_differences_support  <- sign(estimated_differences)
  
  real_differences_support <- real_differences_support[upper.tri(real_differences_support)]
  estimated_differences_support <- estimated_differences_support[upper.tri(estimated_differences_support)]
  
  # sign evaluator
  ST <- sum(real_differences_support == estimated_differences_support)
  STN <- sum(real_differences_support == estimated_differences_support & real_differences_support == 0)
  STP <- ST - STN
  
  SF <- sum(real_differences_support != estimated_differences_support)
  SFP <- sum(real_differences_support != estimated_differences_support & real_differences_support == 0)
  SFN <- SF - SFP
  
  STPR <- STP/(STP+SFN)
  SFPR <- SFP/(STN+SFP)
  
  SACC <- ST/(ST+SF)
  
  SPrecision <- STP/(STP+SFP)
  SRecall <- STP/(STP+SFN)
  
  # support evaluator
  NT <- sum(abs(real_differences_support) == abs(estimated_differences_support))
  NTN <- sum(abs(real_differences_support) == abs(estimated_differences_support) & real_differences_support == 0)
  NTP <- NT - NTN
  
  NF <- sum(abs(real_differences_support) != abs(estimated_differences_support))
  NFP <- sum(abs(real_differences_support) != abs(estimated_differences_support) & real_differences_support == 0)
  NFN <- NF - NFP
  
  NTPR <- NTP/(NTP+NFN)
  NFPR <- NFP/(NTN+NFP)
  
  NACC <- NT/(NT+NF)
  
  NPrecision <- NTP/(NTP+NFP)
  NRecall <- NTP/(NTP+NFN)
  
  data.frame(STP=STP, STN=STN, SFP=SFP, SFN=SFN, STPR=STPR, SFPR=SFPR, SACC=SACC, SPrecision=SPrecision, SRecall=SRecall,
             NTP=NTP, NTN=NTN, NFP=NFP, NFN=NFN, NTPR=NTPR, NFPR=NFPR, NACC=NACC, NPrecision=NPrecision, NRecall=NRecall)
}


solution_path_performance_evaluator <- function(solution_path, differential_structure) {
  performances <- NULL
  for(i in 1:length(solution_path$solution_path)) {
    if(i == 1) {
      estimated_differences <- rep(0, length = nrow(differential_structure) * ncol(differential_structure))
    } else {
      estimated_differences <- rep(0, length = nrow(differential_structure) * ncol(differential_structure))
      estimated_differences[1+solution_path$solution_path[[i-1]]$active_set] <- solution_path$solution_path[[i]]$active_set_values
    }
    estimated_differences <- matrix(estimated_differences, nrow=nrow(differential_structure))
    performances <- rbind(performances, data.frame(lambda = solution_path$solution_path[[i]]$knots_lambdas,
                                                   result_evaluator(real_differences = differential_structure,
                                                                   estimated_differences = estimated_differences)))
  }
  performances
}


aggregate_dtrace_solution_path_resuts <- function(number_of_repetition,
                                              dtrace_solution_path_performances){
  # mean of dtrace solution path results
  set_index <- rep(1, number_of_repetition)
  mean_dtrace_solution_path_performances <- NULL;
  remaining_states <- sum(sapply(dtrace_solution_path_performances, nrow))
  
  number_of_results = number_of_repetition
  while(remaining_states > 0 && number_of_results == number_of_repetition){
    max_index = NA;
    max_lambda_value = 0
    temp_performance_record <- NULL;
    number_of_results <- 0;
    for(i in 1:number_of_repetition){
      if(set_index[i] <= nrow(dtrace_solution_path_performances[[i]])){
        temp_performance_record <- rbind(temp_performance_record,
                                         dtrace_solution_path_performances[[i]][set_index[i],])
        number_of_results = number_of_results+1
        
        if(max_lambda_value < dtrace_solution_path_performances[[i]][set_index[i],1]){
          max_index = i
          max_lambda_value = dtrace_solution_path_performances[[i]][set_index[i],1]
        }
      }
    }
    
    if(max_lambda_value == 0){
      break;
    }
    
    temp_performance_record <- apply(temp_performance_record,2,mean, na.rm=T);
    temp_performance_record[1] <- max_lambda_value
    mean_dtrace_solution_path_performances <- rbind(mean_dtrace_solution_path_performances,
                                                  temp_performance_record)
    set_index[max_index] = set_index[max_index]+1;
    remaining_states <- remaining_states - 1;
  }
  colnames(mean_dtrace_solution_path_performances) <- colnames(dtrace_solution_path_performances[[1]])
  mean_dtrace_solution_path_performances <- head(mean_dtrace_solution_path_performances,-1)
  mean_dtrace_solution_path_performances
}

base_differential_network <- function(numberOfNodes, numberOfSamples=1,
                                      decay_ratio = 0.3,
                                      alternative_value = 0.3,
                                      numberOfChanges = 2) {
  SB = SA = matrix(0, nrow = numberOfNodes, ncol = numberOfNodes)
  for (i in 1:numberOfNodes) {
    for (j in 1:numberOfNodes) {
      SB[i,j] = SA[i,j] = decay_ratio^(abs(i-j))
    }
  }
    for(i in 1:numberOfChanges)
      B[i,numberOfNodes-numberOfChanges+i] = B[numberOfNodes-numberOfChanges+i, i] = alternative_value
    covarianceMatrixA <- solve(SA);
    XA <- rmvnorm(numberOfSamples, mean=rep(0,numberOfNodes), sigma=covarianceMatrixA);

    covarianceMatrixB <- solve(SB);
    XB <- rmvnorm(numberOfSamples, mean=rep(0,numberOfNodes), sigma=covarianceMatrixB);
    list(samplesA = XA, precisionMatrixA = SA, covarianceMatrixA = covarianceMatrixA,
         samplesB = XB, precisionMatrixB = SB, covarianceMatrixB = covarianceMatrixB);
}

