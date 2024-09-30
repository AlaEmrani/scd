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
generate_reference_models <- function(number_of_nodes, number_of_samples,
                                      type="ScaleFree", density_of_graph = 0.2,
                                      power = 1, number_of_changes, mult=1)
{
  #######################################
  ########### Generate Model ############
  #######################################
  set.seed(0)
  
  # zero matrix generation
  A <- matrix(rep(0,number_of_nodes*number_of_nodes), nrow=number_of_nodes)
  B <- matrix(rep(0,number_of_nodes*number_of_nodes), nrow=number_of_nodes)
  total_possible_edges <- choose(number_of_nodes, 2)
  
  # create common model
  random_weights <- rnorm(total_possible_edges)
  random_weights <- random_weights + mult*sign(random_weights)
  
  
  # define base of precision matrix
  A[lower.tri(A)] <- random_weights
  A <- t(A)
  A[lower.tri(A)] <- random_weights
  
  B <- A
  
  
  change_mask <- matrix(data = 0, nrow = number_of_nodes,
                        ncol = number_of_nodes)
  indices <- 1:choose(number_of_nodes,2)
  change_mask[lower.tri(change_mask)] <- indices
  change_mask <- t(change_mask)
  change_mask[lower.tri(change_mask)] <- indices
  mask <- sample(indices,number_of_changes)
  change_mask[change_mask %in% mask] <- -1
  change_mask[change_mask != -1] <- 0
  # create structure
  if(type == "ScaleFree") {
    g <- barabasi.game(number_of_nodes, power, directed = F);
    graph_matrix <- as_adjacency_matrix(g, sparse = F);
    A <- A*graph_matrix;
    B <- B*(graph_matrix+change_mask);
  } else {
    errorCondition(message = "The type value is not standard!")
  }
  
  
  # making positive definite
  minimum_eigen_value <- min(c(eigen(A)$values, eigen(B)$values))
  
  if(minimum_eigen_value < 1)
  {
    A <- A + diag(1-minimum_eigen_value, nrow = number_of_nodes,
                  ncol = number_of_nodes)
    B <- B + diag(1-minimum_eigen_value, nrow = number_of_nodes,
                  ncol = number_of_nodes)
  }
  
  
  #########################################
  ########### Generate Samples ############
  #########################################
  covariance_matrix_A <- solve(A)
  samples_A <- rmvnorm(number_of_samples, mean=rep(0,number_of_nodes), sigma=covariance_matrix_A)
  covariance_matrix_B <- solve(B)
  samples_B <- rmvnorm(number_of_samples, mean=rep(0,number_of_nodes), sigma=covariance_matrix_B)
  
  list(samples_A = samples_A, precision_matrix_A = A, covariance_matrix_A = covariance_matrix_A,
       samples_B = samples_B, precision_matrix_B = B, covariance_matrix_B = covariance_matrix_B)
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
                                      alternate_diag_index = 2) {
  SB = SA = matrix(0, nrow = numberOfNodes, ncol = numberOfNodes)
  for (i in 1:numberOfNodes) {
    for (j in 1:numberOfNodes) {
      SB[i,j] = SA[i,j] = decay_ratio^(abs(i-j))
      if(floor(numberOfNodes/alternate_diag_index) == abs(i-j))
        SB[i,j] = alternative_value
    }
  }
    covarianceMatrixA <- solve(SA);
    XA <- rmvnorm(numberOfSamples, mean=rep(0,numberOfNodes), sigma=covarianceMatrixA);

    covarianceMatrixB <- solve(SB);
    XB <- rmvnorm(numberOfSamples, mean=rep(0,numberOfNodes), sigma=covarianceMatrixB);
    list(samplesA = XA, precisionMatrixA = SA, covarianceMatrixA = covarianceMatrixA,
         samplesB = XB, precisionMatrixB = SB, covarianceMatrixB = covarianceMatrixB);
}

checker <- function(XA, XB, same_indices, change_indices) {
  same_indices = same_indices + 1
  change_indices = change_indices + 1
  n1 <- nrow(XA)
  n2 <- nrow(XB)
  d <- ncol(XA)

  s1 <- t(XA) %*% XA
  s2 <- t(XB) %*% XB

  s <- s1 + s2

  s1_c <- s1[same_indices, same_indices] - s1[same_indices, change_indices] %*%
    solve(s1[change_indices, change_indices]) %*% s1[change_indices, same_indices]
  s2_c <- s2[same_indices, same_indices] - s2[same_indices, change_indices] %*%
    solve(s2[change_indices, change_indices]) %*% s2[change_indices, same_indices]
  s_c <- s1_c + s2_c

  n1_c <- n1 - d + length(same_indices)
  n2_c <- n2 - d + length(same_indices)
  n_c <- n1_c + n2_c

  a21 <- (n1_c^2) * (sum(diag(s1_c %*% s1_c)) - (sum(diag(s1_c))^2)/n1_c) / ( (n1_c+2) * (n1_c-1) )
  a22 <- (n2_c^2) * (sum(diag(s2_c %*% s2_c)) - (sum(diag(s2_c))^2)/n2_c) / ( (n2_c+2) * (n2_c-1) )
  a2 <- (n_c^2) * (sum(diag(s_c %*% s_c)) - (sum(diag(s_c))^2)/n_c) / ( (n_c+2) * (n_c-1) )

  tnm <- ( n1_c*n2_c*(a21 + a22 - 2*sum(diag(s1_c %*% s2_c))) ) / ( 2 * (n1_c + n2_c) * a2 )

  2*pnorm(-abs(tnm), mean = 0, sd = 1)
}

checker_v2 <- function(SA, SB, same_indices, change_indices, repetition=1) {
  same_indices = same_indices + 1
  change_indices = change_indices + 1
  
  nA <- nrow(SA)
  nB <- nrow(SB)
  d <- ncol(SA)

  same_length <- length(same_indices)
  repetition <- repetition*same_length

  p <- NULL
  for(i in 1:repetition){
    c <- matrix(rnorm(same_length), ncol = 1)

    data <- rbind(data.frame(y = SA[, same_indices] %*% c,
                             X = SA[, change_indices]),
                  data.frame(y = SB[, same_indices] %*% c,
                             X = SB[, change_indices]))

    # Load the strucchange package
    library(strucchange)
    library(poolr)
    # Specify the break point
    breakpoint <- nA
    # Perform the Chow test
    chow_test <- sctest(y ~ ., type = "Chow", point = breakpoint, data = data)
    p <- c(p, chow_test$p.value)
  }
  # meta analysis methods with p-values
  list(tippet = tippett(p)$p, fisher = fisher(p)$p,
       invchisq = invchisq(p)$p, stouffer = stouffer(p)$p)
}

checker_v3 <- function(SA, SB, same_indices, change_indices, repetition=1) {
  same_indices = same_indices + 1
  change_indices = change_indices + 1

  nA <- nrow(SA)
  nB <- nrow(SB)
  d <- ncol(SA)

  same_length <- length(same_indices)
  change_length <- length(change_indices)

  p.c <- NULL
  p.f <- NULL
  for(i in same_indices){
    data <- rbind(data.frame(y = SA[, i],
                             X = SA[, change_indices]),
                  data.frame(y = SB[, i],
                             X = SB[, change_indices]))

    # Load the strucchange package
    library(strucchange)
    library(poolr)
    # Specify the break point
    breakpoint <- nA

    formul <- reformulate(response="y", termlabels=c(colnames(data)[-1], "0"))

    # Perform the Chow test
    chow_test <- sctest(formul, type = "Chow", point = breakpoint, data = data)

    OLS_A_residuals <- lm(formula = formul, data = data[1:nA,])$residuals
    OLS_B_residuals <- lm(formula = formul, data = data[(nA+1):(nA+nB),])$residuals
    num = (t(OLS_A_residuals) %*% OLS_A_residuals)/(nA-change_length)
    denom = (t(OLS_B_residuals) %*% OLS_B_residuals)/(nB-change_length)
    f = num/denom
    p.v = 1
    if( f > 1) {
      p.v = pf(f, df = nA-change_length, df2 = nB-change_length, lower.tail = FALSE) +
        pf(1/f, df = nA-change_length, df2 = nB-change_length, lower.tail = TRUE)
    } else {
      p.v = pf(f, df = nA-change_length, df2 = nB-change_length, lower.tail = TRUE) +
        pf(1/f, df = nA-change_length, df2 = nB-change_length, lower.tail = FALSE)
    }

    p.c <- c(p.c, chow_test$p.value)
    p.f <- c(p.f, p.v)
  }

  # meta analysis methods with p-values
  list(tippet = min(tippett(p.c)$p, tippett(p.f)$p),
       fisher = min(fisher(p.c)$p, fisher(p.f)$p),
       invchisq = min(invchisq(p.c)$p, invchisq(p.f)$p),
       stouffer = min(stouffer(p.c)$p, stouffer(p.f)$p))
}
