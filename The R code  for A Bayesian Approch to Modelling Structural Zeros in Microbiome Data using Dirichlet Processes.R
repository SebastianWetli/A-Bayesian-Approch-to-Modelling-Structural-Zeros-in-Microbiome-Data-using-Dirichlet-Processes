library(MCMCpack)
library(dirichletprocess)
library(purrr)
library(plyr)
library(prodlim)
library(mgcv)
library(tibble)
library(ggplot2)
library(profvis)
library(microbenchmark)
library(parallel)
library(Rcpp)
library(multicool)
library(gridExtra)
library(grid)
library(Rcpp)
library(latex2exp)

##### Parameters #####
n <- 200
d <- 5
k <- 10
N <- 30

##### default hyperparameters #####
alpha <- 10
p <- 0.3
lambda <- 1 

### let m be 4. ###
m <- 4 

##### sample true theta and initiation theta #####
start <- function(n, alpha, p, lambda){
  G_0 <- function(reps){
    return(t(sapply(1:reps, function(i)
      rdirichlet(n=1, lambda*c(1, rbinom(d-1, 1, p))))))
  }
  # Simulate data from model
  M <- 1000
  b <- rbeta(M, 1, alpha)
  pivec <- numeric(M)
  pivec[1] <- b[1]
  pivec[2:M] <- sapply(2:M, function(i) b[i] * prod(1 - b[1:(i-1)]))
  theta_star <- G_0(M)
  theta <- theta_star[sample(1:M, size=n, prob = pivec, replace = TRUE),]
}

######################################
######################################
##### estimation of theta solely #####
######################################
######################################

algo8_c_s <- function(n,theta,Y, alpha, p, lambda){
  unitheta <- uniquecombs(theta)
  freq <- c()
  for (i in 1:dim(matrix(unitheta, ncol = d, byrow = TRUE))[1]){
    freq[i] <-   sum(attributes(unitheta)$index==i)
  }
  for (j in 1:n){
    chosen_row <- row.match(theta[j,], matrix(unitheta, ncol = d))
    freq[chosen_row] <- freq[chosen_row] - 1
    kminus <- dim(matrix(unitheta, ncol = d))[1] 
    p <- c() # these are the probabilities to choose a specific value for c_i
    h <- kminus + m # h is found
    for (i in 1:kminus){
      p <- append(p, freq[i]/(n-1+alpha) *
                    dmultinom(Y[j,], prob = matrix(unitheta, ncol = d)[i,1:d]))
    }
    if (freq[chosen_row]>0){
      auxillary <- matrix(NA, nrow = m, ncol = d)
      for (i in 1:m){
        auxillary[i,] <- rdirichlet(n = 1, alpha = c(1,rbinom(n = d-1, size = 1, p = 0.3)))
      }
      phi <- rbind(matrix(unitheta, ncol = d), auxillary)
    } else {
      auxillary <- matrix(NA, nrow = m, ncol = d)
      auxillary[1,] <- theta[j,]
      for (i in 2:m){
        auxillary[i,] <- rdirichlet(n = 1, alpha = c(1,rbinom(n = d-1, size = 1, p = 0.3)))
      }
      phi <- rbind(matrix(unitheta[,1:d], ncol = d), auxillary)
    }
    for (l in 1:m){
      p <- append(p,(alpha/m)/(n - 1 + alpha) * dmultinom(Y[j,], prob = auxillary[l,]))
    }
    p <- p/sum(p)
    sampled_row_number <- sample(1:h, size = 1, prob = p)
    if (sampled_row_number <= dim(matrix(unitheta, ncol = d))[1]){
      freq[sampled_row_number] <- freq[sampled_row_number] + 1
    } else {
      freq <- append(freq,1)
      unitheta <- rbind(matrix(unitheta, ncol = d), phi[sampled_row_number])
    }
    if (freq[chosen_row] == 0){
      freq <- freq[-chosen_row]
      unitheta <- matrix((matrix(unitheta, ncol = d)[-chosen_row,]), ncol = d)
    }
    theta[j,] <- phi[sampled_row_number,]
    }
  theta
}

extend_gamma <- function(x){
  y <- 0
  if (x>0){
    y <- lgamma (x)
  }
  y
}

extend_gamma <- Vectorize(extend_gamma)

phi_sampler <- function(n, theta,Y, alpha, p, lambda){
  c_nums <- dim(uniquecombs(rbind(theta, matrix(0, ncol = d, nrow = 1))))[1] - 1
  for (i in 1:c_nums){
    indices <- uniquecombs(theta)
    rele_Y <- rbind(Y[attributes(indices)$index==i,], matrix(0, nrow = 1, ncol =d))
    zeropos <- colSums(rele_Y)==0
    zeropos[1] <- FALSE
    onepos <- colSums(rele_Y)>0
    onepos <- which(onepos)
    onepos <- unique(c(1,onepos))
    numberzero <- sum(zeropos)
    v_pos <-expand.grid(rep(list(c(0,1)),numberzero))
    if (length(v_pos) > 0){
      for (j in onepos){
        v_pos <- add_column(v_pos, rep(1,2^numberzero),.after = j-1)
      }
      summed_y_s <- colSums(rele_Y)
      log_a_v <- c()
      log_a_v <- apply(v_pos, 1, function(x) (sum(x) - 1) * log(p) + (d - sum(x)) * log(1-p) +
                         extend_gamma(lambda * sum(x)) + sum(extend_gamma(lambda * x + summed_y_s)) -
                         sum(extend_gamma(lambda * x)) - lgamma(sum(lambda * x + summed_y_s)))
      probs <-c()
      for (k in 1:(2^numberzero)){
        #probs[k] <- 1/(sum(exp(log_a_v[k] - log_a_v)))
        probs[k] <- 1/(sum(exp(log_a_v - log_a_v[k])))
      }
      select_v <- v_pos[sample(seq(1:2^numberzero), size = 1, replace = TRUE, prob = probs),]
    } else {
      select_v <- c(rep(1,d))
    }
    alpha_star <- c(lambda * select_v + colSums(rele_Y))
    select_phi <- rdirichlet(1, as.numeric(alpha_star))
    theta <- t(theta)
    theta[,attributes(indices)$index==i] <- select_phi
    theta <- t(theta) 
    
  }
  theta
}

Neal_8_algo <- function(n, theta, Y, k, alpha, p, lambda){
  thetaposterior <- vector("list",k)
  for (j in 1:k){
    theta <- phi_sampler(n,theta,Y, alpha, p, lambda)
    theta <- algo8_c_s(n,theta,Y, alpha, p, lambda)
    thetaposterior[[j]] <- theta 
  }
  return(thetaposterior)
}

sample_path_calculator <- function(n, k, r, alpha, p, lambda){
  sample_path_list <- vector("list", r)
  theta_list <- vector("list", r)
  for (i in 1:r){
    print(i)
    theta <- start(n, alpha, p, lambda)
    Y <- matrix(nrow = n, ncol = d)
    for (j in 1:n){
      Y[j,] <- rmultinom(n=1, size = N, prob = theta[j,])
    }
    start_theta <- start(n, alpha, p, lambda)
    theta_list[[i]] <- theta
    sample_path_list[[i]] <- Neal_8_algo(n, start_theta, Y, k, alpha, p, lambda)
  }
  list("sample_path_list" = sample_path_list, "theta_list" = theta_list)
}

set.seed(100)

sample_path_and_theta <- sample_path_calculator(n = 200, k = 1020, r = 100, alpha = 10, p = 0.3, lambda = 1)

#save(sample_path_and_theta, file = "sample_path_and_theta7")

#load("sample_path_and_theta7")

different_sample_paths <- sample_path_and_theta$sample_path_list

theta_list <- sample_path_and_theta$theta_list

##########################################################
##########################################################
##### Examination of the convergence of theta solely #####
##########################################################
##########################################################

##### Examining the convergence of theta( burn and lag for loss and burn, lag for accuracy), #####
##### where all hyperparameters are known ##### 


burn_theta <- function(n, different_sample_paths, theta_list, number_lags, r){
  loss_matrix <- matrix(NA, nrow = r, ncol = 50)
  loss_theta_burn <- function(number_lags = number_lags, burn_in, num_values = 40, sample_path, theta){
    mean_theta <- Reduce('+',sample_path[seq(burn_in, (num_values - 1) * number_lags +
                                               burn_in, number_lags)])/num_values
    distance <- sum((theta-mean_theta)^2)/n
    distance
  }
  for (i in 1:r){
    sample_path <- different_sample_paths[[i]] 
    theta <- theta_list[[i]]
    for (j in 1:50){
      loss_matrix[i,j] <- loss_theta_burn(number_lags = number_lags, burn_in = 10 * j, sample_path = sample_path, theta = theta)
    }
  }
  loss_matrix
}

burn_loss_matrix <- burn_theta(n = 200, different_sample_paths = different_sample_paths, theta_list = theta_list,
                               number_lags = 7, r = 100)


burn_theta_plot_maker <- function(burn_loss_matrix){
  
  mean_burn_theta <- colMeans(burn_loss_matrix)
  
  lowerquantile_burn_theta <- c()
  upperquantile_burn_theta <- c()
  for (i in 1:25){
    lowerquantile_burn_theta[i] <- quantile(burn_loss_matrix[,i])[2]
    upperquantile_burn_theta[i] <- quantile(burn_loss_matrix[,i])[4]
  }
  
  data_lowerquantile_burn_theta <- data.frame("burn_in" = seq(10, 500, 10), "loss" = lowerquantile_burn_theta)
  data_mean_burn_theta <- data.frame("burn_in" = seq(10, 500, 10), "loss" = mean_burn_theta)
  data_upperquantile_burn_theta <- data.frame("burn_in" = seq(10, 500, 10), "loss" = upperquantile_burn_theta)
  
  p <- ggplot(data = data_mean_burn_theta) + geom_line(aes(burn_in, loss), color = "red", size = 2) +
    geom_line(data = data_lowerquantile_burn_theta, aes(burn_in, loss), color = "red", size = 1.5, alpha = 0.3) +
    geom_line(data = data_upperquantile_burn_theta, aes(burn_in, loss), color = "red", size = 1.5, alpha = 0.3) +
    ggtitle("Loss at different burn-ins") + xlab(TeX("burn-in period"))
  p
}

burn_theta_plot <- burn_theta_plot_maker(burn_loss_matrix)

burn_theta_plot

lag_theta <- function(n, different_sample_paths, theta_list, burn_in, r){
  loss_matrix <- matrix(NA, nrow = r, ncol = 50)
  loss_theta_lag <- function(burn_in = burn_in, number_lags, num_values = 20, sample_path, theta){
    mean_theta <- Reduce('+',sample_path[seq(burn_in, (num_values - 1) * number_lags +
                                               burn_in, number_lags)])/num_values
    distance <- sum((theta-mean_theta)^2)/n
    distance
  }
  for (i in 1:r){
    sample_path <- different_sample_paths[[i]] 
    theta <- theta_list[[i]]
    for (j in 1:50){
      loss_matrix[i,j] <- loss_theta_lag(burn_in = burn_in, number_lags = j, sample_path = sample_path, theta = theta) 
    }
  }
  loss_matrix
}

lag_loss_matrix <- lag_theta(n = 200, different_sample_paths = different_sample_paths, theta_list = theta_list,
                             burn_in = 10, r = 100)

lag_theta_plot_maker <- function(lag_loss_matrix){
  
  mean_lag_theta <- colMeans(lag_loss_matrix)
  
  lowerquantile_lag_theta <- c()
  upperquantile_lag_theta <- c()
  for (i in 1:50){
    lowerquantile_lag_theta[i] <- quantile(lag_loss_matrix[,i])[2]
    upperquantile_lag_theta[i] <- quantile(lag_loss_matrix[,i])[4]
  }
  
  data_lowerquantile_lag_theta <- data.frame("lag" = seq(1, 50, 1), "loss" = lowerquantile_lag_theta)
  data_mean_lag_theta <- data.frame("lag" = seq(1, 50, 1), "loss" = mean_lag_theta)
  data_upperquantile_lag_theta <- data.frame("lag" = seq(1, 50, 1), "loss" = upperquantile_lag_theta)
  
  p <- ggplot(data = data_mean_lag_theta) + geom_line(aes(lag, loss), color = "red", size = 2) +
    geom_line(data = data_lowerquantile_lag_theta, aes(lag, loss), color = "red", size = 1.5, alpha = 0.3) +
    geom_line(data = data_upperquantile_lag_theta, aes(lag, loss), color = "red", size = 1.5, alpha = 0.3) +
    ggtitle("Loss at different lags")
  p
}

lag_theta_plot_maker(lag_loss_matrix)

difference_burn_lag_plot_maker <- function(different_sample_paths, theta_list, burn_in, number_lags, num_values = 20){
  #average_true_theta <- Reduce('+',theta_list[seq(1, 100, 1)])/100
  mean_thetas <- vector("list", 100)
  for (i in 1:100){
    sample_path <- different_sample_paths[[i]]
    mean_thetas[[i]]<- Reduce('+',sample_path[seq(burn_in, (num_values - 1) * number_lags +
                                                    burn_in, number_lags)])/num_values
  }
  abs_difference_list <- vector("list", 100)
  difference_list <- vector("list", 100)
  for (i in 1:100){
  abs_difference_list[[i]] <- abs(theta_list[[i]] - mean_thetas[[i]])
  difference_list[[i]] <- theta_list[[i]] - mean_thetas[[i]]
  }
  list("difference_list" = difference_list, "abs_difference_list" = abs_difference_list)
  #estimate_theta <- Reduce('+',mean_thetas[seq(1,100,1)])/100
  # list("estimate_theta" = estimate_theta, "average_true_theta" = average_true_theta)
}

average_true_and_estimate_theta <- difference_burn_lag_plot_maker(different_sample_paths = different_sample_paths, theta_list = theta_list, burn_in = 10, 
                                                       number_lags = 10)


filled.contour(t(Reduce('+', average_true_and_estimate_theta$abs_difference_list)/100), plot.title =
                 title(main = "absolute difference between\ntrue theta and\n\ estimated theta"))

par(mfrow=c(1,2))

image(t(theta_list[[1]]), xaxt='n', yaxt= 'n') + title("The first true theta")

image(t(different_sample_paths[[1]][[1020]]), xaxt='n', yaxt= 'n') + title("the 1020'th sample")

set.seed(101)

threshold <- c()

for (j in 1:1000){
  print(j)
  first_component <- c()
  for (i in 1:10000){
    alpha2 <- lambda * c(1,rbernoulli(n = d-1 , p = p))
    sample <- rdirichlet(n = 1, alpha = alpha2)
    first_component[i] <- sample[1]
    first_component <- first_component[first_component>0]
  }
  threshold[j] <- unname(quantile(first_component,0.025))
}
threshold <- mean(threshold)

burn_accuracy <- function(n, different_sample_paths, theta_list, number_lags, r, threshold){
  loss_matrix <- matrix(NA, nrow = r, ncol = 50)
  loss_structure_theta_burn <- function(number_lags, burn_in, num_values = 40, sample_path, 
                                        theta , threshold){
    mean_theta <- Reduce('+',sample_path[seq(burn_in, (num_values - 1) * number_lags +
                                               burn_in, number_lags)])/num_values
    loss <- sum((mean_theta[,2:d] < threshold) != (theta[,2:d] ==0))
    loss
  }
  for (i in 1:r){
    theta <- theta_list[[i]]
    sample_path <- different_sample_paths[[i]]
    for (j in 1:50){
      loss_matrix[i,j] <- loss_structure_theta_burn(burn_in = 10 * j, sample_path = sample_path,
                                                    theta = theta, number_lags = number_lags, threshold = threshold) 
    }
  }
  ((n * (d-1) - loss_matrix)/(n * (d-1))) * 100
}

burn_accuracy_matrix <- burn_accuracy(n = 200, different_sample_paths = different_sample_paths, 
                                      theta_list = theta_list, threshold = threshold, number_lags = 7, r = 100)  

burn_accuracy_plot_maker <- function(burn_accuracy_matrix){
  
  mean_burn_accuracy <- colMeans(burn_accuracy_matrix)
  
  lowerquantile_burn_accuracy <- c()
  upperquantile_burn_accuracy <- c()
  for (i in 1:25){
    lowerquantile_burn_accuracy[i] <- quantile(burn_accuracy_matrix[,i])[2]
    upperquantile_burn_accuracy[i] <- quantile(burn_accuracy_matrix[,i])[4]
  }
  
  data_lowerquantile_lag_theta <- data.frame("burn_in" = seq(10, 500, 10), "accuracy" = lowerquantile_burn_accuracy)
  data_mean_lag_theta <- data.frame("burn_in" = seq(10, 500, 10), "accuracy" = mean_burn_accuracy)
  data_upperquantile_lag_theta <- data.frame("burn_in" = seq(10, 500, 10), "accuracy" = upperquantile_burn_accuracy)
  
  p <- ggplot(data = data_mean_lag_theta) + geom_line(aes(burn_in, accuracy), color = "red", size = 2) +
    geom_line(data = data_lowerquantile_lag_theta, aes(burn_in, accuracy), color = "red", size = 1.5, alpha = 0.3) +
    geom_line(data = data_upperquantile_lag_theta, aes(burn_in, accuracy), color = "red", size = 1.5, alpha = 0.3) +
    ggtitle("Accuracy at different burn-ins") + xlab(TeX("burn-in period"))
  p
  
}

burn_zero_one_plot <- burn_accuracy_plot_maker(burn_accuracy_matrix = burn_accuracy_matrix)

lag_accuracy <- function(n,  burn_in, different_sample_paths, theta_list, r, threshold){
  loss_matrix <- matrix(NA, nrow = r, ncol = 50)
  loss_structure_theta_burn <- function(number_lags, burn_in, num_values = 20, sample_path, 
                                        theta , threshold){
    mean_theta <- Reduce('+',sample_path[seq(burn_in, (num_values - 1) * number_lags +
                                               burn_in, number_lags)])/num_values
    loss <- sum((mean_theta[,2:d] < threshold) != (theta[,2:d] ==0))
    loss
  }
  for (i in 1:r){
    theta <- theta_list[[i]]
    sample_path <- different_sample_paths[[i]]
    for (j in 1:50){
      loss_matrix[i,j] <- loss_structure_theta_burn(burn_in = burn_in, sample_path = sample_path,
                                                    theta = theta, number_lags = j, threshold = threshold) 
    }
  }
  ((n * (d-1) - loss_matrix)/(n * (d-1))) * 100
}

lag_accuracy_matrix <- lag_accuracy(n = 200, burn_in = 10, different_sample_paths = different_sample_paths, 
                                    theta_list = theta_list, r = 100, threshold = threshold) 

lag_accuracy_plot_maker <- function(lag_accuracy_matrix){
  
  
  mean_lag_accuracy <- colMeans(lag_accuracy_matrix)
  
  lowerquantile_lag_accuracy <- c()
  upperquantile_lag_accuracy <- c()
  for (i in 1:50){
    lowerquantile_lag_accuracy[i] <- quantile(lag_accuracy_matrix[,i])[2]
    upperquantile_lag_accuracy[i] <- quantile(lag_accuracy_matrix[,i])[4]
  }
  
  data_lowerquantile_lag_theta <- data.frame("lag" = seq(1, 50, 1), "accuracy" = lowerquantile_lag_accuracy)
  data_mean_lag_theta <- data.frame("lag" = seq(1, 50, 1), "accuracy" = mean_lag_accuracy)
  data_upperquantile_lag_theta <- data.frame("lag" = seq(1, 50, 1), "accuracy" = upperquantile_lag_accuracy)
  
  p <- ggplot(data = data_mean_lag_theta) + geom_line(aes(lag, accuracy), color = "red", size = 2) +
    geom_line(data = data_lowerquantile_lag_theta, aes(lag, accuracy), color = "red", size = 1.5, alpha = 0.3) +
    geom_line(data = data_upperquantile_lag_theta, aes(lag, accuracy), color = "red", size = 1.5, alpha = 0.3) +
    ggtitle("Accuracy at different lags")
  p
  
}

lag_accuracy_plot <- lag_accuracy_plot_maker(lag_accuracy_matrix = lag_accuracy_matrix)

optimal_accuracy_calculator <- function(n,  burn_in, number_lags, different_sample_paths, theta_list, r, threshold){
  loss_vector <- matrix(NA, nrow = r, ncol = 1)
  loss_structure_theta_burn <- function(number_lags, burn_in, num_values = 40, sample_path, 
                                        theta , threshold){
    mean_theta <- Reduce('+',sample_path[seq(burn_in, (num_values - 1) * number_lags +
                                               burn_in, number_lags)])/num_values
    loss <- sum((mean_theta[,2:d] < threshold) != (theta[,2:d] ==0))
  }
  for (i in 1:r){
    theta <- theta_list[[i]]
    sample_path <- different_sample_paths[[i]]
    loss_vector[i] <- loss_structure_theta_burn(burn_in = burn_in, sample_path = sample_path,
                                                theta = theta, number_lags = number_lags , threshold = threshold) 
  }
  ((n * (d - 1) - loss_vector)/(n * (d-1))) * 100
}

optimal_accuracy <- colMeans(optimal_accuracy_calculator(n = 200, burn_in = 10, number_lags = 7, 
                                                         different_sample_paths= different_sample_paths, 
                                                         theta_list = theta_list, r = 100, threshold =threshold))

##### Examining the consistency of the convergence of theta solely #####

set.seed(102)

consistency <- function(i){
  loss_vector <- c()
  loss_at_hyperparameters <- function(number_lags = 7, burn_in = 10, num_values = 20, sample_path){
    mean_theta <- Reduce('+',sample_path[seq(burn_in, (num_values - 1) * number_lags +
                                               burn_in, number_lags)])/num_values
    distance <- sum((theta-mean_theta)^2)/n
    distance
  }
  n <- i
  for (j in 1:100){
    print(i*100 +j)
    theta <- start(n, alpha, p, lambda)
    Y <- matrix(nrow = n, ncol = d)
    for (l in 1:n){
      Y[l,] <- rmultinom(n=1, size = N, prob = theta[l,])
    }
    start_theta <- start(n, alpha, p, lambda)
    sample_path <- Neal_8_algo(n = n,start_theta,Y = Y, k = 160, alpha = alpha, p = p, lambda = lambda)
    loss_vector[j] <- loss_at_hyperparameters(sample_path = sample_path)
  }
  loss_vector 
}

consistency_value1 <- consistency(32)
consistency_value2 <- consistency(64)
consistency_value3 <- consistency(128)
consistency_value4 <- consistency(256)
consistency_value5 <- consistency(512)
consistency_value6 <- consistency(1024)


consistency_vector <- c(log(mean(consistency_value1)), log(mean(consistency_value2)), log(mean(consistency_value3)),
                        log(mean(consistency_value4)), log(mean(consistency_value5)), log(mean(consistency_value6)))

#save(consistency_vector, file = "consistency_vector")


consistency_data <- data.frame("logarithm_2_n" = log(c(32, 64, 128, 256, 512, 1024)), "log_e_loss" = consistency_vector)

intercept <- lm(as.matrix(consistency_vector) ~ as.matrix(log(c(32, 64, 128, 256, 512, 1024))))$coefficients[1]
slope <- lm(as.matrix(consistency_vector) ~ as.matrix(log(c(32, 64, 128, 256, 512, 1024))))$coefficients[2]

ggplot(data = consistency_data) + geom_point(aes(logarithm_2_n, log_e_loss), size = 4 ,color = "red") + 
  ggtitle("Consistency/Development of loss as n increases") + geom_abline(intercept = intercept , 
                                                                          slope = slope) + xlab(TeX("log of n")) + ylab("log of loss")

##### Examining the correlation and independence for the convergence of theta #####

set.seed(100)

entrywisecorrelation_calculator_lag <- function(n, different_sample_paths, lag_number, r, k = 1020){
  correlation_list <- vector("list", r)
  reshufled_correlation_list <- vector("list", r)
  for (i in 1:r){
    sample_path <- different_sample_paths[[i]]
    correlation_matrix <- matrix(NA, nrow = n, ncol = d)
    reshufled_correlation_matrix <- matrix(NA, nrow = n, ncol = d)
    for (j in 1:(n*d)){
      time_series <- c()
      for (l in 1:k){
        time_series[l] <- sample_path[[l]][j] 
      }
      if (length(unique(time_series)) > 1){
      correlation_matrix[j] <- abs(acf(time_series, plot= FALSE, lag.max = 70)$acf[lag_number])
      reshufled_correlation_matrix[j] <- abs(acf(sample(time_series, size = k, replace = FALSE), 
                                                 plot= FALSE, lag.max = 70)$acf[lag_number])
      } else if (unique(time_series) == 1 | unique(time_series) == 0) {
        correlation_matrix[j] <- 1
        reshufled_correlation_matrix[j] <- 1
      }
    }
    correlation_list[[i]] <- correlation_matrix
    reshufled_correlation_list[[i]] <- reshufled_correlation_matrix 
  }
  list("correlation_list" = correlation_list, "reshufled_correlation_list" = reshufled_correlation_list)
}

entrywisecorrelation <- function(n, different_sample_paths, r, k = 1020,s){
  correlation_lag_lower <- c()
  correlation_lag_median <- c()
  correlation_lag_ninetyfive <- c()
  correlation_lag_max <- c()
  reshufled_lag_lower <- c()
  reshufled_lag_median <- c()
  reshufled_lag_ninetyfive <- c()
  reshufled_lag_max <- c()
  for (i in 1:s){
    print(i)
    lag_entry_reshufled <- entrywisecorrelation_calculator_lag(n = n, different_sample_paths = different_sample_paths, 
                                                               lag_number = i+1, r = r, k = k)
    lag_entry <- lag_entry_reshufled$correlation_list
    lag_reshufled <- lag_entry_reshufled$reshufled_correlation_list
    mean_lag_entry <- as.vector(Reduce('+',lag_entry[seq(1, r, 1)])/r)
    mean_lag_reshufled <- as.vector(Reduce('+',lag_reshufled[seq(1, r, 1)])/r)
    
    correlation_lag_lower[i] <- quantile(mean_lag_entry, 0.25)
    correlation_lag_median[i] <- median(mean_lag_entry)
    correlation_lag_ninetyfive[i]<- quantile(mean_lag_entry, 0.95)
    correlation_lag_max[i] <- max(mean_lag_entry)
    reshufled_lag_lower[i] <- quantile(mean_lag_reshufled, 0.25)
    reshufled_lag_median[i] <- median(mean_lag_reshufled)
    reshufled_lag_ninetyfive[i] <- quantile(mean_lag_reshufled, 0.95)
    reshufled_lag_max[i] <- max(mean_lag_reshufled)
  }
  list("correlation_lag_lower" = correlation_lag_lower,
       "correlation_lag_median" = correlation_lag_median, "correlation_lag_ninetyfive" = correlation_lag_ninetyfive,
       "correlation_lag_max" = correlation_lag_max, "reshufled_lag_lower" = reshufled_lag_lower,
       "reshufled_lag_median" = reshufled_lag_median,
       "reshufled_lag_ninetyfive" = reshufled_lag_ninetyfive, "reshufled_lag_max" = reshufled_lag_max )
}

entrywisewisecorrelation_percentiles <- entrywisecorrelation(n = 200, different_sample_paths = different_sample_paths,
                                                             r = 100, k = 1020, s = 30)

data_correlation_lag_lower <- data.frame("lag" = seq(1,30,1), "correlation" = 
                                           entrywisewisecorrelation_percentiles$correlation_lag_lower)
data_correlation_lag_median <- data.frame("lag" = seq(1,30,1), "correlation" = 
                                            entrywisewisecorrelation_percentiles$correlation_lag_median)
data_correlation_lag_ninetyfive <- data.frame("lag" = seq(1,30,1), "correlation" = 
                                                entrywisewisecorrelation_percentiles$correlation_lag_ninety)
data_correlation_lag_max <- data.frame("lag" = seq(1,30,1), "correlation" = 
                                         entrywisewisecorrelation_percentiles$correlation_lag_max)
data_reshufled_lag_lower <- mean(entrywisewisecorrelation_percentiles$reshufled_lag_lower)
data_reshufled_lag_median <- mean(entrywisewisecorrelation_percentiles$reshufled_lag_median)
data_reshufled_lag_ninetyfive <- mean(entrywisewisecorrelation_percentiles$reshufled_lag_ninetyfive)
data_reshufled_lag_max <- mean(entrywisewisecorrelation_percentiles$reshufled_lag_max)

ggplot() + 
  geom_line(data = data_correlation_lag_lower, aes(lag,correlation) , size = 1, color = "orange") +
  geom_hline(yintercept = data_reshufled_lag_lower, linetype = "dashed", size = 1, color = "orange", alpha = 1) +
  geom_line(data = data_correlation_lag_median, aes(lag,correlation) , size = 1, color = "red") +
  geom_hline(yintercept = data_reshufled_lag_median, linetype = "dashed", size = 1, color = "red", alpha = 1) +
  geom_line(data = data_correlation_lag_ninetyfive, aes(lag,correlation) , size = 1, color = "blue") + 
  geom_hline(yintercept = data_reshufled_lag_ninetyfive , linetype = "dashed", size = 1, color = "blue", alpha = 1) + 
  geom_line(data = data_correlation_lag_max, aes(lag,correlation) , size = 1, color = "green") + 
  geom_hline(yintercept= data_reshufled_lag_max , size = 1, linetype = "dashed", color = "green", alpha = 1) +
  ggtitle("Development of entrywisecorrelation")

set.seed(103)

correlation_across <- function(n, different_sample_paths,r,k = 1020, ma){
  correlation_list <- vector("list", r)
  correlation_across_matrix <- matrix(NA, nrow = n, ncol = ma)
  reshufled_correlation_list <- vector("list", r)
  reshufled_correlation_across_matrix <- matrix(NA, nrow = n, ncol = ma)
  for (i in 1:r){
    print(i)
    timeseries_across <- matrix(NA, nrow = n, ncol = k)  
    sample_path <- different_sample_paths[[i]]
    N_d <- rnorm(d, 0, 1)
    for (j in 1:k){
      timeseries_across[,j] <- sample_path[[j]] %*% N_d
    }
    for (l in 1:n){
      if (length(unique(timeseries_across[l,])) == 1){
        correlation_across_matrix[l,] <- rep(1,m)
        reshufled_correlation_across_matrix[l,] <- rep(1,m)
      } else {
      for (m in 1:ma){
        
        correlation_across_matrix[l,m] <- abs(acf(timeseries_across[l,], 
                                                  plot= FALSE, lag.max = ma+5)$acf[m+1])
        reshufled_correlation_across_matrix[l,m] <- abs(acf(sample(timeseries_across[l,], size = k, replace = FALSE), 
                                                            plot= FALSE, lag.max = ma+5)$acf[m+1])
      }
    }
    }
    correlation_list[[i]] <- correlation_across_matrix
    reshufled_correlation_list[[i]] <- reshufled_correlation_across_matrix
  }
  acrosscorrelation <- Reduce('+', correlation_list[seq(1, r, 1)])/r
  reshufledacrosscorrelation <- Reduce('+', reshufled_correlation_list[seq(1, r, 1)])/r
  list("acrosscorrelation" = acrosscorrelation, "reshufledacrosscorrelation" = reshufledacrosscorrelation)
}

acrosscorrelation <- function(n, different_sample_paths, r, k = 1020, ma){
  correlation_lag_lower <- c()
  correlation_lag_median <- c()
  correlation_lag_ninetyfive <- c()
  correlation_lag_max <- c()
  reshufled_lag_lower <- c()
  reshufled_lag_median <- c()
  reshufled_lag_ninetyfive <- c()
  reshufled_lag_max <- c()
  lag_across_reshufled <- correlation_across(n = n, different_sample_paths = different_sample_paths, r = r, k = k, ma = ma)
  browser()
  for (i in 1:ma){
    lag_across <- as.vector(lag_across_reshufled$acrosscorrelation[,i])
    lag_reshufled <- as.vector(lag_across_reshufled$reshufledacrosscorrelation[,i])
    correlation_lag_lower[i] <- quantile(lag_across, 0.25)
    correlation_lag_median[i] <- median(lag_across)
    correlation_lag_ninetyfive[i]<- quantile(lag_across, 0.95)
    correlation_lag_max[i] <- max(lag_across)
    reshufled_lag_lower[i] <- quantile(lag_reshufled, 0.25)
    reshufled_lag_median[i] <- median(lag_reshufled)
    reshufled_lag_ninetyfive[i] <- quantile(lag_reshufled, 0.95)
    reshufled_lag_max[i] <- max(lag_reshufled)
  }
  list("correlation_lag_lower" = correlation_lag_lower,
       "correlation_lag_median" = correlation_lag_median, "correlation_lag_ninetyfive" = correlation_lag_ninetyfive,
       "correlation_lag_max" = correlation_lag_max, "reshufled_lag_lower" = reshufled_lag_lower,
       "reshufled_lag_median" = reshufled_lag_median,
       "reshufled_lag_ninetyfive" = reshufled_lag_ninetyfive, "reshufled_lag_max" = reshufled_lag_max )
}

acrosscorrelation_percentiles <- acrosscorrelation(n = 200, different_sample_paths = different_sample_paths,
                                                   r = 100, ma = 65)

a_data_correlation_lag_lower <- data.frame("lag" = seq(2,66,1), "correlation" = 
                                             acrosscorrelation_percentiles$correlation_lag_lower)
a_data_correlation_lag_median <- data.frame("lag" = seq(2,66,1), "correlation" = 
                                              acrosscorrelation_percentiles$correlation_lag_median)
a_data_correlation_lag_ninetyfive <- data.frame("lag" = seq(2,66,1), "correlation" = 
                                                  acrosscorrelation_percentiles$correlation_lag_ninety)
a_data_correlation_lag_max <- data.frame("lag" = seq(2,66,1), "correlation" = 
                                           acrosscorrelation_percentiles$correlation_lag_max)
a_data_reshufled_lag_lower <- mean(acrosscorrelation_percentiles$reshufled_lag_lower)
a_data_reshufled_lag_median <- mean(acrosscorrelation_percentiles$reshufled_lag_median)
a_data_reshufled_lag_ninetyfive <- mean(acrosscorrelation_percentiles$reshufled_lag_ninetyfive)
a_data_reshufled_lag_max <- mean( acrosscorrelation_percentiles$reshufled_lag_max)

ggplot() + 
  geom_line(data = a_data_correlation_lag_lower, aes(lag,correlation) , size = 1, color = "orange") +
  geom_hline(yintercept = a_data_reshufled_lag_lower, linetype = "dashed", size = 1, color = "orange", alpha = 1) +
  geom_line(data = a_data_correlation_lag_median, aes(lag,correlation) , size = 1, color = "red") +
  geom_hline(yintercept = a_data_reshufled_lag_median, linetype = "dashed", size = 1, color = "red", alpha = 1) +
  geom_line(data = a_data_correlation_lag_ninetyfive, aes(lag,correlation) , size = 1, color = "blue") + 
  geom_hline(yintercept = a_data_reshufled_lag_ninetyfive, linetype = "dashed", size = 1, color = "blue", alpha = 1) + 
  geom_line(data = a_data_correlation_lag_max, aes(lag,correlation) , size = 1, color = "green") + 
  geom_hline(yintercept = a_data_reshufled_lag_max, linetype = "dashed", size = 1, color = "green", alpha = 1) +
  ggtitle("Development of the correlation across")

######################################
######################################
##### estimation of alpha solely #####
######################################
######################################

cppFunction('double other_theta_one(double alpha, int d, double p, int k){
   double log_prod = 0; 
   for (int i = 0; i < k - 1 ; ++i){
      log_prod += log(alpha * pow(1-p, d-1) + (i + 1));
   }
   return log_prod;
            }')

cppFunction('double denominator(int n, double alpha){
   double sum = 0;
   for (int i = 0; i < n-1 ; ++i){
      sum += log(alpha + (i + 1));
   }
   double total = -sum;
   return total;
}')

alpha_zero_finder <- function(n){
  f <- function(x){
    x * log(1 + n/x) - n/4
  }
  uniroot(f, interval = c(0.05,n/4))$root
}

log_p_of_alpha_condition_c_p <- function(n, x, theta, p, d){
  output <- (d - 1) * log(1 - p)
  unique_thetas <- uniquecombs(theta)
  truth <- row.match(c(1,rep(0, d - 1)), unique_thetas)
  if (!is.na(truth)){
     num_theta_1 <- sum((colSums(t(theta) == c(1, 0, 0, 0, 0)) == d) == TRUE)
     output <- output + other_theta_one(alpha = x, d = d, p = p, k = num_theta_1)
     l_minus_one <- nrow(unique_thetas) - 1
     output <- output + l_minus_one * log(x * (1 - (1-p)^(d-1)))
     output <- output + denominator(n = n, alpha = x)
  } else {
    l_minus_two <- nrow(unique_thetas) - 1
    output <- output + l_minus_two * log(x * (1 - (1 - p )^(d - 1)))
    output <- output + denominator(n = n, alpha = x)
  }
  output
}

accept_probability_alpha <- function(n, x, y, theta, p, d){
  alpha_zero <- alpha_zero_finder(n = n)
  accept_2 <- exp( log_p_of_alpha_condition_c_p(n = n, x = y, theta = theta, p = p, d = d) -
         log_p_of_alpha_condition_c_p(n = n, x = x, theta = theta, p = p, d = d)
  ) * exp((log(20)/alpha_zero) * (x-y))
  accept <- min(1, accept_2)
  accept
}

set.seed(100)

alpha_sample_path <- function(n, theta, p, d){
  alpha_path <- c()
  alpha <- 5
  for (i in 1:1500){
  y <- abs(alpha + 1  * rnorm(1, mean = 0, sd = 1))  
  U <- runif(1, min = 0, max = 1)
  accept <- accept_probability_alpha(n = n, x = alpha, y = y, theta = theta, p = p, d = d)
  if (U<= accept){
    alpha <- y
  }
  alpha_path[i] <- alpha
  }
  alpha_path
}

##### Examination of the convergence of alpha solely #####

set.seed(100)

means_alpha <- c()
truealpha_vector <- c()

for (j in 1:1000){
   print(j)
   truealpha <- 5 + runif(1, min = -4, max = 15)
   truealpha_vector[j] <- truealpha
   theta_1 <- start(n = 200, alpha = truealpha, p = 0.3, lambda = 1)
   alpha_path <- alpha_sample_path(n = 200, theta = theta_1, p = 0.3, d = 5)
   means_alpha[j] <- mean(alpha_path[1000:1500])
}

data_alpha <- data.frame("true_alpha" = truealpha_vector, "means_alpha" = means_alpha)
save(data_alpha, file = "data_alpha")
load("data_alpha")

set.seed(100)

means_alpha2 <- c()
truealpha_vector2 <- c()

for (j in 1:1000){
  print(j)
  truealpha <- 5 + runif(1, min = -4, max = 15)
  truealpha_vector2[j] <- truealpha
  theta_1 <- start(n = 1000, alpha = truealpha, p = 0.3, lambda = 1)
  alpha_path <- alpha_sample_path(n = 1000, theta = theta_1, p = 0.3, d = 5)
  means_alpha2[j] <- mean(alpha_path[1000:1500])
}

data_alpha2 <- data.frame("true_alpha" = truealpha_vector2, "means_alpha" = means_alpha2)
save(data_alpha2, file = "data_alpha2")
load("data_alpha2")

mean(abs(data_alpha$true_alpha - data_alpha$means_alpha))
mean(abs(data_alpha2$true_alpha - data_alpha2$means_alpha))

var(abs(data_alpha$true_alpha - data_alpha$means_alpha))
var(abs(data_alpha2$true_alpha - data_alpha2$means_alpha))

data_differences <- data.frame("difference" = data_alpha$true_alpha - data_alpha$means_alpha)
data_differences2 <- data.frame("difference" = data_alpha2$true_alpha - data_alpha2$means_alpha)

ggplot(data = data_differences, aes(x = difference)) + geom_histogram(aes(y=..count../sum(..count..))) + 
  ggtitle("The differences between true values\nand estimated values") + ylab("proportion") + geom_vline(xintercept = 0, col = "red", size = 2)

ggplot(data = data_differences) + geom_density(aes(difference), col = "red", fill = "red", alpha = 0.1) + ggtitle("Densities over the difference between\ntrue values and estimated values") +
  xlab(TeX("differences")) + ylab("proportion") + geom_density(data = data_differences2, aes(difference), col = "blue", fill = "blue", alpha = 0.1)

alphaplotdata <- data.frame("true_alphas" = data_alpha$true_alpha, "estimated_alphas" = data_alpha$means_alpha)

ggplot(data = alphaplotdata) + geom_point(aes(true_alpha, means_alpha)) + geom_abline(slope = 1, intercept = 0, col = "red" , size = 2) +
  ggtitle("True alpha values versus estimated values of alpha") + xlab(TeX("True values")) + ylab(TeX("Estimated values"))

rela_dif_alpha <- (data_alpha$means_alpha - data_alpha$true_alpha) * 100 / data_alpha$true_alpha

data_alpha_rela <- data.frame("true_alpha" = data_alpha$true_alpha, "rela_dif_alpha" = rela_dif_alpha )

ggplot(data = data_alpha_rela) + geom_point(aes(true_alpha, rela_dif_alpha)) + ggtitle("relative differences between true alphas and\nestimated alphas") +
  xlab(TeX("True alpha")) + ylab("The relative differnece")
  
##################################
##################################                                                                                     
##### estimation of p solely #####
##################################
##################################

p_estimator <- function(n, alpha, p, lambda, N){
  theta <- start(n = n, alpha = alpha, p = p, lambda = lambda)
  Y <- matrix(nrow = n, ncol = d)
  for (i in 1:n){
    Y[i,] <- rmultinom(n=1, size = N, prob = theta[i,])
  }
  #misclassified<- sum(Y[,1]==0)
  number_zeros <- sum(Y[,-1] == 0)
  number_nonzeros <- sum(Y[,-1] > 0)
  number_nonzeros/(n * (d-1))
}

##### Examination of the convergence of p #####

set.seed(100)

true_prob <- c()
estimated_prob <- c()
for (i in 1:1000){
  print(i)
  true_prob[i] <- runif(n = 1, min = 0, max = 1)
  estimated_prob[i]<- p_estimator(n = 200, alpha = 10, p = true_prob[i], lambda = 1, N = 30)

}

difference_prob_data <- data.frame("differences_prob" = estimated_prob - true_prob)
ggplot(data = difference_prob_data) + geom_histogram(aes(differences_prob)) + xlab(TeX("Differences")) + ylab("Count") + 
  ggtitle("Differences between\ntrue and estimated values of p")

true_est_prob_data <- data.frame("true_prob" = true_prob, "estimated_prob" = estimated_prob) 

ggplot(data = true_est_prob_data) + geom_point(aes(true_prob, estimated_prob)) + geom_abline(intercept = 0, slope = 1, col = "red", size = 2) +
  ggtitle("true p values against estimated p values") + xlab(TeX("true p values")) + ylab(TeX("estimated p values"))

set.seed(101)

differences_N <- matrix(NA, ncol = 15, nrow = 1000) 

for (i in 1:15){
  N <- 15 * 2^i
  for (j in 1:1000){
  print((i - 1) * 1000 + j)
  true_p <- runif(n = 1, min = 0, max = 1)
  estimated_p <-  p_estimator(n = 200, alpha = 10, p = true_p, lambda = 1, N = N)
  differences_N[j,i] <- estimated_p - true_p
  }
}

data_dif_N <- data.frame("log_N" = log(15 * 2^seq(1,15,1)/d), "mean_abs_dif" = colMeans(abs(differences_N))) 
  
ggplot(data = data_dif_N) + geom_line(aes(log_N, log(mean_abs_dif))) + xlab(TeX("log of N")) + ylab(TeX("log of mean of absolute differences")) +
  ggtitle("The mean of absolute differences as N increases")

set.seed(102)

differences_n <- matrix(NA, ncol = 10, nrow = 1000)

for (i in 1:10){
  n <- 15 * 2^i
  for (j in 1:1000){
    print((i - 1) * 1000 + j)
    true_p <- runif(n = 1, min = 0, max = 1)
    estimated_p <- p_estimator(n = n, alpha = 10, p = true_p, lambda = 1, N = 30)
    differences_n[j,i] <- estimated_p - true_p
  }
}

n <- 200

data_dif_n <- data.frame("log_n" = log((15 * 2^seq(1,10,1)), "mean_abs_dif" = colMeans(abs(differences_n))))

ggplot(data = data_dif_n) + geom_line(aes(log_n, log(mean_abs_dif)))  + xlab(TeX("log of n")) + ylab("log of mean of absolute differences") + 
  ggtitle("Development of difference as n increases")

################################
################################
##### estimation of lambda #####
################################
################################

accept_probability_lambda <- function(theta, x, y){
  if (0.01 < y & y < 100){
    unique_theta <- uniquecombs(theta)
    index <- row.match(c(1,rep(0,d-1)), unique_theta)
    if (!is.na(index)){
      unique_theta <- unique_theta[-index,]
    }
    num_unique <- nrow(unique_theta)
    one_structure <- unique_theta > 0
    number_ones <- rowSums(one_structure) 
    total <- (y-x) * sum(log(unique_theta[unique_theta > 0]))
    for (i in 1:num_unique){
      total <- total + lgamma(number_ones[i] * y) - number_ones[i] * lgamma(y) - 
        lgamma(number_ones[i] * x) + number_ones[i] * lgamma(x)
    }
    total <- total + log(x) - log(y)
    total <- exp(total)
    accept <- min(1, total)
  } else {
    accept <- 0
  }
  accept
}

lambda_sample_path_calculator <- function(theta, initlambda, k){
  lambda <- initlambda
  lambda_path <- c()
  for (i in 1:k){
    y <- abs(lambda + rnorm(1, mean = 0, sd = 1))
    U <- runif(1, min = 0, max = 1)
    accept <- accept_probability_lambda(theta = theta, x = lambda, y = y)
    if (U <= accept){
      lambda <- y
    }
    lambda_path[i] <- lambda
  }
  lambda_path
}

##### Examination of the convergence of lambda #####

set.seed(100)

means_lambda <- c()
truelambda_vector <- c()

set.seed(100)

for (j in 1:500){
  print(j)
  truelambda1 <- runif(1, min = 1, max = 10)
  truelambda2 <- 1/truelambda1
  truelambda_vector[2 * j - 1] <- truelambda1 
  truelambda_vector[2 * j] <- truelambda2
  theta_1 <- start(n = 200, alpha = 10, p = 0.3, lambda = truelambda1)
  theta_2 <- start(n = 200, alpha = 10, p = 0.3, lambda = truelambda2)
  lambda_path1 <- lambda_sample_path_calculator(theta = theta_1, initlambda = 1, k = 1500)
  lambda_path2 <- lambda_sample_path_calculator(theta = theta_2, initlambda = 1, k = 1500)
  means_lambda[2 * j -1] <- mean(lambda_path1[1000:1500])
  means_lambda[2 * j] <- mean(lambda_path2[1000:1500])
}

mean(abs(means_lambda - truelambda_vector))

data_differences_lambda <- data.frame("differences" = means_lambda - truelambda_vector)

data_lambda <- data.frame("truelambda" = truelambda_vector, "estimatedlambda" = means_lambda)

save(data_lambda, file = "data_lambda")

load("data_differences_lambda")

ggplot(data = data_differences_lambda) + geom_histogram(aes(differences)) + ggtitle("The differences between\nthe true values and estimated values") +
  xlab(TeX("Difference")) + ylab("Count")

ggplot(data = data_lambda) +  geom_point(aes(truelambda, estimatedlambda)) + geom_abline(intercept = 0, slope = 1, size = 2, col = "red") + 
  xlab(TeX("True lambda")) + ylab("Estimated lambda") + ggtitle("True lambda values versus estimated values of lambda")

rela_lambda <- ((means_lambda - truelambda_vector) / truelambda_vector)* 100

data_rela_lambda <- data.frame("truelambda" = truelambda_vector, "rela_lambda" = rela_lambda)

ggplot(data = data_rela_lambda) + geom_point(aes(truelambda, rela_lambda)) + xlab(TeX("True lambda")) + ylab(TeX("Relative error")) + 
  ggtitle("Relative differences between true lambdas\nand estimated lambdas")

means_lambda2 <- c()
truelambda_vector2 <- c()

set.seed(100)

for (j in 1:500){
  print(j)
  truelambda1 <- runif(1, min = 1, max = 10)
  truelambda2 <- 1/truelambda1
  truelambda_vector2[2 * j - 1] <- truelambda1 
  truelambda_vector2[2 * j] <- truelambda2
  theta_1 <- start(n = 1000, alpha = 10, p = 0.3, lambda = truelambda1)
  theta_2 <- start(n = 1000, alpha = 10, p = 0.3, lambda = truelambda2)
  lambda_path1 <- lambda_sample_path_calculator(theta = theta_1, initlambda = 1, k = 1500)
  lambda_path2 <- lambda_sample_path_calculator(theta = theta_2, initlambda = 1, k = 1500)
  means_lambda2[2 * j -1] <- mean(lambda_path1[1000:1500])
  means_lambda2[2 * j] <- mean(lambda_path2[1000:1500])
}

data_lambda_higher_n <- data.frame("difference" = truelambda_vector2 - means_lambda2)  

save(data_lambda_higher_n, file = "data_lambda_higher_n")

ggplot(data =  data_lambda_higher_n) + geom_density(aes(difference), col = "red") + geom_density(data = data_differences_lambda,aes(differences), col = "blue") + 
  xlim(c(-0.5,0.5)) + ggtitle("Densities for the difference\nfor n=200 and n=1000")

########################################
########################################
##### estimation of the full model #####
########################################
########################################

Neal_8_algo_extended <- function(n, inittheta, initalpha, initlambda, Y, k, p){
  thetaposterior <- vector("list",k)
  alphaposterior <- c()
  lambdaposterior <- c()
  number_nonzeros <- sum(Y[,-1] > 0)
  alpha_sample_path2 <- function(n, theta, p, d, alpha){
    alpha_path <- c()
    for (i in 1:1){
      y <- abs(alpha +  rnorm(1, mean = 0, sd = 1))  
      U <- runif(1, min = 0, max = 1)
      accept <- accept_probability_alpha(n = n, x = alpha, y = y, theta = theta, p = p, d = d)
      if (U<= accept){
        alpha <- y
      }
      alpha_path[i] <- alpha
    }
    alpha_path
  }
  lambda_sample_path_calculator <- function(theta, lambda){
    lambda_path <- c()
    for (i in 1:1){
      y <- abs(lambda +  rnorm(1, mean = 0, sd = 1))
      if (0.01<y & y<100){
      U <- runif(1, min = 0, max = 1)
      accept <- accept_probability_lambda(theta = theta, x = lambda, y = y)
      if (U <= accept){
        lambda <- y
      }
      lambda_path[i] <- lambda
      } else{
        lambda_path[i] <- lambda
      }
    }
    lambda_path
  }
  alpha <- initalpha
  lambda <- initlambda
  theta <- inittheta
  for (j in 1:k){
    print(j)
    theta <- phi_sampler(n = n,theta = theta, Y = Y, alpha = alpha, lambda = lambda, p = p)
    theta <- algo8_c_s(n = n,theta = theta,Y = Y, alpha = alpha, lambda = lambda, p = p)
    thetaposterior[[j]] <- theta 
    alpha <- alpha_sample_path2(alpha, n = n, theta = theta, p = p, d = d)
    alphaposterior[j] <- alpha
    lambda <- lambda_sample_path_calculator(theta = theta, lambda = lambda)
    lambdaposterior[j] <- lambda
  }
  list("thetaposterior" = thetaposterior, "alphaposterior" = alphaposterior, "lambdaposterior" = lambdaposterior,
       "p" = p)
}

set.seed(102)

true_theta_samples_list <- vector("list", 100)
theta_samples_list <- vector("list", 100)
alpha_samples_list <- vector("list", 100)
lambda_samples_list <- vector("list" , 100)
p_samples_list <- vector("list", 100)

for (i in 1:100){
  print(i*1000000)
  truetheta <- start(n = 200, alpha = 15, p = 0.4, lambda = 4)
  true_theta_samples_list[[i]] <- truetheta
  Y <- matrix(NA, nrow = n, ncol = d)
  for (j in 1:n){
    Y[j,] <- rmultinom(n = 1, size = N, prob = truetheta[j,])
  }
  p <- sum(Y[,-1] > 0) /(n * (d - 1)) * (n/ (n - sum(Y[,1] == 0)))
  inittheta <- start(n = 200, alpha = 10, p = p, lambda = 1)
  all_sample_paths <- Neal_8_algo_extended(n = 200, inittheta = inittheta, initalpha = 10, initlambda = 1, Y = Y, k = 1000, p = p)
  
  theta_samples_list[[i]] <- all_sample_paths$thetaposterior
  alpha_samples_list[[i]] <- all_sample_paths$alphaposterior
  lambda_samples_list[[i]] <- all_sample_paths$lambdaposterior
  p_samples_list[[i]] <- all_sample_paths$p
}

fullest <- data.frame("theta_samples_list" = theta_samples_list, "alpha_samples_list" = alpha_samples_list, "lambda_samples_list" = lambda_samples_list,
           "p_samples_list" =  p_samples_list)

save(theta_samples_list, file = "theta_samples_list5")
save(true_theta_samples_list, file = "true_theta_samples_list5")
save(alpha_samples_list, file = "alpha_samples_list5")
save(lambda_samples_list, file = "lambda_samples_list5")
save(p_samples_list, file = "p_samples_list5")


load("theta_samples_list5")
load("true_theta_samples_list5")
load("alpha_samples_list5")
load("lambda_samples_list5")
load("p_samples_list5")

##### Examination of the convergence of the estimation of the full model #####

full_burn_loss_matrix <- burn_theta(n = 200, different_sample_paths = theta_samples_list, theta_list = true_theta_samples_list,
                                    number_lags = 7, r = 100)

burn_theta_plot_maker(full_burn_loss_matrix)

full_lag_loss_matrix <- lag_theta(n = 200, different_sample_paths = theta_samples_list, theta_list = true_theta_samples_list,
                                  burn_in = 10, r = 100)

lag_theta_plot_maker(full_lag_loss_matrix)

set.seed(100)
threshold <- c()
p <- 0.4
lambda <- 4

for (j in 1:1000){
  print(j)
  first_component <- c()
  for (i in 1:10000){
    alpha2 <- lambda * c(1,rbernoulli(n = d-1 , p = p))
    sample <- rdirichlet(n = 1, alpha = alpha2)
    first_component[i] <- sample[1]
    first_component <- first_component[first_component>0]
  }
  threshold[j] <- unname(quantile(first_component,0.025))
}
threshold <- mean(threshold)

full_burn_accuracy_matrix <- burn_accuracy(n = 200, different_sample_paths = theta_samples_list, 
                                           theta_list = true_theta_samples_list, threshold = threshold, number_lags = 7, r = 100)  

burn_accuracy_plot_maker(burn_accuracy_matrix = full_burn_accuracy_matrix)

full_lag_accuracy_matrix <- lag_accuracy(n = 200, burn_in = 10, different_sample_paths = theta_samples_list, 
                                         theta_list = true_theta_samples_list, r = 100, threshold = threshold) 

lag_accuracy_plot_maker(lag_accuracy_matrix = full_lag_accuracy_matrix)

average_true_and_estimate_theta2 <- difference_burn_lag_plot_maker(different_sample_paths = theta_samples_list, theta_list = true_theta_samples_list, burn_in = 10, 
                                                                  number_lags = 10)


filled.contour(t(Reduce('+', average_true_and_estimate_theta2$abs_difference_list)/100), plot.title =
                 title(main = "mean absolute difference between\ntrue thetas and\n\ estimated thetas"))

par(mfrow=c(1,2))

image(t(true_theta_samples_list[[1]]), xaxt='n', yaxt= 'n') + title("The first true theta")

image(t(theta_samples_list[[1]][[1000]]), xaxt='n', yaxt= 'n') + title("the 1000'th sample")

Reduce('+', p_samples_list[seq(1, 100, 1)])/100

mean_alpha_sample_path <- Reduce('+', alpha_samples_list[seq(1,100,1)])/100

data_mean_alpha_sample_path <- data.frame("draw_number" = seq(1,1000,1), "alpha_value" = mean_alpha_sample_path)

ggplot(data = data_mean_alpha_sample_path) + geom_line(aes(draw_number, alpha_value), col = "red") + xlab(TeX("Draw number")) + 
  ylab("Value of alpha") + ggtitle("Mean value of alpha") + geom_abline(slope = 0, intercept = 15, col = "blue", size = 1, linetype = "dashed") +
  ylim(c(9,17))

mean_lambda_sample_path <- Reduce('+', lambda_samples_list[seq(1,100,1)])/100

data_mean_lambda_sample_path <- data.frame("draw_number" = seq(1,1000,1), "lambda_value" = mean_lambda_sample_path)

ggplot(data = data_mean_lambda_sample_path) + geom_line(aes(draw_number, lambda_value), col = "red") + xlab(TeX("Draw number")) + 
  ylab("Value of lambda") + ggtitle("mean value of lambda") + geom_abline(slope = 0, intercept = 4, col = "blue", size = 1, linetype = "dashed")

numclusterstrue <- 0
for (i in 1:100){
  numclusterstrue <- numclusterstrue +  nrow(uniquecombs(true_theta_samples_list[[i]]))
}
numclusterstrue/100

numclusters_sample_path <- c()

for (i in 1:1000){
  print(i)
  numclusters <- 0
  for (j in 1:100){
    numclusters <- numclusters + nrow(uniquecombs(theta_samples_list[[j]][[i]]))
  }
  numclusters_sample_path[i] <- numclusters/100
}

data_numclusters_sample_path <- data.frame("draws" = seq(1,1000,1), "clusters" = numclusters_sample_path)

ggplot(data = data_numclusters_sample_path) + geom_line(aes(draws, clusters), col = "red") + 
  geom_abline(intercept = numclusterstrue/100, slope = 0, col = "blue" , size = 1, linetype = "dashed") + ggtitle("Sample path of clusters against true mean number\nof clusters") +
  xlab(TeX("Draw number"))


set.seed(104)

full_entrywisewisecorrelation_percentiles <- entrywisecorrelation(n = 200, different_sample_paths = theta_samples_list,
                                                                  r = 100, k = 1000, s = 60)

full_data_correlation_lag_lower <- data.frame("lag" = seq(1,60,1), "correlation" = 
                                                full_entrywisewisecorrelation_percentiles$correlation_lag_lower)
full_data_correlation_lag_median <- data.frame("lag" = seq(1,60,1), "correlation" = 
                                                 full_entrywisewisecorrelation_percentiles$correlation_lag_median)
full_data_correlation_lag_ninetyfive <- data.frame("lag" = seq(1,60,1), "correlation" = 
                                                     full_entrywisewisecorrelation_percentiles$correlation_lag_ninety)
full_data_correlation_lag_max <- data.frame("lag" = seq(1,60,1), "correlation" = 
                                              full_entrywisewisecorrelation_percentiles$correlation_lag_max)
full_data_reshufled_lag_lower <- mean(full_entrywisewisecorrelation_percentiles$reshufled_lag_lower)
full_data_reshufled_lag_median <- mean(full_entrywisewisecorrelation_percentiles$reshufled_lag_median)
full_data_reshufled_lag_ninetyfive <- mean(full_entrywisewisecorrelation_percentiles$reshufled_lag_ninetyfive)
full_data_reshufled_lag_max <- mean(full_entrywisewisecorrelation_percentiles$reshufled_lag_max)

ggplot() + 
  geom_line(data = full_data_correlation_lag_lower, aes(lag,correlation) , size = 1, color = "orange") +
  geom_hline(yintercept = full_data_reshufled_lag_lower, linetype = "dashed", size = 1, color = "orange", alpha = 1) +
  geom_line(data = full_data_correlation_lag_median, aes(lag,correlation) , size = 1, color = "red") +
  geom_hline(yintercept = full_data_reshufled_lag_median, linetype = "dashed", size = 1, color = "red", alpha = 1) +
  geom_line(data = full_data_correlation_lag_ninetyfive, aes(lag,correlation) , size = 1, color = "blue") + 
  geom_hline(yintercept = full_data_reshufled_lag_ninetyfive , linetype = "dashed", size = 1, color = "blue", alpha = 1) + 
  geom_line(data = full_data_correlation_lag_max, aes(lag,correlation) , size = 1, color = "green") + 
  geom_hline(yintercept= full_data_reshufled_lag_max , size = 1, linetype = "dashed", color = "green", alpha = 1) +
  ggtitle("Development of entrywisecorrelation")

set.seed(105)

full_acrosscorrelation_percentiles <- acrosscorrelation(n = 200, different_sample_paths = theta_samples_list,
                                                        r = 100, k = 1000, ma = 65)

full_a_data_correlation_lag_lower <- data.frame("lag" = seq(2,66,1), "correlation" = 
                                                  full_acrosscorrelation_percentiles$correlation_lag_lower)
full_a_data_correlation_lag_median <- data.frame("lag" = seq(2,66,1), "correlation" = 
                                                   full_acrosscorrelation_percentiles$correlation_lag_median)
full_a_data_correlation_lag_ninetyfive <- data.frame("lag" = seq(2,66,1), "correlation" = 
                                                       full_acrosscorrelation_percentiles$correlation_lag_ninety)
full_a_data_correlation_lag_max <- data.frame("lag" = seq(2,66,1), "correlation" = 
                                                full_acrosscorrelation_percentiles$correlation_lag_max)
full_a_data_reshufled_lag_lower <- mean(full_acrosscorrelation_percentiles$reshufled_lag_lower)
full_a_data_reshufled_lag_median <- mean(full_acrosscorrelation_percentiles$reshufled_lag_median)
full_a_data_reshufled_lag_ninetyfive <- mean(full_acrosscorrelation_percentiles$reshufled_lag_ninetyfive)
full_a_data_reshufled_lag_max <- mean(full_acrosscorrelation_percentiles$reshufled_lag_max)

ggplot() + 
  geom_line(data = full_a_data_correlation_lag_lower, aes(lag,correlation) , size = 1, color = "orange") +
  geom_hline(yintercept = full_a_data_reshufled_lag_lower, linetype = "dashed", size = 1, color = "orange", alpha = 1) +
  geom_line(data = full_a_data_correlation_lag_median, aes(lag,correlation) , size = 1, color = "red") +
  geom_hline(yintercept = full_a_data_reshufled_lag_median, linetype = "dashed", size = 1, color = "red", alpha = 1) +
  geom_line(data = full_a_data_correlation_lag_ninetyfive, aes(lag,correlation) , size = 1, color = "blue") + 
  geom_hline(yintercept = full_a_data_reshufled_lag_ninetyfive, linetype = "dashed", size = 1, color = "blue", alpha = 1) + 
  geom_line(data = full_a_data_correlation_lag_max, aes(lag,correlation) , size = 1, color = "green") + 
  geom_hline(yintercept = full_a_data_reshufled_lag_max, linetype = "dashed", size = 1, color = "green", alpha = 1) +
  ggtitle("Development of the correlation across")

set.seed(102)

full_consistency <- function(i){
  loss_vector <- c()
  loss_at_hyperparameters <- function(number_lags = 10, burn_in = 10, num_values = 20, sample_path){
    mean_theta <- Reduce('+',sample_path[seq(burn_in, (num_values - 1) * number_lags +
                                               burn_in, number_lags)])/num_values
    distance <- sum((theta-mean_theta)^2)/n
    distance
  }
  n <- i
  for (j in 1:100){
    print(i*100 +j)
    theta <- start(n = n, alpha = 15, p = 0.4, lambda = 4)
    Y <- matrix(nrow = n, ncol = d)
    for (l in 1:n){
      Y[l,] <- rmultinom(n=1, size = N, prob = theta[l,])
    }
    p <- sum(Y[,-1] > 0) /(n * (d - 1)) * (n/ (n - sum(Y[,1] == 0)))
    init_theta <- start(n = n, alpha = 10, p = p, lambda = 1)
    sample_path <- Neal_8_algo_extended(n = n, inittheta = init_theta , initalpha = 10, initlambda = 1, Y, k = 310)$thetaposterior
    loss_vector[j] <- loss_at_hyperparameters(sample_path = sample_path)
  }
  loss_vector 
}

full_consistency_value2 <- full_consistency(64)
full_consistency_value3 <- full_consistency(128)
full_consistency_value4 <- full_consistency(256)
full_consistency_value5 <- full_consistency(512)
full_consistency_value6 <- full_consistency(1024)


full_consistency_vector <- c(mean(full_consistency_value2), mean(full_consistency_value3), mean(full_consistency_value4),
                             mean(full_consistency_value5), mean(full_consistency_value6))

load("full_consistency_vector")

data_full_consistency <- data.frame("log_n" = log(2^seq(6,10,1)), "loss" = log(full_consistency_vector))

consismodel <- lm( log(full_consistency_vector) ~ log(2^seq(6,10,1)))

ggplot(data = data_full_consistency) + geom_point(aes(log_n, loss)) + geom_abline(intercept = consismodel$coefficients[1], slope = 
                                                                                    consismodel$coefficients[2], col = "red", size = 2) +
  ggtitle("Consistency/Development of loss as n increases") + xlab(TeX("log of n")) + ylab(TeX("log of loss"))

#####################################################
#####################################################
##### Estimation of theta solely by algorithm 2 #####
#####################################################
#####################################################
 
p <- 0.3

H_sampling <- function(theta,Y){
  onepos <- which(Y>0)
  onepos <- unique(c(1,onepos))
  zeropos <- Y==0
  zeropos[1] <- FALSE
  numberzero <- sum(zeropos)
  v_pos <-expand.grid(rep(list(c(0,1)),numberzero))
  if (length(v_pos) > 0){
    for (j in onepos){
      v_pos <- add_column(v_pos, rep(1,2^numberzero),.after = j-1)
    }
    log_a_v <- c()
    log_a_v <- apply(v_pos, 1, function(x) (sum(x) - 1) * log(p) + (d - sum(x)) * log(1-p) +
                       extend_gamma(lambda * sum(x)) + sum(extend_gamma(lambda * x + Y)) -
                       sum(extend_gamma(lambda * x)) - lgamma(sum(lambda * x + Y)))
    probs <-c()
    for (k in 1:(2^numberzero)){
      probs[k] <- 1/(sum(exp(log_a_v[k] - log_a_v)))
    }
    select_v <- v_pos[sample(seq(1:2^numberzero), size = 1, replace = TRUE, prob = probs),]
  } else {
    select_v <- c(rep(1,d))
  }
  alpha_star <- c(lambda * select_v + Y)
  select_phi <- rdirichlet(1, as.numeric(alpha_star))
  select_phi
}

Neal_algo_2_first_part <- function(theta, Y){
  unitheta <- count(theta)
  freq <- unitheta$freq
  unitheta <- as.matrix(unname(unitheta[,1:d]))
  n <- dim(theta)[1]
  for (i in 1:n){
    matched_row <- which(apply(unitheta, 1, function(x) identical(theta[i,], x)))
    freq[matched_row] <- freq[matched_row] - 1
    cs <- length(freq)
    firstprobs <- c()
    for (j in 1:cs){
       firstprobs[j] <- (freq[j]/(n - 1 + alpha)) * dmultinom(Y[i,], prob = unitheta[j,])
    }
    posones <- which(Y[i,] + c(1,rep(0,d-1)) >0)
    numberzero <- d - length(posones)
    v_pos <- rep(NA,d)
    v_pos[Y[i,]>0] <- 1
    v_pos[Y[i,]==0] <- 0
    integral <- p^(sum(v_pos) - 1) * (1 - p)^(d - sum(v_pos)) * 
                       exp(extend_gamma(sum(lambda * v_pos)) - sum(extend_gamma(lambda * v_pos)) -
                             extend_gamma(N + sum(lambda*v_pos))  + sum(extend_gamma(lambda * v_pos + Y[i,])))

    posiint <- alpha/(n - 1 + alpha) * integral
    prob <- append(firstprobs, posiint)
    prob <- prob/sum(prob)
    index_sampled <-sample(1:(cs+1), size = 1, prob = prob)
    if (index_sampled <= cs){
      theta[i,] <- unitheta[index_sampled,]
      freq[index_sampled] <- freq[index_sampled] + 1
    } else {
      theta[i,] <- H_sampling(theta[i], Y[i,])
      freq <- append(freq,1)
      unitheta <- rbind(unitheta, theta[i,])
    } 
    if (freq[matched_row]==0){
      freq <- freq[-matched_row]
      unitheta <- unitheta[-matched_row,]
    }
  }
  theta
}

Neal_algo_2 <- function(n,theta, Y, k){
  thetaposterior <- vector("list",k)
  for (j in 1:k){
    print(j)
    theta <- phi_sampler(n = n, theta = theta,Y = Y, alpha = alpha, p = p, lambda = lambda)
    theta <- Neal_algo_2_first_part(theta,Y)
    thetaposterior[[j]] <- theta 
  }
  return(thetaposterior)
}

set.seed(105)

true_theta <- start(n = n, alpha = alpha, p = p, lambda = lambda)
Y <- matrix(nrow = n, ncol = d)
for (i in 1:n){
  Y[i,] <- rmultinom(n=1, size = N, prob = true_theta[i,])
}

starttheta <- start(n = n, alpha = 5, p = 0.4, lambda = 4)
theta_posterior <- Neal_algo_2(theta = starttheta, Y = Y, k = 1000)
image(t(theta_posterior[[1000]]))
image(t(true_theta))
