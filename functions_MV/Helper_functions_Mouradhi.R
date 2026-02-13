rejection_sampling<- function(x){
  #Numerator
  #Divided by l(1), which is the sum of exp(Uj)
  return( exp(max(x))/(sum(exp(x))) )
}
rejection_sampling_modif<- function(x,RiskF,power_g){
  #Numerator ??
  #Divided by ell(1), which is the sum of exp(Uj).
  Numerator<-exp(power_g*RiskF(x))
  print(paste0("Numerator",Numerator))
  Denominator<-(sum(exp(power_g*x)))
  print(paste0("Denominator",Denominator))
  return(Numerator/Denominator)
}
transformed_logsitic_generator<-function(i, sign, column, alpha,r){

  # Sample from q(a,b) ------------------------------------------------------
  # Use the expression depending on a with  ---------------------------------
  #the condition which(sign==i) on the exp(tj) and ajk. 
  ## Simulation of U is not mev, 
  ## that's why the rejection sampling 
  ## is different from the mev-based one...
  T<- -(log(rexp(length(sign),1)))/(1/alpha)+log((column*r)/(gamma(1-(alpha))))
  T[which(sign==i)]<- -(alpha)*(log(rgamma(1,shape=1-alpha,rate=1)))+log((column[which(sign==i)]*r)/(gamma(1-alpha)))
  return(T)
}

matrix_transformation <- function(Sigma) {
  gamma <- (as.matrix(rep(1, nrow(Sigma))) %*% diag(Sigma)) / 2 + (as.matrix(diag(Sigma)) %*% rep(1, nrow(Sigma))) / 2 - Sigma
  diag(gamma) <- diag(Sigma) / 2
  return(gamma)
}

calculus_of_eta <- function(i, mu, gamma) {
  eta <- (gamma[i, -i] / 2)^(1 / 2) + (mu[i] + diag(gamma)[i] - mu[-i] - diag(gamma)[-i]) / (2 * gamma[i, -i])^(1 / 2)
  return(eta)
}

calculus_of_R <- function(i, gamma) {
  R <- (as.matrix(gamma[i, -i]) %*% rep(1, nrow(gamma) - 1) + t(as.matrix(gamma[i, -i]) %*% rep(1, nrow(gamma) - 1)) - gamma[-i, -i]) / (2 * (as.matrix(gamma[i, -i]) %*% rep(1, nrow(gamma) - 1) * t(as.matrix(gamma[i, -i]) %*% rep(1, nrow(gamma) - 1)))^(1 / 2))
  diag(R) <- 1
  return(R)
}

case_d_2 <- function(mu, Gamma) {
  expec <- exp(mu[1] + diag(Gamma)[1]) * pnorm(((2 * Gamma[1, 2])^(1 / 2)) / (2) - ((mu[2] - mu[1] + (diag(Gamma)[2] - diag(Gamma)[1])) / ((2 * Gamma[1, 2])^(1 / 2))), mean = 0, sd = 1) + exp(mu[2] + diag(Gamma)[2]) * pnorm(((2 * Gamma[1, 2])^(1 / 2)) / (2) - ((mu[1] - mu[2] + (diag(Gamma)[1] - diag(Gamma)[2])) / ((2 * Gamma[1, 2])^(1 / 2))), mean = 0, sd = 1)
  return(expec)
}

calculus_d <- function(d, Sigma, mu) {
  if (d == 1) {
    return(exp(mu + (Sigma / 2)))
  } else {
    expectedvalue <- 0
    gamma <- matrix_transformation(Sigma)
    
    if (d == 2) {
      return(case_d_2(mu, gamma))
    } else {
      for (i in 1:d){
        # Computing eta
        eta <- calculus_of_eta(i, mu, gamma)
        # Computing R
        R <- calculus_of_R(i, gamma)
        expectedvalue <- expectedvalue + exp(mu[i] + gamma[i, i]) * pmvnorm(mean = rep(0, d - 1), corr = R, 
                                                                            lower = rep(-Inf, d - 1), 
                                                                            upper = eta)
      }
      return(expectedvalue)
    }
  }
}

mass_of_scenario <- function(d, r, list_of_k_matrices, A) # A is a (d x k) matrix and Sigma is a (dxd) matrix
{
  pi<-numeric(r)
  for (k in 1:r) {
    dimension <- length(which(A[, k] > 0))
    mu <- log(r * A[which(A[, k] > 0), k]) - (1 / 2) * diag(list_of_k_matrices[[k]])[which(A[, k] > 0)]
    Sigma <- list_of_k_matrices[[k]][which(A[, k] > 0), which(A[, k] > 0)]
    pi[k] <- calculus_d(dimension, Sigma, mu)
  }
  pi <- pi / sum(pi)
  return(pi)
}
 
HR_generator<-function(i,j, column, sign, covariance_matrix, r){
  # Formula with j= (i of the article).
  # covariance matrix= (sigma tilde of the article)
  # HR gives the (-(1/2)*diag(covariance matrix))
  Mu<-log(r * column)-(1 / 2) * diag(covariance_matrix)[sign]  + covariance_matrix[sign, i]
  T<-mvrnorm(1, mu = Mu, 
             Sigma = covariance_matrix[sign, sign])
  return(T)
} 
normalize_group <- function(v, group_size) {
  v_normalized <- sapply(seq(1, length(v), by = group_size),
                         function(i) {
                           group <- v[i:(i + group_size - 1)]
                           group_sum <- sum(group)
                           
                           if (group_sum == 0) {
                             return(group)  # Return the original group
                             #if all values are zero
                           } else {
                             return(group / group_sum)  # Otherwise, normalize
                           }
                         })
  
  return(as.vector(v_normalized))
}
W_calculus<-function(k, num_class, X, grid, q){
  W_train<-list()
  W_test<-list()
  N<-nrow(X)
  for(class_k in 1: num_class) {
    low <- ((class_k-1)*(N/num_class)+1)
    up <- (class_k*(N/num_class))
    X_train=X[-c(low:up),]
    X_test=X[c(low:up),]
    R_train<-apply(X_train,2,rank)
    R_test<-apply(X_test,2,rank)
    W_train[[class_k]] <- sapply(c(1:q), function(m)
      stdfEmp(R_train, k = (1-1/num_class)*k, grid[m,]))
    W_test[[class_k]]  <- sapply(c(1:q), function(m)
      stdfEmp(R_test, k = (1/num_class)*k, grid[m,]))
  }
  return(list("train"=W_train, "test"=W_test))
}
shuffleCols <- function(start, A){
  r <- ncol(A)
  k <- ncol(start)
  perms <- permutations(n = k, r = r)
  temp <- apply(perms, 1, function(j) sum( (start[,j] - A)^(2)  ))
  indx <- which(temp == min(temp))
  firstcols <- perms[indx,]
  lastcols <- setdiff(c(1:k),firstcols)
  return(start[,c(firstcols,lastcols)])
}
construct_symmetric_matrix <- function(vec) {
  # Determine d from the length of vec
  d <- (1 + sqrt(1 + 8 * length(vec))) / 2
  if (d != floor(d)) stop("Vector length is incorrect
                          for a symmetric matrix")
  d <- as.integer(d)
  
  # Initialize d x d matrix with 1s on the diagonal
  mat <- diag(1, d, d)
  
  # Fill the upper triangle row by row
  index <- 1
  for (i in 1:(d-1)) {
    for (j in (i+1):d) {
      mat[i, j] <- vec[index]
      index <- index + 1
    }
  }
  
  # Make the matrix symmetric
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  return(mat)
}
cross_validation<-function(d , r, grid, lambda, num_col = NULL, start , type = c("SSR_row_HR", "SSR_row_log"), p , w , num_class=10){
  print(lambda)
  # Algorithm 2 in the paper.
  
  # Return cross-validation score -----------------------------------------------
  
  #start <- c(start, 0.5)
  w_train <- w$train
  w_test <- w$test
  q <- nrow(grid)
  if(is.null(num_col)){num_col <- r} #if not specified,
  #use the correct number of columns
  l <- d * num_col
  CV <- vector(length = num_class)
  for (class_k in 1:num_class){
    
    # Calls optimisation in each class. ---------------------------------------
    optimizer_minus_class_k <- param_estim(d = d , r = r ,
                                           grid = grid , lambda = lambda ,
                                           num_col = num_col , start = start  , type = type ,
                                           p = p  ,  w = w_train[[class_k]] )
    
    v_A <- as.vector( t(optimizer_minus_class_k$pls_matrix) )
    v_alpha <- optimizer_minus_class_k$pls_dep
    if(type == "SSR_row_log"){
      
      # Function loss of the model -------------------------------------------------------
      CV[class_k] <- .C( type , as.double(p) , as.double(0) ,
                         as.double(v_A) , as.integer(d),
                         as.integer(num_col) , as.integer(q) ,
                         as.double(rep(v_alpha, num_col)) ,
                         as.double(w_test[[class_k]]) ,
                         as.double(c(t(grid))) , R = double(1))$R
    }
    if(type == "SSR_row_HR"){
      CV[class_k] <- .C( type , as.double(p) , as.double(0) ,
                         as.double(v_A) , as.integer(d),
                         as.integer(num_col) , as.integer(q) ,
                         as.double(v_alpha) ,
                         as.double(w_test[[class_k]]) ,
                         as.double(c(t(grid))) , R = double(1))$R
    }
  }
  return(mean(CV))
}

