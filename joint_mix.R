library('mvtnorm')
library('statmod')
library('flexmix')

# this is a modification of the blockGLasso function in the BayesianGLasso package
blockGLasso_mix <- function (X, z, mu, mu0, K, 
		lambdaPriora = 1, 	#r	these are the notations in Perrakis
    		lambdaPriorb = 1/10,    #s	these are the notations in Perrakis
    		Sigma,			# current value of variance matrices
    		verbose = FALSE
    		){	
#	X: 	n x p data matrix
#	z: 	n vector of allocation variables
#	mu:	p x K matrix of means
#	Sigma:  p x p x K array of variance-covariance per component
#	mu0:	p vector of prior mean
#	K:	number of components
#	lambdaPriora: Shrinkage hyperparameter (lambda) gamma distribution
#          shape
#
#	lambdaPriorb: Shrinkage hyperparameter (lambda) gamma distribution
#          scale
	Psi <- numeric(K)
	p <- dim(X)[2]
	indMat <- matrix(1:p^2, ncol = p, nrow = p)
	perms <- matrix(NA, nrow = p - 1, ncol = p)
	permInt <- 1:p
	for (i in 1:ncol(perms)) {
		perms[, i] <- permInt[-i]
	}
	Sigma_New <- Sigma
	Omega_New <- Sigma
	for(k in 1:K){
		ind <- which(z == k)
		n_k <- length(ind)
		S_k <- (mu[,k] - mu0)%*% t(mu[,k] - mu0)
		if(n_k > 0){
			X_k <- X[ind, ] - matrix(mu[,k], n_k, p, byrow = TRUE)		
			S_k <- S_k + t(X_k) %*% X_k
		}
#		Omega_k <- MASS::ginv(Sigma_New[,,k])
		tmp <- chol(Sigma_New[,,k])
		Omega_k <- chol2inv(tmp)
#		print(det(Omega_k))
		# step 13 **************************************************************
		#	lambda corresponds to \psi_k
		lambdaPosta <- (lambdaPriora + (p * (p + 1)/2))	
	        lambdaPostb <- (lambdaPriorb + sum(abs(c(Omega_k)))/2)
		lambda <- stats::rgamma(1, shape = lambdaPosta, scale = 1/lambdaPostb)	        
		Psi[k] <- lambda
		#***********************************************************************
		
#		 step 12 **************************************************************		
		Omega_kTemp <- Omega_k[lower.tri(Omega_k)]
		rinvgaussFun <- function(x) {
		    x <- ifelse(x < 1e-12, 1e-12, x)
		    return(statmod::rinvgauss(n = 1, mean = x, shape = lambda^2))
		}
		tau <- matrix(NA, nrow = p, ncol = p)		
		tau[lower.tri(tau)] <- 1/sapply(sqrt(lambda^2/(Omega_kTemp^2)), 
		    rinvgaussFun)
		tau[upper.tri(tau)] <- t(tau)[upper.tri(t(tau))]
		#***********************************************************************		
		for (i in 1:p) {
			tauI <- as.matrix(tau[perms[, i], i])
			Sigma11 <- Sigma_New[perms[, i], perms[, i], k]
			Sigma12 <- Sigma_New[perms[, i], i, k]
			S21 <- S_k[i, perms[, i]]
			Omega_k11inv <- Sigma11 - Sigma12 %*% t(Sigma12)/Sigma_New[i, i, k]
			Ci <- (S_k[i, i] + lambda) * Omega_k11inv + diag(1/tauI)
#			print(det(Ci))
			CiChol <- 1
			try( CiChol <- chol(Ci) )
			if(length(CiChol) > 1){
				mui <- solve(-Ci, S_k[perms[, i], i])
				beta <- mui + solve(CiChol, stats::rnorm(p - 1))
				Omega_k[perms[, i], i] <- beta
				Omega_k[i, perms[, i]] <- beta
				gamm <- stats::rgamma(n = 1, shape = n_k/2 + 1, rate = (S_k[1, 1] + lambda)/2)
				Omega_k[i, i] <- gamm + t(beta) %*% Omega_k11inv %*% beta
				Omega_kInvTemp <- Omega_k11inv %*% beta
				Sigma_New[perms[, i], perms[, i], k] <- Omega_k11inv + (Omega_kInvTemp %*% t(Omega_kInvTemp))/gamm
				Sigma_New[perms[, i], i, k] <- Sigma_New[i, perms[, i], k] <- (-Omega_kInvTemp/gamm)
				Sigma_New[i, i, k] <- 1/gamm
				Omega_New[,,k] <- Omega_k
			}else{
				Sigma_New[i, i, k] <- 1e+6 + 1
			}
		}
		if(max(Sigma_New[,,k]) > 1e+6 ){
#			if(verbose){
#				cat(paste0('[WARNING] generated Sigma matrix for component: ', k, ' too large.'),'\n')
#				cat(paste0('        the number of assigned observations is: ', length(ind), '.'),'\n')
#			}
			Sigma_New[,,k] = matrix(0,p,p)
			diag(Sigma_New[,,k]) <- runif(1)*diag(var(X))
		}
		tmp <- chol(Sigma_New[,,k])		#these 
		Omega_New[,,k] <- chol2inv(tmp)	# are added by me in order to reduce numerical errors
		
	} 
	results <- vector('list', length = 3)
	results[[1]] <- Sigma_New
	results[[2]] <- Omega_New
	results[[3]] <- Psi
	names(results) <- c('Sigma', 'Omega', 'Psi')
	#Sigma is the variance-covariance per component
	#Omega are the corresponding inverses
	return(results)   		
    		
 }





# gibbs step for Z
update_z <- function (pr, mu, alpha, beta, sigma2, Omega, x_data, y_data){
	#	pr: current value of mixing proportions (K-vector)
	#	mu: current value of mu (p x p x K)
	#	alpha: currengt value of alpha (K-vector)
	#	beta: current value of beta (p x K)
	#	sigma2: current value of sigma2 (K-vector)
	#	Omega: current value of inv Sigma (p x p x K)
	#	x_data: X data matrix (n x p)
	#	y_data: y data (n-vector)
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	K <- length(pr)
	probs <- array(data = 0, dim = c(n, K))
	logL <- 0
	for (k in 1:K) {
#		print(det(Omega[,,k]))
		chol_Omega <- chol(Omega[,,k])
		log_det_term <- sum(log(diag(chol_Omega))) # note that 2 cancels out with 0.5
		center_x <- x_data - matrix(mu[,k], nrow = n, ncol = p,  byrow = TRUE)
		probs[, k] <- log(pr[k]) - 0.5 * apply(center_x, 1, 
			function(tmp) {
			  return(as.numeric(t(tmp) %*% Omega[,,k] %*% tmp))
			}) + log_det_term # the last one corresponds to 0.5 * log(det(Omega[,,k])) 
		probs[, k] <- probs[, k] - 0.5 * log(sigma2[k]) - 0.5 * ((y_data - alpha[k] - t(beta[,k]) %*% t(x_data))^2)/sigma2[k]
		logL <- logL + exp(probs[,k])
	}
	logL <- -n*(p+1)*log(2*pi)/2 + sum(log(logL))
	probs <- array(t(apply(probs, 1, function(tmp) {
		return(exp(tmp - max(tmp)))
		})), dim = c(n, K))
	z <- apply(probs, 1, function(tmp) {
		if (anyNA(tmp)) {
		    tmp <- rep(1, K)
		}
		return(sample(K, 1, prob = tmp))
	})
	results <- vector("list", length = 3)
	names(results) <- c("probs", "z", "logL")
	results[[1]] <- probs
	results[[2]] <- z
	results[[3]] <- logL
	return(results)
}



# proposal M1
prop_M1 <- function (pr, mu, alpha, beta, sigma2, Omega, x_data, y_data, z){
	#	pr: current value of mixing proportions (K-vector)
	#	mu: current value of mu (p x p x K)
	#	alpha: currengt value of alpha (K-vector)
	#	beta: current value of beta (p x K)
	#	sigma2: current value of sigma2 (K-vector)
	#	Omega: current value of inv Sigma (p x p x K)
	#	x_data: X data matrix (n x p)
	#	y_data: y data (n-vector)
	n <- dim(x_data)[1]
	p <- dim(x_data)[2]
	K <- length(pr)

	n_j <- table(z)
	comps <- as.numeric(sample(names(n_j), 2))
	p1 <- rbeta(1, 1, 1)
	z_prop <- z
	ind <- vector('list', length = 2)
	ind[[1]] <- which(z == comps[1])
	ind[[2]] <- which(z == comps[2])	
	z_prop[c(ind[[1]], ind[[2]])] <- sample(comps, length(ind[[1]]) + length(ind[[2]]), replace = TRUE, prob = c(p1, 1-p1))

	ind_prop <- vector('list', length = 2)
	ind_prop[[1]] <- which(z_prop == comps[1])
	ind_prop[[2]] <- which(z_prop == comps[2])	

	probs <- array(data = 0, dim = c(n, K))
	ar <- 0
	j <- 0
	for (k in comps) {
		j <- j+1
		chol_Omega <- chol(Omega[,,k])	
		log_det_term <- sum(log(diag(chol_Omega))) # note that 2 cancels out with 0.5			
		center_x <- x_data[ind[[j]], ] - matrix(mu[,k], nrow = length(ind[[j]]), ncol = p,  byrow = TRUE)
		probs1 <- - 0.5 * apply(center_x, 1, 
			function(tmp) {
			  return(as.numeric(t(tmp) %*% Omega[,,k] %*% tmp))
			}) + log_det_term # the last one corresponds to 0.5 * log(det(Omega[,,k])) 
		mm <- matrix(x_data[ind[[j]],], nrow = length(ind[[j]]), ncol = p) 			
		probs1 <- probs1 - 0.5 * log(sigma2[k]) - 0.5 * ((y_data[ind[[j]]] - alpha[k] - t(beta[,k]) %*% t(mm))^2)/sigma2[k]
		ar <- ar - sum(probs1)

		if(length(ind_prop[[j]]) > 0){
		center_x_new <- x_data[ind_prop[[j]], ] - matrix(mu[,k], nrow = length(ind_prop[[j]]), ncol = p,  byrow = TRUE)
		probs2 <- - 0.5 * apply(center_x_new, 1, 
			function(tmp) {
			  return(as.numeric(t(tmp) %*% Omega[,,k] %*% tmp))
			}) + log_det_term # the last one corresponds to 0.5 * log(det(Omega[,,k]))
		mm <- matrix(x_data[ind_prop[[j]],], nrow = length(ind_prop[[j]]), ncol = p) 
		probs2 <- probs2 - 0.5 * log(sigma2[k]) - 0.5 * ((y_data[ind_prop[[j]]] - alpha[k] - t(beta[,k]) %*% t(mm))^2)/sigma2[k]
		}else{probs2 = 0}
		ar <- ar + sum(probs2)
	}

	results <- vector("list", length = 2)
	names(results) <- c("ar", "z_prop")
	results[[1]] <- ar
	results[[2]] <- z_prop
	return(results)
}


myDirichlet <- function (alpha){
    k <- length(alpha)
    theta <- rgamma(k, shape = alpha, rate = 1)
    return(theta/sum(theta))
}


gibbs_lasso_shrinkage_sampler <- function(x_data, y_data, mcmc_iter = 2000, thin = 10, K,
			alpha_dirichlet = 1, # dirichlet prior concentration parameter (common for all clusters)
			r = 1,
			s = 10^{-1},
			sigma2_alpha = 10^3,
			a = 2.01,
			b = 1.01,
			verbose = FALSE,
			initial_values = NULL,
			overfitting_start = FALSE,
			warm_up = 500,
			initial_warm_up = 500
		){
	if(mcmc_iter %% thin != 0){stop("`mcmc_iter` should be a multiple of `thin`.")}
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	# fixed prior parameters
	m0 <- numeric(p)
	alpha_dirichlet <- rep(alpha_dirichlet, K)
	if(overfitting_start){
		mcmc_iter = warm_up + mcmc_iter
		alpha_dirichlet_overfitting = rep(100, K)
	}else{
		warm_up = 0
	}	
	# define objects that store mcmc output
	pr_sim <- matrix(NA, mcmc_iter/thin, K)
	mu_sim <- array(NA, c(mcmc_iter/thin, p, K))
	alpha_sim <- matrix(NA, mcmc_iter/thin, K)
	beta_sim <- array(NA, c(mcmc_iter/thin, p, K))
	sigma2_sim <- matrix(NA, mcmc_iter/thin, K)
	Sigma_sim <- array(data = NA, dim = c(mcmc_iter/thin, p, p, K))
	delta_sim <- matrix(NA, mcmc_iter/thin, K)
	tau2_sim <- array(NA, c(mcmc_iter/thin, p, K))
	lambda_sim <- matrix(NA, mcmc_iter/thin, K)
	z_sim <- matrix(NA, mcmc_iter/thin, n)
	K_sim <- numeric(mcmc_iter/thin)
	logL <- numeric(mcmc_iter/thin)
	#=======================================================
	# initialization
	#=======================================================
	# main parameters (random start)
	pr <- runif(K)
	pr <- pr/sum(pr)
	mu <- matrix(rnorm(p*K), p, K)
	alpha <- rnorm(K)
	beta <- matrix(rnorm(p*K), p, K)
	sigma2 <- rgamma(K, shape = 1, rate = 0.5)
	Omega <- Sigma <- array(data = 0, dim = c(p,p,K))
	for(k in 1:K){
		diag(Sigma[,,k]) <- diag(var(x_data))  #cov(x_data)
		Omega[,,k] <- solve(Sigma[,,k])
	}
	# hyper-parameters (random start)
	delta <- rgamma(K, shape = 0.5, rate = 0.5)
	tau2 <- matrix(rgamma(p*K, shape = 1, rate = 1), p, K)
	lambda <- rexp(K, rate = 1)
	n_k <- numeric(K)
	Z <- matrix(0, n, K)
	diagZ <- array(data = 0, dim = c(n,n,K))
	# Step 1. update latent allocations via the `update_z` function
	up_z <- update_z(pr = pr, mu = mu, alpha = alpha, 
		beta = beta, sigma2 = sigma2, Omega = Omega, 
		x_data = x_data, y_data = y_data)
	z <- up_z$z
	
	if( is.null(initial_values$z)==FALSE ){
	#=======================================================
		 #in case a vector of initial clusters is provided
		for(iter in 1:initial_warm_up){
			z <- initial_values$z	
			K_plus <- length(table(z))
			n_k <- numeric(K)
			Z <- matrix(0, n, K)
			diagZ <- array(data = 0, dim = c(n,n,K))
			for(k in 1:K){	
				ind <- which(z == k)	
				n_k[k] <- length(ind)
				Z[ind,k] <- 1
				diagZ[,,k] <- diag(Z[,k]) 			
				if(iter < initial_warm_up/3){
					if(length(ind) > p-1){
						Sigma[,,k] <- var(x_data[ind,])
					}else{
						Sigma[,,k] <- matrix(0, p, p)
						diag(Sigma[,,k]) <- diag(var(x_data))  #cov(x_data)
					}
					Omega[,,k] <- solve(Sigma[,,k])
				}
			# Step 9. update mu_k
				mu_mean <- (t(x_data) %*% Z[,k] + m0)/(n_k[k] + 1)
				mu_var <- Sigma[,,k]/(n_k[k] + 1)
				if(all(is.finite(mu_var))){
					mu[,k] <- rmvnorm(1, mean = mu_mean, sigma = mu_var)
				}
			# Step 4. update beta_k	
				TAU_k <- diag(tau2[,k])
				TAU_k_inv <- diag(1/tau2[,k])
				A_k <- t(x_data) %*% diagZ[,,k] %*% x_data 
				A_k <- A_k + TAU_k_inv
				chol_A_k <- chol(A_k)
				A_k_inv <- chol2inv(chol_A_k)
				beta_var <- sigma2[k] * A_k_inv
				beta_mean <- A_k_inv %*% t(x_data) %*% diagZ[,,k] %*% (y_data - alpha[k])
				if(all(is.finite(beta_var))){
					beta[,k] <- rmvnorm(1, mean = beta_mean, sigma = beta_var)
				}
			# Step 3. update alpha_k
				if(n_k[k] > 0){
					w <- sigma2_alpha/(sigma2_alpha + sigma2[k]/n_k[k])
					alpha_mean <- w * t(y_data - x_data %*% beta[,k]) %*% Z[,k]/n_k[k]
					alpha_var <- w * sigma2[k]/n_k[k]
				}else{
					alpha_mean <- 0
					alpha_var <- sigma2_alpha
				}
				alpha[k] <- rnorm(1, alpha_mean, sqrt(alpha_var))
			# Step 5. update sigma2_k
				e_k <- y_data - alpha[k] - x_data %*% beta[,k]
				a_new <- a + 0.5 * (n_k[k] + p)
				b_new <- b + 0.5 * (t(e_k) %*% diagZ[,,k] %*% e_k + t(beta[,k]) %*% TAU_k_inv %*% beta[,k])
				sigma2[k] <- 1/rgamma(1, shape = a_new, rate = b_new)
			# Step 6. update tau2_k
				inv_g_mean <- sqrt(sigma2[k]) * lambda[k]/abs(beta[,k])
				inv_g_mean <- ifelse(inv_g_mean < 1e-12, 1e-12, inv_g_mean)
				tau2[,k] <- 1/rinvgauss(n = p, mean = inv_g_mean, shape = lambda[k]^2)
			# Step 7. update lambda_k
				lambda[k] <- sqrt(rgamma(1, shape = p + 0.5, rate = 0.5*( sum(tau2[,k]) + delta[k] )))	
			# Step 8. update delta_k	
				delta[k] <- rgamma(1, shape = 1, rate = 0.5*(lambda[k]^2 + 1))
				
			}
			if(iter > initial_warm_up/3){
				block_lasso_step <- blockGLasso_mix(X = x_data, z = z, mu = mu, mu0 = m0, K = K, 
						lambdaPriora = r, 
				    		lambdaPriorb = s,
				    		Sigma = Sigma,
				    		verbose = verbose
				    		)
				Sigma <-  block_lasso_step$Sigma   		
				Omega <-  block_lasso_step$Omega
			}	
		}


	}
	
	m1_ar <- 0
#---------------------------------------------------------------------------------	
	# the main loop start here
#---------------------------------------------------------------------------------	
	for (iter in 1:mcmc_iter){
		n_k <- numeric(K)
		Z <- matrix(0, n, K)
		diagZ <- array(data = 0, dim = c(n,n,K))

		for(k in 1:K){
		# Step 2. define allocation related quantities
			ind <- which(z == k)
			n_k[k] <- length(ind)
			Z[ind,k] <- 1
			diagZ[,,k] <- diag(Z[,k]) 
		# Step 3. update alpha_k
			if(n_k[k] > 0){
				w <- sigma2_alpha/(sigma2_alpha + sigma2[k]/n_k[k])
				alpha_mean <- w * t(y_data - x_data %*% beta[,k]) %*% Z[,k]/n_k[k]
				alpha_var <- w * sigma2[k]/n_k[k]
			}else{
				alpha_mean <- 0
				alpha_var <- sigma2_alpha
			}
			alpha[k] <- rnorm(1, alpha_mean, sqrt(alpha_var))
		# Step 4. update beta_k	
			TAU_k <- diag(tau2[,k])
			TAU_k_inv <- diag(1/tau2[,k])
			A_k <- t(x_data) %*% diagZ[,,k] %*% x_data 
			A_k <- A_k + TAU_k_inv
	#		A_k_inv <- solve(A_k)
			chol_A_k <- chol(A_k)
			A_k_inv <- chol2inv(chol_A_k)
			beta_var <- sigma2[k] * A_k_inv
			beta_mean <- A_k_inv %*% t(x_data) %*% diagZ[,,k] %*% (y_data - alpha[k])
			if(all(is.finite(beta_var))){
				beta[,k] <- rmvnorm(1, mean = beta_mean, sigma = beta_var)
			}

		# Step 5. update sigma2_k
			e_k <- y_data - alpha[k] - x_data %*% beta[,k]
			a_new <- a + 0.5 * (n_k[k] + p)
			b_new <- b + 0.5 * (t(e_k) %*% diagZ[,,k] %*% e_k + t(beta[,k]) %*% TAU_k_inv %*% beta[,k])
			sigma2[k] <- 1/rgamma(1, shape = a_new, rate = b_new)
		# Step 6. update tau2_k
			inv_g_mean <- sqrt(sigma2[k]) * lambda[k]/abs(beta[,k])
			inv_g_mean <- ifelse(inv_g_mean < 1e-12, 1e-12, inv_g_mean)
			tau2[,k] <- 1/rinvgauss(n = p, mean = inv_g_mean, shape = lambda[k]^2)
		# Step 7. update lambda_k
			lambda[k] <- sqrt(rgamma(1, shape = p + 0.5, rate = 0.5*( sum(tau2[,k]) + delta[k] )))	
		# Step 8. update delta_k	
			delta[k] <- rgamma(1, shape = 1, rate = 0.5*(lambda[k]^2 + 1))
		# Step 9. update mu_k
			mu_mean <- (t(x_data) %*% Z[,k] + m0)/(n_k[k] + 1)
			mu_var <- Sigma[,,k]/(n_k[k] + 1)
			if(all(is.finite(mu_var))){
				mu[,k] <- rmvnorm(1, mean = mu_mean, sigma = mu_var)
			}

		}

		# Steps 10, 11, 12, 13. update Omega via the `blockGLasso_mix` function


		block_lasso_step <- blockGLasso_mix(X = x_data, z = z, mu = mu, mu0 = m0, K = K, 
				lambdaPriora = r, 
		    		lambdaPriorb = s,
		    		Sigma = Sigma,
		    		verbose = verbose
		    		)
		Sigma <-  block_lasso_step$Sigma   		
		Omega <-  block_lasso_step$Omega
		  
		# Step 14: update mixing proportions
		if(iter > warm_up){
			pr <- myDirichlet(alpha_dirichlet + n_k)
		}else{
			pr <- myDirichlet(alpha_dirichlet_overfitting + n_k)
		}

		# Step 1. update latent allocations via the `update_z` function
		up_z <- update_z(pr = pr, mu = mu, alpha = alpha, 
			beta = beta, sigma2 = sigma2, Omega = Omega, 
			x_data = x_data, y_data = y_data)
		z <- up_z$z


		# try this
#		m1 <- prop_M1(pr = pr, mu = mu, alpha = alpha, 
#			beta = beta, sigma2 = sigma2, Omega = Omega, 
#			x_data = x_data, y_data = y_data, z = z)
#		u <- log(runif(1))
#		if(u < m1$ar){
#			z <- m1$z_prop
#			m1_ar <- m1_ar + 1
#		}




		
		# store output every thin iterations
		if(iter %% thin == 0){
			pr_sim[iter/thin, ] <- pr
			mu_sim[iter/thin,,] <- mu
			alpha_sim[iter/thin,] <- alpha
			beta_sim[iter/thin, , ] <- beta
			sigma2_sim[iter/thin, ] <- sigma2
			Sigma_sim[iter/thin,,,] <- Sigma
			delta_sim[iter/thin,] <- delta
			tau2_sim[iter/thin,,] <- tau2
			lambda_sim[iter/thin,] <- lambda
			z_sim[iter/thin,] <- z
			K_sim[iter/thin] <- length(table(z))
			logL[iter/thin] <- up_z$logL
		}	

		if(iter %% 100 == 0){
#			if(iter == 100){
#				matplot(mu, type = 'l', lwd = 2)
#			}else{
#				matplot(mu, type = 'l', lwd = 2, add = TRUE)
#			}
#			pairs(cbind(y_data, x_data), col = z)\
			mymain = paste0('clusters at iteration: ', iter)
			par(mfrow = c(1,3))
			plot(logL[1:(iter/thin)], type = 'l', xlab = 'iteration (thinned)', ylab ='log-likelihood' )
			plot(K_sim[1:(iter/thin)], type = 'l', xlab = 'iteration (thinned)', ylab ='number of clusters' )
			matplot(t(cbind(y_data,x_data)),col = z, type = 'l', xaxt = 'n', main = mymain)
			axis(side = 1, at = 1:(p+1), labels = c('y',paste0('x',1:p)))
			legend('bottomleft', 
				paste0('component: ', names(table(z)), ', n', names(table(z)), ' = ', table(z)), 
				col = as.numeric(names(table(z))), lty = 1)
			cat(paste('iteration: ', iter), '\n')
			print(table(z))
			cat(paste0('Move 1: ',100*m1_ar/iter),'\n')			
		}
	}
	n_parameters <- K *(2 + 2*p + p*(p+1)/2) + K - 1 #+ K*(2+p)
	BIC <- -2*max(logL) + n_parameters * log(n)
	AIC <- -2*max(logL) + n_parameters * 2
	burn <- floor(0.1*mcmc_iter/thin)
	DIC_2 <-  -4 * mean(logL[-(1:burn)]) + 2 * max(logL[-(1:burn)]) 
#	compute ICL
	ICL_BIC <- BIC
	if(K > 1){
		ind_max <- which.max(logL)
		for(k in 1:K){
			tmp <- chol(Sigma_sim[ind_max,,,k])
			Omega[,,k] <- chol2inv(tmp)
		}
		up_z <- update_z(pr = pr_sim[ind_max,], mu = mu_sim[ind_max,,], alpha = alpha_sim[ind_max,], 
			beta = beta_sim[ind_max,,], sigma2 = sigma2_sim[ind_max,], Omega = Omega, 
			x_data = x_data, y_data = y_data)
		nz <- up_z$probs
		thresh <- -744
		exp_thresh <- exp(thresh)
		ind <- 1:K
		for (i in 1:n) {
			index <- ind[nz[i, ] < exp_thresh]
			if(length(index) > 0){
				nz[i, index] <- rep(exp_thresh, length(index))
				nz[i, ] <- nz[i, ]/sum(nz[i, ])
			}
		}
		ICL_BIC <- BIC - 2 * sum(nz * log(nz))
	}
#################################################3
	results <- vector('list', length = 15)
	results[[1]] <- pr_sim
	results[[2]] <- mu_sim
	results[[3]] <- alpha_sim
	results[[4]] <- beta_sim
	results[[5]] <- sigma2_sim
	results[[6]] <- Sigma_sim
	results[[7]] <- delta_sim
	results[[8]] <- tau2_sim
	results[[9]] <- lambda_sim		
	results[[10]] <- z_sim			
	results[[11]] <- logL		
	results[[12]] <- BIC
	results[[13]] <- ICL_BIC
	results[[14]] <- AIC
	results[[15]] <- DIC_2	
	names(results) <- c('pr', 'mu', 'alpha', 'beta', 'sigma2', 'SIGMA', 
	'delta', 'tau2', 'lambda', 'z', 'logL', 'BIC', 'ICL', 'AIC', 'DIC_2')
	return(results)
}


bnb_prior <- function(K, a_lambda = 1, a_pi = 4, b_pi = 3){
	return(lgamma(a_lambda + K - 1) + lbeta(a_lambda + a_pi, K - 1 + b_pi) - 
		- lgamma(a_lambda) - lgamma(K) - lbeta(a_pi, b_pi))
	
}

K_conditional <- function(K_max = 100, log_prior_K, z, a){
	# z should be named with entries 1, ..., K_plus!
	myNames <- as.numeric(names(table(z)))
	K_plus <- length(table(z))
	z_check <- sum(myNames > K_plus)
	if(z_check > 0){
		stop('not valid z')
	}
	K = K_plus:K_max
	n_k <- table(z)
	log_prob <- log_prior_K[K_plus:K_max] + K_plus * log(a) + lgamma(K+1) - K_plus * log(K) - lgamma(K - K_plus + 1)
	for(k in 1:K_plus){
		log_prob <- log_prob + lgamma(n_k[k] + a/K) - lgamma(1 + a/K)
	}
	full_prob <- numeric(K_max)
	for(k in K){
		full_prob[k] <- 1/sum(exp(log_prob - log_prob[k - K_plus + 1]))
	}
	return(full_prob)	
}

metropolis_move_for_a <- function(z, K, a_old, s2_a, nu_l = 6, nu_r = 3){
	myNames <- as.numeric(names(table(z)))
	K_plus <- length(table(z))
	n <- length(z)
	z_check <- sum(myNames > K_plus)
	if(z_check > 0){
		stop('not valid z')
	}
	n_k <- numeric(K)
	for(k in 1:K){
		ind <- which(z == k)
		n_k[k] <- length(ind)
	}
	mh <- 0
	a_new <- rlnorm(1, log(a_old), s2_a)
	ar <- df(a_new, df1 = nu_l, df2 = nu_r, log = TRUE) - df(a_old, df1 = nu_l, df2 = nu_r, log = TRUE) + 
		K_plus * (log(a_new) - log(a_old)) + lgamma(a_new) - lgamma(a_old) - lgamma(n+a_new) + lgamma(n+a_old) 
	for(k in 1:K_plus){
		ar <- ar + lgamma(n_k[k] + a_new/K) - lgamma(n_k[k] + a_old/K) - lgamma(1+a_new/K) + lgamma(1+a_old/K)
	}
	ar <- ar + log(a_new) - log(a_old)
	if(is.na(ar)){ar <- -70000; a_old = 1*runif(1)}
	if(log(runif(1)) < ar){
		mh <- 1
	}else{
		a_new <- a_old
	}
	results <- vector('list', length = 2)
	results[[1]] <- a_new
	results[[2]] <- mh
	names(results) <- c('a', 'accepted')
	return(results)
}


log_dirichlet_pdf <- function(alpha, pr){
	K <- length(pr)
	norm_const <- K*lgamma(alpha) - lgamma(K*alpha)
	lpd <- (alpha - 1)*sum(log(pr)) - norm_const
	return(lpd)
}



initialize_sampler <- function(x_data, y_data, mcmc_iter = 30, thin = 1, K_max,
			r = 1,			
			s = 10^{-2},
			sigma2_alpha = 10^3,
			a = 0.01,
			b = 0.01,
			verbose = FALSE,
			initial_values = NULL,
			warm_up = 5,
			initial_warm_up = 2,
			nu_l = 6, nu_r = 3,	#dirichlet alpha prior parameters
			a_lambda = 1, a_pi = 4, b_pi = 3, #bnb prior parameters
			s2_a = 5, nRuns = 10, burn = 20, flexmix_runs = 30){

	myData <- as.data.frame(cbind(y_data, x_data))
	colnames(myData)[1] <- 'y'
	colnames(myData)[-1] <- paste0('x',1:p)

	kvals <- max(c(floor(K_max/2 - 4),2)):floor(K_max/2 + 2)
	k = sample(kvals, 1, replace = TRUE)
	
	init <- tryCatch({
		pissa <- initFlexmix(y ~ ., data = myData, 
			k = k, model = FLXMRglm(family='gaussian'), 
			nrep = flexmix_runs, control = list(minprior = 0.05))
		}, error = function(pissa) {
		pissa <- initFlexmix(y ~ ., data = myData, 
			k = k, model = FLXMRglm(family='gaussian'), 
			nrep = flexmix_runs, control = list(minprior = 0.1))
		return(pissa)
		}, finally = pissa)
	initial_values <- vector('list', length = 1)

	initial_values$z <- init@cluster
	
        myRuns <- telescoping_sampler(x_data = x_data, y_data = y_data, mcmc_iter = mcmc_iter, thin = thin, K_max = K_max,
                        r = r,                  
                        s = s,
                        sigma2_alpha = sigma2_alpha,
                        a = a,
                        b = b,
                        initial_values = initial_values,
                        nu_l = nu_l, nu_r = nu_r,     #dirichlet alpha prior parameters
                        a_lambda = a_lambda, a_pi = a_pi, b_pi = b_pi, #bnb prior parameters
                        s2_a = s2_a, # scale for mh proposal
                        warm_up = warm_up, 
                        initial_warm_up = initial_warm_up
                )
#	logP <- mean(myRuns$logP[-(1:burn)])
	logP <- myRuns$logP[mcmc_iter]
	z <- myRuns$z[mcmc_iter,]
	logP_best <- logP
	cat('iter: ', 1, ', logP = ', logP_best,'\n')
	print(table(z))

	for(iter in 2:nRuns){
		k = sample(kvals, 1, replace = TRUE)
		init <- tryCatch({
			pissa <- initFlexmix(y ~ ., data = myData, 
				k = k, model = FLXMRglm(family='gaussian'), 
				nrep = flexmix_runs, control = list(minprior = 0.05))
			}, error = function(pissa) {
			pissa <- initFlexmix(y ~ ., data = myData, 
				k = k, model = FLXMRglm(family='gaussian'), 
				nrep = flexmix_runs, control = list(minprior = 0.1))
			return(pissa)
			}, finally = pissa)
		initial_values <- vector('list', length = 1)

		initial_values$z <- init@cluster
		myRuns <- telescoping_sampler(x_data = x_data, y_data = y_data, mcmc_iter = mcmc_iter, thin = thin, K_max = K_max,
			        r = r,                  
			        s = s,
			        sigma2_alpha = sigma2_alpha,
			        a = a,
			        b = b,
			        initial_values = initial_values,
			        nu_l = nu_l, nu_r = nu_r,     #dirichlet alpha prior parameters
			        a_lambda = a_lambda, a_pi = a_pi, b_pi = b_pi, #bnb prior parameters
			        s2_a = s2_a, # scale for mh proposal
			        warm_up = warm_up, 
			        initial_warm_up = initial_warm_up
			)
		logP <- myRuns$logP[mcmc_iter]
		z_overfitting <- myRuns$z[mcmc_iter,]


		if(logP > logP_best){
			z <- z_overfitting
			logP_best <- logP
		}
		cat('iter: ', iter, ', logP = ', logP_best,'\n')
		print(table(z))

	}
	results <- vector('list', length = 2)
	results[[1]] <- logP
	results[[2]] <- z
	names(results) <- c('logP', 'z')
	return(results)		
}

telescoping_sampler <- function(x_data, y_data, mcmc_iter = 2000, thin = 10, K_max,
			r = 1,			
			s = 10^{-1},
			sigma2_alpha = 10^3,
			a = 2.01,
			b = 1.01,
			verbose = FALSE,
			initial_values = NULL,
			warm_up = 500,
			initial_warm_up = 1000,
#######################################################################################3
			nu_l = 6, nu_r = 3,	#dirichlet alpha prior parameters
			a_lambda = 1, a_pi = 4, b_pi = 3, #bnb prior parameters
			s2_a = 5 # scale for mh proposal
		){
		
#################################################3
	if(mcmc_iter %% thin != 0){stop("`mcmc_iter` should be a multiple of `thin`.")}
	p <- dim(x_data)[2]
	n <- dim(x_data)[1]
	# fixed prior parameters
	m0 <- numeric(p)
	# define objects that store mcmc output
	alpha_dirichlet_sim <- numeric(mcmc_iter/thin)
	pr_sim <- matrix(NA, mcmc_iter/thin, K_max)
	mu_sim <- array(NA, c(mcmc_iter/thin, p, K_max))
	alpha_sim <- matrix(NA, mcmc_iter/thin, K_max)
	beta_sim <- array(NA, c(mcmc_iter/thin, p, K_max))
	sigma2_sim <- matrix(NA, mcmc_iter/thin, K_max)
	Sigma_sim <- array(data = NA, dim = c(mcmc_iter/thin, p, p, K_max))
	delta_sim <- matrix(NA, mcmc_iter/thin, K_max)
	tau2_sim <- array(NA, c(mcmc_iter/thin, p, K_max))
	lambda_sim <- matrix(NA, mcmc_iter/thin, K_max)
	z_sim <- matrix(NA, mcmc_iter/thin, n)
	K_sim <- numeric(mcmc_iter/thin)
	K_plus_sim <- numeric(mcmc_iter/thin)	
	logL <- logP <- numeric(mcmc_iter/thin)
	mh_rate <- 0
	#=======================================================
	# initialization
	#=======================================================
	dvX <- var(x_data)
	# main parameters (random start)
	K = K_max
	pr <-  runif(K)
	pr <- pr/sum(pr)
	mu <- matrix(rnorm(p*K), p, K)
	alpha <- rnorm(K)
	beta <- matrix(rnorm(p*K), p, K)
	sigma2 <- rgamma(K, shape = 1, rate = 0.5)
	Omega <- Sigma <- array(data = 0, dim = c(p,p,K))
	for(k in 1:K){
		diag(Sigma[,,k]) <- diag(var(x_data))  #cov(x_data)
		Omega[,,k] <- solve(Sigma[,,k])  
	}
	# hyper-parameters (random start)
	alpha_dirichlet <- 10*runif(1)
	delta <- rgamma(K, shape = 0.5, rate = 0.5)
	tau2 <- matrix(rgamma(p*K, shape = 1, rate = 1), p, K)
	lambda <- rexp(K, rate = 1)
	n_k <- numeric(K)
	Z <- matrix(0, n, K)
	diagZ <- array(data = 0, dim = c(n,n,K))
	# Step 1. update latent allocations via the `update_z` function
	up_z <- update_z(pr = pr, mu = mu, alpha = alpha, 
		beta = beta, sigma2 = sigma2, Omega = Omega, 
		x_data = x_data, y_data = y_data)
	z <- up_z$z
	
	# Psi stems from blockLasso
	log_posterior <- function(alpha_dirichlet, pr, mu, alpha, beta, sigma2, Omega, Sigma, lambda, logL, z, K, Psi){
		K_plus <- length(table(z))
		p <- dim(beta)[1]
		#inv gamma
		log_sigma2_prior <- sum(a * log(b) - (a+1)*log(sigma2)	- b/sigma2 - lgamma(a))
		# half-cauchy(0,1)
		log_lambda_prior <- -sum(log(1+lambda^2)) + K*log(2/pi)
		#normal(0, sigma2_alpha)
		log_alpha_prior <- sum(dnorm(alpha,0, sqrt(sigma2_alpha), log = TRUE))
		#laplace (0, sigma_k/lambda_k)
		log_beta_prior <- 0
		for(k in 1:K){
			log_beta_prior <- log_beta_prior - p*log(2*sqrt(sigma2[k])/lambda[k]) - sum(abs(beta[,k]))*lambda[k]/sqrt(sigma2[k])
		}
		# Psi prior: 	
		log_psi_prior <- sum(dgamma(Psi, shape = r, rate = s, log = TRUE))
		# (Omega, mu) prior
		log_omega_prior <- 0
		log_mu_prior <- 0
		for(k in 1:K){
			my_mat <- Omega[,,k]
			low_part <- my_mat[lower.tri(my_mat)]
			#diagonal entries
			log_omega_prior <- log_omega_prior + sum(dexp(diag(my_mat), rate = Psi[k]/2, log = TRUE))
			# upper-diagonal entries
			log_omega_prior <- log_omega_prior + 0.5*p*(p-1)*log(Psi[k]/2) - sum(abs(low_part))*Psi[k]
			log_mu_prior <- log_mu_prior + dmvnorm(t(mu[,k]), mean = m0, sigma = Sigma[,,k], log = TRUE)
		}
		logP = logL + df(alpha_dirichlet, nu_l, nu_r, log = TRUE) + 
			bnb_prior(K, a_lambda, a_pi, b_pi) + 
			log_dirichlet_pdf(alpha_dirichlet/K, pr) + 
			log_sigma2_prior + log_lambda_prior + log_beta_prior + log_alpha_prior + 
			log_psi_prior + log_omega_prior + log_mu_prior
		return(logP)
	}

	
	
	if( is.null(initial_values$z)==FALSE ){
	#=======================================================
	 #in case a vector of initial clusters is provided
	 #in case a vector of initial clusters is provided
		z <- initial_values$z	
		K_plus <- length(table(z))
		n_k <- numeric(K)
		Z <- matrix(0, n, K)
		diagZ <- array(data = 0, dim = c(n,n,K))
		for(k in 1:K){	
			ind <- which(z == k)	
			n_k[k] <- length(ind)
			Z[ind,k] <- 1
			diagZ[,,k] <- diag(Z[,k]) 			
			if(length(ind) > p-1){
				Sigma[,,k] <- var(x_data[ind,])
			}else{
				Sigma[,,k] <- matrix(0, p, p)
				diag(Sigma[,,k]) <- diag(var(x_data))  #cov(x_data)
			}
			Omega[,,k] <- solve(Sigma[,,k])
		# Step 9. update mu_k
			mu_mean <- (t(x_data) %*% Z[,k] + m0)/(n_k[k] + 1)
			mu_var <- Sigma[,,k]/(n_k[k] + 1)
			if(all(is.finite(mu_var))){
				mu[,k] <- rmvnorm(1, mean = mu_mean, sigma = mu_var)
			}

		# Step 4. update beta_k	
			TAU_k <- diag(tau2[,k])
			TAU_k_inv <- diag(1/tau2[,k])
			A_k <- t(x_data) %*% diagZ[,,k] %*% x_data 
			A_k <- A_k + TAU_k_inv
			chol_A_k <- chol(A_k)
			A_k_inv <- chol2inv(chol_A_k)
			beta_var <- sigma2[k] * A_k_inv
			beta_mean <- A_k_inv %*% t(x_data) %*% diagZ[,,k] %*% (y_data - alpha[k])
			if(all(is.finite(beta_var))){
				beta[,k] <- rmvnorm(1, mean = beta_mean, sigma = beta_var)
			}

		# Step 3. update alpha_k
			if(n_k[k] > 0){
				w <- sigma2_alpha/(sigma2_alpha + sigma2[k]/n_k[k])
				alpha_mean <- w * t(y_data - x_data %*% beta[,k]) %*% Z[,k]/n_k[k]
				alpha_var <- w * sigma2[k]/n_k[k]
			}else{
				alpha_mean <- 0
				alpha_var <- sigma2_alpha
			}
			alpha[k] <- rnorm(1, alpha_mean, sqrt(alpha_var))
		# Step 5. update sigma2_k
			e_k <- y_data - alpha[k] - x_data %*% beta[,k]
			a_new <- a + 0.5 * (n_k[k] + p)
			b_new <- b + 0.5 * (t(e_k) %*% diagZ[,,k] %*% e_k + t(beta[,k]) %*% TAU_k_inv %*% beta[,k])
			sigma2[k] <- 1/rgamma(1, shape = a_new, rate = b_new)
		}
		# prior for K
		log_prior_K <- bnb_prior(1:K_max, a_lambda = a_lambda, a_pi = a_pi, b_pi = b_pi)
		for(iter in 1:initial_warm_up){
			K_plus <- length(table(z))
			z_labels <- as.numeric(names(table(z)))
	###################################################################
	# 	relabeling step so the first K_+ components are not empty
			if(K_plus < K){
				if(max(z) > K_plus){
					z_relabeled <- as.numeric(as.factor(z))
					new_z_labels <- as.numeric(names(table(z_relabeled)))		

					perm <- c(z_labels, (1:K_max)[-c(z_labels)])
					pr <- pr[perm]
					mu <- mu[,perm]
					alpha <- alpha[perm]
					beta <- beta[,perm]
					sigma2 <- sigma2[perm]
					Omega <- Omega[,,perm]
					Sigma <- Sigma[,,perm]
					delta <- delta[perm]
					tau2 <- tau2[,perm]
					lambda <- lambda[perm]
					z <- z_relabeled
				}
			}
	#		print(K_plus)	
	#####################################################################
		
	#		step 3 new values of K and alpha_dirichlet

			k_cond_distribution <- K_conditional(K_max = K_max, log_prior_K = log_prior_K, z = z, a = alpha_dirichlet)
			K <- sample(1:K_max, 1, prob = k_cond_distribution)
	#		print(alpha_dirichlet)

	#		print(s2_a)
			mh_step <- metropolis_move_for_a(z = z, K = K, a_old = alpha_dirichlet, s2_a = s2_a, nu_l = nu_l, nu_r = nu_r)
	#		print(mh_step)
			alpha_dirichlet <- mh_step$a
			n_k <- numeric(K)
			Z <- matrix(0, n, K)
			diagZ <- array(data = 0, dim = c(n,n,K))
	 
			for(k in 1:K_plus){
			# Step 2. define allocation related quantities
				ind <- which(z == k)
				n_k[k] <- length(ind)
				Z[ind,k] <- 1
				diagZ[,,k] <- diag(Z[,k]) 
			# Step 3. update alpha_k
				w <- sigma2_alpha/(sigma2_alpha + sigma2[k]/n_k[k])
				alpha_mean <- w * t(y_data - x_data %*% beta[,k]) %*% Z[,k]/n_k[k]
				alpha_var <- w * sigma2[k]/n_k[k]
				alpha[k] <- rnorm(1, alpha_mean, sqrt(alpha_var))
			# Step 4. update beta_k	
				TAU_k <- diag(tau2[,k])
				TAU_k_inv <- diag(1/tau2[,k])
				A_k <- t(x_data) %*% diagZ[,,k] %*% x_data 
				A_k <- A_k + TAU_k_inv
		#		A_k_inv <- solve(A_k)
				chol_A_k <- chol(A_k)
				A_k_inv <- chol2inv(chol_A_k)
				beta_var <- sigma2[k] * A_k_inv
				beta_mean <- A_k_inv %*% t(x_data) %*% diagZ[,,k] %*% (y_data - alpha[k])
				if(all(is.finite(beta_var))){
					beta[,k] <- rmvnorm(1, mean = beta_mean, sigma = beta_var)
				}
			# Step 5. update sigma2_k
				e_k <- y_data - alpha[k] - x_data %*% beta[,k]
				a_new <- a + 0.5 * (n_k[k] + p)
				b_new <- b + 0.5 * (t(e_k) %*% diagZ[,,k] %*% e_k + t(beta[,k]) %*% TAU_k_inv %*% beta[,k])
				sigma2[k] <- 1/rgamma(1, shape = a_new, rate = b_new)
			# Step 6. update tau2_k
				inv_g_mean <- sqrt(sigma2[k]) * lambda[k]/abs(beta[,k])
				inv_g_mean <- ifelse(inv_g_mean < 1e-12, 1e-12, inv_g_mean)
				tau2[,k] <- 1/rinvgauss(n = p, mean = inv_g_mean, shape = lambda[k]^2)
			# Step 7. update lambda_k
				lambda[k] <- sqrt(rgamma(1, shape = p + 0.5, rate = 0.5*( sum(tau2[,k]) + delta[k] )))	
			# Step 8. update delta_k	
				delta[k] <- rgamma(1, shape = 1, rate = 0.5*(lambda[k]^2 + 1))
			# Step 9. update mu_k
				mu_mean <- (t(x_data) %*% Z[,k] + m0)/(n_k[k] + 1)
				mu_var <- Sigma[,,k]/(n_k[k] + 1)
				if(all(is.finite(mu_var))){
					mu[,k] <- rmvnorm(1, mean = mu_mean, sigma = mu_var)
				}

			}
			if(K_plus < K){
			for(k in (K_plus+1):K){
				alpha_mean <- 0
				alpha_var <- sigma2_alpha
				alpha[k] <- rnorm(1, alpha_mean, sqrt(alpha_var))
				sigma2[k] <- 1/rgamma(1, shape = a, rate = b)
				if(!is.finite(sigma2[k])){sigma2[k] = 1/rgamma(1, shape = 1, rate = 0.1)}			
				if(sigma2[k] > 100){sigma2[k] = 1/rgamma(1, shape = 1, rate = 1)}							
				lambda[k] <- rt(1, df = 1)
				while(lambda[k] < 0){
					lambda[k] <- rt(1, df = 1)
				}
	#			cat(paste0('k = ', k),'\n')
				tau2[,k] <- rexp(p, rate = lambda[k]^2)
				TAU_k <- diag(tau2[,k])
				zar <- sigma2[k] * TAU_k
				zar[!is.finite(zar)] <- 1/rgamma(sum(!is.finite(zar)*1), shape = 1, rate = 0.1)	
				beta[,k] <- rmvnorm(1, mean = rep(0,p), sigma = zar)
				delta[k] <- rgamma(1, shape = 0.5, rate = 0.5)	
				Sigma[,,k] <- 1*runif(1)*dvX	
				mu[,k] <- rmvnorm(1, mean = m0, sigma = Sigma[,,k])
			}
			}


			# Steps 10, 11, 12, 13. update Omega via the `blockGLasso_mix` function



			block_lasso_step <- blockGLasso_mix(X = x_data, z = z, mu = as.matrix(mu[,1:K]), mu0 = m0, K = K, 
					lambdaPriora = r, 
			    		lambdaPriorb = s,
			    		Sigma = array(Sigma[,,1:K],dim=c(p,p,K)),
			    		verbose = verbose
			    		)
			Sigma[,,1:K] <-  block_lasso_step$Sigma   		
			Omega[,,1:K] <-  block_lasso_step$Omega
			  
			# Step 14: update mixing proportions
			pr[1:K] <- myDirichlet(alpha_dirichlet/K + n_k)
		
		}

	}


##################################################################
	#MAIL LOOP STARTS HERE
	# prior for K
	log_prior_K <- bnb_prior(1:K_max, a_lambda = a_lambda, a_pi = a_pi, b_pi = b_pi)
	Psi <- block_lasso_step$Psi  
	up_z <- update_z(pr = pr[1:K], mu = as.matrix(mu[,1:K]), alpha = alpha[1:K], 
		beta = as.matrix(beta[,1:K]), sigma2 = sigma2[1:K], Omega = array(Omega[,,1:K],dim=c(p,p,K)), 
		x_data = x_data, y_data = y_data)
	z <- up_z$z

	for (iter in 1:mcmc_iter){
#		print(iter)

		K_plus <- length(table(z))
		z_labels <- as.numeric(names(table(z)))
###################################################################
# 	relabeling step so the first K_+ components are not empty
		if(K_plus < K){
			if(max(z) > K_plus){
				z_relabeled <- as.numeric(as.factor(z))
				new_z_labels <- as.numeric(names(table(z_relabeled)))		

				perm <- c(z_labels, (1:K_max)[-c(z_labels)])
				pr <- pr[perm]
				mu <- mu[,perm]
				alpha <- alpha[perm]
				beta <- beta[,perm]
				sigma2 <- sigma2[perm]
				Omega <- Omega[,,perm]
				Sigma <- Sigma[,,perm]
				delta <- delta[perm]
				tau2 <- tau2[,perm]
				lambda <- lambda[perm]
				z <- z_relabeled
				Psi <- Psi[perm]
			}
		}


		log_post <- log_posterior(alpha_dirichlet = alpha_dirichlet, pr = pr[1:K], 
			mu = as.matrix(mu[,1:K]), alpha = alpha[1:K], beta = as.matrix(beta[,1:K]), sigma2 = sigma2[1:K], 
			Omega = array(Omega[,,1:K], dim = c(p,p,K)), Sigma = array(Sigma[,,1:K], dim = c(p,p,K)), 
			lambda = lambda[1:K], logL = up_z$logL, 
			z = z, K = K, Psi = Psi) 
		
#		print(K_plus)	
#####################################################################
	
#		step 3 new values of K and alpha_dirichlet

		k_cond_distribution <- K_conditional(K_max = K_max, log_prior_K = log_prior_K, z = z, a = alpha_dirichlet)
		K <- sample(1:K_max, 1, prob = k_cond_distribution)
#		print(alpha_dirichlet)

#		print(s2_a)
		mh_step <- metropolis_move_for_a(z = z, K = K, a_old = alpha_dirichlet, s2_a = s2_a, nu_l = nu_l, nu_r = nu_r)
#		print(mh_step)
		alpha_dirichlet <- mh_step$a
		mh_rate <- mh_rate + mh_step$accepted
		n_k <- numeric(K)
		Z <- matrix(0, n, K)
		diagZ <- array(data = 0, dim = c(n,n,K))
 
		for(k in 1:K_plus){
		# Step 2. define allocation related quantities
			ind <- which(z == k)
			n_k[k] <- length(ind)
			Z[ind,k] <- 1
			diagZ[,,k] <- diag(Z[,k]) 
		# Step 3. update alpha_k
			w <- sigma2_alpha/(sigma2_alpha + sigma2[k]/n_k[k])
			alpha_mean <- w * t(y_data - x_data %*% beta[,k]) %*% Z[,k]/n_k[k]
			alpha_var <- w * sigma2[k]/n_k[k]
			alpha[k] <- rnorm(1, alpha_mean, sqrt(alpha_var))
		# Step 4. update beta_k	
			TAU_k <- diag(tau2[,k])
			TAU_k_inv <- diag(1/tau2[,k])
			A_k <- t(x_data) %*% diagZ[,,k] %*% x_data 
			A_k <- A_k + TAU_k_inv
	#		A_k_inv <- solve(A_k)
			chol_A_k <- chol(A_k)
			A_k_inv <- chol2inv(chol_A_k)
			beta_var <- sigma2[k] * A_k_inv
			beta_mean <- A_k_inv %*% t(x_data) %*% diagZ[,,k] %*% (y_data - alpha[k])
			if(all(is.finite(beta_var))){
				beta[,k] <- rmvnorm(1, mean = beta_mean, sigma = beta_var)
			}
		# Step 5. update sigma2_k
			e_k <- y_data - alpha[k] - x_data %*% beta[,k]
			a_new <- a + 0.5 * (n_k[k] + p)
			b_new <- b + 0.5 * (t(e_k) %*% diagZ[,,k] %*% e_k + t(beta[,k]) %*% TAU_k_inv %*% beta[,k])
			sigma2[k] <- 1/rgamma(1, shape = a_new, rate = b_new)
		# Step 6. update tau2_k
			inv_g_mean <- sqrt(sigma2[k]) * lambda[k]/abs(beta[,k])
			inv_g_mean <- ifelse(inv_g_mean < 1e-12, 1e-12, inv_g_mean)
			tau2[,k] <- 1/rinvgauss(n = p, mean = inv_g_mean, shape = lambda[k]^2)
		# Step 7. update lambda_k
			lambda[k] <- sqrt(rgamma(1, shape = p + 0.5, rate = 0.5*( sum(tau2[,k]) + delta[k] )))	
		# Step 8. update delta_k	
			delta[k] <- rgamma(1, shape = 1, rate = 0.5*(lambda[k]^2 + 1))
		# Step 9. update mu_k
			mu_mean <- (t(x_data) %*% Z[,k] + m0)/(n_k[k] + 1)
			mu_var <- Sigma[,,k]/(n_k[k] + 1)
			if(all(is.finite(mu_var))){
				mu[,k] <- rmvnorm(1, mean = mu_mean, sigma = mu_var)
			}

		}
		if(K_plus < K){
		for(k in (K_plus+1):K){
			alpha_mean <- 0
			alpha_var <- sigma2_alpha
			alpha[k] <- rnorm(1, alpha_mean, sqrt(alpha_var))
			sigma2[k] <- 1/rgamma(1, shape = a, rate = b)
			if(!is.finite(sigma2[k])){sigma2[k] = 1/rgamma(1, shape = 1, rate = 0.1)}			
			if(sigma2[k] > 100){sigma2[k] = 1/rgamma(1, shape = 1, rate = 1)}										
			lambda[k] <- rt(1, df = 1)
			while(lambda[k] < 0){
				lambda[k] <- rt(1, df = 1)
			}
#			cat(paste0('k = ', k),'\n')
			tau2[,k] <- rexp(p, rate = lambda[k]^2)
			TAU_k <- diag(tau2[,k])
			zar <- sigma2[k] * TAU_k
			zar[!is.finite(zar)] <- 1/rgamma(sum(!is.finite(zar)*1), shape = 1, rate = 1)	
			beta[,k] <- rmvnorm(1, mean = rep(0,p), sigma = zar)
			delta[k] <- rgamma(1, shape = 0.5, rate = 0.5)	
			Sigma[,,k] <- 1*runif(1)*dvX	
			mu[,k] <- rmvnorm(1, mean = m0, sigma = Sigma[,,k])
		}
		}


		# Steps 10, 11, 12, 13. update Omega via the `blockGLasso_mix` function



		block_lasso_step <- blockGLasso_mix(X = x_data, z = z, mu = as.matrix(mu[,1:K]), mu0 = m0, K = K, 
				lambdaPriora = r, 
		    		lambdaPriorb = s,
		    		Sigma = array(Sigma[,,1:K],dim=c(p,p,K)),
		    		verbose = verbose
		    		)
		Sigma[,,1:K] <-  block_lasso_step$Sigma   		
		Omega[,,1:K] <-  block_lasso_step$Omega
		Psi <- block_lasso_step$Psi  
		# Step 14: update mixing proportions
		pr[1:K] <- myDirichlet(alpha_dirichlet/K + n_k)
#		print(length(pr))
		# Step 1. update latent allocations via the `update_z` function
		up_z <- update_z(pr = pr[1:K], mu = as.matrix(mu[,1:K]), alpha = alpha[1:K], 
			beta = as.matrix(beta[,1:K]), sigma2 = sigma2[1:K], Omega = array(Omega[,,1:K],dim=c(p,p,K)), 
			x_data = x_data, y_data = y_data)
		z <- up_z$z
#		log_post <- log_posterior(alpha_dirichlet = alpha_dirichlet, pr = pr[1:K], 
#			mu = mu[,1:K], alpha = alpha[1:K], beta = beta[,1:K], sigma2 = sigma2[1:K], 
#			Omega = Omega[,,1:K], Sigma = Sigma[,,1:K], lambda = lambda[1:K], logL = up_z$logL, 
#			z = z, K = K, Psi = Psi) 
		# store output every thin iterations
		if(iter %% thin == 0){
			pr_sim[iter/thin, ] <- pr
			mu_sim[iter/thin,,] <- mu
			alpha_sim[iter/thin,] <- alpha
			beta_sim[iter/thin, , ] <- beta
			sigma2_sim[iter/thin, ] <- sigma2
			Sigma_sim[iter/thin,,,] <- Sigma
			delta_sim[iter/thin,] <- delta
			tau2_sim[iter/thin,,] <- tau2
			lambda_sim[iter/thin,] <- lambda
			z_sim[iter/thin,] <- z
			K_sim[iter/thin] <- K
			K_plus_sim[iter/thin] <- length(table(z))
			logL[iter/thin] <- up_z$logL
			logP[iter/thin] <- log_post
			alpha_dirichlet_sim[iter/thin] <- alpha_dirichlet
		}	
		if(verbose){
		if(iter %% 100 == 0){
#			if(iter == 100){
#				matplot(mu, type = 'l', lwd = 2)
#			}else{
#				matplot(mu, type = 'l', lwd = 2, add = TRUE)
#			}
#			pairs(cbind(y_data, x_data), col = z)\
			mymain = paste0('clusters at iteration: ', iter)
			par(mfrow = c(2,2))
#			myYlim <- c(quantile(logP[1:iter], 0.05), max(logP[1:iter]))
			plot(logP[1:(iter/thin)], type = 'l', xlab = 'iteration (thinned)', ylab ='log-posterior' )
			plot(alpha_dirichlet_sim[1:(iter/thin)], type = 'l', xlab = 'iteration (thinned)', ylab ='alpha-dirichlet' )
			matplot(cbind(K_sim[1:(iter/thin)],K_plus_sim[1:(iter/thin)]), type = 'l', col = 1:2, lty = 1)
			matplot(t(cbind(y_data,x_data)),col = z, type = 'l', xaxt = 'n', main = mymain)
			axis(side = 1, at = 1:(p+1), labels = c('y',paste0('x',1:p)))
			legend('bottomleft', 
				paste0('component: ', names(table(z)), ', n', names(table(z)), ' = ', table(z)), 
				col = as.numeric(names(table(z))), lty = 1)
			cat(paste('iteration: ', iter), '\n')
			print(table(z))
			cat(paste0('MH rate: ', round(100*mh_rate/iter, 2), '%.'), '\n')
		}
		}
	}

	mh_rate <- 100*mh_rate/iter



	results <- vector('list', length = 16)
	results[[1]] <- pr_sim
	results[[2]] <- mu_sim
	results[[3]] <- alpha_sim
	results[[4]] <- beta_sim
	results[[5]] <- sigma2_sim
	results[[6]] <- Sigma_sim
	results[[7]] <- delta_sim
	results[[8]] <- tau2_sim
	results[[9]] <- lambda_sim		
	results[[10]] <- z_sim			
	results[[11]] <- logL		
	results[[12]] <- K_sim
	results[[13]] <- K_plus_sim				
	results[[14]] <- alpha_dirichlet_sim
	results[[15]] <- mh_rate
	results[[16]] <- logP
	names(results)[1:16] <- c('pr', 'mu', 'alpha', 'beta', 'sigma2', 'SIGMA', 
	'delta', 'tau2', 'lambda', 'z', 'logL', 'nComponents', 'nClusters', 'alpha_dirichlet', 'MH_acceptance_rate', 'logP')
	return(results)		
	}


# fit: object returned by telescoping sampler
# newdata: matrix-like object with same columns as X 
predict.bgcwm <- function(fit, newdata, alpha = 0.05, burn = 0.1){

	nIter <- length(fit$nClusters)
	burn <- floor(length(fit$nClusters)/burn^{-1})
	iter <- burn + 1
	res = c()
	for(j in 1:dim(newdata)[1]){
		y_new <- numeric(nIter - burn)
		z_new <- numeric(nIter - burn)
		for(iter in (burn+1):nIter){
			K <- fit$nClusters[iter]
			w <- fit$pr[iter,1:K]
			mu <- fit$mu[iter, ,1:K]
			SIGMA <- fit$SIGMA[iter,,,1:K]
	#		sample z = k|x,...
			probs <- numeric(K)
			for(k in 1:K){
				probs[k] = w[k] * mvtnorm::dmvnorm(x = newdata[j,], mean = mu[,k], sigma = SIGMA[,,k])
			}
			k <- sample(1:K, 1, prob = probs)
			z_new[iter - burn] <- k
	#		sample y|x, z = k,...		
			s2 <- fit$sigma[iter, k]
			b <- fit$beta[iter, , k]
			a <- fit$alpha[iter, k]
			cond.mean <- a + t(b) %*% newdata[j,]
			y_new[iter - burn] <- rnorm(1, mean = cond.mean, sd = sqrt(s2))		
		}
		
		y.pred <- mean(y_new)	
		y.low <- quantile(y_new, probs = alpha/2)
		y.up <- quantile(y_new, probs = 1 - alpha/2)
		
		d <- density(y_new)
		modes_idx <- which(diff(sign(diff(d$y))) == -2) + 1
		modes <- d$x[modes_idx]
		main_mode <- modes[which.max(d$y[modes_idx])]
#		map1 <- as.numeric(which.max(table(z_new)))
#		map1 <- modes[map1]
#		k <- sample(K, 1, prob =as.numeric(table(z_new)))
#		est1 <- modes[k]
		res <- rbind(res, c(y.pred, y.low, y.up, main_mode))
		cat(paste0("Calculating model predictions for new observation: ", j,'.'),'                 \r')
	}
	cat('\n')
	colnames(res)[c(1,4)] <- c('posterior_mean', 'map')
	return(res)
}



