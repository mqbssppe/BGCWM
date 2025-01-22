source('joint_mix.R')
library('label.switching')
library('flexmix')
library('mclust')

mcmc_iter <- 10000
df <- read.table(file = "real_dataset.txt")
# response variable
y <- as.numeric(df$y)
# ground-truth classification
z <- df$z
# covariates
X <- as.matrix(df[,2:16])

p <- dim(X)[2]

#	INITIALIZATION using flexmix
set.seed(9)
       init <- initialize_sampler(x_data = X, y_data = y, mcmc_iter = 2, thin = 1, K_max = 30,
                        r = 1,                  
                        s = 10^{-2},
                        sigma2_alpha = 10^3,
                        a = 0.01,
                        b = 0.01,
                        verbose = FALSE,
                        initial_values = NULL,
                        warm_up = 2,
                        initial_warm_up = 1,
                        nu_l = 6, nu_r = 3,     #dirichlet alpha prior parameters
                        a_lambda = 1, a_pi = 4, b_pi = 3, #bnb prior parameters
                        s2_a = 5, nRuns = 10, burn = 20, flexmix_runs = 30)
	initial_values <- vector('list', length = 1)
	initial_values$z <- init$z
#	The following is the main function for the telescoping sampler of the paper
#		here we will only run a single chain for 5000 MCMC iterations (for illustration purposes)	
        fit <- telescoping_sampler(x_data = X, y_data = y, mcmc_iter = 5000, thin = 1, K_max = 30,
                        r = 1,                  
                        s = 10^{-2},
                        sigma2_alpha = 10^3,
                        a = 0.01,
                        b = 0.01,
                        verbose = TRUE,
                        initial_values = initial_values,
                        nu_l = 6, nu_r = 3,     #dirichlet alpha prior parameters
                        a_lambda = 1, a_pi = 4, b_pi = 3, #bnb prior parameters
                        s2_a = 5, # scale for mh proposal
			warm_up = 50,
			initial_warm_up = 1  
                )
#	plot of the sampled values of the number of clusters
	plot(fit$nClusters)
#	burn the first 1000 iterations
	burn <- 1000
#	retained MCMC draws of the number of clusters
	nClusters <- fit$nClusters[-(1:burn)]
	table(nClusters)
#	retained MCMC draws of allocation variables	
	z_sim <- fit$z[-(1:burn),]	
	K_hat =  as.numeric(names(which.max(table(nClusters))))
#	log-posterior density (up to normalizing const)	
	logP <- fit$logP[-(1:burn)]
#	undo label switching conditional on K_hat = 4
	z_values <- z_sim[which(nClusters == K_hat),]
	l_values <- logP[which(nClusters == K_hat)]       
	map_index <- which.max(l_values)
	zpivot = z_values[map_index,]
	myZ <- t(apply(z_values, 1, function(x)as.numeric(as.character(factor(x, labels = 1:K_hat)))  ))
	ls <- label.switching(method = 'ECR', zpivot = zpivot, z = myZ, groundTruth = z)
#	estimated clustering after undoing label-switching	
	z_hat <- ls$clusters[1,]
#	confusion matrix between estimated and true clusters	
	table(z, z_hat)	
#	adjusted Rand Index
	adjustedRandIndex(z, z_hat)	


#########################################################################################################################
#	NOTE: the results in the paper correspond to a much larger number of MCMC iterations and multiple MCMC chains!
#########################################################################################################################	
	
		




