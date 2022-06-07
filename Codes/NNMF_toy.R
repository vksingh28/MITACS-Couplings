set.seed(1)
v1 = c(1,2,3)
v2 = c(5,8,10)
M = matrix(c(v1,v1,v1,v1,v2,v2,v2,v2), nrow = 3, ncol = 8)
M

# loss function / penalty function
penalty = function(C, R, M) {
	# trying to minimize l-0 norm and penalizing negative terms
	loss = 10*sum((M-C%*%R) != 0) + 100*(sum(C<0) + sum(R<0)) + sum(abs(M-C%*%R))
	print(paste("loss:", loss))
	return (loss)
}

# probability distribution (beta is a hyperparameter) to draw C and R
prob_dist = function(CR, M) {
	beta = 0.1
	C = matrix(CR[1:3], nrow = 3, ncol = 1)
	R = matrix(CR[4:11], nrow = 1, ncol = 8)
	# discussed probability distribution using the penalty function
	value = exp(-beta*penalty(C, R, M))
	print(value)
	return (value)
}

# using gibbs sampling with normal proposal
reps = 1e4
mcmc_nnmf = function(M) {
	CR = list()
	CR[[1]] = c(v1, 1,1,1,1,3,3,3,3)
	for(i in 2:reps){
		current_x = CR[[i-1]]
		sd = 1
		CR[[i]] = numeric(11)
		proposed_x = current_x + rnorm(11, 0, sd)
		mh_ratio = (prob_dist(proposed_x, M)/prob_dist(current_x, M))
		if(runif(1) < min(1, mh_ratio)){
			CR[[i]] = proposed_x
		} else {
			CR[[i]] = current_x
		}
	}
	return (CR)
}

CR = list()
CR = mcmc_nnmf(M)
C = matrix(CR[[reps]][1:3], nrow = 3, ncol = 1)
R = matrix(CR[[reps]][4:11], nrow = 1, ncol = 8)
C
R
M-C%*%R
