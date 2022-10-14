set.seed(1)
# v1
# v2
# v3
# M - mxn
# m,n,k

energy = function(C, R, M){
	loss = 1000*sum(abs(M-C%*%R)^0.001) + sum((abs(M-C%*%R))^0.5)^2 + 200*(sum(C<0) + sum(R<0))
	print(paste("loss: ", loss))
	return (loss)
}

target_density = function(CR, M){
	beta = 0.01
	# C
	# R
	prob = exp(-beta*energy(C, R, M))
	print("density: ", prob)
	return(prob)
}

reps = 1e3
# todo; implement randomization in selecting the index later, learn how to return accept
mcmc = function(M){
	CR = list()
	accept = 0
	# CR[[1]] =
	for(i in 2:reps) {
		print(paste("iteration: ", i))
		current = CR[[i-1]]
		for(j in c()){
			proposed_x = current_x
			proposed_x[j] = current_x[j] + rnorm(1,0,.05)
			mh_ratio = target_density(proposed_x, M)/target_density(current_x, M)
			print(paste("mh ratio: ", mh_ratio)
			if(runif(1) < min(1, mh_ratio)){
				current_x = proposedd_x
				accept = accept+1
			}
		}
		CR[[i]] = current_x
	}
	CR[[reps+1]] = accept
	return (CR)
}

CR = list()
CR = mcmc(M, k)
accept = CR[[reps+1]]
# C
# R
# M-C%*%R