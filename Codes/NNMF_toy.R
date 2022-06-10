set.seed(1)
v1 = c(1,2,3)
v2 = c(5,8,10)
v3 = c(13,15,26)
M = matrix(c(v1,v1,v1,v2,v2,v2,v3,v3,v3), nrow = 3, ncol = 9)
M

# loss function / penalty function
penalty = function(C, R, M) {
	# trying to minimize l-0 norm and penalizing negative terms
	loss = 1000*sum(abs(M-C%*%R)^0.001) + 10*sum((abs(M-C%*%R))^0.5)^2 + 200*(sum(C<0) + sum(R<0))
	# loss = 1000*sum((M-C%*%R)!=0) + sum((abs(M-C%*%R))^0.5)^2 + 200*(sum(C<0) + sum(R<0)) + sum(abs(M-C%*%R))
	print(paste("loss:", loss))
	return (loss)
}

# probability distribution (beta is a hyperparameter) to draw C and R
prob_dist = function(CR, M) {
	beta = .01
	C = matrix(CR[1:6], nrow = 3, ncol = 2)
	R = matrix(CR[7:24], nrow = 2, ncol = 9)
	# discussed probability distribution using the penalty function
	value = exp(-beta*penalty(C, R, M))
	print(value)
	return (value)
}

# using metropolis with single particle moves
reps = 1e4
accept = 0
gibbs_nnmf = function(M){
	CR = list()
	# CR[[1]] = c(v2,v1,0,0,.9,.9,1,1,.9,.9,0,0,2,2)
	v = runif(3, min=1, max=2)*(sqrt(sum(v3^2)))
	v
	CR[[1]] = c(v3,v1,0,1,0,1,0,1,0,0,0,0,0,0,1,0,1,0,1,0)
	for(i in 2:reps){
		print(i)
		current_x = CR[[i-1]]
		for(j in c(13,14,15,16,17,18)){
			proposed_x = current_x
			proposed_x[j] = current_x[j] + rnorm(1, 0, .05)
			mh_ratio = prob_dist(proposed_x, M)/prob_dist(current_x, M)
			print(paste("mh: ", mh_ratio))
			if(runif(1) < min(1, mh_ratio)){
				current_x = proposed_x
				accept = accept + 1
				print(accept)
			}
		}
		CR[[i]] = current_x
	}
	return (CR)
}

CR = list()
CR = gibbs_nnmf(M)
C = matrix(CR[[reps]][1:6], nrow = 3, ncol = 2)
R = matrix(CR[[reps]][7:24], nrow = 2, ncol = 9)
C
R
M-C%*%R
C%*%R
M
