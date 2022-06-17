library("msm")
set.seed(1)
# v1 = c(1,2,3)
# v2 = c(5,8,10)
# v3 = c(13,15,26)
m = 5
n = 5*2
k = 3
v = list()
for(i in 1:m){
	v[[i]] = numeric(m)
	v[[i]][i] = 1
}
# v[[1]] = rtnorm(m, mean = 0, sd = 10, lower=0)
# v[[2]] = rtnorm(m, mean = 5, sd = 2, lower=0)
# v[[3]] = rtnorm(m, mean = 3, sd = 6, lower=0)
# v[[4]] = rtnorm(m, mean = 20, sd = 1, lower=0)
# v[[5]] = rtnorm(m, mean = 15, sd = 0.1, lower=0)
# v[[6]] = rtnorm(m, mean = 1, sd = 5, lower=0)
M = matrix(c(v[[1]],v[[1]],v[[2]],v[[2]],v[[3]],v[[3]],v[[4]],v[[4]],v[[5]],v[[5]]), nrow = m, ncol = n)
M

# loss function / penalty function
penalty = function(C, R, M) {
	# trying to minimize l-0 norm
	loss = 0.01*sum(R*(R>1)) + 0.01*sum((abs(M-C%*%R))^2)^0.5 + sum((abs(M-C%*%R))^0.5)^2
	# loss = 1000*sum((M-C%*%R)!=0) + sum((abs(M-C%*%R))^0.5)^2 + 200*(sum(C<0) + sum(R<0)) + sum(abs(M-C%*%R))
	return (loss)
}

# probability distribution (beta is a hyperparameter) to draw C and R
prob_dist = function(CR, M) {
	beta = 10
	x = m*k
	y = m*k+1
	z = (m+n)*k
	C = matrix(CR[1:x], nrow = m, ncol = k)
	R = matrix(CR[y:z], nrow = k, ncol = n)
	# discussed probability distribution using the penalty function
	# value = exp(-beta*penalty(C, R, M))
	value = -beta*penalty(C,R,M)
	return (value)
}

# using metropolis with single particle moves
reps = 1e3
accept = 0
e = list()
for(i in 1:k){
	e[[i]] = numeric(k)
	e[[i]][i] = 1
}

gibbs_nnmf = function(M){
	CR = list()
	acceptance_prob = c()
	prob_density = c()
	CR[[1]] = c(v[[3]], v[[2]], v[[5]], e[[1]],e[[1]], e[[2]],e[[2]], numeric(3), numeric(3), numeric(3), numeric(3), e[[3]],e[[3]]) + c(rtnorm(15, mean=0, sd=1, lower=0), rtnorm(30, mean=0, sd=1, lower=0))
	for(i in 2:reps){
		current_x = CR[[i-1]]
		# possible_indices = (m*k+1):(m*k+k*n) - ()
		# possible_indices = c(16:18,19:21,34:36,37:39, 40:42, 43:45)
		for(j in 1:45){
			proposed_x = current_x
			proposed_x[j] = current_x[j] + rnorm(1, 0, .1)
			prob_density = c(prob_density, prob_dist(proposed_x, M))
			# mh_ratio = prob_dist(proposed_x, M)/prob_dist(current_x, M)
			mh_ratio = (prob_dist(proposed_x, M)-prob_dist(current_x, M))
			if(proposed_x[j] >= 0){
				acceptance_prob = c(acceptance_prob, min(1,exp(mh_ratio)))
			}
			print(paste("mh: ", exp(mh_ratio)))
			if((runif(1) < exp(mh_ratio)) && (proposed_x[j] >= 0)){
				current_x = proposed_x
				accept = accept + 1
				print(accept)
			}
		}
		CR[[i]] = current_x
	}

	hist(prob_density)
	hist(acceptance_prob)
	return (CR)
}

CR = list()
CR = gibbs_nnmf(M)
x = m*k
y = m*k+1
z = (m+n)*k
C = matrix(CR[[reps]][1:x], nrow = m, ncol = k)
R = matrix(CR[[reps]][y:z], nrow = k, ncol = n)
C1 = matrix(CR[[1]][1:x], nrow = m, ncol = k)
R1 = matrix(CR[[1]][y:z], nrow = k, ncol = n)
C1
C

R1
R

M-C%*%R
M-C1%*%R1
# C%*%R
# M

#           [,1]       [,3]      [,5]      [,7]
# [1,]  7.551818   6.476649  8.663017  20.59390
# [2,] 11.816428   6.151563  7.927327  20.91898
# [3,]  1.457067   4.389223  9.749586  20.78214
# [4,]  1.397953   8.023562  2.730398  20.07456
# [5,]  4.360686   5.779686  2.902858  18.01065
#          [,9]    [,10]      [,12]
# [1,]  15.06198  3.0897078
# [2,]  14.99439  7.7933978
# [3,]  14.98442  0.4860614
# [4,]  14.85292  2.9383581
# [5,]  14.95218  0.7309748

# v2,v3,v6
