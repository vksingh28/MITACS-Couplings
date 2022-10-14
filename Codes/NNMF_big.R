library("msm")
m = 10
n = 4*3
k = 2
v = list()
for(i in 1:m){
	v[[i]] = numeric(m)
	v[[i]][i] = 1
}
M = matrix(c(v[[1]],v[[1]],v[[1]],v[[6]],v[[6]],v[[6]],v[[5]],v[[5]],v[[5]],v[[8]], v[[8]], v[[8]]), nrow = m, ncol = n)
M

# loss function / penalty function
penalty = function(C, R, M) {
	# trying to minimize l-0 norm
	# 0<p<1
	p = 0.5
	q = 1/p
	loss = sum((abs(M-C%*%R))^p)^(q) + 0.01*sum((abs(M-C%*%R))^2)^0.5 + 0.01*sum((R > 1) && (R < 0.001))
	loss
	return (loss)
}

# probability distribution (beta is a hyperparameter) to draw C and R
prob_dist = function(CR, M) {
	beta = 5
	x = m*k
	y = m*k+1
	z = (m+n)*k
	C = matrix(CR[1:x], nrow = m, ncol = k)
	R = matrix(CR[y:z], nrow = k, ncol = n)
	value = -beta*penalty(C,R,M)
	return (value)
}

# using metropolis with single particle moves
reps = 44*1e4
accept = 0
e = list()
for(i in 1:k){
	e[[i]] = numeric(k)
	e[[i]][i] = 1
}

initial_node = c(v[[1]], v[[5]], e[[1]], e[[1]], e[[1]], numeric(k), numeric(k), numeric(k), e[[2]],e[[2]], e[[2]], numeric(k), numeric(k), numeric(k))
# initial_node = rtnorm(44, mean=0, sd = 10, lower = 0)
gibbs_nnmf = function(M){
	CR = list()
	acceptance_prob = c()
	prob_density = c()
	CR[[1]] = initial_node + c(rtnorm(20, mean=0, sd=5, lower=0), rtnorm(24, mean=0, sd=1, lower=0))
	# indices = c(1:15,16,19,22,25,28,31,34,37,40,43)
	indices = c(1:44)
	for(i in 2:reps){
		current_x = CR[[i-1]]
		j = sample(indices, 1)
		proposed_x = current_x
		proposed_x[j] = current_x[j] + rnorm(1, 0, .1)
		prob_density = c(prob_density, prob_dist(proposed_x, M))
		mh_ratio = (prob_dist(proposed_x, M)-prob_dist(current_x, M))
		if(proposed_x[j] >= 0){
			acceptance_prob = c(acceptance_prob, min(1,exp(mh_ratio)))
		}
		print(paste("mh: ", mh_ratio))
		if((log(runif(1)) < mh_ratio) && (proposed_x[j] >= 0)){
			current_x = proposed_x
			accept = accept + 1
		}
		CR[[i]] = current_x
	}
	print(accept)
	print(accept/(reps))
	# hist(prob_density)
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
C1 = matrix(initial_node[1:x], nrow = m, ncol = k)
R1 = matrix(initial_node[y:z], nrow = k, ncol = n)
C1
C
penalty(C, R, M)
penalty(C1, R1, M)
R1
R

M-C%*%R
x = matrix(M-C%*%R, nrow=1)
hist(x)
M-C1%*%R1
# C%*%R
# M
CR[[reps]][16:18]
# v3,v2,v5
# v3,v1,v5

# for(k in c(16,19,22,25,28,31,34,37,40,43)){
# 	# r = k+ceiling(runif(2, 0, 3))-1
# 	for(rep in 1:3){
# 		r = k+c((rep-1)%%3, (rep)%%3)
# 		print(r)
# 	}
# }
