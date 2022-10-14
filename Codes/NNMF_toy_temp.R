library("msm")
library("bayesplot")
m = 5
n = 5*2
k = 3
x = m*k
y = m*k+1
z = (m+n)*k
v = list()
for(i in 1:m){
	v[[i]] = numeric(m)
	v[[i]][i] = 1
}
M = matrix(c(v[[1]],v[[1]],v[[2]],v[[2]],v[[3]],v[[3]],v[[4]],v[[4]],v[[5]],v[[5]]), nrow = m, ncol = n)
M

# loss function / penalty function
penalty = function(C, R, M) {
	# trying to minimize l-0 norm
	# 0<p<1
	p = 0.5
	q = 1/p
	loss = sum((abs(M-C%*%R))^p)^(q) + 0.01*sum((abs(M-C%*%R))^2)^0.5
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
reps = 1e3*45
accept = 0
e = list()
for(i in 1:k){
	e[[i]] = numeric(k)
	e[[i]][i] = 1
}

initial_node = c(v[[3]], v[[2]], v[[5]], numeric(3), numeric(3), e[[2]],e[[2]], e[[1]],e[[1]], numeric(3), numeric(3), e[[3]],e[[3]])
# initial_node = rtnorm(45, mean=0, sd = 10, lower = 0)
gibbs_nnmf = function(M){
	CR = list()
	acceptance_prob = c()
	prob_density = c()
	penalty_vector = c()
	CR[[1]] = initial_node + c(rtnorm(15, mean=0, sd=2, lower=0), rtnorm(30, mean=0, sd=1, lower=0))
	C_temp = matrix(CR[[1]][1:x], nrow = m, ncol = k)
	R_temp = matrix(CR[[1]][y:z], nrow = k, ncol = n)
	penalty_vector = c(penalty_vector, penalty(C_temp, R_temp, M))
	for(i in c(16,19,22,25,28,31,34,37,40,43)){
		j = i+2
		CR[[1]][i:j] = CR[[1]][i:j]/(sum(CR[[1]][i:j]))
	}
	indices = c(1:15,16,19,22,25,28,31,34,37,40,43)
	for(i in 2:reps){
		current_x = CR[[i-1]]
			j = sample(indices, 1)
			if(j <= 15){
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
			} else {
					proposed_x = current_x
					r = sample(c(j, j+1, j+2), 2, replace=FALSE)
					rw = rnorm(1, 0, .1)
					proposed_x[r[1]] = current_x[r[1]] + rw
					proposed_x[r[2]] = current_x[r[2]] - rw
					prob_density = c(prob_density, prob_dist(proposed_x, M))
					mh_ratio = (prob_dist(proposed_x, M)-prob_dist(current_x, M))
					if((proposed_x[r[1]] >= 0) && (proposed_x[r[2]] >= 0)){
						acceptance_prob = c(acceptance_prob, min(1,exp(mh_ratio)))
					}
					print(paste("mh: ", mh_ratio))
					if((log(runif(1)) < mh_ratio) && (proposed_x[r[1]] >= 0) && (proposed_x[r[2]] >= 0)){
						current_x = proposed_x
						accept = accept + 1
					}
				# }
			}
		# }
		CR[[i]] = current_x
		C_temp = matrix(CR[[i]][1:x], nrow = m, ncol = k)
		R_temp = matrix(CR[[i]][y:z], nrow = k, ncol = n)
		penalty_vector = c(penalty_vector, penalty(C_temp, R_temp, M))
	}
	print(accept)
	print(accept/(reps))
	CR_data = data.frame(t(sapply(CR, c)))
	# plot(CR_data, penalty_vector)
	# hist(prob_density)
	# hist(acceptance_prob)
	# plot(penalty_vector[1e4:reps], prob_density[1e4:reps], type="p", lwd=0.05)
	# mcmc_acf(CR)
	
	# acf(CR)
	return (CR)
}
CR = list()
CR = gibbs_nnmf(M)
CR_data = data.frame(t(sapply(CR, c)))
plot(CR_data)
acf(CR_data)
 # x = m*k
# y = m*k+1
# z = (m+n)*k
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