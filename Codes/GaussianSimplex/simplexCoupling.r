n = 5
reps = 1e1

gibbs_simplex_distributions = function(n, reps) {
	CR = list()
	starting_point = (numeric(n) + 1)/n
	CR[[1]] = starting_point
	lambda = runif(reps)
	for(i in 2:reps){
		current_x = CR[[i-1]]
		proposed_x = current_x
		indices = sample(1:n, 2)
		proposed_x[min(indices)] = lambda[i]*(current_x[indices[1]] + current_x[indices[2]])
		proposed_x[max(indices)] = (1-lambda[i])*(current_x[indices[1]] + current_x[indices[2]])
		current_x = proposed_x
		CR[[i]] = current_x
		print(CR[[i]])
	}
	return(CR)
}

gibbs_simplex_coupling = function(n, reps) {
	X = list()
	Y = list()
	starting_x = c(0.1092629, 0.0457090, 0.1062776, 0.5549756, 0.1837749)
	starting_y = c(0.35648087, 0.08650361, 0.02071886, 0.05238302, 0.48391364)
	X[[1]] = starting_x
	Y[[1]] = starting_y
	lambda = runif(reps)
	for(i in 2:reps){
		current_x = X[[i-1]]
		current_y = Y[[i-1]]

		proposed_x = current_x
		proposed_y = current_y

		indices = sample(1:n, 2)
		proposed_x[min(indices)] = lambda[i]*(current_x[indices[1]] + current_x[indices[2]])
		proposed_x[max(indices)] = (1-lambda[i])*(current_x[indices[1]] + current_x[indices[2]])
		current_x = proposed_x

		proposed_y[min(indices)] = lambda[i]*(current_y[indices[1]] + current_y[indices[2]])
		proposed_y[max(indices)] = (1-lambda[i])*(current_y[indices[1]] + current_y[indices[2]])
		current_y = proposed_y

		X[[i]] = current_x
		Y[[i]] = current_y
	}
	X_Y = list("X" = X, "Y" = Y)
	return(X_Y)
}

CR = list()
CR = gibbs_simplex_distributions(n, reps)
CR

X_Y = list()
X_Y = gibbs_simplex_coupling(n, reps)
X = list()
X = X_Y$X
Y = list()
Y = X_Y$Y
X[[reps]]-Y[[reps]]
# 0.1092629 0.0457090 0.1062776 0.5549756 0.1837749
# 0.35648087 0.08650361 0.02071886 0.05238302 0.48391364