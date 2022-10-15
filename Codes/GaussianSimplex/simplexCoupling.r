n = 5
reps = 1e4

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
		CR[[i]] = current_x
		print(CR[[i]])
	}
	return(CR)
}

CR = list()
CR = gibbs_simplex_distributions(n, reps)
CR
