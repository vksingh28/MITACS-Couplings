




n = 4
reps = 1e2
gibbs_simplex_distributions = function(n) {
	CR = list()
	starting_point = (numeric(n) + 1)/n
	CR[[1]] = starting_point
	for(i in 2:reps){
		current_x = CR[[i-1]]
		proposed_x = current_x
		indices = sample(1:n, 2)
		lambda = runif(1)
		proposed_x[min(indices)] = lambda*(current_x[indices[1]] + current_x[indices[2]])
		proposed_x[max(indices)] = (1-lambda)*(current_x[indices[1]] + current_x[indices[2]])

		current_x = proposed_x
		CR[[i]] = current_x
		print(i)
		print(indices)
		print(lambda)
		print(CR[[i]])
	}
	return(CR)
}

CR = list()
CR = gibbs_simplex_distributions(4)
CR
