library('msm')
# N = #number of games
# n = #number of teams
# Choose constants

N = 20 # Number of teams
n = 1000 # Number of games
# a = (1:N)/(N+3) # True skill
a = sort(rtnorm(N, 0, sd=100, lower=0, upper=10000))

###########################
# Define helper functions #
###########################


# This makes you win if your skill is at all better.

WinThresh = function(a1,a2) {
	if(a1 > a2){
		return(1)
	}
	return(-1)
}

# This makes you win with some probability proportional to your skill

WinProp = function(a1,a2) {
	if(runif(1,0,1) < a1/(a1+a2)){
		return(1)
	}
	return(-1)
}

#########################
# Main Simulation Stuff #
#########################

# Simulates data. The WinFct is the function used by the BT model (that is: P[team 1 wins] = WinFct(team1, team2)). By default the better team wins with probability 1, which is stupid but easy.

BTSim = function(N,n,a, WinFct=WinThresh) {
	res = matrix(0,nrow=n,ncol=3)
	for(i in 1:n) {
		teams = sample(N,2)
		res[i,1] = teams[1]
		res[i,2] = teams[2]
		res[i,3] = WinFct(a[teams[1]], a[teams[2]])
	}
	return(res)
}
outcomes = BTSim(N,n,a)
outcomes[1, 1]


penalty = function(n, order) {
	loss = 0
	for(i in 1:n){
		home = outcomes[i, 1]
		away = outcomes[i, 2]
		loss = loss + (order[home] > order[away])*(outcomes[i, 3] == 1) + (order[home] < order[away])*(outcomes[i, 3] == -1)
	}
	return (loss)
}

prob_dist = function(n, order) {
	beta = 1
	value = -beta*penalty(n, order)
	# print(prob_dist)
	return (value)
}

reps = 1e3
gibbs_rankings = function(n) {
	CR = list()
	acceptance_prob = c()
	initial_order = c(1:N)
	prob_density = c(prob_dist(n, initial_order))
	CR[[1]] = initial_order
	accept = 0
	for(i in 2:reps){
		current_x = CR[[i-1]]
		proposed_x = current_x
		indices = sample(initial_order, 2)
		temp = proposed_x[indices[1]]
		proposed_x[indices[1]] = proposed_x[indices[2]]
		proposed_x[indices[2]] = temp
		mh_ratio = prob_dist(n, proposed_x) - prob_dist(n, current_x)
		prob_density = c(prob_density, prob_dist(n, proposed_x))
		print(prob_dist(n, proposed_x))
		acceptance_prob = c(acceptance_prob, min(1, exp(mh_ratio)))
		if((log(runif(1)) < mh_ratio)){
			current_x = proposed_x
			accept = accept+1
		}
		CR[[i]] = current_x
	}
	print("acceptance")
	print(accept/(reps-1))

	return (CR)
}

CR = list()
CR = gibbs_rankings(n)
print(CR[[reps]])