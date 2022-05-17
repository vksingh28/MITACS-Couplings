set.seed(1)

# Function to draw from \Pi_{X, Y}(.)
# Rejection sampling using standard cauchy distribution
Pi_X.Y = function(x, y) {
	# unnormalised pdf of Pi_{X, Y}
	Pi = function(t){
		min(dnorm(t, mean = x, sd = 1), dnorm(t, mean = y, sd = 1))
	}

	# interval to calculate max of f/g
	interval = c(min(x, y)-10, max(x, y)+10)
	M = max(Pi(interval)/dcauchy(interval))
	accept = 0
	while(accept == 0){
		U = runif(1)
		prop = rcauchy(1)
		ratio = Pi(prop)/(M*dcauchy(prop))
		if(U <= ratio){
			accept = 1
			return(prop)
		}
	}
}

# Function to draw from r_X(.)
# Rejection sampling using standard cauchy distribution
Remainder_X = function(x, y){
	# unnormalised pdf of r_X(.)
	r_X = function(t){
		dnorm(t, mean = x, sd = 1) - min(dnorm(t, mean = x, sd = 1), dnorm(t, mean = y, sd = 1))
	}
	# interval used to calculate max of f/g
	interval = c(min(x, y)-10, max(x, y)+10)
	M = max(r_X(interval)/dcauchy(interval))
	accept = 0
	while(accept == 0){
		U = runif(1)
		prop = rcauchy(1)
		ratio = r_X(prop)/(M*dcauchy(prop))
		if(U <= ratio){
			accept = 1
			return(prop)
		}
	}
}

# Function to draw from r_Y(.)
Remainder_Y = function(x, y){
	# unnormalised pdf of r_X(.)
	r_Y = function(t){
		dnorm(t, mean = y, sd = 1) - min(dnorm(t, mean = x, sd = 1), dnorm(t, mean = y, sd = 1))
	}
	# interval used to calculate max of f/g
	interval = c(min(x, y)-10, max(x, y)+10)
	M = max(r_Y(interval)/dcauchy(interval))
	accept = 0
	while(accept == 0){
		U = runif(1)
		prop = rcauchy(1)
		ratio = r_Y(prop)/(M*dcauchy(prop))
		if(U <= ratio){
			accept = 1
			return(prop)
		}
	}
}

# Function that returns normalizing constant of \Pi_{X, Y}(.)
Normalize_Pi = function(x, y){
	integrand = function(t){
		min(dnorm(t, mean = x, sd = 1), dnorm(t, mean = y, sd = 1))
	}
	return(integrate(Vectorize(integrand), lower = -Inf, upper = Inf))
}

# Number of MC draws
N = 1e2

# Parallel MCs
X = numeric(N)
Y = numeric(N)

# Initial draws for both the chains
X[1] = rnorm(1)
Y[1] = rnorm(1)

# Unif r.v.s
U = runif(N, 0, 1)

# Final algorithm
for(i in 2:N){
	if(U[i] < Normalize_Pi(X[i-1], Y[i-1])$value){
		X[i] = Y[i] = Pi_X.Y(X[i-1], Y[i-1])
	} else {
		X[i] = Remainder_X(X[i-1], Y[i-1])
		Y[i] = Remainder_Y(X[i-1], Y[i-1])
	}
}
print(head(X-Y, n = 20))
