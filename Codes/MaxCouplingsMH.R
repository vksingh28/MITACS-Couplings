# Function to draw from \Pi_{X, Y}(.)
Pi_X.Y = function(x, y) {

}

# Function to draw from \r_X(.)
Remainder_X = function(x, y){

}

# Function to draw from \r_Y(.)
Remainder_Y = function(x, y){

}

# Function that returns normalizing constant of \Pi_{X, Y}(.)
Normalize_Pi = function(x, y){

}

# Number of MC draws
N = 1e4

# Parallel MCs
X = numeric(N)
Y = numeric(N)

# Initial draws for both the chains
X[1] = rnorm(1, 0, 1)
Y[1] = rnorm(1, 0, 1)

# Unif r.v.s
U = runif(N, 0, 1)

for(i in 2:N){
	if(U[i] < Normalize_Pi(X[i-1], Y[i-1])){
		X[i] = Y[i] = Pi_X.Y(X[i-1], Y[i-1])
	} else {
		X[i] = Remainder_X(X[i-1], Y[i-1])
		Y[i] = Remainder_Y(X[i-1], Y[i-1])
	}
}
