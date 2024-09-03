
set.seed(123)

#####################################################################################
# Libraries
# install.packages("rBeta2009")
library(rBeta2009)
# install.packages("Rmpfr")
library(Rmpfr)
#####################################################################################
# GMM Sampler

# Define parameters for the Gaussian Mixture Model
d = 100	# Number of families (dimension of the simplex)
means = runif(d, -25, 25);
variances = runif(d, 5, 10);
proportions = as.vector(rdirichlet(1, rep(1, d)))  # Proportions for each component

# Number of samples to generate
num_samples = 10000

# Generate data from the Gaussian Mixture Model
generate_gmm_sample = function(n, means, variances, proportions) {
  component = sample.int(length(means), n, replace = TRUE, prob = proportions)
  samples = numeric(n)
  families = numeric(n)
  for (i in 1:n) {
    samples[i] = rnorm(1, mean = means[component[i]], sd = sqrt(variances[component[i]]))
    families[i] = component[i]
  }
  return(list(samples = samples, families = families))
}

generated_data = generate_gmm_sample(num_samples, means, variances, proportions)

# Access generated samples and their corresponding families
generated_samples = generated_data$samples
component_families = generated_data$families
total_families = as.vector(table(factor(component_families, levels = 1:d)))

# Plot the generated data
hist(generated_samples, breaks = 30, freq = FALSE, main = "Generated Data from GMM", xlab = "Value")
lines(density(generated_samples), col = "blue", lwd = 2)


#####################################################################################

# gibbs_MH_simplex = function(c, d, total_families) {
# 	c_new = c
# 	lambda = runif(1)
# 	indices = sample(1:d, 2, replace = FALSE)
# 	i = indices[1]; j = indices[2];
# 	c_new[i] = lambda * (c[i] + c[j])
# 	c_new[j] = (1-lambda) * (c[i] + c[j])
# 	log_rho = sum(total_families)
# 	# MH step
# 	log_rho = sum(total_families*log(c_new)) - sum(total_families*log(c))
# 	log_unif = log(runif(1))

# 	if(log_unif < min(0, log_rho)){
# 		return(c_new)
# 	}
# 	return(c)
# }
gibbs_MH_simplex = function(c, d, total_families, delta) {
	c_new = c
	lambda = runif(1,-1,1)
	indices = sample(1:d, 2, replace = FALSE)
	i = indices[1]; j = indices[2];
	c_new[i] = c[i] + lambda*delta
	c_new[j] =  c[j] - lambda*delta
	if ((min(c_new[i], c_new[j]) < 0) | (max(c_new[i], c_new[j]) > 1)){
		return(c)
	}
	log_rho = sum(total_families)
	# MH step
	log_rho = sum(total_families*log(c_new)) - sum(total_families*log(c))
	log_unif = log(runif(1))

	if(log_unif < min(0, log_rho)){
		return(c_new)
	}
	return(c)
}

# Gibbs sampler + MH
R = 1e5
burn_in = 1e4
C = matrix(NA, nrow = d, ncol = R)
c = as.vector(rdirichlet(1, rep(1, d)))
r = 1
delta = 0.01	# step size for proposal
for(r in 1:(R + burn_in)){
	if(r%%1000 == 0 & r <= burn_in){
      print(paste("Warm-up: iteration:", r))
    } else if (r%%1000 == 0 & r > burn_in){
      print(paste("Sampling: iteration:", r-burn_in))
    }

	# Proposal from uniform d-simplex using Gibbs (not anymore)
	# c = gibbs_MH_simplex(c, d, total_families)
	c = gibbs_MH_simplex(c, d, total_families, delta)

	# Store values after burnin done
	if(r > burn_in){
		C[, r - burn_in] = c
	}
}
# pretty stable estimates
cat("Row means =", rowMeans(C), "\nrow variances =", apply(C, 1, var), "\n")
# proportions
cat("actual props =", proportions)
cat("error =", mean((proportions-rowMeans(C))^2), "\n")


# TODO: Need to change the algorithm for proportional coupling everywhere
# TODO: Need to fix the bug as mentioned by aaron in 4.b
#####################################################################################
# MH-coupling
prop_coupling = function(total_families, x, y, d){
	# delta = exp(-iter/10)
	x_new = x; y_new = y;
	lambda = runif(1,0,1)
	indices = sample(1:d, 2, replace = FALSE)
	i = indices[1]; j = indices[2];
# old stuff
	x_new[i] = lambda * (x[i] + x[j])
	x_new[j] = (1-lambda) * (x[i] + x[j])
	y_new[i] = lambda * (y[i] + y[j])
	y_new[j] = (1-lambda) * (y[i] + y[j])
# new stuff
	# flag_x = 0
	# x_new[i] = x[i] + lambda*delta;
	# x_new[j] = x[j] - lambda*delta;
	# if ((min(x_new[i], x_new[j]) < 0) | (max(x_new[i], x_new[j]) > 1)){
	# 	flag_x = 1
	# }
	# flag_y = 0
	# y_new[i] = y[i] + lambda*delta;
	# y_new[j] = y[j] - lambda*delta
	# if ((min(y_new[i], y_new[j]) < 0) | (max(y_new[i], y_new[j]) > 1)){
	# 	flag_y = 1
	# }
	# MH Step
	log_rho_x = sum(total_families*log(x_new)) - sum(total_families*log(x))
	log_rho_y = sum(total_families*log(y_new)) - sum(total_families*log(y))

	# can use other couplings, rn using the same thing
	log_unif_x = log_unif_y = log(runif(1))

	# if(!flag_x & log_unif_x < min(0, log_rho_x)){
	# 	x = x_new
	# }
	# if(!flag_y & log_unif_y < min(0, log_rho_y)){
	# 	y = y_new
	# }

	if(log_unif_x < log_rho_x){
		# cat("x, iter ", r, "\n")
		x = x_new
	}
	if(log_unif_y < log_rho_y){
		# cat("y, iter ", r, "\n")
		y = y_new
	}
	# x = x_new; y = y_new;
	# cat(mean((x-y)^2), "\n")
	return(list("x" = x, "y" = y))
}
#####################################################################################
#####################################################################################
# First Coupling Stage
#####################################################################################

# Proportional Coupling
exp_bound = 5
R_prop = ceiling(3/2*exp_bound*d*log(d))*200 # Trying out the Theorem of Simplex Mixing Time
X = matrix(NA, nrow = d, ncol = R_prop)
# X = double(matrix(NA, nrow = d, ncol = R_prop))	# Initialize chain randomly
Y = matrix(NA, nrow = d, ncol = R_prop)	# Initialize chain from the actual distribution using the gibbs sample
# Y = double(matrix(NA, nrow = d, ncol = R_prop))
X[, 1] = as.vector(rdirichlet(1, rep(1, d)));
x = X[, 1]
# x = X[, 1]
Y[, 1] = C[, R];
y = Y[, 1]
# y = Y[, 1]
for(r in 2:R_prop){
	if(r %% 1000 == 0){
		print(paste("iteration: ", r))
	}
	coupled_val = prop_coupling(total_families, x, y, d)
	x = coupled_val$x; y = coupled_val$y
	X[, r] = x; Y[, r] = y
	if(r %% 1000 == 0){
		print(mean((x-y)^2))
	}
}
print(mean((X[, R_prop]-Y[, R_prop])^2))
print(2*d^(-exp_bound))

#####################################################################################
# Current goal right now and some thoughts
#####################################################################################
# we shrinked X and Y very close to each other
# now we need to create the sequence of graphs/partitions
# first we sample a shit ton of index pairs uniformly
# then we construct the partition pairs
# then we do the coupling steps which is if the time agrees, do subset, otherwise to proportional
# a little caveat here, when doing subset coupling, we'll also check for the MH step as we're not sampling
# from the uniform distribution, we have our own posterior. so our algorithm might take more time to couple,
# but probably only by a factor of something and won't change the complexity is what I feel!
#####################################################################################


# Second coupling stage
#####################################################################################
eps = 1;
T = ceiling((0.5 + eps)*d*log(d));

# matrix with a lot if index pairs
indices = matrix(NA, nrow = T+1, ncol = 2);
for(i in 2:(T+1)){
	indices[i, ] = sort(sample(1:d, 2, replace = FALSE));
}

# vector containing all the marked times
marked_times = c();
S.t.1 = list();
S.t.2 = list();
partitions = list();
partitions[[(T+1)]] = lapply(1:d, function(x) c(x));
for(t in (T+1):2){
	partition_curr = partitions[[t]];
	partition_next = list();
	i = indices[t-1, 1]; j = indices[t-1, 2];
	for(set in partition_curr){
		# S_i is the set that contains i, S_j is the set that contains j
		if(i %in% set){
			S_i = set;
		}
		if(j %in% set){
			S_j = set;
		}
	}
	if(!identical(S_i, S_j)){
		S.t.1[[t]] = sort(S_i);
		S.t.2[[t]] = sort(S_j);
		marked_times = c(marked_times, t);
		S_new = sort(c(S_i, S_j));
		partition_next[[length(partition_next)+1]] = S_new;
		for(set in partition_curr){
			if(!identical(set, S_i) & !identical(set, S_j)){
				partition_next[[length(partition_next)+1]] = sort(set);
			}
		}
	} else {
		partition_next = partition_curr;
	}
	partitions[[t-1]] = partition_next;
}
partitions[[1]]


# new markov chains
X_new = matrix(NA, nrow=d, ncol=T+1);
Y_new = matrix(NA, nrow=d, ncol=T+1);

has_coupled = numeric(T+1);
slope_side = numeric(T+1);
x = X_new[, 1] = X[, R_prop];
y = Y_new[, 1] = Y[, R_prop];
# Subset coupling + proportional coupling
t = 1
for(t in 1:(T)){
	i = indices[T+1-t, 1]; j = indices[T+1-t, 2];
	x_new = x; y_new = y;
	if((T+1-t+1) %in% marked_times){
		# subset coupling
		slope = (y[i]+y[j])/(x[i]+x[j]);
		if(slope > 1){
			slope_side[t] = 1;
			diff_vec = y-x;
			lambda_y = runif(1);
			set_diff = S.t.1[[T+1-t+1]]; set_diff = set_diff[set_diff != i];
			lambda_x = lambda_y*slope + 1/(x[i]+x[j])*sum(diff_vec[set_diff]);
			if(lambda_x < 1 & lambda_x > 0){
				has_coupled[t] = 1;
			} else {
				lambda_x = runif(1);
			}
		} else {
			slope_side[t] = 2;
			diff_vec = x-y;
			lambda_x = runif(1);
			set_diff = S.t.1[[T+1-t+1]]; set_diff = set_diff[set_diff != i];
			lambda_y = lambda_x/slope + 1/(y[i]+y[j])*sum(diff_vec[set_diff]);
			if(lambda_y < 1 & lambda_y > 0){
				has_coupled[t] = 1;
			} else {
				lambda_y = runif(1);
			}
		}
		if(has_coupled[t]){
			if(slope_side[t] == 1){
				x_new[i] = y_new[i] = lambda_y * (y[i] + y[j]);
				x_new[j] = y_new[j] = (1-lambda_y) * (y[i] + y[j]);
			} else if (slope_side[t] == 2) {
			   	y_new[i] = x_new[i] = lambda_x * (x[i] + x[j]);
				y_new[j] = x_new[j] = (1-lambda_x) * (x[i] + x[j]);
			}
		} else {
			x_new[i] = lambda_x * (x[i] + x[j]);
			x_new[j] = (1-lambda_x) * (x[i] + x[j]);
			y_new[i] = lambda_y * (y[i] + y[j]);
			y_new[j] = (1-lambda_y) * (y[i] + y[j]);
		}
		cat(t, has_coupled[t], x_new[i], x_new[j], y_new[i], y_new[j], "\n");

		# MH Step
		log_rho_x = sum(total_families*log(x_new)) - sum(total_families*log(x))
		log_rho_y = sum(total_families*log(y_new)) - sum(total_families*log(y))

		# can use other couplings, rn using the same thing
		log_unif_x = log_unif_y = log(runif(1))

		if(log_unif_x < min(0, log_rho_x)){
			x = x_new
		}
		if(log_unif_y < min(0, log_rho_y)){
			y = y_new
		}
	}
	else {
		# proportional coupling
		lambda = runif(1);
		x_new[i] = lambda * (x[i] + x[j])
		x_new[j] = (1-lambda) * (x[i] + x[j])
		y_new[i] = lambda * (y[i] + y[j])
		y_new[j] = (1-lambda) * (y[i] + y[j])
		# MH Step
		log_rho_x = sum(total_families*log(x_new)) - sum(total_families*log(x))
		log_rho_y = sum(total_families*log(y_new)) - sum(total_families*log(y))

		# can use other couplings, rn using the same thing
		log_unif_x = log_unif_y = log(runif(1));

		if(log_unif_x < min(0, log_rho_x)){
			x = x_new
		}
		if(log_unif_y < min(0, log_rho_y)){
			y = y_new
		}
	}
	X_new[, 1+t] = x; Y_new[, 1+t] = y;
}

range(X_new[, 1+T] - Y_new[, 1+T])
mean(X_new[, 1+T] - Y_new[, 1+T] == 0)
