
set.seed(123)

#####################################################################################
# Libraries
# install.packages("rBeta2009")
library(rBeta2009)

#####################################################################################
# GMM Sampler

# Define parameters for the Gaussian Mixture Model
d = 100 	# Number of families (dimension of the simplex)
# means = c(1, 4, 7, 10, 13, 16, 19, 22, 25, 28)  # Mean values for each component
# variances = numeric(d)+1  # Variances for each component
means = runif(d, -50, 50);
variances = runif(d, 1, 50);
proportions = as.vector(rdirichlet(1, rep(1, d)))  # Proportions for each component

# Number of samples to generate
num_samples = 100000

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
total_families = table(component_families)

# Plot the generated data
hist(generated_samples, breaks = 30, freq = FALSE, main = "Generated Data from GMM", xlab = "Value")
lines(density(generated_samples), col = "blue", lwd = 2)

#####################################################################################

# Gibbs sampler + MH
R = 1e4
burn_in = 1e3
C = matrix(NA, nrow = d, ncol = R)
c = as.vector(rdirichlet(1, rep(1, d)))
for(r in 1:(R + burn_in)){
	if(r%%1000 == 0 & r <= burn_in){
      print(paste("Warm-up: iteration:", r))
    } else if (r%%1000 == 0 & r > burn_in){
      print(paste("Sampling: iteration:", r-burn_in))
    }

	# Proposal from uniform d-simplex using Gibbs
	# TODO: make this into a function
	c_new = c
	lambda = runif(1)
	indices = sample(1:d, 2, replace = FALSE)
	i = indices[1]
	j = indices[2]
	c_new[i] = lambda * (c[i] + c[j])
	c_new[j] = (1-lambda) * (c[i] + c[j])

	# MH step
	log_rho = sum(total_families*log(c_new)) - sum(total_families*log(c))
	log_unif = log(runif(1))

	if(log_unif < min(0, log_rho)){
		c = c_new
	}

	# Store values after burnin done
	if(r > burn_in){
		C[, r - burn_in] = c
	}
}
# pretty stable estimates
cat("Row means =", rowMeans(C), "\nrow variances =", apply(C, 1, var), "\n")
proportions
#####################################################################################
# First Coupling Stage
#####################################################################################
# Proportional Coupling
# R_prop = 1e3
exp_bound = 2
R_prop = 2*ceiling(3/2*exp_bound*d*log(d)) # Trying out the Theorem of Simplex Mixing Time

X = matrix(NA, nrow = d, ncol = R_prop)	# Initialize chain randomly
Y = matrix(NA, nrow = d, ncol = R_prop)	# Initialize chain from the actual distribution using the gibbs sample
X[, 1] = as.vector(rdirichlet(1, rep(1, d))); x = X[, 1]
Y[, 1] = C[, R]; y = Y[, 1]
for(r in 2:R_prop){
	if(r %% 1000 == 0){
		print(paste("iteration: ", r))
	}
	x_new = x; y_new = y
	lambda = runif(1)
	indices = sample(1:d, 2, replace = FALSE)
	i = indices[1]; j = indices[2]
	x_new[i] = lambda * (x[i] + x[j])
	x_new[j] = (1-lambda) * (x[i] + x[j])
	y_new[i] = lambda * (y[i] + y[j])
	y_new[j] = (1-lambda) * (y[i] + y[j])

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
	X[, r] = x; Y[, r] = y
}

print(sum((X[, R_prop]-Y[, R_prop])^2))
print(d^(-exp_bound))

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


# Second coupling stage
#####################################################################################
#TODO: we need a function to generate a shit ton of indices
#TODO: we need a function to perform proportional coupling
#TODO: we need a function to perform subset coupling
#####################################################################################
eps = 1;
T = ceiling((0.5 + eps)*d*log(d));
indices = matrix(NA, nrow = T+1, ncol = 2);
for(i in 2:(T+1)){
	indices[i, ] = sort(sample(1:d, 2, replace = FALSE));
}
# we need a vector called marked_times to store the marked times
# we need a two lists to store the subsets S(t, 1) and S(t, 2)
# we need a list to store the list of partitions at each step
# TODO: this shit's gonna take too much memory, need to reduce it's complexity once it starts working.
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
		cat(t, !identical(S_i, S_j), "\n");
		S.t.1[[t]] = sort(S_i);
		S.t.2[[t]] = sort(S_j);
		marked_times = c(marked_times, t);
		S_new = sort(c(S_i, S_j));
		# partition_next = c(partition_next, c(S_i, S_j));
		partition_next[[length(partition_next)+1]] = S_new;
		for(set in partition_curr){
			if(!identical(set, S_i) & !identical(set, S_j)){
				# partition_next = c(partition_next, set);
				partition_next[[length(partition_next)+1]] = sort(set);
			}
		}
	} else {
		partition_next = partition_curr;
	}
	partitions[[t-1]] = partition_next;
}

# new markov chains
X_new = matrix(NA, nrow=d, ncol=T+1);
Y_new = matrix(NA, nrow=d, ncol=T+1);
x = X_new[, 1] = X[, R_prop];
y = Y_new[, 1] = Y[, R_prop];
# Subset coupling + proportional coupling
t = 1
for(t in 1:(T+1)){
	i = indices[T+1-t, 1]; j = indices[T+1-t, 2];
	x_new = x; y_new = y;
	if((T+1-t+1) %in% marked_times){
		print(t)
		# subset coupling
		slope = (y[i]+y[j])/(x[i]+x[j]);
		if(slope > 1){
			diff_vec = y-x;
			lambda_y = runif(1);
			set_diff = S.t.1[[T+1-t+1]]; set_diff = set_diff[set_diff != i];
			lambda_x = lambda_y*slope + 1/(x[i]+x[j])*sum(diff_vec[set_diff]);
			print(sum(diff_vec[set_diff]))
			cat("x", lambda_x, "\n");
			lambda_x = ifelse(lambda_x < 1 & lambda_x > 0, lambda_x, runif(1));
		} else {
			diff_vec = x-y;
			lambda_x = runif(1);
			set_diff = S.t.1[[T+1-t+1]]; set_diff = set_diff[set_diff != i];
			lambda_y = lambda_x/slope + 1/(y[i]+y[j])*sum(diff_vec[set_diff]);
			cat("y", lambda_y, "\n");
			lambda_y = ifelse(lambda_y < 1 & lambda_y > 0, lambda_y, runif(1));
		}
		x_new[i] = lambda_x * (x[i] + x[j]);
		x_new[j] = (1-lambda_x) * (x[i] + x[j]);
		y_new[i] = lambda_y * (y[i] + y[j]);
		y_new[j] = (1-lambda_y) * (y[i] + y[j]);
		cat(t, x_new[i], x_new[j], y_new[i], y_new[j], "\n");

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
		# x = x_new; y = y_new;
	}
	else {
		# proportional coupling
		lambda = runif(1);
		x_new[i] = lambda * (x[i] + x[j])
		x_new[j] = (1-lambda) * (x[i] + x[j])
		y_new[i] = lambda * (y[i] + y[j])
		y_new[j] = (1-lambda) * (y[i] + y[j])
		cat("prop", t, all(x_new > 0), all(y_new > 0), "\n");
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
		# x = x_new; y = y_new;
	}
	X[, 1+t] = x; Y[, 1+t] = y;
}

(X[, 1+T] - Y[, 1+T])
