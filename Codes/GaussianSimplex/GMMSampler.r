
set.seed(123)

#####################################################################################
# GMM Sampler

# Define parameters for the Gaussian Mixture Model
d = 3 	# Number of families (dimension of the simplex)
means = c(1, 4, 8)  # Mean values for each component
variances = c(1, 1, 1)  # Variances for each component
proportions = c(0.3, 0.4, 0.3)  # Proportions for each component

# Number of samples to generate
num_samples = 1000

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
#TODO: change this step to be a uniform sample from the simplex using dirichlet distribution
c = c(0.1, 0.2, 0.7)
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

#####################################################################################
# Proportional Coupling
