# R program to calculate maximum likelihood estimate of parameter
# beta in negative exponential model of mating probability as a
# function of population density (product binomial likelihood).

# Data from Nishigaki, J.  1963.  The effect of low population
# density on the mating chance and the fecundity of the azuki bean
# weevil, Callosobruchus chinensis L.  Japanese Journal of Ecology
# 13:178-184).

# Model and analysis from Dennis, B.  1989.  Allee effects: population
# growth, critical density, and the chance of extinction.  Natural
# Resource Modeling 3:481-538.

# Data and optimization settings.
k=c(20,20,21,20,20,24,21,24,27,20)           # Total number females.
n=c(1,2,3,4,5,6,7,8,9,10)                    # Number males per 100 cm^2.
y=c(7,6,11,12,12,19,16,19,25,19)             # Number females mated.

betalo=.1    # Boundaries of interval on which objective
betahi=.3    # function is to be optimized.

# Define the product binomial log likelihood function.
lnlike=function(beta,K,N,Y)
  sum(lfactorial(K)-lfactorial(Y)-lfactorial(K-Y)+
        Y*log(1-exp(-beta*N))-(K-Y)*beta*N)

# R function optimize():  one-dimensional optimization function.
ML=optimize(f=lnlike,K=k,N=n,Y=y,lower=betalo,upper=betahi,
            maximum=TRUE)
beta=ML$maximum               #  ML estimate of beta.
log.likelihood=ML$objective   #  Maximized log-likelihood.

# Print the results.
beta
log.likelihood

# Calculate mating probability curve as well as observed mating
# probabilities for plotting.
nn=seq(0,11,.01)
prob.mating=1-exp(-beta*nn)
obs.probs=y/k

# Plot the results.
# First plot is the model 9does Response Model)
# the points command adds the observed data in the form of proportions
plot(nn,prob.mating,type="l",pch=1,cex=1.5,xlab="male density",
     ylab="female mating probability")		#pch is the plotting character
points(n,obs.probs,type="p")		#points command adds points to an open graph