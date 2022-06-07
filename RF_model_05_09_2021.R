#Final Project Reed Frost Model of Epidemics
#chain binomial

z=.03		#probability one susceptible will get infected by contact with one infected
		#or probability of "adequate contact".
nsim=20	#number of population trajectories to simulate
T=13		#end time of simulations
inf=matrix(0,T+1,nsim)		#inf starts as matrix of 0s with
					#T+1 rows and nsim columns.  Each
					#column of inf will become a population trajectory.
					#each value in the entire first row
					#of inf is set to 19.  Initial condition.

sus=inf			#sets up a second matrix for the susceptibles.
inf[1,]=2			#initial condition of matrix of infectives.  each value in the first row is set to 2
				#These are the number of infectives at the beginning of each time step.
sus[1,]=100			#initial condition of matrix of susceptibles.  First row is set
				#to 100 infectives.

for (t in 1:T){
	inf[t+1,]=rbinom(nsim,sus[t,],1-(1-z)^inf[t])	#each subsequent row of x 
									#generated from current row.
sus[t+1,]=sus[t,]-inf[t+1,]			#Each new value(s) of susceptibles.
}


pop=cbind(sus,inf)
time=0:T
matplot(time,sus,type="l", main=" Reed Frost Model", col=4, xlab="time step",
ylab="population percentage")				#Matrix plot
								#Plots columns of x versus time.

matlines(time,inf,type="l", col=2)
legend("topright", legend=c("susceptible individuals", "Infected individual"), fill=c(4,2))

