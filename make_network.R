# Written by Matt Hitchings, 11th Jan 2017
## Updated by Rebecca Kahn May 2021
# Make social network for disease simulation

make_network <- function(ave_community_size, community_size_range, num_communities, rate_within, rate_between, prop_hr) {
  # Function to make a network of closely-connected scommunities that are more sparsely
  # connected with each other. I use a stochastic block model.
  # Inputs:
  # Ave_community_size is the average number of people in a community
  # Community_size_range is the range of community sizes. Currently size of a community
  # is being chosen from a uniform distribution on (ave-range/2, ave+range/2)
  # Num_communities is the number of communities in the study population
  # Note that as the code stands the study population size is not fixed. 
  # rate_within is the probability of an edge between any two nodes within the same community
  # rate_between is the probability of an edge between any two nodes in different communities
  # prop_hr is proportion high risk

  require(NetSurv)
  require(Matrix)
  require(Rlab)
  require(igraph)
  require(deSolve)
  require(reshape2)
  require(ggplot2)

  # Create the network, and assign all members a community number
  community_sizes <- ave_community_size + round(runif(num_communities,-community_size_range/2,community_size_range/2))
  studypop_size <- sum(community_sizes)
  # Currently all communities have the same connectedness, and all communities are equally
  # connected to each other
  within_rates <- diag(nrow=num_communities,ncol=num_communities,x=rate_within)
  between_rates <- matrix(rate_between,nrow=num_communities,ncol=num_communities) -
    diag(nrow=num_communities,ncol=num_communities,x=rate_between)
  rates<-within_rates+between_rates
  
  g <- sample_sbm(studypop_size,rates,community_sizes)
  # Give the nodes a name so that igraph remembers them
  V(g)$name<-1:studypop_size
  V(g)$community<-rep(1:num_communities,community_sizes)
  # Trial status will track whether a node is not in the trial (NA), in the control arm (0) or
  # in the vaccine arm (1)
  ## Enrollment day is the day a node is enrolled into the trial
  ## Symptomatic will track whether they are symptomatic or asymptomatically infected
  V(g)$trialstatus<-NA
  V(g)$enrollmentday<-NA
  V(g)$symptomatic<-NA
  V(g)$risk <- rbinom(studypop_size,1,prop_hr) # risk group
  
  return(g)
  
}
