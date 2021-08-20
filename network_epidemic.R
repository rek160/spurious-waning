# Written Matt Hitchings, 20th Dec 2016
## Edited by Rebecca Kahn, May 2021
# Simulate disease on a social network which includes a certain proportion of asymptomatic infection

# Parameters
library(tidyverse)  
library(deSolve)

setwd("/n/home00/rkahn/VE")
#setwd("/Users/rebeccakahn/Dropbox/Git_repo/Asymptomatic-Infection-Vaccine-Efficacy v2")

args=(commandArgs(TRUE))

j<-1
ave_community_size<-20000#4000
community_size_range<-0 #400
num_communities<-1 #5
rate_within<- 0.0018174 # 0.002182, 0.0018174, 0.001455
rate_between<-0 
## betaS is force of infection in symptomatics and betaAS in asymptomatics -- simplifying assumption is that they are the same
betaS<-0.005
betaAS<-0.005
## pS is proportion of the infected in the control group that becomes symptomatic and vpS is the proportion in the vaccine group
pS<-0.5
vpS<-0.5
VT <- "leaky"
## introductions at beginning
num_introductions<-10
direct_VE<-0.6
rel_VE_PI <- 1 # relative VE of those previously infected
incperiod_shape<-3000 
incperiod_rate<-600
infperiod_shape<-4200
infperiod_rate<-600
trial_startday<-1
trial_length<-365 # if assume protection only lasts 90 days
enrollment_period<-1
num_ind_enrolled_per_day<-num_communities/enrollment_period ## number of clusters enrolled per day since using stratified indivdual randomization
sample_percent<-0.1
cluster_coverage<-0.25
ave_inc_period<-ceiling(incperiod_shape/incperiod_rate)
ave_inf_period<-ceiling(infperiod_shape/infperiod_rate)
prop_rec <- 0.33 # proportion previously infected

prop_hr <- 0.5 # proportion high risk
hr <- 1.5 # increase in beta for high risk
lr <- 0.5 # decrease in beta for low risk
days_reinf <- 90 # number of days before have chance of reinfection
reinf <- 1 # reduction in probability of reinfection 

R0 <- (1 - (infperiod_rate/(infperiod_rate+betaS))^infperiod_shape) *
  (((ave_community_size-1)*(1-rate_within)*rate_within +
      (num_communities-1)*ave_community_size*(1-rate_between)*rate_between +
      ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)^2)/
     ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)-1)
R0

for (i in 1:length(args)) { 
  eval (parse (text = args[[i]] )) 
} 

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

source("make_network.R")

network_epidemic<-function(g,betaS,betaAS,num_introductions,VE,rel_VE_PI,hr,lr,days_reinf,reinf,
                           incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,pS,vpS,prop_rec,
                           bTrial,bCluster,VT,
                           trial_startday,trial_length,
                           num_ind_enrolled_per_day,enrollment_period,cluster_coverage) {
  # Inputs:
  # g - the graph to run the epidemic on
  ## betaS - In the current model, every infectious individual contacts all their neighbours in a time step
  ## and infects each susceptible with probability betaS if they are symptomatic. So betaS represents the probability 
  ## of infection from one contact between a symptomatic infectious and a susceptible.
  ## betaAS - In the current model, every infectious individual contacts all their neighbours in a time step
  # and infects each susceptible with probability betaAS if they are asymptomatic. So betaAS represents the probability 
  ## of infection from one contact between an asymptomatic infectious and a susceptible.
  # num_introductions - how many separate introductions we expect on average from the main epidemic. This is used to calibrate
  # the external force of infection
  # VE - direct leaky efficacy of the vaccine
  # inc/inf period shape/rate - parameters for the incubation period and the gamma distributed infectious period
  ## pS - In the current model, every infected individual has a probability of becoming infectious and symptomatic 
  ## and a probability of becoming infectious and asymptomatic. pS represents the probability an infected person will 
  ## develop symptoms after the incubation period (assumed to be same for symptomatic and asymptomatic).
  ## Note: 1-pS is the probability an infected person will become infectious and asymptomatic
  ## vpS - In the current model, every vaccinated infected individual has a probability of becoming infectious and symptomatic 
  ## and a probability of becoming infectious and asymptomatic.vpS represents the probability a vaccinated person who gets 
  ## infected will develop symptoms after the incubation period (assumed to be same for symptomatic and asymptomatic).
  ## Note: 1-vpS is the probability a vaccinated person who gets infected will be asymptomatic. 
  # bTrial - whether we are running a trial or not
  # bCluster - indicator of whether we are running the cRCT (1) or the iRCT (0) -- in this scenario it is always set to 0 as we are always
  # VT is vaccine type - leaky or all or nothing
  ## running an iRCT
  # trial_startday - first day of trial enrollment
  # trial_length - end of follow-up of trial partcipants, counting from the first day of enrollment
  # num_enrolled_per_day - number of clusters enrolled / individually randomized per day
  # enrollment_period - length of enrollment period
  # cluster_coverage -  The proportion of each cluster we expect to enroll -- note we are using stratified individual randomization so 
  ## essentially the same number of people from each community are enrolled and invidiually randomized
  
  require(NetSurv)
  require(Matrix)
  require(Rlab)
  require(igraph)
  require(deSolve)
  require(reshape2)
  require(ggplot2)
  require(caTools)
  
  list <- structure(NA,class="result")
  "[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
      a <- args[[i]]
      if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
  }
  
  # Recover and spread functions
  recover<-function(e_nodes,ev_nodes,iS_nodes,iAS_nodes,iS_nodes_susc,iS_nodes_vacc,iAS_nodes_susc,iAS_nodes_vacc,
                    r_nodes,infperiod_shape,infperiod_rate,pS,vpS,t) {
    # Input is a list of the exposed nodes that were not vaccinated (e_nodes) and exposed nodes that were vaccinated (ev_nodes), 
    # with number of days since infection and total incubation/latent
    # period, and equivalently for the infectious nodes.
    ## iS_nodes are the symptomatic infectious nodes; iAS nodes are the asymptomatic infectious nodes;
    ## 'susc' means the nodes were not vaccinated; 'vacc' means the nodes were vaccinated.
    ## Therefore, iS_nodes_vacc are the nodes that were vaccinated, infected and became symptomatic
    # For each of the iS and iAS nodes, we will add it to recovered (symptomatic or asymptomatic) if the number of days infected has
    # reached the total length of the infectious period
    ## For each of the e_nodes and ev_nodes, we will advance to infectious compartments when their number of days exposed has 
    ## reached the total length of the incubation period
    
    # Advance infectious nodes first, otherwise you will doubly advance any nodes switching from
    # exposed to infectious at this time step
    indices_to_removeS <- iS_nodes[2,]>=iS_nodes[3,]
    indices_to_removeAS <- iAS_nodes[2,]>=iAS_nodes[3,]
    newremoved<-c(as.vector(iS_nodes[1,indices_to_removeS]),as.vector(iAS_nodes[1,indices_to_removeAS]))
    
    # Add one day to the length of each infected individual's time infected
    iS_nodes[2,] <- iS_nodes[2,]+1
    iAS_nodes[2,] <- iAS_nodes[2,]+1
    
    # Remove any recovered from iS_nodes and iAS_nodes and add to r_nodes
    ## Also remove recovered nodes from iS and iAS susc and vacc lists
    iS_nodes <- iS_nodes[,!(iS_nodes[1,] %in% newremoved),drop=FALSE]
    iAS_nodes <- iAS_nodes[,!(iAS_nodes[1,] %in% newremoved),drop=FALSE]
    iS_nodes_susc <- setdiff(iS_nodes_susc, newremoved)
    iS_nodes_vacc <- setdiff(iS_nodes_vacc, newremoved)
    iAS_nodes_susc <- setdiff(iAS_nodes_susc, newremoved)
    iAS_nodes_vacc <- setdiff(iAS_nodes_vacc, newremoved)
    r_nodes <- cbind(r_nodes,rbind(newremoved,rep(t,length(newremoved))))
    
    # Now advance exposed unvaccinated nodes
    ## Determine if exposed nodes will become symptomatic (with probability pS) or asymptomatic
    newinfectious_S <- NULL
    newinfectious_AS <- NULL
    newinfectious <- NULL
    if (length(e_nodes)>0) {
      indices_to_remove <- e_nodes[2,]>=e_nodes[3,]
      newinfectious<-as.vector(e_nodes[1,indices_to_remove])
      if (length(newinfectious)>0) {
        num_symptomatic <-rbinom(1,length(newinfectious),pS)
        newinfectious_S<-newinfectious[sample.int(length(newinfectious),num_symptomatic)]
        newinfectious_AS<-setdiff(newinfectious,newinfectious_S)
      }
    }
    
    ## Same above for some vaccinated (probability of being symptomatic is vpS)
    ## Because vaccination status affects probability of becoming symptomatic or asymptomatic, 
    ## I have created ev_nodes for vaccinated and infected nodes
    newinfectious_S_v <- NULL
    newinfectious_AS_v <- NULL
    newinfectious_v <- NULL
    if (length(ev_nodes)>0) {
      indices_to_remove_v <- ev_nodes[2,]>=ev_nodes[3,]
      newinfectious_v<-as.vector(ev_nodes[1,indices_to_remove_v])
      if (length(newinfectious_v)>0) {
        num_symptomatic_v<-rbinom(1,length(newinfectious_v),vpS)
        newinfectious_S_v<-newinfectious_v[sample.int(length(newinfectious_v),num_symptomatic_v)]
        newinfectious_AS_v<-setdiff(newinfectious_v,newinfectious_S_v)
      }
    }  
    
    ## list of newly infectious symptomatic and asymptomatic nodes (combining vaccinated and unvaccinated)
    newinfectious_S_total <- c(newinfectious_S,newinfectious_S_v)
    newinfectious_AS_total <- c(newinfectious_AS,newinfectious_AS_v)
    newinfectious_total <- c(newinfectious_S_total,newinfectious_AS_total)
    
    # Add one day to the length of each infected individual's time infected
    e_nodes[2,] <- e_nodes[2,]+1
    ev_nodes[2,] <- ev_nodes[2,]+1
    
    # Remove any progressing from e_nodes and ev_nodes and add to iS_nodes or iAS_nodes
    e_nodes <- e_nodes[,!(e_nodes[1,] %in% newinfectious),drop=FALSE]
    ev_nodes <- ev_nodes[,!(ev_nodes[1,] %in% newinfectious_v),drop=FALSE]
    # Give each newly infected node an infectious period
    inf_periods_S <- rgamma(length(newinfectious_S_total),infperiod_shape,infperiod_rate)
    inf_periods_AS <- rgamma(length(newinfectious_AS_total),infperiod_shape,infperiod_rate)
    iS_nodes <- cbind(iS_nodes,rbind(newinfectious_S_total,rep(0,length(newinfectious_S_total)),inf_periods_S))
    iAS_nodes <- cbind(iAS_nodes,rbind(newinfectious_AS_total,rep(0,length(newinfectious_AS_total)),inf_periods_AS))
    
    ## Also update the lists distinguishing susceptible and vaccinated iS and iAS nodes
    iS_nodes_susc <- c(iS_nodes_susc,newinfectious_S)
    iS_nodes_vacc <- c(iS_nodes_vacc,newinfectious_S_v)
    iAS_nodes_susc <- c(iAS_nodes_susc,newinfectious_AS)
    iAS_nodes_vacc <- c(iAS_nodes_vacc,newinfectious_AS_v)
    
    i_nodes<-cbind(iS_nodes,iAS_nodes)
    
    list(e_nodes, ev_nodes,iS_nodes,iAS_nodes,i_nodes,iS_nodes_susc,iS_nodes_vacc,iAS_nodes_susc,iAS_nodes_vacc,
         r_nodes,newinfectious_S_total,newinfectious_AS_total,sort(newinfectious_total))
  }
  
  spread<-function(g,s_nodes,nv_nodes,v_nodes,e_nodes,ev_nodes,iS_nodes,iAS_nodes, 
                   iS_nodes_susc,iS_nodes_vacc,iAS_nodes_susc,iAS_nodes_vacc,r_nodes,
                   betaS,betaAS,VE,rel_VE_PI,hr,lr,days_reinf,reinf,
                   incperiod_shape,incperiod_rate,
                   connected_nodes,external_inf_F,source_num_inf,t_step){
   
    # Spread will create new infected nodes from two sources: infectious nodes within the study
    # population, and external pressure from the source population
    # Inputs:
    # g is the graph, used to find neighbours of infected nodes
    ## s_nodes are all susceptible nodes in population
    ## nv_nodes are susceptible, non-vaccinated nodes enrolled in the trial
    ## v_nodes are vaccinated nodes enrolled in the trial
    ## e_nodes and ev_nodes are exposed susceptible and exposed vaccinated nodes
    ## iS_nodes and iAS_nodes are infected symptomatic and infected asymptomatic nodes, respectively
    ## betaS is the hazard of infection for one contact with symptomatic infectious
    ## betaAS is the hazard of infection for one contact with asymptomatic infectious
    # VE is direct VE 
    # incperiod_shape and rate are used to assign each newly exposed node a latent/incubation period
    # length, currently drawn from a gamma distribution
    # connected_nodes is a list of nodes that are connected the the source population (under the current
    # model this is just all nodes)
    # external_inf_F is a constant of proportionality that defines infectious pressure from source pop
    # source_num_inf is the number of infectious individuals in the source population
    # t_step is the time since the simulation start
    
    # Process: go through list of iS_nodes and iAS_nodes, and choose a random number of their susceptible neighbours to
    # be infected, according to betaS and betaAS and choose a random number of their susceptible vaccinated neighbours to
    # be infected, according to betaS(1-VE) and betaAS(1-VE)
    # Then go through list of nodes that are connected to source population and infect each susceptible
    # one with probability 1-exp(-FI), where I is the number/proportion of infectious, and F is a constant
    betaS_v <- betaS*(1-VE)
    betaAS_v <- betaAS*(1-VE)
    highrisk <- V(g)[V(g)$risk==1]$name
    lowrisk <- V(g)[V(g)$risk==0]$name
    infectees_susc_S<-c()
    infectees_susc_AS<-c()
    infectees_vacc_S<-c()
    infectees_vacc_AS<-c()
    
    # which recovered nodes eligible for reinfection 
    r_nodes_reinf <- t > (r_nodes[2,] + days_reinf) 
    r_nodes_s <-as.vector(r_nodes[1,r_nodes_reinf]) # those susceptible again
    
    if (VT == "aon"){
      r_nodes_nv <- r_nodes_s
      r_nodes_v <- data.frame()
    } else if (VT == "leaky"){ 
      r_nodes_nv <- intersect(r_nodes_s,V(g)[which(V(g)$trialstatus==0)]$name)
      r_nodes_v <- intersect(r_nodes_s,V(g)[which(V(g)$trialstatus==1)]$name)
    }
    
    if (ncol(iS_nodes)>0) {
      # Get a list of all neighbours of all infected symptomatic nodes
      potential_contacts_S<-lapply(iS_nodes[1,],function(x) neighbors(g,x))
      susc_contacts_S<-lapply(potential_contacts_S,function(x,susceptibles) intersect(x,susceptibles),susceptibles=s_nodes)
      susc_contacts_S_high <- lapply(susc_contacts_S,function(x,highrisk) intersect(x,highrisk),highrisk=highrisk)
      susc_contacts_S_low <- lapply(susc_contacts_S,function(x,lowrisk) intersect(x,lowrisk),lowrisk=lowrisk)
      num_neighbours_susc_S_high<-rapply(susc_contacts_S_high,length)
      num_neighbours_susc_S_low<-rapply(susc_contacts_S_low,length)
      # Sample from each group of neighbours in turn
      # First choose how many high risk neighbours each node infects
      num_contacts_susc_S_high<-rbinom(length(num_neighbours_susc_S_high),num_neighbours_susc_S_high,1-exp(-betaS*hr))
      # Then sample from the neighbours
      # If one node gets picked twice by different nodes, just discard the duplicate.
      # In the rare case that each i_nodes makes a number of new infectees equal to the number
      # of infectious nodes, mapply will make a matrix and unlist won't work. Therefore, the c() around
      # it ensures we turn the matrix into a vector. Unique then removes duplicates.
      infectees_susc_S_high<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=susc_contacts_S_high,y=num_contacts_susc_S_high)))
      infectees_susc_S_high<-unique(infectees_susc_S_high)
      
      # now choose low
      num_contacts_susc_S_low<-rbinom(length(num_neighbours_susc_S_low),num_neighbours_susc_S_low,1-exp(-betaS*lr))
      infectees_susc_S_low<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=susc_contacts_S_low,y=num_contacts_susc_S_low)))
      infectees_susc_S_low<-unique(infectees_susc_S_low)
      
      # combine
      infectees_susc_S <- c(infectees_susc_S_high,infectees_susc_S_low)
      
      if (length(v_nodes)>0) {
        # Same as above but with vaccinated susceptible nodes
        vacc_contacts_S<-lapply(potential_contacts_S,function(x,vacc) intersect(x,vacc),vacc=v_nodes)
        
        vacc_contacts_S_high <- lapply(vacc_contacts_S,function(x,highrisk) intersect(x,highrisk),highrisk=highrisk)
        vacc_contacts_S_low <- lapply(vacc_contacts_S,function(x,lowrisk) intersect(x,lowrisk),lowrisk=lowrisk)
        num_neighbours_vacc_S_high<-rapply(vacc_contacts_S_high,length)
        num_neighbours_vacc_S_low<-rapply(vacc_contacts_S_low,length)
        num_contacts_vacc_S_high<-rbinom(length(num_neighbours_vacc_S_high),num_neighbours_vacc_S_high,1-exp(-betaS_v*hr))
        infectees_vacc_S_high<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=vacc_contacts_S_high,y=num_contacts_vacc_S_high)))
        infectees_vacc_S_high<-unique(infectees_vacc_S_high)
        num_contacts_vacc_S_low<-rbinom(length(num_neighbours_vacc_S_low),num_neighbours_vacc_S_low,1-exp(-betaS_v*lr))
        infectees_vacc_S_low<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=vacc_contacts_S_low,y=num_contacts_vacc_S_low)))
        infectees_vacc_S_low<-unique(infectees_vacc_S_low)
        infectees_vacc_S <- c(infectees_vacc_S_high,infectees_vacc_S_low)
        
        
      } else {
        infectees_vacc_S<-c()
      }
      

      
      if (length(r_nodes_nv)>0) {
        
        # Same as above but with vaccinated susceptible nodes
        rec_contacts_S_nv<-lapply(potential_contacts_S,function(x,rec) intersect(x,rec),rec=r_nodes_nv)
        
        rec_contacts_S_nv_high <- lapply(rec_contacts_S_nv,function(x,highrisk) intersect(x,highrisk),highrisk=highrisk)
        rec_contacts_S_nv_low <- lapply(rec_contacts_S_nv,function(x,lowrisk) intersect(x,lowrisk),lowrisk=lowrisk)
        num_neighbours_rec_S_nv_high<-rapply(rec_contacts_S_nv_high,length)
        num_neighbours_rec_S_nv_low<-rapply(rec_contacts_S_nv_low,length)
        num_contacts_rec_S_nv_high<-rbinom(length(num_neighbours_rec_S_nv_high),num_neighbours_rec_S_nv_high,1-exp(-betaS*hr*reinf))
        infectees_rec_S_nv_high<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=rec_contacts_S_nv_high,y=num_contacts_rec_S_nv_high)))
        infectees_rec_S_nv_high<-unique(infectees_rec_S_nv_high)
        num_contacts_rec_S_nv_low<-rbinom(length(num_neighbours_rec_S_nv_low),num_neighbours_rec_S_nv_low,1-exp(-betaS*lr*reinf))
        infectees_rec_S_nv_low<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=rec_contacts_S_nv_low,y=num_contacts_rec_S_nv_low)))
        infectees_rec_S_nv_low<-unique(infectees_rec_S_nv_low)
        infectees_rec_S_nv <- c(infectees_rec_S_nv_high,infectees_rec_S_nv_low)
        
        
      } else {
        infectees_rec_S_nv<-c()
      }
      
      
      if (length(r_nodes_v)>0) {
        
        # Same as above but with vaccinated susceptible nodes
        rec_contacts_S_v<-lapply(potential_contacts_S,function(x,rec) intersect(x,rec),rec=r_nodes_v)
        
        rec_contacts_S_v_high <- lapply(rec_contacts_S_v,function(x,highrisk) intersect(x,highrisk),highrisk=highrisk)
        rec_contacts_S_v_low <- lapply(rec_contacts_S_v,function(x,lowrisk) intersect(x,lowrisk),lowrisk=lowrisk)
        num_neighbours_rec_S_v_high<-rapply(rec_contacts_S_v_high,length)
        num_neighbours_rec_S_v_low<-rapply(rec_contacts_S_v_low,length)
        num_contacts_rec_S_v_high<-rbinom(length(num_neighbours_rec_S_v_high),num_neighbours_rec_S_v_high,1-exp(-betaS_v*hr*reinf*rel_VE_PI))
        infectees_rec_S_v_high<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=rec_contacts_S_v_high,y=num_contacts_rec_S_v_high)))
        infectees_rec_S_v_high<-unique(infectees_rec_S_v_high)
        num_contacts_rec_S_v_low<-rbinom(length(num_neighbours_rec_S_v_low),num_neighbours_rec_S_v_low,1-exp(-betaS_v*lr*reinf*rel_VE_PI))
        infectees_rec_S_v_low<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=rec_contacts_S_v_low,y=num_contacts_rec_S_v_low)))
        infectees_rec_S_v_low<-unique(infectees_rec_S_v_low)
        infectees_rec_S_v <- c(infectees_rec_S_v_high,infectees_rec_S_v_low)
        
        
      } else {
        infectees_rec_S_v<-c()
      }
      
    } else {
      infectees_susc_S<-c()
      infectees_vacc_S<-c()
      infectees_rec_S_nv<-c()
      infectees_rec_S_v<-c()
    }
    
    if (ncol(iAS_nodes)>0) {  
      ## Repeat above for all infected asymptomatic nodes
      potential_contacts_AS<-lapply(iAS_nodes[1,],function(x) neighbors(g,x))
      susc_contacts_AS<-lapply(potential_contacts_AS,function(x,susceptibles) intersect(x,susceptibles),susceptibles=s_nodes)
      susc_contacts_AS_high <- lapply(susc_contacts_AS,function(x,highrisk) intersect(x,highrisk),highrisk=highrisk)
      susc_contacts_AS_low <- lapply(susc_contacts_AS,function(x,lowrisk) intersect(x,lowrisk),lowrisk=lowrisk)
      num_neighbours_susc_AS_high<-rapply(susc_contacts_AS_high,length)
      num_neighbours_susc_AS_low<-rapply(susc_contacts_AS_low,length)
      # Sample from each group of neighbours in turn
      # First choose how many high risk neighbours each node infects
      num_contacts_susc_AS_high<-rbinom(length(num_neighbours_susc_AS_high),num_neighbours_susc_AS_high,1-exp(-betaAS*hr))
      # Then sample from the neighbours
      # If one node gets picked twice by different nodes, just discard the duplicate.
      # In the rare case that each i_nodes makes a number of new infectees equal to the number
      # of infectious nodes, mapply will make a matrix and unlist won't work. Therefore, the c() around
      # it ensures we turn the matrix into a vector. Unique then removes duplicates.
      infectees_susc_AS_high<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=susc_contacts_AS_high,y=num_contacts_susc_AS_high)))
      infectees_susc_AS_high<-unique(infectees_susc_AS_high)
      
      # now choose low
      num_contacts_susc_AS_low<-rbinom(length(num_neighbours_susc_AS_low),num_neighbours_susc_AS_low,1-exp(-betaAS*lr))
      infectees_susc_AS_low<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=susc_contacts_AS_low,y=num_contacts_susc_AS_low)))
      infectees_susc_AS_low<-unique(infectees_susc_AS_low)
      
      # combine
      infectees_susc_AS <- c(infectees_susc_AS_high,infectees_susc_AS_low)
      
      if (length(v_nodes)>0) {
        ## asymptomatic vaccinated
        vacc_contacts_AS<-lapply(potential_contacts_AS,function(x,vacc) intersect(x,vacc),vacc=v_nodes)
        vacc_contacts_AS_high <- lapply(vacc_contacts_AS,function(x,highrisk) intersect(x,highrisk),highrisk=highrisk)
        vacc_contacts_AS_low <- lapply(vacc_contacts_AS,function(x,lowrisk) intersect(x,lowrisk),lowrisk=lowrisk)
        num_neighbours_vacc_AS_high<-rapply(vacc_contacts_AS_high,length)
        num_neighbours_vacc_AS_low<-rapply(vacc_contacts_AS_low,length)
        # Sample from each group of neighbours in turn
        # First choose how many high risk neighbours each node infects
        num_contacts_vacc_AS_high<-rbinom(length(num_neighbours_vacc_AS_high),num_neighbours_vacc_AS_high,1-exp(-betaAS_v*hr))
        # Then sample from the neighbours
        # If one node gets picked twice by different nodes, just discard the duplicate.
        # In the rare case that each i_nodes makes a number of new infectees equal to the number
        # of infectious nodes, mapply will make a matrix and unlist won't work. Therefore, the c() around
        # it ensures we turn the matrix into a vector. Unique then removes duplicates.
        infectees_vacc_AS_high<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=vacc_contacts_AS_high,y=num_contacts_vacc_AS_high)))
        infectees_vacc_AS_high<-unique(infectees_vacc_AS_high)
        
        # now choose low
        num_contacts_vacc_AS_low<-rbinom(length(num_neighbours_vacc_AS_low),num_neighbours_vacc_AS_low,1-exp(-betaAS_v*lr))
        infectees_vacc_AS_low<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=vacc_contacts_AS_low,y=num_contacts_vacc_AS_low)))
        infectees_vacc_AS_low<-unique(infectees_vacc_AS_low)
        
        # combine
        infectees_vacc_AS <- c(infectees_vacc_AS_high,infectees_vacc_AS_low)
        
      } else {
        infectees_vacc_AS<-c()
      }
      

      if (length(r_nodes_nv)>0) {
        
        # Same as above but with vaccinated susceptible nodes
        rec_contacts_AS_nv<-lapply(potential_contacts_AS,function(x,rec) intersect(x,rec),rec=r_nodes_nv)
        
        rec_contacts_AS_nv_high <- lapply(rec_contacts_AS_nv,function(x,highrisk) intersect(x,highrisk),highrisk=highrisk)
        rec_contacts_AS_nv_low <- lapply(rec_contacts_AS_nv,function(x,lowrisk) intersect(x,lowrisk),lowrisk=lowrisk)
        num_neighbours_rec_AS_nv_high<-rapply(rec_contacts_AS_nv_high,length)
        num_neighbours_rec_AS_nv_low<-rapply(rec_contacts_AS_nv_low,length)
        num_contacts_rec_AS_nv_high<-rbinom(length(num_neighbours_rec_AS_nv_high),num_neighbours_rec_AS_nv_high,1-exp(-betaAS*hr*reinf))
        infectees_rec_AS_nv_high<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=rec_contacts_AS_nv_high,y=num_contacts_rec_AS_nv_high)))
        infectees_rec_AS_nv_high<-unique(infectees_rec_AS_nv_high)
        num_contacts_rec_AS_nv_low<-rbinom(length(num_neighbours_rec_AS_nv_low),num_neighbours_rec_AS_nv_low,1-exp(-betaAS*lr*reinf))
        infectees_rec_AS_nv_low<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=rec_contacts_AS_nv_low,y=num_contacts_rec_AS_nv_low)))
        infectees_rec_AS_nv_low<-unique(infectees_rec_AS_nv_low)
        infectees_rec_AS_nv <- c(infectees_rec_AS_nv_high,infectees_rec_AS_nv_low)
        
        
      } else {
        infectees_rec_AS_nv<-c()
      }
      

      if (length(r_nodes_v)>0) {
        
        # Same as above but with vaccinated susceptible nodes
        rec_contacts_AS_v<-lapply(potential_contacts_AS,function(x,rec) intersect(x,rec),rec=r_nodes_v)
        
        rec_contacts_AS_v_high <- lapply(rec_contacts_AS_v,function(x,highrisk) intersect(x,highrisk),highrisk=highrisk)
        rec_contacts_AS_v_low <- lapply(rec_contacts_AS_v,function(x,lowrisk) intersect(x,lowrisk),lowrisk=lowrisk)
        num_neighbours_rec_AS_v_high<-rapply(rec_contacts_AS_v_high,length)
        num_neighbours_rec_AS_v_low<-rapply(rec_contacts_AS_v_low,length)
        num_contacts_rec_AS_v_high<-rbinom(length(num_neighbours_rec_AS_v_high),num_neighbours_rec_AS_v_high,1-exp(-betaAS_v*hr*reinf*rel_VE_PI))
        infectees_rec_AS_v_high<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=rec_contacts_AS_v_high,y=num_contacts_rec_AS_v_high)))
        infectees_rec_AS_v_high<-unique(infectees_rec_AS_v_high)
        num_contacts_rec_AS_v_low<-rbinom(length(num_neighbours_rec_AS_v_low),num_neighbours_rec_AS_v_low,1-exp(-betaAS_v*lr*reinf*rel_VE_PI))
        infectees_rec_AS_v_low<-c(unlist(mapply(function(x,y) x[sample.int(length(x),y)],x=rec_contacts_AS_v_low,y=num_contacts_rec_AS_v_low)))
        infectees_rec_AS_v_low<-unique(infectees_rec_AS_v_low)
        infectees_rec_AS_v <- c(infectees_rec_AS_v_high,infectees_rec_AS_v_low)
        
        
      } else {
        infectees_rec_AS_v<-c()
      }
      

    } else {
      infectees_susc_AS<-c()
      infectees_vacc_AS<-c()
      infectees_rec_AS_nv<-c()   
      infectees_rec_AS_v<-c()
    }

    infectees_susc<- unique(c(infectees_susc_S,infectees_susc_AS))
    infectees_vacc<- unique(c(infectees_vacc_S,infectees_vacc_AS))
    infectees_rec_nv <- unique(c(infectees_rec_S_nv,infectees_rec_AS_nv))
    infectees_rec_v <- unique(c(infectees_rec_S_v,infectees_rec_AS_v))

    newinfected_susc <- c(infectees_susc)
    newinfected_vacc <- c(infectees_vacc)
    newinfected <- unique(c(newinfected_susc, newinfected_vacc,infectees_rec_nv,infectees_rec_v))
    
    ## Double check no duplicated nodes                                   
    if ((anyDuplicated(newinfected))>0) {
      sink("Error.txt",append=TRUE)
      print(iS_nodes[1,])
      print(iAS_nodes[1,])
      print(e_nodes[1,])
      print(ev_nodes[1,])
      print(s_nodes)
      print(v_nodes)
      print(infectees_susc)
      print(conn_inf_susc)
      print(infectees_vacc)
      print(conn_inf_vacc)
    }
    
    if (length(newinfected)>0) {
      
      # Give each newly exposed node an incubation/latent period
      ## Under current parameters, incubation period will be essentially 6 for all but infectious period will still have distribution
      ## Vacc and non vaccinated have same incubation period here
      inc_periods_susc <- rgamma(length(newinfected_susc),incperiod_shape,incperiod_rate)
      inc_periods_vacc <- rgamma(length(newinfected_vacc),incperiod_shape,incperiod_rate)
      
      inc_periods_rec_nv <- rgamma(length(infectees_rec_nv),incperiod_shape,incperiod_rate)
      inc_periods_rec_v <- rgamma(length(infectees_rec_v),incperiod_shape,incperiod_rate)

      ## Add them to e_nodes and ev_nodes and remove from s_nodes, nv_nodes and v_nodes
      e_nodes <- cbind(e_nodes,rbind(newinfected_susc,rep(0,length(newinfected_susc)),inc_periods_susc))
      e_nodes <- cbind(e_nodes,rbind(infectees_rec_nv,rep(0,length(infectees_rec_nv)),inc_periods_rec_nv))
      
      s_nodes<-setdiff(s_nodes,newinfected)
      nv_nodes<-setdiff(nv_nodes,newinfected)
      v_nodes <- setdiff(v_nodes,newinfected)
      
      r_nodes <- r_nodes[,!(r_nodes[1,] %in% newinfected),drop=FALSE]
      
      
      if (length(newinfected_vacc)>0){
        ev_nodes <- cbind(ev_nodes,rbind(newinfected_vacc,rep(0,length(newinfected_vacc)),inc_periods_vacc))
      }
      
      if (length(infectees_rec_v)>0){
        ev_nodes <- cbind(ev_nodes,rbind(infectees_rec_v,rep(0,length(infectees_rec_v)),inc_periods_rec_v))
      }
    }
    
    list(s_nodes,nv_nodes,v_nodes,e_nodes,ev_nodes,r_nodes)
  }
  
                                     
  studypop_size<-length(V(g))
 
  # Number of timesteps to run the epidemic - only need to go until the end of the trial
  num_timesteps <- trial_startday + trial_length + enrollment_period - 1
  
  if (bTrial) {
    
    # Parameters to do with trial recruitment
    # For iRCT, enrollment per day is number of people enrolled per day
    # For equal assignment to each arm these need to be even numbers.
    enrollment_schedule <- rep(num_ind_enrolled_per_day,enrollment_period)
    enroll_endday <- trial_startday+enrollment_period-1
    
    non_trial_clusters <- 1:max(V(g)$community)
    
  }
  
  ## Initialize the S, E, EV, IS, IAS and R nodes. I seed the epidemic from an SEIR curve in a source population,
  # so initially all nodes in the study population are susceptible
  ## e_nodes,ev_nodes,iS_nodes and iAS_nodes are matrices. The first row is the identity of the node. The second row
  # is the number of days since infection/infectiousness. The third row is the total incubation/infectious period, 
  # drawn from a distribution when it becomes infected/infectious.
  e_nodes<-matrix(nrow=3,ncol=0)
  ev_nodes<-matrix(nrow=3,ncol=0)
  iS_nodes<-matrix(nrow=3,ncol=0)
  iAS_nodes<-matrix(nrow=3,ncol=0)
  i_nodes<-matrix(nrow=3,ncol=0)
  v_nodes<-c()
  s_nodes<-as.vector(V(g))
  nv_nodes<-c()
  r_nodes<-matrix(nrow=2,ncol=0)
  iS_nodes_susc<-c()
  iS_nodes_vacc<-c()
  iAS_nodes_susc<-c() 
  iAS_nodes_vacc<-c()
  

  
  if (prop_rec > 0){ # if no one is recovered to start then initialize introductions only in e_nodes; if not, also add 3 times to iS_nodes/iAS_nodes
  # otherwise do 3 times introductions in iS nodes
    # move initial introductions to e
    init_nodes <- sample(s_nodes,num_introductions*2,replace=FALSE)
    s_nodes <-setdiff(s_nodes,init_nodes)
    e_nodes <- rbind(init_nodes,rep(0,length(init_nodes)),rep(ave_inc_period,length(init_nodes)))
    
    init_nodes_iS <- sample(s_nodes,num_introductions*10*pS,replace=FALSE)
    s_nodes <-setdiff(s_nodes,init_nodes_iS)
    iS_nodes <- rbind(init_nodes_iS,rep(0,length(init_nodes_iS)),rep(ave_inf_period,length(init_nodes_iS)))
    init_nodes_iAS <- sample(s_nodes,num_introductions*10*(1-pS),replace=FALSE)
    s_nodes <-setdiff(s_nodes,init_nodes_iAS)
    iAS_nodes <- rbind(init_nodes_iAS,rep(0,length(init_nodes_iAS)),rep(ave_inf_period,length(init_nodes_iAS)))
  
    # move some to recovered to start but do by risk 
    high_risk <- V(g)[risk==1]$name
    low_risk <- V(g)[risk==0]$name
    
    s_nodes_hr <- intersect(s_nodes,high_risk)
    s_nodes_lr <- intersect(s_nodes,low_risk)
    
    prior_infection_hr <- sample(s_nodes_hr,round(prop_rec*length(s_nodes_hr)*1.5),replace=FALSE)
    prior_infection_lr <- sample(s_nodes_lr,round(prop_rec*length(s_nodes_lr)*0.5),replace=FALSE)
    prior_infection <- c(prior_infection_hr,prior_infection_lr)
    s_nodes <-setdiff(s_nodes,prior_infection)
    r_nodes <- rbind(prior_infection,c(sample(seq(-100,-1,1),length(prior_infection),replace=TRUE)))
  } else{
    prior_infection <- c()
    # move initial introductions to e
    init_nodes <- sample(s_nodes,num_introductions,replace=FALSE)
    s_nodes <-setdiff(s_nodes,init_nodes)
    e_nodes <- rbind(init_nodes,rep(0,length(init_nodes)),rep(ave_inc_period,length(init_nodes)))
  }
  # Initialize results.
  # Results will be, for each newly-infected symptomatic node, the identity of the node, the day it was infected,
  # the community of which it is a member, its trial status at the time of infection 
  ## and whether it's symptomatic or asymptomatic. 
  ## Symptomatic will receive value 1 and Asymptomatic will receive value 0
  # This should be enough information to run a Cox PH.
  # Make a data frame the size of the study pop and fill it in, then trim at the end.
  results<-data.frame("InfectedNode"=rep(NA,studypop_size),"DayInfected"=rep(NA,studypop_size),
                      "Community"=rep(NA,studypop_size),"TrialStatus"=rep(NA,studypop_size),
                      "Risk"=rep(NA,studypop_size), "Symptomatic"=rep(NA,studypop_size))
  popsize<-rep(NA,num_timesteps)
  infperiods<-rep(NA,studypop_size)
  numinfectious<-0
  Epidemic_Curve<-data.frame()
  
  for (t in 1:num_timesteps) {
    #cat(t,length(v_nodes),"\n")
    cat(t,ncol(iS_nodes),ncol(iAS_nodes),"\n")
    
    # I'm recovering first, so I need to ensure that everyone has at least one chance to infect.
    # I do this by initializing an infectious node with 0 days since infection, seeing whether they
    # recover, then advancing them one day along their infectious period.
    
    if (bTrial) {
      
      # Recruit and randomize if during the enrollment period
      if ((t>=trial_startday) && (t<=enroll_endday)) {
        
        num_to_enroll <- enrollment_schedule[t-trial_startday+1]
        
        if (bCluster == 0) {
          # Need to choose from those not already enrolled a number, then randomly assign each one
          # to vaccine or control. If vaccine, change their trialstatus value in the graph to 1, and if
          # they're susceptible move them to v_nodes. If control, just change their trialstatus value
          # in the graph to 0.
          # I sample from those who are non-infectious
          
          ## Stratified randomization
          
          # Need to choose from the clusters not already enrolled
          new_clusters <- sample(non_trial_clusters,num_to_enroll)
          
          # From the chosen clusters, choose a fraction of the non-infectious individual. That fraction is defined in the inputs
          # I will then vaccinate half of each chosen sample
          # For each new cluster, I sample from that cluster a proportion of the whole cluster,
          # but only from the susceptible or exposed individuals. If I'm trying to sample more than are available (because
          # there are lots of infectious/recovered individuals), just sample all of them.
          # These are the people who are recruited - I then assign half and half to vaccine or control
          new_recruits <- lapply(new_clusters,
                                 function(x) sample(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes,r_nodes[1,])),
                                                    min(round(cluster_coverage*length(V(g)[community==x]$name)),
                                                        length(intersect(V(g)[community==x]$name,c(e_nodes[1,],s_nodes,r_nodes[1,]))))))
          
          unlist(new_recruits) -> nr
          intersect(nr,r_nodes[1,]) -> temp
          
          new_vacc <- unlist(lapply(new_recruits,
                                    function(x) sample(x,round(length(x)/2))))
          new_controls <- setdiff(unlist(new_recruits),new_vacc)
          
          
          non_trial_clusters<-setdiff(non_trial_clusters,new_clusters)
          
        } 
        
        V(g)[name %in% new_controls]$trialstatus<-0
        V(g)[name %in% new_vacc]$trialstatus<-1
        V(g)[name %in% new_controls]$enrollmentday<-t
        V(g)[name %in% new_vacc]$enrollmentday<-t
        
        if (VT == "leaky"){
          # Move the vaccinated susceptibles from s_nodes to v_nodes
          vacc_susc <- intersect(s_nodes,new_vacc)
          nonvacc_susc <- intersect(s_nodes,new_controls)
          s_nodes <- setdiff(s_nodes,vacc_susc)
          v_nodes <- c(v_nodes,vacc_susc)
          nv_nodes <- c(nv_nodes,nonvacc_susc)
        } else if (VT == "aon"){
          new_vacc_aon <- sample(new_vacc,length(new_vacc)*direct_VE,replace=FALSE)
          vacc_aon <- intersect(s_nodes,new_vacc_aon)
          s_nodes <- setdiff(s_nodes,vacc_aon)
          r_nodes <- cbind(r_nodes,rbind(vacc_aon,rep(1000,length(vacc_aon)))) # make 1000 so never get infected again
          nonvacc_susc <- intersect(s_nodes,new_controls)
          nv_nodes <- c(nv_nodes,nonvacc_susc)
        }
      }
    }
    
    # Only need to recover if there are any infected or exposed
    if ((ncol(iS_nodes)>0) || (ncol(iAS_nodes)>0) || (ncol(e_nodes)>0) || (ncol(ev_nodes)>0)) {
      list[e_nodes,ev_nodes,iS_nodes,iAS_nodes,i_nodes,iS_nodes_susc,iS_nodes_vacc,iAS_nodes_susc,iAS_nodes_vacc,
           r_nodes,newinfectious_S_total, newinfectious_AS_total,newinfectious]<-
        recover(e_nodes,ev_nodes,iS_nodes,iAS_nodes,iS_nodes_susc,iS_nodes_vacc,iAS_nodes_susc,iAS_nodes_vacc,
                r_nodes,infperiod_shape,infperiod_rate,pS,vpS,t)
      
    } else {
      newinfectious_S_total <- c()
      neinfectious_AS_total <- c()
      newinfectious <- c()
      }
    
    list[s_nodes,nv_nodes,v_nodes,e_nodes,ev_nodes,r_nodes]<-
      spread(g,s_nodes,nv_nodes,v_nodes,e_nodes,ev_nodes,iS_nodes,iAS_nodes,
             iS_nodes_susc,iS_nodes_vacc,iAS_nodes_susc,iAS_nodes_vacc,r_nodes,
             betaS,betaAS,VE,rel_VE_PI,hr,lr,days_reinf,reinf,incperiod_shape,incperiod_rate,
             connected_to_source,extF,out$I[t],t) #out$I1[t]+out$I2[t]+out$I3[t]
    
    ## Get list of infected nodes that were not vaccinated that are in trial (do not need to do for
    ## infected nodes that were vaccinated because only includes trial nodes)
    trial_nodes_t <- V(g)[!is.na(V(g)$trialstatus)]$name
    if(length(trial_nodes_t)!=5000){cat(t,length(trial_nodes_t),"\n")}
    iS_nodes_susc_trial<-setdiff(intersect(trial_nodes_t,iS_nodes[1,]),iS_nodes_vacc)
    iAS_nodes_susc_trial<-setdiff(intersect(trial_nodes_t,iAS_nodes[1,]),iAS_nodes_vacc)
    ## Update epidemic curve for trial population with: time step, # of nonvaccinated susceptible nodes, 
    ## # of vaccinated susceptible nodes, # of symptomatic infected non-vaccinated nodes, 
    ## # of symptomatic infected vaccinated nodes, # of asymptomatic non-vaccinated nodes, # asymptomatic vaccinated nodes
    Epidemic_Curve_List<-c(t,length(nv_nodes),length(v_nodes),
                           length(iS_nodes_susc_trial),length(iS_nodes_vacc),
                           length(iAS_nodes_susc_trial),length(iAS_nodes_vacc))
    Epidemic_Curve<-rbind(Epidemic_Curve,Epidemic_Curve_List)
    colnames(Epidemic_Curve)<-c('t','S','V','iS_s','iS_v','iAS_s','iAS_v')
    
    # Update results
    popsize[t]<-ncol(iS_nodes)+ncol(iAS_nodes)+ncol(e_nodes)+ncol(ev_nodes)+length(s_nodes)+ncol(r_nodes)+length(v_nodes)
    
    # Another error check - the population size should never change                                    
    if (ncol(iS_nodes)+ncol(iAS_nodes)+ncol(e_nodes)+ncol(ev_nodes)+length(s_nodes)+ncol(r_nodes)+length(v_nodes) != studypop_size) {
      sink("Error.txt",append=TRUE)
      print(iS_nodes[1,])
      print(iAS_nodes[1,])
      print(e_nodes[1,])
      print(ev_nodes[1,])
      print(s_nodes)
      print(v_nodes)
      print(r_nodes)
      stop("Unequal lengths")
      
    }
    
    numnewinfectious<-length(newinfectious)
    if (numnewinfectious>0) {
      
      # Update results
      results$InfectedNode[(numinfectious+1):(numinfectious+numnewinfectious)]<-newinfectious
      results$DayInfected[(numinfectious+1):(numinfectious+numnewinfectious)]<-rep(t,numnewinfectious)
      results$Community[(numinfectious+1):(numinfectious+numnewinfectious)]<-V(g)[name %in% newinfectious]$community
      results$TrialStatus[(numinfectious+1):(numinfectious+numnewinfectious)]<-V(g)[name %in% newinfectious]$trialstatus
      results$Risk[(numinfectious+1):(numinfectious+numnewinfectious)]<-V(g)[name %in% newinfectious]$risk
      V(g)[name %in% iS_nodes[1,]]$symptomatic<-1
      V(g)[name %in% iAS_nodes[1,]]$symptomatic<-0
      results$Symptomatic[(numinfectious+1):(numinfectious+numnewinfectious)]<-V(g)[name %in% newinfectious]$symptomatic
      infperiods[(numinfectious+1):(numinfectious+numnewinfectious)]<-i_nodes[3,i_nodes[1,] %in% newinfectious]
      
      numinfectious <- numinfectious+numnewinfectious
      
    }
  }
  
  trial_nodes <- V(g)[!is.na(V(g)$trialstatus)]$name
  
  # Tidy up results
  results<-results[1:numinfectious,]
  infperiods<-infperiods[1:numinfectious]
  
  if (nrow(results)>1){
    
    results$DayEnrolled <- V(g)[results$InfectedNode]$enrollmentday
    results$DayInfected <- ifelse(is.na(results$DayEnrolled),results$DayInfected,results$DayInfected - results$DayEnrolled)
    
    # Get a list of nodes that were enrolled in the trial but never infected
    noninf<-setdiff(trial_nodes,results$InfectedNode)
    # Get list of nodes that became infectious while they were in the trial
    results_analysis<-results[!is.na(results$TrialStatus),]
    # Get a list of nodes who were infected after their follow-up time was over
    # (i.e. those enrolled at the beginning but infected right at the end)
    censored <- results_analysis[results_analysis$DayInfected>trial_length,]
    results_analysis<-results_analysis[results_analysis$DayInfected+results_analysis$DayEnrolled<=trial_length+trial_startday,]
    # Assign them eventstatus=1 for the Cox analysis
    results_analysis$eventstatus<-rep(1,nrow(results_analysis))
    # Make data frame for those who were never infected (i.e. censored by end of study)
    # They will have follow-up time defined to be the end of the follow-up (=trial length)
    ## Symptomatic status will be NA
    noninfdf<-data.frame(InfectedNode=noninf,DayInfected=rep(trial_length,length(noninf)),
                         Community=V(g)$community[noninf],
                         TrialStatus=V(g)$trialstatus[noninf],
                         Risk=V(g)$risk[noninf],
                         eventstatus=rep(0,length(noninf)),
                         DayEnrolled=V(g)[noninf]$enrollmentday,
                         Symptomatic=rep(NA,length(noninf)))
    if (nrow(censored)>0) {
      censored$DayInfected<-trial_length
      censored$eventstatus<-0
      censored$Symptomatic<-NA
    }
    
    results_analysis<-rbind(results_analysis,noninfdf,censored)
    
    # Finally, exclude any cases who were infected during the first n days of follow-up
    # This is to get rid of those who were already latently infected when enrolled
    results_analysis<-results_analysis[results_analysis$DayInfected>ave_inc_period,]
    
    results_analysis %>%
      mutate(PI = ifelse(InfectedNode %in% prior_infection,1,0)) -> results_analysis # prior infection
  
  #return(nrow(results))
  list(results_analysis,results,popsize,infperiods,trial_nodes,g,
       s_nodes,nv_nodes,v_nodes,e_nodes,ev_nodes,i_nodes,
       iS_nodes,iAS_nodes,iS_nodes_susc,iS_nodes_vacc,iAS_nodes_susc,iAS_nodes_vacc,r_nodes,
       Epidemic_Curve,num_timesteps)
  
  }
}
  
# j is the number of simulatiton on -- can make a for loop through the j's or run simultaneously
j<-j
track <- 0

#master
VT <- "leaky"
results_analysis_master <- data.frame()

for (rate_within in c(0.0018174, 0.002182)){ #0.001455,
  for (direct_VE in c(0.7,0.9)){
    for (vpS in c(0.25)){
      for(prop_rec in c(0)){ # 0,0
        for(prop_hr in c(0.5)){ #0,
          for(reinf in c(0.05)){
            for(rel_VE_PI in c(1)){

              pS <- 0.5

              track <- track + 1
              R0 <- round((1 - (infperiod_rate/(infperiod_rate+betaS))^infperiod_shape) *
                            (((ave_community_size-1)*(1-rate_within)*rate_within +
                                (num_communities-1)*ave_community_size*(1-rate_between)*rate_between +
                                ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)^2)/
                               ((ave_community_size-1)*rate_within + (num_communities-1)*ave_community_size*rate_between)-1),2)
              R0


              # Make the network
              g<-make_network(ave_community_size, community_size_range, num_communities, rate_within, rate_between, prop_hr)
              # Run the iRCT (in the inputs below, the 1 means that we are running the trial, and the 0 means that
              # it's an iRCT)
              list[results_analysis,results,popsize,infperiods,trial_nodes,g,
                   s_nodes,nv_nodes,v_nodes,e_nodes,ev_nodes,i_nodes,
                   iS_nodes,iAS_nodes,iS_nodes_susc,iS_nodes_vacc,iAS_nodes_susc,iAS_nodes_vacc,r_nodes,
                   Epidemic_Curve,num_timesteps]<-
                network_epidemic(g,betaS,betaAS,num_introductions,direct_VE,rel_VE_PI,hr,lr,days_reinf,reinf,
                                 incperiod_shape,incperiod_rate,infperiod_shape,infperiod_rate,pS,vpS,prop_rec,1,0,VT,
                                 trial_startday,trial_length,num_ind_enrolled_per_day,enrollment_period,cluster_coverage)

              results_analysis %>%
                cbind(R0=R0,direct_VE=direct_VE,vpS=vpS,prop_hr=prop_hr,prop_rec=prop_rec,reinf=reinf,rel_VE_PI=rel_VE_PI,j=j) %>%
                bind_rows(results_analysis_master) -> results_analysis_master

              write.csv(results_analysis_master, paste0(j,"_results_analysis_master.csv"))
            }
          }
        }
      }
    }
  }
}



  