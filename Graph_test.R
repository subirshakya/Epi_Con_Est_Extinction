## function to compute asymmetry of root node
 library(ips)
 library(ape)
 library(TreeSim) 

# root asymmetry (ra) function
# 'tr' is a object of class 'phylo' from the ape library 
ra<- function(tr)
{
	N<- Ntip(tr)
	rn<- N+1
	ee<- tr$edge
	ees<- ee[ee[,1]==rn, ]	# extract edges that stem from root node
	desc.node<- ees[,2]
	
	if(any(desc.node < N))	{ RS<- 1; LS<- N-1 	# if either descendant of root is a terminal tip
	}	else	{
		RN<- descendants(tr, node=desc.node[1], type="t")	# get descdendants
		LN<- descendants(tr, node=desc.node[2], type="t")
		RS<- length(RN)
		LS<- length(LN)
	}	
	
	ext<- abs(RS - LS) / (2* max(c(RS, LS)))  # compute asymmetry metric
	res<- c(ext, max(c(RS, LS)), min(c(RS, LS)))
	names(res)<- c("asym", "right", "left")

	# side with more tips is considered the "right" side of tree
	return(res)
	
}

nrep<- 1000
layout(1:3)

# no extinction
tr<- sim.bd.taxa(n=32, numbsim=nrep, lambda=1, mu=0)
rav<- t(sapply(tr, ra))
hist(rav[,1], 25, col="pink", main="No extinction")

# low extinction
tr<- sim.bd.taxa(n=32, numbsim=nrep, lambda=1, mu=0.2)
rav<- t(sapply(tr, ra))
hist(rav[,1], 25, col="pink", main="Low extinction, mu=0.2")

# high extinction
tr<- sim.bd.taxa(n=32, numbsim=nrep, lambda=1, mu=0.8)
rav<- t(sapply(tr, ra))
hist(rav[,1], 25, col="pink", main="High extinction, mu=0.8")

