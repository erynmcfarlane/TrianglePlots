#### simulate two different populations without fixed markers

### This is from Alex Buerkle's Population Genetics lab notes
###simulate some data
set.seed(42)
nloci<-10000  ## we will have data from nloci
nind<-2000  ## we will have data from nind that were sampled from the true population
sim.theta<-5  # high parameter means high polymorphism, must be positive and > zero, can be thought of as analygous to nucleotide polymorphism

## simulate allele frequencies at nloci number of loci, by random draws from a beta
### I've made these beta distributions really uneven. But I could make them a lot more even and see what happens.
sim.p<-rbeta(nloci, 1, 5)
hist(sim.p, breaks=seq(0,1,0.1))
sim.q<-rbeta(nloci, 5, 1)
hist(sim.q, breaks=seq(0,1,0.1))
# simulate genotypes in the sample
sim.x <- matrix(rbinom(nloci*nind, 2, prob=sim.p), nrow=nind, ncol=nloci) ### this gives us individuals in rows and snps in columns
str(sim.x)

sim.y <- matrix(rbinom(nloci*nind, 2, prob=sim.q), nrow=nind, ncol=nloci) ### this gives us individuals in rows and snps in columns
str(sim.y)

#### First generation

###make parental generation gametes
gametes.x.parental<-matrix(nrow=nrow(sim.x), ncol=ncol(sim.x))
for(i in 1:nrow(sim.x)){
  for(j in 1:ncol(sim.x)){
    gametes.x.parental[i, j]<-ifelse(sim.x[i, j]==0, 0, 
           ifelse(sim.x[i, j]==1, sample(0:1, 1, replace=FALSE), 1))
  }
}

gametes.y.parental<-matrix(nrow=nrow(sim.y), ncol=ncol(sim.y))
for(i in 1:nrow(sim.y)){
  for(j in 1:ncol(sim.y)){
    gametes.y.parental[i, j]<-ifelse(sim.y[i, j]==0, 0, 
                                     ifelse(sim.y[i, j]==1, sample(0:1, 1, replace=FALSE), 1))
  }
}


###sample parents
### these are the parents who become the parents of F1 hybrids.
### to get individuals in the population who are not hybrids, I can pull from the distributions again

x.parents<-sample(1:2000, 1000, replace=FALSE)
y.parents<-sample(1:2000, 1000, replace=FALSE)

gametes.x.parental[x.parents,]+gametes.y.parental[y.parents,]->F1_genomes

sim.x.gen1 <- matrix(rbinom(nloci*(nind-nrow(F1_genomes))/2, 2, prob=sim.p), nrow=(nind-nrow(F1_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns
sim.y.gen1 <- matrix(rbinom(nloci*(nind-nrow(F1_genomes))/2, 2, prob=sim.q), nrow=(nind-nrow(F1_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns

gen1<-rbind(sim.x.gen1, sim.y.gen1, F1_genomes)
### let's do a triangle plot here?
### need to get heterozygosity and Q12 for each individual



###gen 2 - now that hybridization has started, I'm still going to have half of the individuals be 'pure parental types' and half be the result of this population interbreeding

gametes.gen1<-matrix(nrow=nrow(gen1), ncol=ncol(gen1))
for(i in 1:nrow(gen1)){
  for(j in 1:ncol(gen1)){
    gametes.gen1[i, j]<-ifelse(gen1[i, j]==0, 0, 
                                     ifelse(gen1[i, j]==1, sample(0:1, 1, replace=FALSE), 1))
  }
}
