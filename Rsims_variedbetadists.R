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
Q12<-matrix(nrow=nrow(gen1), ncol=1)
het<-matrix(nrow=nrow(gen1), ncol=ncol(gen1))
for(i in 1:nrow(gen1)){
  for(j in 1:ncol(gen1)){
  het[i,j]<-ifelse(gen1[i, j]==1, 1, 0)
  }
}
Q12<-rowSums(het)/ncol(gen1)

#gen<-rep(1, 2000)
#individual<-(1:2000)
#phenotype<-rep(-9, 2000)
#sex<-sample(1:2, 2000, replace=TRUE)

###change gen1 into plink happy format
ifelse(gen1==0, "A A", ifelse(gen1==1, "A C", "C C"))->gen1_AC
cbind(rep(1, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), gen1_AC)->gen1_ped

chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

gen1_map<-cbind.data.frame(chromosomes, locus_name, rep(0,2000), locus)

write.table(gen1_ped, file="gen1.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen1_map, file="gen1.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen1 --make-bed --out gen1")
system("~/admixture gen1.bed 2")
read.table("gen1.2.Q")->gen1_q

###gen 2 - now that hybridization has started, I'm still going to have half of the individuals be 'pure parental types' and half be the result of this population interbreeding

gametes.gen1<-matrix(nrow=nrow(gen1), ncol=ncol(gen1))
for(i in 1:nrow(gen1)){
  for(j in 1:ncol(gen1)){
    gametes.gen1[i, j]<-ifelse(gen1[i, j]==0, 0, 
                                     ifelse(gen1[i, j]==1, sample(0:1, 1, replace=FALSE), 1))
  }
}

x.parents.g1<-sample(1:2000, 2000, replace=FALSE)

gametes.gen1[x.parents.g1[1:1000],]+gametes.gen1[x.parents.g1[1001:2000],]->F2_genomes

sim.x.gen2 <- matrix(rbinom(nloci*(nind-nrow(F2_genomes))/2, 2, prob=sim.p), nrow=(nind-nrow(F2_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns
sim.y.gen2 <- matrix(rbinom(nloci*(nind-nrow(F2_genomes))/2, 2, prob=sim.q), nrow=(nind-nrow(F2_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns

gen2<-rbind(sim.x.gen2, sim.y.gen2, F2_genomes)
### let's do a triangle plot here?
### need to get heterozygosity and Q12 for each individual
### this takes annoyingly long to run. Don't panic
Q12_gen2<-matrix(nrow=nrow(gen2), ncol=1)
het_gen2<-matrix(nrow=nrow(gen2), ncol=ncol(gen2))
for(i in 1:nrow(gen2)){
  for(j in 1:ncol(gen2)){
    het_gen2[i,j]<-ifelse(gen2[i, j]==1, 1, 0)
  }
}

Q12_gen2<-rowSums(het_gen2)/ncol(gen2)

###change gen2 into plink happy format
ifelse(gen2==0, "A A", ifelse(gen2==1, "A C", "C C"))->gen2_AC
cbind(rep(2, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), gen2_AC)->gen2_ped

chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

gen2_map<-cbind.data.frame(chromosomes, locus_name, rep(0,2000), locus)

write.table(gen2_ped, file="gen2.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen2_map, file="gen2.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen2 --make-bed --out gen2")

### Run admixture to get q
system("~/admixture gen2.bed 2")
read.table("gen2.2.Q")->gen2_q


### gen 3 ###


gametes.gen2<-matrix(nrow=nrow(gen2), ncol=ncol(gen2))
for(i in 1:nrow(gen2)){
  for(j in 1:ncol(gen2)){
    gametes.gen2[i, j]<-ifelse(gen2[i, j]==0, 0, 
                               ifelse(gen2[i, j]==1, sample(0:1, 1, replace=FALSE), 1))
  }
}

x.parents.g2<-sample(1:2000, 2000, replace=FALSE)

gametes.gen2[x.parents.g2[1:1000],]+gametes.gen2[x.parents.g2[1001:2000],]->F3_genomes

sim.x.gen3 <- matrix(rbinom(nloci*(nind-nrow(F3_genomes))/2, 2, prob=sim.p), nrow=(nind-nrow(F3_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns
sim.y.gen3 <- matrix(rbinom(nloci*(nind-nrow(F3_genomes))/2, 2, prob=sim.q), nrow=(nind-nrow(F3_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns

gen3<-rbind(sim.x.gen3, sim.y.gen3, F3_genomes)
### let's do a triangle plot here?
### need to get heterozygosity and Q12 for each individual
Q12_gen3<-matrix(nrow=nrow(gen3), ncol=1)
het_gen3<-matrix(nrow=nrow(gen3), ncol=ncol(gen3))
for(i in 1:nrow(gen3)){
  for(j in 1:ncol(gen3)){
    het_gen3[i,j]<-ifelse(gen3[i, j]==1, 1, 0)
  }
}
Q12_gen3<-rowSums(het_gen3)/ncol(gen3)

###change gen3 into plink happy format
ifelse(gen3==0, "A A", ifelse(gen3==1, "A C", "C C"))->gen3_AC
cbind(rep(3, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), gen3_AC)->gen3_ped

chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

gen3_map<-cbind.data.frame(chromosomes, locus_name, rep(0,2000), locus)

write.table(gen3_ped, file="gen3.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen3_map, file="gen3.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen3 --make-bed --out gen3")

### Run admixture to get q
system("~/admixture gen3.bed 2")
read.table("gen3.2.Q")->gen3_q

### gen 4
### parents gen(i)-1
gametes.gen3<-matrix(nrow=nrow(gen3), ncol=ncol(gen3))
for(i in 1:nrow(gen3)){
  for(j in 1:ncol(gen3)){
    gametes.gen3[i, j]<-ifelse(gen3[i, j]==0, 0, 
                               ifelse(gen3[i, j]==1, sample(0:1, 1, replace=FALSE), 1))
  }
}

x.parents.g3<-sample(1:2000, 2000, replace=FALSE)

gametes.gen3[x.parents.g3[1:1000],]+gametes.gen3[x.parents.g3[1001:2000],]->F4_genomes

sim.x.gen4 <- matrix(rbinom(nloci*(nind-nrow(F4_genomes))/2, 2, prob=sim.p), nrow=(nind-nrow(F4_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns
sim.y.gen4 <- matrix(rbinom(nloci*(nind-nrow(F4_genomes))/2, 2, prob=sim.q), nrow=(nind-nrow(F4_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns

###offspring gen(i)
gen4<-rbind(sim.x.gen4, sim.y.gen4, F4_genomes)
### let's do a triangle plot here?
### need to get heterozygosity and Q12 for each individual
Q12_gen4<-matrix(nrow=nrow(gen4), ncol=1)
het_gen4<-matrix(nrow=nrow(gen4), ncol=ncol(gen4))
for(i in 1:nrow(gen4)){
  for(j in 1:ncol(gen4)){
    het_gen4[i,j]<-ifelse(gen4[i, j]==1, 1, 0)
  }
}
Q12_gen4<-rowSums(het_gen4)/ncol(gen4)

###change gen4 into plink happy format
ifelse(gen4==0, "A A", ifelse(gen4==1, "A C", "C C"))->gen4_AC
cbind(rep(4, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), gen4_AC)->gen4_ped

chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

gen4_map<-cbind.data.frame(chromosomes, locus_name, rep(0,2000), locus)

write.table(gen4_ped, file="gen4.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen4_map, file="gen4.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen4 --make-bed --out gen4")

### Run admixture to get q
system("~/admixture gen4.bed 2", wait=TRUE)
read.table("gen4.2.Q")->gen4_q

### gen 5

###parents gen(i)-1
gametes.gen4<-matrix(nrow=nrow(gen4), ncol=ncol(gen2))
for(i in 1:nrow(gen4)){
  for(j in 1:ncol(gen4)){
    gametes.gen4[i, j]<-ifelse(gen4[i, j]==0, 0, 
                               ifelse(gen4[i, j]==1, sample(0:1, 1, replace=FALSE), 1))
  }
}

x.parents.g4<-sample(1:2000, 2000, replace=FALSE)

gametes.gen4[x.parents.g4[1:1000],]+gametes.gen4[x.parents.g4[1001:2000],]->F5_genomes

###offspring gen(i)
sim.x.gen5 <- matrix(rbinom(nloci*(nind-nrow(F4_genomes))/2, 2, prob=sim.p), nrow=(nind-nrow(F4_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns
sim.y.gen5 <- matrix(rbinom(nloci*(nind-nrow(F4_genomes))/2, 2, prob=sim.q), nrow=(nind-nrow(F4_genomes))/2, ncol=nloci) ### this gives us individuals in rows and snps in columns

gen5<-rbind(sim.x.gen5, sim.y.gen5, F5_genomes)
### let's do a triangle plot here?
### need to get heterozygosity and Q12 for each individual
Q12_gen5<-matrix(nrow=nrow(gen5), ncol=1)
het_gen5<-matrix(nrow=nrow(gen5), ncol=ncol(gen5))
for(i in 1:nrow(gen5)){
  for(j in 1:ncol(gen5)){
    het_gen5[i,j]<-ifelse(gen5[i, j]==1, 1, 0)
  }
}
Q12_gen5<-rowSums(het_gen5)/ncol(gen5)
###change gen5 into plink happy format
ifelse(gen5==0, "A A", ifelse(gen5==1, "A C", "C C"))->gen5_AC
cbind(rep(5, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), gen5_AC)->gen5_ped

chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

gen5_map<-cbind.data.frame(chromosomes, locus_name, rep(0,2000), locus)

write.table(gen5_ped, file="gen5.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen5_map, file="gen5.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen5 --make-bed --out gen5")

### Run admixture to get q
system("~/admixture gen5.bed 2", wait=TRUE)
read.table("gen5.2.Q")->gen5_q