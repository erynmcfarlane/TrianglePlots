#### simulate two different populations without fixed markers


### NEED TO PULL Q12 MODEL FROM ZACH's PAPER
### This is from Alex Buerkle's Population Genetics lab notes
###simulate some data
set.seed(42)
nloci<-10000  ## we will have data from nloci
nind<-2000  ## we will have data from nind that were sampled from the true population
#sim.theta<-5  # high parameter means high polymorphism, must be positive and > zero, can be thought of as analygous to nucleotide polymorphism

#### need to change here when I change nloci
chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

source('plot_triangles.R') ###Josh's triangle plot function

## simulate allele frequencies at nloci number of loci, by random draws from a beta
### I've made these beta distributions really uneven. But I could make them a lot more even and see what happens.
fst_high<-c(3, 1) ### these thetas are way too high
fst_medium<-c(6,4)
fst_low<-c(10,8)

######################## FST_high for 5 generations ######
sim.p<-rbeta(nloci, fst_high[1], fst_high[2])
hist(sim.p, breaks=seq(0,1,0.1))
sim.q<-rbeta(nloci, fst_high[2], fst_high[1])
hist(sim.q, breaks=seq(0,1,0.1))
# simulate genotypes in the sample
sim.x <- matrix(rbinom(nloci*nind, 2, prob=sim.p), nrow=nind, ncol=nloci) ### this gives us individuals in rows and snps in columns
str(sim.x)

sim.y <- matrix(rbinom(nloci*nind, 2, prob=sim.q), nrow=nind, ncol=nloci) ### this gives us individuals in rows and snps in columns
str(sim.y)

### add FST between the two here:

parents<-rbind(sim.x, sim.y)
### make parental populations into plink .bed file
ifelse(parents==0, "A A", ifelse(parents==1, "A C", "C C"))->parents_AC
cbind(rep(1, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), parents_AC)->parents_ped

parents_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(parents_ped, file="parents.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(parents_map, file="parents.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file parents --make-bed --out parents", wait=TRUE)
system("~/admixture parents.bed 2", wait=TRUE)


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
Q12_gen1<-rowSums(het)/ncol(gen1)

#gen<-rep(1, 2000)
#individual<-(1:2000)
#phenotype<-rep(-9, 2000)
#sex<-sample(1:2, 2000, replace=TRUE)

###change gen1 into plink happy format
ifelse(gen1==0, "A A", ifelse(gen1==1, "A C", "C C"))->gen1_AC
cbind(rep(1, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), gen1_AC)->gen1_ped

gen1_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen1_ped, file="gen1.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen1_map, file="gen1.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen1 --make-bed --out gen1", wait=TRUE)
system("~/admixture gen1.bed 2", wait=TRUE)
read.table("gen1.2.Q")->gen1_q
gen1_qQ<-cbind(gen1_q[,1], Q12_gen1)

triangleplot("gen1_fst_high.pdf", gen1_qQ, 1)

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

gen2_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen2_ped, file="gen2.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen2_map, file="gen2.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen2 --make-bed --out gen2", wait=TRUE)

### Run admixture to get q
system("~/admixture gen2.bed 2", wait=TRUE)
read.table("gen2.2.Q")->gen2_q
gen2_qQ<-cbind(gen2_q[,1], Q12_gen2)
triangleplot("gen2_fst_high.pdf", gen2_qQ, 2)

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

gen3_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen3_ped, file="gen3.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen3_map, file="gen3.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen3 --make-bed --out gen3", wait=TRUE)

### Run admixture to get q
system("~/admixture gen3.bed 2", wait=TRUE)
read.table("gen3.2.Q")->gen3_q
gen3_qQ<-cbind(gen3_q[,1], Q12_gen3)

triangleplot("gen3_fst_high.pdf", gen3_qQ, 3)

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

gen4_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen4_ped, file="gen4.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen4_map, file="gen4.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen4 --make-bed --out gen4", wait=TRUE)

### Run admixture to get q
system("~/admixture gen4.bed 2", wait=TRUE)
read.table("gen4.2.Q")->gen4_q
gen4_qQ<-cbind(gen4_q[,1], Q12_gen4)

triangleplot("gen4_fst_high.pdf", gen4_qQ, 4)

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

#### need to change here when I change nloci
chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

gen5_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen5_ped, file="gen5.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen5_map, file="gen5.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen5 --make-bed --out gen5", wait=TRUE)

### Run admixture to get q
system("~/admixture gen5.bed 2", wait=TRUE)
read.table("gen5.2.Q")->gen5_q
gen5_qQ<-cbind(gen5_q[,1], Q12_gen5)

triangleplot("gen5_fst_high.pdf", gen5_qQ, 5)

### what's going to be really interesting is adding a function to get the diagnostic markers out, and then doing the triangle plots for those too

######################## FST_medium for 5 generations ######
sim.p<-rbeta(nloci, fst_medium[1], fst_medium[2])
hist(sim.p, breaks=seq(0,1,0.1))
sim.q<-rbeta(nloci, fst_medium[2], fst_medium[1])
hist(sim.q, breaks=seq(0,1,0.1))
# simulate genotypes in the sample
sim.x <- matrix(rbinom(nloci*nind, 2, prob=sim.p), nrow=nind, ncol=nloci) ### this gives us individuals in rows and snps in columns
str(sim.x)

sim.y <- matrix(rbinom(nloci*nind, 2, prob=sim.q), nrow=nind, ncol=nloci) ### this gives us individuals in rows and snps in columns
str(sim.y)

### add FST between the two here:

parents<-rbind(sim.x, sim.y)
### make parental populations into plink .bed file
ifelse(parents==0, "A A", ifelse(parents==1, "A C", "C C"))->parents_AC
cbind(rep(1, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), parents_AC)->parents_ped

parents_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(parents_ped, file="parents.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(parents_map, file="parents.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file parents --make-bed --out parents", wait=TRUE)
system("~/admixture parents.bed 2", wait=TRUE)


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
Q12_gen1<-rowSums(het)/ncol(gen1)

#gen<-rep(1, 2000)
#individual<-(1:2000)
#phenotype<-rep(-9, 2000)
#sex<-sample(1:2, 2000, replace=TRUE)

###change gen1 into plink happy format
ifelse(gen1==0, "A A", ifelse(gen1==1, "A C", "C C"))->gen1_AC
cbind(rep(1, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), gen1_AC)->gen1_ped

gen1_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen1_ped, file="gen1.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen1_map, file="gen1.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen1 --make-bed --out gen1", wait=TRUE)
system("~/admixture gen1.bed 2", wait=TRUE)
read.table("gen1.2.Q")->gen1_q
gen1_qQ<-cbind(gen1_q[,1], Q12_gen1)

triangleplot("gen1_fst_medium.pdf", gen1_qQ, 1)

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

gen2_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen2_ped, file="gen2.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen2_map, file="gen2.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen2 --make-bed --out gen2", wait=TRUE)

### Run admixture to get q
system("~/admixture gen2.bed 2", wait=TRUE)
read.table("gen2.2.Q")->gen2_q
gen2_qQ<-cbind(gen2_q[,1], Q12_gen2)
triangleplot("gen2_fst_medium.pdf", gen2_qQ, 2)

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

gen3_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen3_ped, file="gen3.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen3_map, file="gen3.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen3 --make-bed --out gen3", wait=TRUE)

### Run admixture to get q
system("~/admixture gen3.bed 2", wait=TRUE)
read.table("gen3.2.Q")->gen3_q
gen3_qQ<-cbind(gen3_q[,1], Q12_gen3)

triangleplot("gen3_fst_medium.pdf", gen3_qQ, 3)

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

gen4_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen4_ped, file="gen4.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen4_map, file="gen4.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen4 --make-bed --out gen4", wait=TRUE)

### Run admixture to get q
system("~/admixture gen4.bed 2", wait=TRUE)
read.table("gen4.2.Q")->gen4_q
gen4_qQ<-cbind(gen4_q[,1], Q12_gen4)

triangleplot("gen4_fst_medium.pdf", gen4_qQ, 4)

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

#### need to change here when I change nloci
chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

gen5_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen5_ped, file="gen5.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen5_map, file="gen5.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen5 --make-bed --out gen5", wait=TRUE)

### Run admixture to get q
system("~/admixture gen5.bed 2", wait=TRUE)
read.table("gen5.2.Q")->gen5_q
gen5_qQ<-cbind(gen5_q[,1], Q12_gen5)

triangleplot("gen5_fst_medium.pdf", gen5_qQ, 5)


######################## fst_low for 5 generations ######
sim.p<-rbeta(nloci, fst_low[1], fst_low[2])
hist(sim.p, breaks=seq(0,1,0.1))
sim.q<-rbeta(nloci, fst_low[2], fst_low[1])
hist(sim.q, breaks=seq(0,1,0.1))
# simulate genotypes in the sample
sim.x <- matrix(rbinom(nloci*nind, 2, prob=sim.p), nrow=nind, ncol=nloci) ### this gives us individuals in rows and snps in columns
str(sim.x)

sim.y <- matrix(rbinom(nloci*nind, 2, prob=sim.q), nrow=nind, ncol=nloci) ### this gives us individuals in rows and snps in columns
str(sim.y)

### add FST between the two here:

parents<-rbind(sim.x, sim.y)
### make parental populations into plink .bed file
ifelse(parents==0, "A A", ifelse(parents==1, "A C", "C C"))->parents_AC
cbind(rep(1, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), parents_AC)->parents_ped

parents_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(parents_ped, file="parents.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(parents_map, file="parents.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file parents --make-bed --out parents", wait=TRUE)
system("~/admixture parents.bed 2", wait=TRUE)


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
Q12_gen1<-rowSums(het)/ncol(gen1)

#gen<-rep(1, 2000)
#individual<-(1:2000)
#phenotype<-rep(-9, 2000)
#sex<-sample(1:2, 2000, replace=TRUE)

###change gen1 into plink happy format
ifelse(gen1==0, "A A", ifelse(gen1==1, "A C", "C C"))->gen1_AC
cbind(rep(1, 2000), (1:2000), sample(1:2, 2000, replace=TRUE),rep(0, 2000), rep(0, 2000), rep(-9, 2000), gen1_AC)->gen1_ped

gen1_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen1_ped, file="gen1.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen1_map, file="gen1.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen1 --make-bed --out gen1", wait=TRUE)
system("~/admixture gen1.bed 2", wait=TRUE)
read.table("gen1.2.Q")->gen1_q
gen1_qQ<-cbind(gen1_q[,1], Q12_gen1)

triangleplot("gen1_fst_low.pdf", gen1_qQ, 1)

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

gen2_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen2_ped, file="gen2.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen2_map, file="gen2.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen2 --make-bed --out gen2", wait=TRUE)

### Run admixture to get q
system("~/admixture gen2.bed 2", wait=TRUE)
read.table("gen2.2.Q")->gen2_q
gen2_qQ<-cbind(gen2_q[,1], Q12_gen2)
triangleplot("gen2_fst_low.pdf", gen2_qQ, 2)

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

gen3_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen3_ped, file="gen3.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen3_map, file="gen3.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen3 --make-bed --out gen3", wait=TRUE)

### Run admixture to get q
system("~/admixture gen3.bed 2", wait=TRUE)
read.table("gen3.2.Q")->gen3_q
gen3_qQ<-cbind(gen3_q[,1], Q12_gen3)

triangleplot("gen3_fst_low.pdf", gen3_qQ, 3)

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

gen4_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen4_ped, file="gen4.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen4_map, file="gen4.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen4 --make-bed --out gen4", wait=TRUE)

### Run admixture to get q
system("~/admixture gen4.bed 2", wait=TRUE)
read.table("gen4.2.Q")->gen4_q
gen4_qQ<-cbind(gen4_q[,1], Q12_gen4)

triangleplot("gen4_fst_low.pdf", gen4_qQ, 4)

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

#### need to change here when I change nloci
chromosomes<-rep(1:20, each=500)
locus<-rep(1:500, 20)
locus_name<-paste0(chromosomes,":", locus)

gen5_map<-cbind.data.frame(chromosomes, locus_name, rep(0,nloci), locus)

write.table(gen5_ped, file="gen5.ped", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(gen5_map, file="gen5.map", quote=FALSE, col.names=FALSE, row.names=FALSE)

system("/Users/Eryn/PATH/plink --file gen5 --make-bed --out gen5", wait=TRUE)

### Run admixture to get q
system("~/admixture gen5.bed 2", wait=TRUE)
read.table("gen5.2.Q")->gen5_q
gen5_qQ<-cbind(gen5_q[,1], Q12_gen5)

triangleplot("gen5_fst_low.pdf", gen5_qQ, 5)

save.image("Triangles161123.RData")
