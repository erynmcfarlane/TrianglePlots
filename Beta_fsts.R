sim.p<-rbeta(nloci, 17, 15)### gives fst 0f 0.558
hist(sim.p, breaks=seq(0,1,0.1))
sim.q<-rbeta(nloci, 15, 17)
hist(sim.q, breaks=seq(0,1,0.1))

sim.p<-rbeta(nloci, 20, 18) ### fst of 0.532
hist(sim.p, breaks=seq(0,1,0.1))
sim.q<-rbeta(nloci, 18, 20)
hist(sim.q, breaks=seq(0,1,0.1))

sim.p<-rbeta(nloci, 25, 24) ### fst = 0.455

sim.p<-rbeta(nloci, 30, 28) ### fst = 0.334

sim.p<-rbeta(nloci, 40, 38) ### fst = 0.320

sim.p<-rbeta(nloci, 80, 70) ### fst = 0.197

sim.p<-rbeta(nloci, 10, 8) ### fst = 0.745


#### MOAR Markers means much lower FSt