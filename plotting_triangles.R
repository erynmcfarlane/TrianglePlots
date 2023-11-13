#### once all of the simulations are done, move into R to do all this ####
library(tidyverse)
library(MASS)
library(data.table)
library(MetBrewer)
library(introgress)

datafiles<-list.files("/gscratch/emcfarl2/TrianglePlots/",  pattern="*u*.main", recursive=TRUE, include.dirs=TRUE)
###this would work if I didn't have some other random stuff in there
basenames<-basename(datafiles)
m<-str_extract(basenames, "(\\d+\\.*\\d*)")
c<-str_match(basenames,"c(\\d+\\.*\\d*)")[,2]
c[is.na(c)]<-0 #### I think this is right, as c is the measure of selection?
#mech<-str_extract(basenames, "^([^_]+_){1}([^_])") 


### I can pull q and Q12 from the data, but how do I get heterozygosity?

alldata<-list()
for(i in 1:length(datafiles)){
  alldata[[i]]<-fread(datafiles[i], sep=",", header=T)
  alldata[[i]]$m<-as.numeric(rep(m[i], nrow(alldata[[i]]))) # just giving all individuals in the sim the same m and c
  alldata[[i]]$c<-as.numeric(rep(c[i], nrow(alldata[[i]])))
  #alldata[[i]]$mech<-as.factor(rep(mech[i], nrow(alldata[[i]])))
}


alldata_df<-do.call(rbind.data.frame, alldata)

rm(alldata)

genotypes<-ifelse(alldata_df[,9:518] == 1, 1, 0) 
### it's still not clear to me how this is different from Q12?
alldata_df$average_het<-rowMeans(genotypes)

alldata_df[which(alldata_df$gen==10 & alldata_df$m==0.01),]->gen_10_0.01
alldata_df[which(alldata_df$gen==10 & alldata_df$m==0.2),]->gen_10_0.2
alldata_df[which(alldata_df$gen==100 & alldata_df$m==0.01 ),]->gen_100_0.01
alldata_df[which(alldata_df$gen==100 & alldata_df$m==0.2 ),]->gen_100_0.2

###plot triangles

library(MetBrewer)
triangle_cols <- met.brewer("OKeeffe1", 20)


#quartz(height=12, width=15)
#par(mar=c(5,5,1,1), mfrow=c(4,5))
#for (i in 1:20)
#{
  #rep_sub <- subset(gen10_deme6, gen10_deme6[,1]==i)
pdf(file="triangles_for_triangles_10_0.01.pdf", width=10, height=10)
  plot(0, type="n", xlab="q", ylab="Q", xlim=c(0, 1), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25)
  segments(0, 0, 0.5, 1, lwd=3)
  segments(1, 0, 0.5, 1, lwd=3)
  for (i in 1:500)
  {
    points(gen_10_0.01[i,5], gen_10_0.01[i,6], pch=21, bg=adjustcolor(triangle_cols[1], alpha.f=0.75), cex=2)
  }
  box(lwd=2)
dev.off()  
  
pdf(file="triangles_for_triangles100.pdf", width=10, height=10)
plot(0, type="n", xlab="q", ylab="Q", xlim=c(0, 1), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25)
segments(0, 0, 0.5, 1, lwd=3)
segments(1, 0, 0.5, 1, lwd=3)
for (i in 1:500)
{
  points(gen_100[i,5], gen_100[i,6], pch=21, bg=adjustcolor(triangle_cols[1], alpha.f=0.75), cex=2)
}
box(lwd=2)
dev.off()  


