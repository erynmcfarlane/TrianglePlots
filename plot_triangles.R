### Triangle plots for all generations before this:
library(MetBrewer)
triangle_cols <- met.brewer("OKeeffe1", 20)

triangleplot<-function(filename, qQ12_data, generation){
  ### file name is a .pdf in quotes
  ### qQ12_data is q in first column, Q12 in second column, and each row is an individual
  ### generation is to set colours to be consistent for each generation
  pdf(file=filename, width=10, height=10)
  plot(0, type="n", xlab="q", ylab="Q", xlim=c(0, 1), ylim=c(0, 1), las=1, cex.lab=1.5, cex.axis=1.25)
  segments(0, 0, 0.5, 1, lwd=3)
  segments(1, 0, 0.5, 1, lwd=3)
  for (i in 1:500)
  {
    points(qQ12_data[i,1], qQ12_data[i,2], pch=21, bg=adjustcolor(triangle_cols[generation], alpha.f=0.75), cex=2)
  }
  box(lwd=2)
  dev.off()  
}