res=read.table("c:/data/prores4.txt", header=T)	
mres=read.table("c:/data/promres4.txt", header=T)


############################################################################################

# Plots used in the paper

par(mfrow=c(2, 2))

trtc1=ifelse(res$trt==0, "Placebo", "Prednisone")
trtc3=ifelse(mres$trt==0, "Placebo", "Prednisone")

boxplot(res ~ trtc1, data = res, col = "lightgray", ylab="Residual", main="A")
boxplot(res ~ trtc3, data = mres, col = "lightgray", ylab="Martigale Residual", main="B")

plot(res$start, res$res, xlab="time", ylab="Residual", col="black", type="p", pch=".", cex=.6, main="C")
lines(lowess(res$start, res$res, iter=0))

############################################################################################
