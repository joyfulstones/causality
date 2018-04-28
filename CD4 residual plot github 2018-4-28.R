# CD4 residual plot for Model IV

res=read.table("c:/data/cd4res4.txt", header=T)	
mres=read.table("c:/data/cd4mres4.txt", header=T)

############################################################################################

# Plots used in the paper

par(mfrow=c(2, 2))

trtc1=ifelse(res$trt==0, "ddI", "ddC")
trtc3=ifelse(mres$trt==0, "ddI", "ddC")

boxplot(res ~ trtc1, data = res, col = "lightgray", ylab="Residual", main="A")
boxplot(res ~ trtc3, data = mres, col = "lightgray", ylab="Martigale Residual", main="B")

plot(res$year, res$res, xlab="year", ylab="Residual", col="black", type="p", pch=".", cex=.6, main="C")
lines(lowess(res$year, res$res, iter=0))

hemobl2=mres$hemobl+12
plot(hemobl2, mres$res, xlab="Base Hgb", ylab="Martigale Residual", col="black", type="p", pch=".", cex=.6, main="D")
lines(lowess(hemobl2, mres$res, iter=0))
