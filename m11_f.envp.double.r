#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO OBTAIN THE 'SENSITIVITY ENVELOP' USING TWO DIFFERENT CALIBRATION RUN

f.envp.double = function(
sp = "HCRE",
envp.1 = hcre.p1$envp, d1 = "C:/Alban/Alban HAWAII/HAKALAU/R code/CALIBRATION/hcre-p1", 
envp.2 = hcre.p2$envp, d2 = "C:/Alban/Alban HAWAII/HAKALAU/R code/CALIBRATION/hcre-p2", 
n.envp = 10, y.obs = y.obs.HCRE, col = "blue") {

#----------

n1 = dim(envp.1)[1]; n2 = dim(envp.2)[1]

wch.envp = rep(1, n1); envp.1 = cbind(wch.envp, envp.1)
wch.envp = rep(2, n2); envp.2 = cbind(wch.envp, envp.2)

envp.1.s = subset(envp.1, select = c("s.ad", "rat.s", "fec", "fec.1", "rat.f", "K.b", "K.nb.1", "Sm.ac", "gamma.mov", "R.ter", "fidelity.ad", "m.natal",  "psi.DD"))
envp.2.s = subset(envp.2, select = c("s.ad", "rat.s", "fec", "fec.1", "rat.f", "K.b", "K.nb.1", "Sm.ac", "gamma.mov", "R.ter", "fidelity.ad", "m.natal",  "psi.DD"))

i.id = j.rm = c()
	for(i in 1 : n1) {
		for(j in 1 : n2) {
			u = unique(as.vector(envp.1.s[i,] == envp.2.s[j,]))
				if( length(u) == 1 & u[1] == TRUE ) { i.id = c(i.id, i); j.rm = c(j.rm, j) ; break }
		}
	}
j.rm = unique(j.rm)
	
if( length(j.rm) > 0 ) envp.2 = envp.2[-j.rm,]

envp = rbind(envp.1, envp.2)
envp = envp[order(envp[,"data.fit"], decreasing = TRUE),]
envp = envp[1:n.envp,]

# plot
cex = 0.7
grad = seq(envp$grad.max[1], envp$grad.min[1], length.out = envp$nr[1])

plot(grad, y.obs, main = sp, ylim = c(0, 1.1*max(y.obs)), xlim = c(1, 2), pch = 3, bty = "n", xlab = "Elevation (Km)", ylab = "# pairs / ha of native forest", col = col, cex = 0.9)

# plot the n.envp best run

for(i in 1 : n.envp) {
	if(envp$wch.envp[i] == 1) setwd(d1) else setwd(d2) 
	add.s(d = paste("s", envp$total.sim[i], ".1.Rdata", sep = ""), col.p = "pink", pch = 15, cex = 1.2, col.l = "pink", lty = 2)
	points(grad, y.obs, pch = 3, col = col, cex = 0.9)
}

envp

} # f.envp.double

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
f.envp.double = cmpfun(f.envp.double)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

