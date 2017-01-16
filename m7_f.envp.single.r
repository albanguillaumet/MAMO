#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO OBTAIN THE 'SENSITIVITY ENVELOP' USING A SINGLE CALIBRATION RUN

f.envp.single = function(
sp = "HCRE", output.dir = "C:/Alban/Alban HAWAII/HAKALAU/R code/CALIBRATION/hcre-p1", y.obs = y.obs.HCRE,
n.envp = 10, grad = seq(1.9, 1.0, length = 10), col = "blue",
d = read.table("C:/Alban/Alban HAWAII/HAKALAU/Data/Model calibration/Parameters/table.calib.txt", header = T, sep = "\t", dec = ".")
){

cex = 0.7
plot(grad, y.obs, main = sp, ylim = c(0, 1.1*max(y.obs)), xlim = c(1, 2), pch = 3, bty = "n", xlab = "Elevation (Km)", ylab = "# pairs / ha of native forest", col = col, cex = 0.9)

#----------

# read the calibration runs
setwd(output.dir)
d.sim = read.table("t.sim.txt", header = T, sep = "\t", dec = ".")
n = max(d.sim$total.sim, na.rm = TRUE)
d.sim = d.sim[1:n,]

#----------

# Identify which K.nb.1 has been used
K.nb.1 = rep(0, n)
if(sp == "IIWI" | sp == "APAP")  {
	v = d.sim$K.nb.1.high
	K.nb.1 = rep(NA, n)
	for(i in 1 : n) {
		if(round(v[i], 5) == 3.23333) K.nb.1[i] = 2003
		if(round(v[i], 5) == 9.58000) K.nb.1[i] = 2004
		if(round(v[i], 5) == 6.40667) K.nb.1[i] = 2003.5
	}
}
d.sim = cbind(d.sim, K.nb.1)

#----------

# extract the n.envp best run
r = d.sim[order(d.sim$data.fit, decreasing = TRUE),]
envp = r[1 : n.envp,]

# plot the n.envp best run
b = envp$total.sim
for(i in 1 :  n.envp) add.s(d = paste("s", b[i], ".1.Rdata", sep = ""), col.p = "pink", pch = 15, cex = 1.2, col.l = "pink", lty = 2)
points(grad, y.obs, pch = 3, col = col, cex = 0.9)

#----------

# export a table with parameter values for 11-13 parameters
stat.envp = as.data.frame(matrix(NA, nr = 5, nc = 13))
colnames(stat.envp) = v = c("s.ad", "rat.s", "fec", "fec.1", "rat.f", "K.b", "K.nb.1", "Sm.ac", "gamma.mov", "R.ter", "fidelity.ad", "m.natal",  "psi.DD")
rownames(stat.envp) = c("mean.envp", "sd.envp", "min.envp", "max.envp", "exp.mean")

for(j in 1 : 13) stat.envp[1, j] = mean( envp[ ,v[j] ], na.rm = TRUE )
for(j in 1 : 13) stat.envp[2, j] = sd( envp[ ,v[j] ], na.rm = TRUE )
for(j in 1 : 13) stat.envp[3, j] = min( envp[ ,v[j] ], na.rm = TRUE )
for(j in 1 : 13) stat.envp[4, j] = max( envp[ ,v[j] ], na.rm = TRUE )
for(j in 1 : 13) stat.envp[5, j] = mean( unique( d.sim[ ,v[j] ] ) )

#----------

# create a modified data set for another round of calibration step, if needed
if(dim(d)[1] == 1) rownames(d) = sp
if(dim(d)[1] == 8) rownames(d) = d$species

i = 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "s.ad_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "s.ad_max"] = stat.envp[5, i]
d[sp, "s.juv_min"] = d[sp, "s.ad_min"] / 2; d[sp, "s.juv_max"] = d[sp, "s.ad_max"] / 2
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "rat.s_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "rat.s_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "fec_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "fec_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "fec.1_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "fec.1_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "rat.f_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "rat.f_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "K.b_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "K.b_max"] = stat.envp[5, i]
i = i + 1
if(sp == "IIWI" | sp == "APAP") 	{ 
	if(stat.envp[1, i] >= stat.envp[5, i]) { d[sp, "K.nb.1_min"] = 2003.5; d[sp, "K.nb.1_max"] = 2004 }
	if(stat.envp[1, i] < stat.envp[5, i]) { d[sp, "K.nb.1_min"] = 2003; d[sp, "K.nb.1_max"] = 2003.5 }
} else { 
d[sp, "K.nb.1_min"] = d[sp, "K.nb.1_max"] = NA
}
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "Sm.ac_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "Sm.ac_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "gamma.mov_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "gamma.mov_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "R.ter_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "R.ter_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "fidelity.ad_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "fidelity.ad_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "m.natal_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "m.natal_max"] = stat.envp[5, i]
i = i + 1
if(stat.envp[1, i] >= stat.envp[5, i]) d[sp, "psi.DD_min"] = stat.envp[5, i]
if(stat.envp[1, i] < stat.envp[5, i]) d[sp, "psi.DD_max"] = stat.envp[5, i]

#----------

list(d.sim = d.sim, envp = envp, stat.envp = round(stat.envp, 5), d = d)

} # f.envp.single

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
f.envp.single = cmpfun(f.envp.single)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
