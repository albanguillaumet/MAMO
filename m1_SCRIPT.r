#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Chapter 2 
#THE BASIC MAMO FUNCTION: mamo

# BLOCK 1 
# SOURCE CODE - Please copy and paste the code below before each use of the software (mandatory)

# Location of the source files
setwd("C:/Programs/MAMO")

# load libraries
source("m2_libraries.r")

# load data sets
source("m3_data.r")

# load the mamo function
source("m4_mamo.r")

# load ancillary functions
source("m5_ancillary.functions.r")

# load the f.calibr function
source("m6_f.calibr.r")

# load the f.envp.single function
source("m7_f.envp.single.r")

# load the f.run function
source("m8_f.run.r")

# load the f.plot.univar function
source("m9_f.plot.univar.r")

# load the f.plot.bivar function
source("m10_f.plot.bivar.r")

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 2
# run an example of MAMO

mamo.ex = mamo(
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, t.b = 242, f.nb.1 = 0.5, min.fledg = 100, peak.fledg = (2/3), SD.fledg = 0,
# Initial conditions
init.1 = 300, init.2 = 300, 
# Survival (ad = annual, juv = from fledging to breeding age)
s.ad = 0.729, rat.s = 0, s.juv = (0.729/2),
# Reproduction and habitat quality
fec = 3, fec.1 = 3, rat.f = 0, K.b = 600, thr.DD = 300, K.nb.1 = K.nb.1.2004, K.nb.2 =  K.nb.2, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, Sm.ac = 0.13, 
# Movements
gamma.mov = 0.541, calc.gamma.d = "fast.risky", n.sim.disp = 10000, R.ter = 0.0234, fidelity.ad = 0.95, m.natal = 0.3, SD.natal = 0.3, psi.DD = 1,
# Other options 
add.cline = TRUE)

# see output

mamo.ex

# Add real data points for comparison

points(seq(1.9, 1.0, length.out = 10), (y.obs.IIWI.1*100), pch = 3, col = "black", cex = 1, lwd = 4)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Chapter 3
#CALIBRATING MAMO

# BLOCK 3
# CALIBRATE MAMO FOR A SINGLE SPECIES: THE ALTITUDINAL MIGRANT IIWI ('IIWI paper')

run_f.calibr = FALSE

#---

if(run_f.calibr == TRUE) {

f.calibr( 
# species
sp = "IIWI", output.dir = "C:/Programs/MAMO/CALIBRATION/run/IIWI.1_1",
y.obs = y.obs.IIWI.1, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.IIWI.1.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 3, input.direct.s.ad = TRUE, s.ad.direct = c(0.729, 0.78, 0.877), n.rat.s = 1,
# Reproduction and habitat quality
n.fec = 3, input.direct.fec = TRUE, fec.direct = c(3, 2.444, 1.5), n.fec.1 = 3, input.direct.fec.1 = TRUE, fec.1.direct = c(3, 2.444, 1.5),
paired.s.ad.fec = TRUE, n.rat.f = 1, n.K.b = 3, K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg), n.K.nb.1 = 3, K.nb.2 = K.nb.2, 
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 3, input.direct.Sm.ac = TRUE, 
Sm.ac.direct = c(0.02, 0.07, 0.13),
# Movements
n.gamma.mov = 4, input.direct.gamma.mov = TRUE, gamma.mov.direct = c(-10, 0.159, 0.541, 10), calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 1, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1)

}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 4
# Get the 'sensitivity envelop' for the 'IIWI paper'

IIWI.1_1.calibr = f.envp.single(sp = "IIWI", output.dir = "C:/Programs/MAMO/CALIBRATION/run/IIWI.1_1", 
y.obs = y.obs.IIWI.1, n.envp = 10, col = "blue",
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.IIWI.1.txt", header = T, sep = "\t", dec = ".") )

#---

envp.IIWI.1_1 = IIWI.1_1.calibr$envp

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Chapter 4
#SIMULATION STUDY

#BLOCK 5
# Run management scenarios based on the 'sensitivity envelop' for the 'IIWI paper'

run_f.run = FALSE

#---

if(run_f.run == TRUE) {

IIWI.1_1.run = f.run( 
# species-specific
sp = "IIWI",
envp = envp.IIWI.1_1,
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1",
y.obs = y.obs.IIWI.1, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.IIWI.1.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 1,
mg.RAT = rep(1, 10),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 3,
mg.K.NB.2 = list( rep(1, 10), c(rep(1.5, 3), rep(1, 7)), c(rep(2, 3), rep(1, 7)) ), 
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = rep(0.5, 10),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 3, evol.AC = c(3, 2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 6
# Fig. 5 of the 'IIWI paper'

f.plot.univar(
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
var.excl = NA, val.excl = NA,
ylab = "# pairs IIWI / ha of native forest", main = NA, col.obs = "red", margin.up = 0,
y.obs = y.obs.IIWI.1,
wch.plot = "RISK", col.cat = c("grey90", "grey70", "pink", "grey40", "black"), 
val.RISK = NA, val.AC = 3, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
leg.x = 1, leg.y = 4.5, leg.text = c("Past", "Pres / 2", "Pres-obs", "Pres-sim", "2100/2", "2100"), 
leg.pch = c(15,15,3,15,15,15), leg.lty = c(2,2,1,2,2,2), leg.col = c("grey90", "grey70", "red", "pink", "grey40", "black")
)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 7
# Fig. 6 of the 'IIWI paper'

z = f.plot.bivar(
	output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
	var.excl = NA, val.excl = NA,
	take.subset = TRUE, val.RISK = NA, val.AC = NA, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
	y = "np.metapop",  ylab = "# pairs IIWI in metapopulation", main = NA,
	x = "AC", xlab = "Malaria mortality", 
	o = "RISK", title.o = "RISK", legend.o = TRUE, l.pos = NA
	)
	
	# add observed data point
	
	np.metapop.obs = sum(y.obs.IIWI.1* 100)  * 2 # *100 :  # pairs / ha -> # pairs / km2; sum: patch -> all gradient (1 column); *2 -> 2 columns
	points(3, np.metapop.obs, col = "red", pch = 3, cex = 1.8, lwd = 2)
	
	# differentiate the different categories of the parameter gamma.mov in the future
	k = 0.2
	di = z$d.subset
	di = di[di$AC == 3 & di$RISK == 5,]
	di.m1 = di[di$gamma.mov == 0.159,]; di.m2 = di[di$gamma.mov == 0.541,]; di.m3 = di[di$gamma.mov == 10,]
	points(di.m1$AC+k, pch = 25, cex = 1.8, di.m1$np.metapop, col = "blue", lwd = 2)
	quantile(di.m1$np.metapop, probs = c(0.025, 0.975))
	points(di.m2$AC+k, pch = 25, cex = 1.8, di.m2$np.metapop, col = "purple", lwd = 2)
	quantile(di.m2$np.metapop, probs = c(0.025, 0.975))
	points(di.m3$AC+k, pch = 25, cex = 1.8, di.m3$np.metapop, col = "red", lwd = 2)
	quantile(di.m3$np.metapop, probs = c(0.025, 0.975))
	# statistical relationship between np.metapop and gamma.mov
	summary(lm(np.metapop ~ gamma.mov, di))
	r.f = di[, "batch"]
	summary(lmer(np.metapop ~ gamma.mov + (1 | r.f), di))
	# relationship between more variables
	d.acp = subset(di, select = c(np.metapop, fec, K.b, m.elev.K.nb.1, gamma.mov, fidelity.ad, m.natal, psi.DD) )
	cor(d.acp)
	acp = dudi.pca(d.acp, scannf = FALSE, nf = 2)
	pc1 = (acp$l1)$RS1; pc2 = (acp$l1)$RS2
	acp; acp$co

	# calculate the fraction of IIWI remaining as compared to the pre-malaria era
	d1 = z$d.subset[z$d.subset$RISK == 1, ]
	d2 = z$d.subset[z$d.subset$AC == 3 & z$d.subset$RISK == 2,]
	d3 = z$d.subset[z$d.subset$AC == 3 & z$d.subset$RISK == 3,]
	d4 = z$d.subset[z$d.subset$AC == 3 & z$d.subset$RISK == 4,]
	d5 = z$d.subset[z$d.subset$AC == 3 & z$d.subset$RISK == 5,]
	#	
	c(mean(d1$np.metapop), mean(d2$np.metapop), mean(d3$np.metapop), mean(d4$np.metapop), mean(d5$np.metapop))
	#
	np.metapop.obs / mean(d1$np.metapop) # currently observed / pre-malaria
	mean(d3$np.metapop)  / mean(d1$np.metapop) # current (simulated) / pre-malaria
	mean(d5$np.metapop) / mean(d1$np.metapop) # future (simulated) / pre-malaria
	quantile(d5$np.metapop, probs = c(0.025, 0.975)) / mean(d1$np.metapop) # future (simulated) / pre-malaria
	
	# Median IIWI elevation in the study area (1000-1900 m)
	summary(d1$m.elev.np) # pre-malaria
	summary(d3$m.elev.np) # present
	summary(d5$m.elev.np) # future
	
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 8
# Fig. 7 of the 'IIWI paper'

layout(matrix(c(1:3), 3, 1)); par(mar = c(3.6,3.6,0,0)+0.5, cex.main = 1)

i = f.plot.bivar(
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
var.excl = NA, val.excl = NA,
take.subset = TRUE, val.RISK = NA, val.AC = NA, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
y = "np.high",  ylab = "# pairs IIWI - 1800m", main = NA,
x = "AC", xlab = "", 
o = "RISK", title.o = "RISK", legend.o = TRUE, l.pos = NA
)

# add observed data point
np.1800.obs = y.obs.IIWI.1[2] * 100 
points(3, np.1800.obs, col = "red", pch = 3, cex = 1.8, lwd = 2)

# differentiate the different categories of the parameter gamma.mov in the future
k = 0.2
di = i$d.subset
di = di[di$AC == 3 & di$RISK == 5,]
di.m1 = di[di$gamma.mov == 0.159,]; di.m2 = di[di$gamma.mov == 0.541,]; di.m3 = di[di$gamma.mov == 10,]
points(di.m1$AC+k, pch = 25, cex = 1.8, di.m1$np.high, col = "blue", lwd = 2)
quantile(di.m1$np.high, probs = c(0.025, 0.975))
points(di.m2$AC+k, pch = 25, cex = 1.8, di.m2$np.high, col = "purple", lwd = 2)
quantile(di.m2$np.high, probs = c(0.025, 0.975))
points(di.m3$AC+k, pch = 25, cex = 1.8, di.m3$np.high, col = "red", lwd = 2)
quantile(di.m3$np.high, probs = c(0.025, 0.975))
# statistical relationship between np.high and gamma.mov
summary(lm(np.high ~ gamma.mov, di))
r.f = di[, "batch"]
summary(lmer(np.high ~ gamma.mov + (1 | r.f), di))
# relationship between more variables
d.acp = subset(di, select = c(np.high, fec, K.b, m.elev.K.nb.1, gamma.mov, fidelity.ad, m.natal, psi.DD) )
cor(d.acp)
acp = dudi.pca(d.acp, scannf = FALSE, nf = 2)
pc1 = (acp$l1)$RS1; pc2 = (acp$l1)$RS2
acp; acp$co

# calculate the fraction of IIWI remaining as compared to the pre-malaria era		
d1 = i$d.subset[i$d.subset$RISK == 1,]
d3 = i$d.subset[i$d.subset$AC == 3 & i$d.subset$RISK == 3,]
d5 = i$d.subset[i$d.subset$AC == 3 & i$d.subset$RISK == 5,]
np.1800.obs / mean(d1$np.high) # currently observed / pre-malaria
mean(d3$np.high) / mean(d1$np.high) # current (simulated) / pre-malaria
mean(d5$np.high) / mean(d1$np.high) # future (simulated) / pre-malaria
quantile(d5$np.high, probs = c(0.025, 0.975)) / mean(d1$np.high) # future (simulated) / pre-malaria

#---

j = f.plot.bivar( 
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
var.excl = NA, val.excl = NA,
take.subset = TRUE, val.RISK = NA, val.AC = NA, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
y = "np.mid",  ylab = "# pairs IIWI - 1500 m", main = NA,
x = "AC", xlab = "", 
o = "RISK", title.o = "RISK", legend.o = TRUE, l.pos = NA
)

# add observed data point
np.1500.obs = y.obs.IIWI.1[5] * 100 
points(3, np.1500.obs, col = "red", pch = 3, cex = 1.8, lwd = 2)

# calculate the fraction of IIWI remaining as compared to the pre-malaria era	
d1 = j$d.subset[j$d.subset$RISK == 1,]
d3 = j$d.subset[j$d.subset$AC == 3 & j$d.subset$RISK == 3,]
d5 = j$d.subset[j$d.subset$AC == 3 & j$d.subset$RISK == 5,]
np.1500.obs / mean(d1$np.mid) # currently observed / pre-malaria
mean(d3$np.mid) / mean(d1$np.mid) # current (simulated) / pre-malaria
mean(d5$np.mid) / mean(d1$np.mid) # future (simulated) / pre-malaria
quantile(d5$np.mid, probs = c(0.025, 0.975)) / mean(d1$np.mid) # future (simulated) / pre-malaria

#---

k = f.plot.bivar( 
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
var.excl = NA, val.excl = NA,
take.subset = TRUE, val.RISK = NA, val.AC = NA, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
y = "np.low",  ylab = "# pairs IIWI - 1200 m", main = NA,
x = "AC", xlab = "Malaria mortality", 
o = "RISK", title.o = "RISK", legend.o = TRUE
)

# add observed data point
np.1200.obs = y.obs.IIWI.1[8] * 100 
points(3, np.1200.obs, col = "red", pch = 3, cex = 1.8, lwd = 2)

# calculate the fraction of IIWI remaining as compared to the pre-malaria era	
d1 = k$d.subset[k$d.subset$RISK == 1,]
d3 = k$d.subset[k$d.subset$AC == 3 & k$d.subset$RISK == 3,]
d5 = k$d.subset[k$d.subset$AC == 3 & k$d.subset$RISK == 5,]
np.1200.obs / mean(d1$np.low) # currently observed / pre-malaria
mean(d3$np.low) / mean(d1$np.low) # current (simulated) / pre-malaria
mean(d5$np.low) / mean(d1$np.low) # future (simulated) / pre-malaria
quantile(d5$np.low, probs = c(0.025, 0.975)) / mean(d1$np.low) # future (simulated) / pre-malaria

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 9
# Fig. 8 of the 'IIWI paper'

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
ylab = "# pairs IIWI / ha of native forest", fig.title = NA, main = NA, margin.up = -0.5,
wch.plot = "RES.2", col.cat = c("black", "grey70", "grey90"), add.cline.cat = rep("TRUE", 3),
val.RISK = 3, val.AC = 3, val.RAT = 1, val.RES.1 = 1, val.RES.2 = NA,
add.leg = TRUE,
leg.x = 1, leg.y = 5,
leg.text = c("Present", "Present (RES×1.5 | high)", "Present (RES×2 | high)", "2100/2", "2100/2 (RES×1.5 | high)", "2100/2 (RES×2 | high)", "2100", "2100 (RES×1.5 | high)", "2100 (RES×2 | high)"),
leg.pch = NA, leg.lty = c(rep(1,3), rep(2,3),rep(3,3)), leg.col = c(rep(c("black", "grey70", "grey90"), 3)),
create.plot = TRUE, lty.mean = 1
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
ylab = "# pairs IIWI / ha of native forest", fig.title = NA, main = NA, margin.up = -0.5,
wch.plot = "RES.2", col.cat = c("black", "grey70", "grey90"), add.cline.cat = rep("TRUE", 3),
val.RISK = 4, val.AC = 3, val.RAT = 1, val.RES.1 = 1, val.RES.2 = NA,
add.leg = FALSE,
create.plot = FALSE, lty.mean = 2
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
ylab = "# pairs IIWI / ha of native forest", fig.title = NA, main = NA, margin.up = -0.5,
wch.plot = "RES.2", col.cat = c("black", "grey70", "grey90"), add.cline.cat = rep("TRUE", 3),
val.RISK = 5, val.AC = 3, val.RAT = 1, val.RES.1 = 1, val.RES.2 = NA,
add.leg = FALSE,
create.plot = FALSE, lty.mean = 3
)

abline(v = 1.7, lty = 4)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 10
# Fig. S1 in Appendix S4 of the 'IIWI paper'

setwd("C:/Programs/MAMO/CALIBRATION/run/IIWI.1_1")
d = read.table("t.sim.txt", header = T, sep = "\t", dec = ".")[1:5184,]
d = d[d$Sm.ac == 0.13,]

d1 = d[d$m.natal == 0.30,]
d2 = d[d$m.natal == 0.95,]
wilcox.test(d1$np.low, d2$np.low)
quantile(d1$np.low, probs = c(0.025, 0.975))
quantile(d2$np.low, probs = c(0.025, 0.975))
y.obs.IIWI.1[8]

ggplot(d, aes(x = np.low, fill = as.factor(m.natal))) +
    geom_histogram(binwidth = 1, position = "dodge") +
	 guides(fill = guide_legend(title = "m.natal")) +
	 geom_vline(aes(xintercept = y.obs.IIWI.1[8]*100), color = "red", linetype = "dashed", size = 1) +
	 theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
	 ylab("Frequency") +
	 xlab("# pairs IIWI / patch at 1200 m")

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
