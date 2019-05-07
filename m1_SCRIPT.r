#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

# load the f.envp.double function
source("m11_f.envp.double.r")

# load the f.envp.double function
source("m12_f.plot.composite.r")

# load the f.data.species function
source("m13_f.data.species.r")

# load the f.data.species function
source("m14_f.data.community.r")

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

#edit(IIWI.1_1.calibr$d.sim)
#edit(IIWI.1_1.calibr$envp)
#edit(IIWI.1_1.calibr$stat.envp)
#edit(IIWI.1_1.calibr$d)

#---

envp.IIWI.1_1 = IIWI.1_1.calibr$envp

#setwd("C:/Programs/MAMO/SAUV/code Fig IIWI paper")
#write.table(envp.IIWI.1_1, "envp.IIWI.1_1.txt", sep = "\t", dec = ".", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
# Color version of Fig. 6 of the 'IIWI paper'

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
# Color version of Fig. 7 of the 'IIWI paper'

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
leg.text = c("Present", "Present (RESÃ1.5 | high)", "Present (RESÃ2 | high)", "2100/2", "2100/2 (RESÃ1.5 | high)", "2100/2 (RESÃ2 | high)", "2100", "2100 (RESÃ1.5 | high)", "2100 (RESÃ2 | high)"),
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
# Appendix S4 Fig. S1 of the 'IIWI paper'

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

#MAMO AT THE COMMUNITY LEVEL

# BLOCK 11
# First run of MAMO calibration for each species separately (8 native Hakalau species) using f.calibr

run_f.calibr = FALSE

#---

if(run_f.calibr == TRUE) {

f.calibr( 
# species
sp = "AKIP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/akip-p1",
y.obs = y.obs.AKIP,
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "AKEP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/akep-p1",
y.obs = y.obs.AKEP,
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "HCRE", output.dir = "C:/Programs/MAMO/CALIBRATION/run/hcre-p1",
y.obs = y.obs.HCRE,
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "IIWI", output.dir = "C:/Programs/MAMO/CALIBRATION/run/iiwi-p1",
y.obs = y.obs.IIWI,
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg), n.K.nb.1 = 3, K.nb.2 = K.nb.2, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 3, input.direct.gamma.mov = TRUE, gamma.mov.direct = c(0.159, 0.541, 10), calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "ELEP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/elep-p1",
y.obs = y.obs.ELEP,
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "OMAO", output.dir = "C:/Programs/MAMO/CALIBRATION/run/omao-p1",
y.obs = y.obs.OMAO,
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "APAP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/apap-p1",
y.obs = y.obs.APAP,
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg), n.K.nb.1 = 3, K.nb.2 = K.nb.2, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 3, input.direct.gamma.mov = TRUE, gamma.mov.direct = c(0.159, 0.541, 10), calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "HAAM", output.dir = "C:/Programs/MAMO/CALIBRATION/run/haam-p1",
y.obs = y.obs.HAAM,
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."), 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 12
# Obtain the input for the second calibration run using f.envp.single

layout(matrix(c(1:8), 4, 2)); par(mar = c(2.2,2.2,0,0)+0.5, cex.main = 1)

akip.p1 = f.envp.single(sp = "AKIP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/akip-p1",
y.obs = y.obs.AKIP, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

akep.p1 = f.envp.single(sp = "AKEP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/akep-p1",
y.obs = y.obs.AKEP, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

hcre.p1 = f.envp.single(sp = "HCRE", output.dir = "C:/Programs/MAMO/CALIBRATION/run/hcre-p1",
y.obs = y.obs.HCRE, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

iiwi.p1 = f.envp.single(sp = "IIWI", output.dir = "C:/Programs/MAMO/CALIBRATION/run/iiwi-p1",
y.obs = y.obs.IIWI, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

elep.p1 = f.envp.single(sp = "ELEP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/elep-p1",
y.obs = y.obs.ELEP, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

omao.p1 = f.envp.single(sp = "OMAO", output.dir = "C:/Programs/MAMO/CALIBRATION/run/omao-p1",
y.obs = y.obs.OMAO, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

apap.p1 = f.envp.single(sp = "APAP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/apap-p1",
y.obs = y.obs.APAP, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

haam.p1 = f.envp.single(sp = "HAAM", output.dir = "C:/Programs/MAMO/CALIBRATION/run/haam-p1",
y.obs = y.obs.HAAM, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#BLOCK 13
# Second round of MAMO calibration using f.calibr

run_f.calibr = FALSE

#---

if(run_f.calibr == TRUE) {

f.calibr( 
# species
sp = "AKIP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/akip-p2",
y.obs = y.obs.AKIP,
d = akip.p1$d, 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "AKEP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/akep-p2",
y.obs = y.obs.AKEP,
d = akep.p1$d, 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "HCRE", output.dir = "C:/Programs/MAMO/CALIBRATION/run/hcre-p2",
y.obs = y.obs.HCRE,
d = hcre.p1$d, 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "IIWI", output.dir = "C:/Programs/MAMO/CALIBRATION/run/iiwi-p2",
y.obs = y.obs.IIWI,
d = iiwi.p1$d, 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg), n.K.nb.1 = 2, K.nb.2 = K.nb.2, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 2, input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "ELEP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/elep-p2",
y.obs = y.obs.ELEP,
d = elep.p1$d, 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "OMAO", output.dir = "C:/Programs/MAMO/CALIBRATION/run/omao-p2",
y.obs = y.obs.OMAO,
d = omao.p1$d, 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "APAP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/apap-p2",
y.obs = y.obs.APAP,
d = apap.p1$d, 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg), n.K.nb.1 = 2, K.nb.2 = K.nb.2, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 2, input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

f.calibr( 
# species
sp = "HAAM", output.dir = "C:/Programs/MAMO/CALIBRATION/run/haam-p2",
y.obs = y.obs.HAAM,
d = haam.p1$d, 
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 2, input.direct.s.ad = FALSE, n.rat.s = 2,
# Reproduction and habitat quality
n.fec = 2, input.direct.fec = FALSE, n.fec.1 = 2, input.direct.fec.1 = FALSE, paired.s.ad.fec = FALSE, 
n.rat.f = 2, n.K.b = 2, K.nb.1 = NA, n.K.nb.1 = 1, K.nb.2 = NA, reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1, n.Sm.ac = 2, input.direct.Sm.ac = FALSE, 
# Movements
n.gamma.mov = 1,  input.direct.gamma.mov = FALSE, calc.gamma.d = "fast.risky",
n.sim.disp = 10000, n.R.ter = 2, n.fidelity.ad = 2, n.m.natal = 2, n.psi.DD = 2,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 1, design = "simple", batch = 1
) 

}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#BLOCK 14
# get the 'sensitivity envelop' for each species based on the second round of calibration

layout(matrix(c(1:8), 4, 2)); par(mar = c(2.2,2.2,0,0)+0.5, cex.main = 1)

akip.p2 = f.envp.single(sp = "AKIP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/akip-p2",
y.obs = y.obs.AKIP, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

akep.p2 = f.envp.single(sp = "AKEP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/akep-p2",
y.obs = y.obs.AKEP, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

hcre.p2 = f.envp.single(sp = "HCRE", output.dir = "C:/Programs/MAMO/CALIBRATION/run/hcre-p2",
y.obs = y.obs.HCRE, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

iiwi.p2 = f.envp.single(sp = "IIWI", output.dir = "C:/Programs/MAMO/CALIBRATION/run/iiwi-p2",
y.obs = y.obs.IIWI, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

elep.p2 = f.envp.single(sp = "ELEP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/elep-p2",
y.obs = y.obs.ELEP, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

omao.p2 = f.envp.single(sp = "OMAO", output.dir = "C:/Programs/MAMO/CALIBRATION/run/omao-p2",
y.obs = y.obs.OMAO, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

apap.p2 = f.envp.single(sp = "APAP", output.dir = "C:/Programs/MAMO/CALIBRATION/run/apap-p2",
y.obs = y.obs.APAP, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

haam.p2 = f.envp.single(sp = "HAAM", output.dir = "C:/Programs/MAMO/CALIBRATION/run/haam-p2",
y.obs = y.obs.HAAM, n.envp = 10, col = "blue", 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 15
# Calibration: get the final 'sensitivity envelop' for each species
# Color version of Fig. 2 of the community-level paper


layout(matrix(c(1:8), 4, 2)); par(mar = c(2.2,2.2,0,0)+0.5, cex.main = 1)

envp.akip = f.envp.double(
sp = "AKIP",
envp.1 = akip.p1$envp, d1 = "C:/Programs/MAMO/CALIBRATION/run/akip-p1", 
envp.2 = akip.p2$envp, d2 = "C:/Programs/MAMO/CALIBRATION/run/akip-p2", 
n.envp = 10, y.obs = y.obs.AKIP, col = "blue")

envp.akep = f.envp.double(
sp = "AKEP",
envp.1 = akep.p1$envp, d1 = "C:/Programs/MAMO/CALIBRATION/run/akep-p1", 
envp.2 = akep.p2$envp, d2 = "C:/Programs/MAMO/CALIBRATION/run/akep-p2", 
n.envp = 10, y.obs = y.obs.AKEP, col = "blue")

envp.hcre = f.envp.double(
sp = "HCRE",
envp.1 = hcre.p1$envp, d1 = "C:/Programs/MAMO/CALIBRATION/run/hcre-p1", 
envp.2 = hcre.p2$envp, d2 = "C:/Programs/MAMO/CALIBRATION/run/hcre-p2", 
n.envp = 10, y.obs = y.obs.HCRE, col = "blue")

envp.iiwi = f.envp.double(
sp = "IIWI",
envp.1 = iiwi.p1$envp, d1 = "C:/Programs/MAMO/CALIBRATION/run/iiwi-p1", 
envp.2 = iiwi.p2$envp, d2 = "C:/Programs/MAMO/CALIBRATION/run/iiwi-p2", 
n.envp = 10, y.obs = y.obs.IIWI, col = "blue")

envp.elep = f.envp.double(
sp = "ELEP",
envp.1 = elep.p1$envp, d1 = "C:/Programs/MAMO/CALIBRATION/run/elep-p1", 
envp.2 = elep.p2$envp, d2 = "C:/Programs/MAMO/CALIBRATION/run/elep-p2", 
n.envp = 10, y.obs = y.obs.ELEP, col = "blue")

envp.omao = f.envp.double(
sp = "OMAO",
envp.1 = omao.p1$envp, d1 = "C:/Programs/MAMO/CALIBRATION/run/omao-p1", 
envp.2 = omao.p2$envp, d2 = "C:/Programs/MAMO/CALIBRATION/run/omao-p2", 
n.envp = 10, y.obs = y.obs.OMAO, col = "blue")

envp.apap = f.envp.double(
sp = "APAP",
envp.1 = apap.p1$envp, d1 = "C:/Programs/MAMO/CALIBRATION/run/apap-p1", 
envp.2 = apap.p2$envp, d2 = "C:/Programs/MAMO/CALIBRATION/run/apap-p2", 
n.envp = 10, y.obs = y.obs.APAP, col = "blue")

envp.haam = f.envp.double(
sp = "HAAM",
envp.1 = haam.p1$envp, d1 = "C:/Programs/MAMO/CALIBRATION/run/haam-p1", 
envp.2 = haam.p2$envp, d2 = "C:/Programs/MAMO/CALIBRATION/run/haam-p2", 
n.envp = 10, y.obs = y.obs.HAAM, col = "blue")

#setwd("C:/Programs/MAMO/SAUV")
#write.table(envp.haam, "envp.haam.txt", sep = "\t", dec = ".", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 16
# Run management scenarios based on the 'sensitivity envelop' for each species independently

run_f.run = FALSE

#---

if(run_f.run == TRUE) {

akip.run.10 = f.run( 
# species-specific
sp = "AKIP",
envp = envp.akip,
output.dir = "C:/Programs/MAMO/RUN/akip.10",
y.obs = y.obs.AKIP, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 10), c(rep(0.3, 4), rep(1, 6)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 10), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 6)),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

akip.run.11 = f.run( 
# species-specific
sp = "AKIP",
envp = envp.akip,
output.dir = "C:/Programs/MAMO/RUN/akip.11",
y.obs = c(NA, y.obs.AKIP), 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 11, nc = 2, grad = c(2000, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 11), c(rep(0.3, 4), rep(1, 7)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 11), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 11), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 7)),
alpha.b = alpha.b_nc.11, alpha.nb.1 = alpha.nb.1_nc.11, alpha.nb.2 = alpha.nb.2_nc.11, alpha.1 = alpha.1_nc.11,
alpha.b.2100 = alpha.b.2100_nc.11, alpha.nb.1.2100 = alpha.nb.1.2100_nc.11, alpha.nb.2.2100 = alpha.nb.2.2100_nc.11, alpha.1.2100 = alpha.1.2100_nc.11,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

#---

akep.run.10 = f.run( 
# species-specific
sp = "AKEP",
envp = envp.akep,
output.dir = "C:/Programs/MAMO/RUN/akep.10",
y.obs = y.obs.AKEP, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 10), c(rep(0.3, 4), rep(1, 6)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 10), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 6)),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

akep.run.11 = f.run( 
# species-specific
sp = "AKEP",
envp = envp.akep,
output.dir = "C:/Programs/MAMO/RUN/akep.11",
y.obs = c(NA, y.obs.AKEP), 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 11, nc = 2, grad = c(2000, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 11), c(rep(0.3, 4), rep(1, 7)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 11), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 11), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 7)),
alpha.b = alpha.b_nc.11, alpha.nb.1 = alpha.nb.1_nc.11, alpha.nb.2 = alpha.nb.2_nc.11, alpha.1 = alpha.1_nc.11,
alpha.b.2100 = alpha.b.2100_nc.11, alpha.nb.1.2100 = alpha.nb.1.2100_nc.11, alpha.nb.2.2100 = alpha.nb.2.2100_nc.11, alpha.1.2100 = alpha.1.2100_nc.11,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

#---

hcre.run.10 = f.run( 
# species-specific
sp = "HCRE",
envp = envp.hcre,
output.dir = "C:/Programs/MAMO/RUN/hcre.10",
y.obs = y.obs.HCRE, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 10), c(rep(0.3, 4), rep(1, 6)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 10), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 6)),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

hcre.run.11 = f.run( 
# species-specific
sp = "HCRE",
envp = envp.hcre,
output.dir = "C:/Programs/MAMO/RUN/hcre.11",
y.obs = c(NA, y.obs.HCRE), 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 11, nc = 2, grad = c(2000, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 11), c(rep(0.3, 4), rep(1, 7)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 11), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 11), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 7)),
alpha.b = alpha.b_nc.11, alpha.nb.1 = alpha.nb.1_nc.11, alpha.nb.2 = alpha.nb.2_nc.11, alpha.1 = alpha.1_nc.11,
alpha.b.2100 = alpha.b.2100_nc.11, alpha.nb.1.2100 = alpha.nb.1.2100_nc.11, alpha.nb.2.2100 = alpha.nb.2.2100_nc.11, alpha.1.2100 = alpha.1.2100_nc.11,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

#---

iiwi.run.10 = f.run( 
# species-specific
sp = "IIWI",
envp = envp.iiwi,
output.dir = "C:/Programs/MAMO/RUN/iiwi.10",
y.obs = y.obs.IIWI, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 10), c(rep(0.3, 4), rep(1, 6)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 10), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 6)),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

iiwi.run.11 = f.run( 
# species-specific
sp = "IIWI",
envp = envp.iiwi,
output.dir = "C:/Programs/MAMO/RUN/iiwi.11",
y.obs = c(NA, y.obs.IIWI), 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 11, nc = 2, grad = c(2000, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 11), c(rep(0.3, 4), rep(1, 7)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 11), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 11), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 7)),
alpha.b = alpha.b_nc.11, alpha.nb.1 = alpha.nb.1_nc.11, alpha.nb.2 = alpha.nb.2_nc.11, alpha.1 = alpha.1_nc.11,
alpha.b.2100 = alpha.b.2100_nc.11, alpha.nb.1.2100 = alpha.nb.1.2100_nc.11, alpha.nb.2.2100 = alpha.nb.2.2100_nc.11, alpha.1.2100 = alpha.1.2100_nc.11,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

#---

elep.run.10 = f.run( 
# species-specific
sp = "ELEP",
envp = envp.elep,
output.dir = "C:/Programs/MAMO/RUN/elep.10",
y.obs = y.obs.ELEP, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 10), c(rep(0.3, 4), rep(1, 6)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 10), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 6)),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

elep.run.11 = f.run( 
# species-specific
sp = "ELEP",
envp = envp.elep,
output.dir = "C:/Programs/MAMO/RUN/elep.11",
y.obs = c(NA, y.obs.ELEP), 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 11, nc = 2, grad = c(2000, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 11), c(rep(0.3, 4), rep(1, 7)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 11), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 11), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 7)),
alpha.b = alpha.b_nc.11, alpha.nb.1 = alpha.nb.1_nc.11, alpha.nb.2 = alpha.nb.2_nc.11, alpha.1 = alpha.1_nc.11,
alpha.b.2100 = alpha.b.2100_nc.11, alpha.nb.1.2100 = alpha.nb.1.2100_nc.11, alpha.nb.2.2100 = alpha.nb.2.2100_nc.11, alpha.1.2100 = alpha.1.2100_nc.11,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

#---

omao.run.10 = f.run( 
# species-specific
sp = "OMAO",
envp = envp.omao,
output.dir = "C:/Programs/MAMO/RUN/omao.10",
y.obs = y.obs.OMAO, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 10), c(rep(0.3, 4), rep(1, 6)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 10), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 6)),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

omao.run.11 = f.run( 
# species-specific
sp = "OMAO",
envp = envp.omao,
output.dir = "C:/Programs/MAMO/RUN/omao.11",
y.obs = c(NA, y.obs.OMAO), 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 11, nc = 2, grad = c(2000, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 11), c(rep(0.3, 4), rep(1, 7)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 11), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 11), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 7)),
alpha.b = alpha.b_nc.11, alpha.nb.1 = alpha.nb.1_nc.11, alpha.nb.2 = alpha.nb.2_nc.11, alpha.1 = alpha.1_nc.11,
alpha.b.2100 = alpha.b.2100_nc.11, alpha.nb.1.2100 = alpha.nb.1.2100_nc.11, alpha.nb.2.2100 = alpha.nb.2.2100_nc.11, alpha.1.2100 = alpha.1.2100_nc.11,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

#---

apap.run.10 = f.run( 
# species-specific
sp = "APAP",
envp = envp.apap,
output.dir = "C:/Programs/MAMO/RUN/apap.10",
y.obs = y.obs.APAP, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 10), c(rep(0.3, 4), rep(1, 6)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 10), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 6)),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

apap.run.11 = f.run( 
# species-specific
sp = "APAP",
envp = envp.apap,
output.dir = "C:/Programs/MAMO/RUN/apap.11",
y.obs = c(NA, y.obs.APAP), 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 11, nc = 2, grad = c(2000, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 11), c(rep(0.3, 4), rep(1, 7)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 11), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 11), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 7)),
alpha.b = alpha.b_nc.11, alpha.nb.1 = alpha.nb.1_nc.11, alpha.nb.2 = alpha.nb.2_nc.11, alpha.1 = alpha.1_nc.11,
alpha.b.2100 = alpha.b.2100_nc.11, alpha.nb.1.2100 = alpha.nb.1.2100_nc.11, alpha.nb.2.2100 = alpha.nb.2.2100_nc.11, alpha.1.2100 = alpha.1.2100_nc.11,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

#---

haam.run.10 = f.run( 
# species-specific
sp = "HAAM",
envp = envp.haam,
output.dir = "C:/Programs/MAMO/RUN/haam.10",
y.obs = y.obs.HAAM, 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 10), c(rep(0.3, 4), rep(1, 6)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 10), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 10), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 6)),
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
alpha.b.2100 = alpha.b.2100, alpha.nb.1.2100 = alpha.nb.1.2100, alpha.nb.2.2100 = alpha.nb.2.2100, alpha.1.2100 = alpha.1.2100,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

haam.run.11 = f.run( 
# species-specific
sp = "HAAM",
envp = envp.haam,
output.dir = "C:/Programs/MAMO/RUN/haam.11",
y.obs = c(NA, y.obs.HAAM), 
d = read.table("C:/Programs/MAMO/CALIBRATION/Starting parameters/param_calib.HAKALAU.txt", header = T, sep = "\t", dec = "."),
# Spatial structure
nr = 11, nc = 2, grad = c(2000, 1000), unit = 1,
# Time frame
T = 60, Tm = 5, SD.fledg = 0,
# Survival 
n.RAT = 2,
mg.RAT = list( rep(1, 11), c(rep(0.3, 4), rep(1, 7)) ),
# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 1,
mg.K.NB.1 = rep(1, 11), # no management
K.nb.2 = K.nb.2,
n.K.NB.2 = 1,
mg.K.NB.2 = rep(1, 11), # no management
reproduction.malaria = "simple",
# Malaria parameters (daily except Sm.ac)
n.RISK = 5, mg.RISK = c(rep(0.5, 4), rep(1, 7)),
alpha.b = alpha.b_nc.11, alpha.nb.1 = alpha.nb.1_nc.11, alpha.nb.2 = alpha.nb.2_nc.11, alpha.1 = alpha.1_nc.11,
alpha.b.2100 = alpha.b.2100_nc.11, alpha.nb.1.2100 = alpha.nb.1.2100_nc.11, alpha.nb.2.2100 = alpha.nb.2.2100_nc.11, alpha.1.2100 = alpha.1.2100_nc.11,
n.AC = 2, evol.AC = c(2, 1), 
# Movements
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
# Other options 
add.cline = FALSE,
# Simulations 
n.sim = 2, design = "simple", batch = 1
)

}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 17
# CREATE DATA SET FOR COMMUNITY-LEVEL ANALYSIS
# Step # 1 - create species-specific data sets

# weights to be applied to each species in the community-level average

w = as.data.frame(matrix(NA, nr = 8, nc = 3))
rownames(w) = c("AKIP", "AKEP", "HCRE", "IIWI", "ELEP", "OMAO", "APAP", "HAAM")
colnames(w) = c("w1", "w2", "w3")
w[,1] = 1
w[,2] = c(3, 3, 3, 2, 2, 2, 1, 1)
w[,3] = 1 / c(1900, 12000, 14000, 360000, 190000, 170000, 1300000, 810000)

# info needed to tease apart the separate effects of malaria transmission risk and management

table.risk = as.data.frame(matrix(c(1:5, c(0, 1, 1, 2, 2), c(0, 1, 0, 1, 0)), nc = 3))
colnames(table.risk) = c("RISK_run", "RISK_analysis", "mg.RISK")

# species-specific data sets

AKIP = f.data.species(sp = "AKIP", w = w, 
directory = c("C:/programs/MAMO/RUN/akip.10", "C:/programs/MAMO/RUN/akip.11"),
vec.hab = c(0, 1), table.risk = table.risk)

AKEP = f.data.species(sp = "AKEP", w = w, 
directory = c("C:/programs/MAMO/RUN/akep.10", "C:/programs/MAMO/RUN/akep.11"),
vec.hab = c(0, 1), table.risk = table.risk)

HCRE = f.data.species(sp = "HCRE", w = w, 
directory = c("C:/programs/MAMO/RUN/hcre.10", "C:/programs/MAMO/RUN/hcre.11"),
vec.hab = c(0, 1), table.risk = table.risk)

IIWI = f.data.species(sp = "IIWI", w = w, 
directory = c("C:/programs/MAMO/RUN/iiwi.10", "C:/programs/MAMO/RUN/iiwi.11"),
vec.hab = c(0, 1), table.risk = table.risk)

ELEP = f.data.species(sp = "ELEP", w = w, 
directory = c("C:/programs/MAMO/RUN/elep.10", "C:/programs/MAMO/RUN/elep.11"),
vec.hab = c(0, 1), table.risk = table.risk)

OMAO = f.data.species(sp = "OMAO", w = w, 
directory = c("C:/programs/MAMO/RUN/omao.10", "C:/programs/MAMO/RUN/omao.11"),
vec.hab = c(0, 1), table.risk = table.risk)

APAP = f.data.species(sp = "APAP", w = w, 
directory = c("C:/programs/MAMO/RUN/apap.10", "C:/programs/MAMO/RUN/apap.11"),
vec.hab = c(0, 1), table.risk = table.risk)

HAAM = f.data.species(sp = "HAAM", w = w, 
directory = c("C:/programs/MAMO/RUN/haam.10", "C:/programs/MAMO/RUN/haam.11"),
vec.hab = c(0, 1), table.risk = table.risk)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 18
# CREATE DATA SET FOR COMMUNITY-LEVEL ANALYSIS
# Step # 2 - create a community-level data set

d.c = f.data.community(
AKIP = AKIP, AKEP = AKEP, HCRE = HCRE, IIWI = IIWI, 
ELEP = ELEP, OMAO = OMAO, APAP = APAP, HAAM = HAAM,
var_ = c("mg.hab", "risk", "mg.risk", "ac", "mg.rat", "mg.res.1", "mg.res.2", "envp", "sim")
)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 19
# Statistical models: mixed-effect linear models

# Present

# take the subset of data corresponding to present (both for malaria mortality and transmission risk)
d.p = d.c[d.c$ac == 2 & d.c$risk == 1,]

# mixed linear model
m1.p = lmer(n.w1 ~ mg.risk + mg.hab + mg.rat + (1 | envp) + (1 | sim), data = d.p)
summary(m1.p)
shapiro.test(residuals(m1.p))

m2.p = lmer(n.w2 ~ mg.risk + mg.hab + mg.rat + (1 | envp) + (1 | sim), data = d.p)
summary(m2.p)
shapiro.test(residuals(m2.p))

m3.p = lmer(n.w3 ~ mg.risk + mg.hab + mg.rat + (1 | envp) + (1 | sim), data = d.p)
summary(m3.p)
shapiro.test(residuals(m3.p))

# Future

# take the subset of data corresponding to future malaria transmission risk and no evolution of malaria mortality
d.f = d.c[d.c$ac == 2 & d.c$risk == 2,]

# mixed linear model
m1.f = lmer(n.w1 ~ mg.risk + mg.hab + mg.rat + (1 | envp) + (1 | sim), data = d.f)
summary(m1.f)
shapiro.test(residuals(m1.f))

m2.f = lmer(n.w2 ~ mg.risk + mg.hab + mg.rat + (1 | envp) + (1 | sim), data = d.f)
summary(m2.f)
shapiro.test(residuals(m2.f))

m3.f = lmer(n.w3 ~ mg.risk + mg.hab + mg.rat + (1 | envp) + (1 | sim), data = d.f)
summary(m3.f)
shapiro.test(residuals(m3.f))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 20
# Effect size

f.effect.size = function(){

# f.delta function based on info provided in http://www.ats.ucla.edu/stat/r/faq/deltamethod.htm and http://www.inside-r.org/packages/cran/msm/docs/deltamethod
f.delta = function(m = m1.p, b = 2){ # b: risk = 2, hab = 3, rat = 4
	a = 1 # intercept
	cf = c(fixef(m)[a], fixef(m)[b])
	#cf
	k = matrix(NA, 2, 2)
	v = vcov(m)
	k[1,1] = v[a,a]
	k[1,2] = k[2,1] = v[a,b]
	k[2,2] = v[b,b]
	deltamethod(~ (x1+x2) / x1, cf, k)
}

#---

m = list(m1.p, m2.p, m3.p, m1.f, m2.f, m3.f)

est_ = list(); est_[[1]] = summary(m1.p)$coef[,1]; est_[[2]] = summary(m2.p)$coef[,1]
est_[[3]] = summary(m3.p)$coef[,1]; est_[[4]] = summary(m1.f)$coef[,1]
est_[[5]] = summary(m2.f)$coef[,1]; est_[[6]] = summary(m3.f)$coef[,1]

se_ = list(); se_[[1]] = summary(m1.p)$coef[,2]; se_[[2]] = summary(m2.p)$coef[,2]
se_[[3]] = summary(m3.p)$coef[,2]; se_[[4]] = summary(m1.f)$coef[,2]
se_[[5]] = summary(m2.f)$coef[,2]; se_[[6]] = summary(m3.f)$coef[,2]

fixef_ = list(); fixef_[[1]] = fixef(m1.p); fixef_[[2]] = fixef(m2.p); fixef_[[3]] = fixef(m3.p)
fixef_[[4]] = fixef(m1.f); fixef_[[5]] = fixef(m2.f); fixef_[[6]] = fixef(m3.f)

w = rep(NA, 6); for(i in 1 : 6) w[i] = unique(est_[[i]] == fixef_[[i]]) # True

es = as.data.frame(matrix(NA, nr = 6, nc = 22))
colnames(es) = c(
"Time", "response", 
"est.int", "est.risk", "est.hab", "est.rat", 
"se.int","se.risk", "se.hab", "se.rat", 
"eff.risk", "eff.hab", "eff.rat",
"se.eff.risk", "se.eff.hab", "se.eff.rat",
"CI1.eff.risk", "CI2.eff.risk",
"CI1.eff.hab", "CI2.eff.hab",
"CI1.eff.rat", "CI2.eff.rat"
)
es$Time = c(rep("Present", 3), rep("Future", 3))
es$response = rep(c("n.w1", "n.w2", "n.w3"), 2)
for(i in 1: 6) {
	es$est.int[i] = fixef_[[i]][1]
	es$est.risk[i] = fixef_[[i]][2]
	es$est.hab[i] = fixef_[[i]][3]
	es$est.rat[i] = fixef_[[i]][4]
	es$se.int[i] = se_[[i]][1]
	es$se.risk[i] = se_[[i]][2]
	es$se.hab[i] = se_[[i]][3]
	es$se.rat[i] = se_[[i]][4]
	es$eff.risk[i] = ( fixef_[[i]][2] + fixef_[[i]][1] ) / fixef_[[i]][1]
	es$eff.hab[i] = ( fixef_[[i]][3] + fixef_[[i]][1] ) / fixef_[[i]][1]
	es$eff.rat[i] = ( fixef_[[i]][4] + fixef_[[i]][1] ) / fixef_[[i]][1]
	es$se.eff.risk[i] = f.delta(m[[i]], 2)
	es$se.eff.hab[i] = f.delta(m[[i]], 3)
	es$se.eff.rat[i] = f.delta(m[[i]], 4)
}

for(i in 1: 6) {
	es$CI1.eff.risk[i] = es$eff.risk[i] - 2 * es$se.eff.risk[i]
	es$CI2.eff.risk[i] = es$eff.risk[i] + 2 * es$se.eff.risk[i]
	es$CI1.eff.hab[i] = es$eff.hab[i] - 2 * es$se.eff.hab[i]
	es$CI2.eff.hab[i] = es$eff.hab[i] + 2 * es$se.eff.hab[i]
	es$CI1.eff.rat[i] = es$eff.rat[i] - 2 * es$se.eff.rat[i]
	es$CI2.eff.rat[i] = es$eff.rat[i] + 2 * es$se.eff.rat[i]
}


nd = 3
for(i in 3: 22) es[,i] = round(es[,i], nd)

es

}

#---

es = f.effect.size()
es

#setwd("C:/programs/MAMO/SAUV")
#write.table(es, "effect.size.txt", sep = "\t", dec = ".", row.names = F)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 21
# Plot and compare alternative scenarios

f.plot.scenarios = function(d.c = d.c, es = es, Time = "Present",
val.mg.hab = NA, val.risk = NA, val.mg.risk = NA, val.ac = NA, val.mg.rat = NA, val.mg.res.1 = NA, val.mg.res.2 = NA,
y = "n.w1", ylim = c(0, 1), main = "w1 - present", add.x.labels = TRUE){

d.s = d.c

# Data subset

if(!is.na(val.mg.hab)) d.s = d.s[d.s$mg.hab == val.mg.hab,]
if(!is.na(val.risk)) d.s = d.s[d.s$risk == val.risk,]
if(!is.na(val.mg.risk)) d.s = d.s[d.s$mg.risk == val.mg.risk,]
if(!is.na(val.ac)) d.s = d.s[d.s$ac == val.ac,]
if(!is.na(val.mg.rat)) d.s = d.s[d.s$mg.rat == val.mg.rat,]
if(!is.na(val.mg.res.1)) d.s = d.s[d.s$mg.res.1 == val.mg.res.1,]
if(!is.na(val.mg.res.2)) d.s = d.s[d.s$mg.res.2 == val.mg.res.2,]

#---

n.s = dim(d.s)[1]
scenarios = rep(NA, n.s)

for(i in 1 : n.s) {

		if(d.s$mg.hab[i] == 0 & d.s$mg.risk[i] == 0 & d.s$mg.rat[i] == 0) scenarios[i] = 1
		if(d.s$mg.hab[i] == 0 & d.s$mg.risk[i] == 1 & d.s$mg.rat[i] == 0) scenarios[i] = 2
		if(d.s$mg.hab[i] == 0 & d.s$mg.risk[i] == 0 & d.s$mg.rat[i] == 1) scenarios[i] = 3
		if(d.s$mg.hab[i] == 1 & d.s$mg.risk[i] == 0 & d.s$mg.rat[i] == 0) scenarios[i] = 4
		if(d.s$mg.hab[i] == 0 & d.s$mg.risk[i] == 1 & d.s$mg.rat[i] == 1) scenarios[i] = 5
		if(d.s$mg.hab[i] == 1 & d.s$mg.risk[i] == 1 & d.s$mg.rat[i] == 0) scenarios[i] = 6
		if(d.s$mg.hab[i] == 1 & d.s$mg.risk[i] == 0 & d.s$mg.rat[i] == 1) scenarios[i] = 7
		if(d.s$mg.hab[i] == 1 & d.s$mg.risk[i] == 1 & d.s$mg.rat[i] == 1) scenarios[i] = 8
		
}

d.s = cbind(scenarios, d.s)

d1 = d.s[d.s$scenarios == 1,]; d2 = d.s[d.s$scenarios == 2,]; d3 = d.s[d.s$scenarios == 3,]; d4 = d.s[d.s$scenarios == 4,]
d5 = d.s[d.s$scenarios == 5,]; d6 = d.s[d.s$scenarios == 6,]; d7 = d.s[d.s$scenarios == 7,]; d8 = d.s[d.s$scenarios == 8,]

plot(d.s$scenarios, d.s[, y], ylim = ylim, xaxt = "n", bty = "n", pch = 20, col = "black", xlab = "Scenarios", ylab = "Community response",
main = main, cex.main = 0.7)
n.scenarios = max(unique(scenarios))
if(add.x.labels == TRUE) {
axis(1, at = 1 : n.scenarios, labels = c("000", "m00", "0r0", "00h", "mr0", "m0h", "0rh", "mrh")) } else {
axis(1, at = 1 : n.scenarios, labels = rep("", n.scenarios)) 
}
points(1, mean(d1[, y]), pch = 4, cex = 2, lwd = 2, col = "blue")
points(2, mean(d2[, y]), pch = 4, cex = 2, lwd = 2, col = "blue")
points(3, mean(d3[, y]), pch = 4, cex = 2, lwd = 2, col = "blue")
points(4, mean(d4[, y]), pch = 4, cex = 2, lwd = 2, col = "blue")
points(5, mean(d5[, y]), pch = 4, cex = 2, lwd = 2, col = "blue")
points(6, mean(d6[, y]), pch = 4, cex = 2, lwd = 2, col = "blue")
points(7, mean(d7[, y]), pch = 4, cex = 2, lwd = 2, col = "blue")
points(8, mean(d8[, y]), pch = 4, cex = 2, lwd = 2, col = "blue")

mean.y = c(mean(d1[,y]), mean(d2[,y]), mean(d3[,y]), mean(d4[,y]), 
mean(d5[,y]), mean(d6[,y]), mean(d7[,y]), mean(d8[,y]))

if(y == "n.w1" & Time == "Present") i = 1
if(y == "n.w2" & Time == "Present") i = 2
if(y == "n.w3" & Time == "Present") i = 3
if(y == "n.w1" & Time == "Future") i = 4
if(y == "n.w2" & Time == "Future") i = 5
if(y == "n.w3" & Time == "Future") i = 6

v = rep(NA, n.scenarios)
v[1] = es$est.int[i]
v[2] = v[1] + es$est.risk[i]
v[3] = v[1] + es$est.rat[i]
v[4] = v[1] + es$est.hab[i]
v[5] = v[1] + es$est.risk[i] + es$est.rat[i]
v[6] = v[1] + es$est.risk[i] + es$est.hab[i]
v[7] = v[1] + es$est.rat[i] + es$est.hab[i]
v[8] = v[1] + es$est.risk[i] + es$est.rat[i] + es$est.hab[i]

points(1, v[1], pch = 3, cex = 2, lwd = 2, col = "red")
points(2, v[2], pch = 3, cex = 2, lwd = 2, col = "red")
points(3, v[3], pch = 3, cex = 2, lwd = 2, col = "red")
points(4, v[4], pch = 3, cex = 2, lwd = 2, col = "red")
points(5, v[5], pch = 3, cex = 2, lwd = 2, col = "red")
points(6, v[6], pch = 3, cex = 2, lwd = 2, col = "red")
points(7, v[7], pch = 3, cex = 2, lwd = 2, col = "red")
points(8, v[8], pch = 3, cex = 2, lwd = 2, col = "red")

list(d.s = d.s, mean.y = mean.y)

} # function

#---

# Appendix S2 Fig. S1 of the community-level paper

layout(matrix(c(1:6), 3, 2)); par(mar = c(2.3,2.3,0,0)+0.5, cex.main = 1)
f.plot.scenarios(d.c = d.c, es = es, Time = "Present", val.risk = 1, val.ac = 2, y = "n.w1", ylim = c(0, 0.6), main = "n.w1 - present", add.x.labels = FALSE)
f.plot.scenarios(d.c = d.c, es = es, Time = "Present", val.risk = 1, val.ac = 2, y = "n.w2", ylim = c(0, 0.6), main = "n.w2 - present", add.x.labels = FALSE)
f.plot.scenarios(d.c = d.c, es = es, Time = "Present", val.risk = 1, val.ac = 2, y = "n.w3", ylim = c(0, 0.6), main = "n.w3 - present", add.x.labels = TRUE)
f.plot.scenarios(d.c = d.c, es = es, Time = "Future", val.risk = 2, val.ac = 2, y = "n.w1", ylim = c(0, 0.6), main = "n.w1 - future", add.x.labels = FALSE)
f.plot.scenarios(d.c = d.c, es = es, Time = "Future", val.risk = 2, val.ac = 2, y = "n.w2", ylim = c(0, 0.6), main = "n.w2 - future", add.x.labels = FALSE)
f.plot.scenarios(d.c = d.c, es = es, Time = "Future", val.risk = 2, val.ac = 2, y = "n.w3", ylim = c(0, 0.6), main = "n.w3 - future", add.x.labels = TRUE)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 22
# Comparative response to management practices as indicated by a principal component analysis

f.pca = function(d.c = d.c, val.risk = 1, val.ac = 2, nf = 5){

d.pca = d.c[d.c$risk == val.risk & d.c$ac == val.ac, ]

d.pca = subset(d.pca, select = c(mg.hab, mg.risk, mg.rat, mg.res.1, mg.res.2, AKIP, AKEP, HCRE, IIWI, ELEP, OMAO, APAP, HAAM) )
#cor(d.pca)

pca = dudi.pca(d.pca, scannf = FALSE, nf = nf)
#pc1 = (pca$l1)$RS1; pc2 = (pca$l1)$RS2; ...

pca 

}

# Conduct PCAs
pca_p = f.pca(d.c = d.c, val.risk = 1, val.ac = 2, nf = 3)
pca_f = f.pca(d.c = d.c, val.risk = 2, val.ac = 2, nf = 3)

# Appendix S3 Fig. S1 of the community-level paper
# note that in R the four sections (A-D) have to be obtained separately
fviz_pca_var(pca_p, axes = c(1, 2)) # section labelled 'A'
fviz_pca_var(pca_p, axes = c(1, 3)) # section labelled 'B'
fviz_pca_var(pca_f, axes = c(1, 2)) # section labelled 'C'
fviz_pca_var(pca_f, axes = c(1, 3)) # section labelled 'D'

# Appendix S3 Table S1 of the community-level paper
t2 = cbind(round(pca_p$co, 2), round(pca_f$co, 2))
t2 = as.data.frame(t2)
colnames(t2) = c("PC1-pres", "PC2-pres", "PC3-pres", "PC1-futur", "PC2-futur", "PC3-futur")
t2
#setwd("C:/programs/MAMO/SAUV")
#write.table(t2, "table.pca.txt", sep = "\t", dec = ".", row.names = T)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 23
# Appendix S3 Fig. S2; note that this version also includes ELEP AND OMAO:
# the algorithm used failed to reliably estimate the parameters of the cline for these two species (clines were likely too shallow)

layout(matrix(c(1:8), 4, 2)); par(mar = c(2.2,2.2,0,0)+0.5, cex.main = 1)

col.cat = c("grey90", "grey70", "grey50", "grey30", "black")
add.cline.cat = c("FALSE", rep("TRUE", 4))

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/akip.10", 
ylab = "# pairs AKIP / ha", fig.title = "AKIP", main = NA, margin.up = -0.9,
wch.plot = "RISK", col.cat = col.cat, add.cline.cat = add.cline.cat,
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = FALSE, create.plot = TRUE, lty.mean = 1
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/akep.10", 
ylab = "# pairs AKEP / ha", fig.title = "AKEP", main = NA, margin.up = 0,
wch.plot = "RISK", col.cat = col.cat, add.cline.cat = add.cline.cat,
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = FALSE, create.plot = TRUE, lty.mean = 1
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/hcre.10", 
ylab = "# pairs HCRE / ha", fig.title = "HCRE", main = NA, margin.up = -0.6,
wch.plot = "RISK", col.cat = col.cat, add.cline.cat = add.cline.cat,
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = FALSE, create.plot = TRUE, lty.mean = 1
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/iiwi.10", 
ylab = "# pairs IIWI / ha", fig.title = "IIWI", main = NA, margin.up = 0,
wch.plot = "RISK", col.cat = col.cat, add.cline.cat = add.cline.cat,
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = FALSE, create.plot = TRUE, lty.mean = 1
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/elep.10", 
ylab = "# pairs ELEP / ha", fig.title = "ELEP", main = NA, margin.up = -0.3,
wch.plot = "RISK", col.cat = col.cat, add.cline.cat = add.cline.cat,
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = FALSE, create.plot = TRUE, lty.mean = 1
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/omao.10", 
ylab = "# pairs OMAO / ha", fig.title = "OMAO", main = NA, margin.up = -0.8,
wch.plot = "RISK", col.cat = col.cat, add.cline.cat = add.cline.cat,
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = FALSE, create.plot = TRUE, lty.mean = 1
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/apap.10", 
ylab = "# pairs APAP / ha", fig.title = "APAP", main = NA, margin.up = -0.8,
wch.plot = "RISK", col.cat = col.cat, add.cline.cat = add.cline.cat,
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = FALSE, create.plot = TRUE, lty.mean = 1
)

f.plot.composite(
output.dir = "C:/Programs/MAMO/RUN/haam.10", 
ylab = "# pairs HAAM / ha", fig.title = "HAAM", main = NA, margin.up = -0.8,
wch.plot = "RISK", col.cat = col.cat, add.cline.cat = add.cline.cat,
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = FALSE, create.plot = TRUE, lty.mean = 1
)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 24
# Appendix S3 Fig. S3; note that this version also includes ELEP AND OMAO:
# the algorithm used failed to reliably estimate the parameters of the cline for these two species (clines were likely too shallow)

layout(matrix(c(1:8), 4, 2)); par(mar = c(2.2,2.2,0,0)+0.5, cex.main = 1)

# ALL SPECIES
var.excl = rep("RISK", 3); val.excl = c(1, 2, 4)
main = NA
wch.plot = "RISK"; col.cat = c("grey50", "black")
add.leg = FALSE

#AKIP
output.dir = "C:/Programs/MAMO/RUN/akip.10"
ylab = "# pairs AKIP / ha"; fig.title = "AKIP"
margin.up = -0.93

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = TRUE, lty.mean = 3)

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = FALSE, lty.mean = 1)

abline(v = 1.6, lty = 4)

#AKEP
output.dir = "C:/Programs/MAMO/RUN/akep.10"
ylab = "# pairs AKEP / ha"; fig.title = "AKEP"
margin.up = -1.0

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = TRUE, lty.mean = 3)

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = FALSE, lty.mean = 1)

abline(v = 1.6, lty = 4)

#HCRE
output.dir = "C:/Programs/MAMO/RUN/hcre.10"
ylab = "# pairs HCRE / ha"; fig.title = "HCRE"
margin.up = -0.6

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = TRUE, lty.mean = 3)

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = FALSE, lty.mean = 1)

abline(v = 1.6, lty = 4)

#IIWI
output.dir = "C:/Programs/MAMO/RUN/iiwi.10"
ylab = "# pairs IIWI / ha"; fig.title = "IIWI"
margin.up = -0.3

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = TRUE, lty.mean = 3)

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = FALSE, lty.mean = 1)

abline(v = 1.6, lty = 4)

#ELEP
output.dir = "C:/Programs/MAMO/RUN/elep.10"
ylab = "# pairs ELEP / ha"; fig.title = "ELEP"
margin.up = -0.3

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = TRUE, lty.mean = 3)

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = FALSE, lty.mean = 1)

abline(v = 1.6, lty = 4)

#OMAO
output.dir = "C:/Programs/MAMO/RUN/omao.10"
ylab = "# pairs OMAO / ha"; fig.title = "OMAO"
margin.up = -1

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = TRUE, lty.mean = 3)

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = FALSE, lty.mean = 1)

abline(v = 1.6, lty = 4)

#APAP
output.dir = "C:/Programs/MAMO/RUN/apap.10"
ylab = "# pairs APAP / ha"; fig.title = "APAP"
margin.up = -0.8

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = TRUE, lty.mean = 3)

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = FALSE, lty.mean = 1)

abline(v = 1.6, lty = 4)

#HCRE
output.dir = "C:/Programs/MAMO/RUN/haam.10"
ylab = "# pairs HAAM / ha"; fig.title = "HAAM"
margin.up = -0.8

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = TRUE, lty.mean = 3)

f.plot.composite(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
add.leg = add.leg, create.plot = FALSE, lty.mean = 1)

abline(v = 1.6, lty = 4)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 25
# Appendix S1 Table S2 of the community-level paper: average parameter value in the 'sensitivity envelop'

sp = c("AKIP", "AKEP", "HCRE", "IIWI", "ELEP", "OMAO", "APAP", "HAAM")
e = list(); e[[1]] = envp.akip; e[[2]] = envp.akep; e[[3]] = envp.hcre; e[[4]] = envp.iiwi
e[[5]] = envp.elep; e[[6]] = envp.omao; e[[7]] = envp.apap; e[[8]] = envp.haam

m_s.ad = m_fec = m_fec.1 = m_rat.f = m_rat.s = m_K.b = m_K.nb.1 = 
m_gamma.mov = m_Sm.ac = m_R.ter = m_fidelity.ad = m_m.natal = m_psi.DD = rep(NA, 8)

for(i in 1 : 8) {
	m_s.ad[i] = mean(e[[i]]$s.ad)
	m_fec[i] = mean(e[[i]]$fec)
	m_fec.1[i] = mean(e[[i]]$fec.1)
   m_rat.f[i] = mean(e[[i]]$rat.f)
   m_rat.s[i] = mean(e[[i]]$rat.s)
   m_K.b[i] = mean(e[[i]]$K.b)
   m_K.nb.1[i] = mean(e[[i]]$K.nb.1)
   m_gamma.mov[i] = mean(e[[i]]$gamma.mov)
   m_Sm.ac[i] = mean(e[[i]]$Sm.ac)
   m_R.ter[i] = mean(e[[i]]$R.ter)
   m_fidelity.ad[i] = mean(e[[i]]$fidelity.ad)   
   m_m.natal[i] = mean(e[[i]]$m.natal)
   m_psi.DD[i] = mean(e[[i]]$psi.DD)
}

 app = data.frame(cbind(sp, m_s.ad, m_fec, m_fec.1, m_rat.f, m_rat.s, m_K.b, m_K.nb.1,
 m_gamma.mov, m_Sm.ac, m_R.ter, m_fidelity.ad, m_m.natal, m_psi.DD))
 
 for(i in 1 : 8) {
	if(app$m_gamma.mov[i] == -10) app$m_gamma.mov[i] = NA
	if(app$m_K.nb.1[i] == 0) app$m_K.nb.1[i] = NA
}

#setwd("C:/programs/MAMO/SAUV")
#write.table(app, "appendix.A4.txt", sep = "\t", dec = ".", row.names = F)
 
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 26
# Fig. 3  of the community-level paper

# f.plot.composite modified: cex.main = 1.2, xlab = "Elevation (km)"

f.plot.composite.2 = function(
 output.dir = "C:/Programs/MAMO/RUN/IIWI.1_3", 
 #var.excl = c("alpha.b.low", "alpha.b.low"), val.excl = c(0.00449374, 0.00898748),
 var.excl = NA, val.excl = NA,
 ylab = "# pairs Iiwi / ha of native forest", fig.title = NA, main = "", margin.up = 0.5,
 wch.plot = "RISK", col.cat = c("grey90", "grey70", "pink", "grey40", "black"), add.cline.cat = c("FALSE", rep("TRUE", 4)),
 val.RISK = NA, val.AC = 3, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
 add.leg = FALSE, leg.x = 1, leg.y = 4.5, leg.text = c("Past", "Pres / 2", "Pres-obs", "Pres-sim", "2100/2", "2100"), 
 leg.pch = c(15,15,3,15,15,15), leg.lty = c(2,2,1,2,2,2), leg.col = c("grey90", "grey70", "red", "pink", "grey40", "black"),
 create.plot = TRUE, lty.mean = 2
){
 
 # Data set
 
 setwd(output.dir)
 d.sim = read.table("t.sim.txt", header = T, sep = "\t", dec = ".")
 d.factors = read.table("factors.txt", header = T, sep = "\t", dec = ".")
 
 d = d.sim
 d = d[!is.na(d[,1]),]
 n = dim(d)[1]
 
 if(dim(d.factors)[1] != n) stop("d.sim and d.factors have different lengths")
 
 n.sim = max(d.factors$SIM)
 
 d.f = cbind(d, d.factors)
 
 if(!is.na(val.excl[1])) for(i in 1 : length(val.excl)) d.f = d.f[ d.f[, var.excl[i] ] != val.excl[i], ]
 
 #----------
 
 # Data subset
 
 if(is.na(val.RISK)) d.subset = d.f[d.f$AC == val.AC & d.f$RAT.S == val.RAT & d.f$RES.1 == val.RES.1 & d.f$RES.2 == val.RES.2,]
 if(is.na(val.AC)) d.subset = d.f[d.f$RISK == val.RISK & d.f$RAT.S == val.RAT & d.f$RES.1 == val.RES.1 & d.f$RES.2 == val.RES.2,] 
 if(is.na(val.RAT)) d.subset = d.f[d.f$RISK == val.RISK & d.f$AC == val.AC & d.f$RES.1 == val.RES.1 & d.f$RES.2 == val.RES.2,] 
 if(is.na(val.RES.1)) d.subset = d.f[d.f$RISK == val.RISK & d.f$AC == val.AC & d.f$RAT.S == val.RAT & d.f$RES.2 == val.RES.2,] 
 if(is.na(val.RES.2)) d.subset = d.f[d.f$RISK == val.RISK & d.f$AC == val.AC & d.f$RAT.S == val.RAT & d.f$RES.1 == val.RES.1,] 
 
 #----------
 
 # Figure title
 
 if(!is.na(fig.title)) main = fig.title else {
   
   if(!is.na(main) & is.na(val.RISK)) main = paste(main, "|",  "AC", "=", val.AC, ",", "RAT", "=", val.RAT, ",", "RES.1", "=", val.RES.1, ",", "RES.2", "=", val.RES.2)
   if(!is.na(main) & is.na(val.AC)) main = paste(main, "|", "RISK", "=", val.RISK, ",", "RAT", "=", val.RAT, ",", "RES.1", "=", val.RES.1, ",", "RES.2", "=", val.RES.2)
   if(!is.na(main) & is.na(val.RAT)) main = paste(main, "|", "RISK", "=", val.RISK, ",", "AC", "=", val.AC, ",",  "RES.1", "=", val.RES.1, ",","RES.2", "=", val.RES.2)
   if(!is.na(main) & is.na(val.RES.1)) main = paste(main, "|", "RISK", "=", val.RISK, ",", "AC", "=", val.AC, ",", "RAT", "=", val.RAT, ",", "RES.2", "=", val.RES.2)
   if(!is.na(main) & is.na(val.RES.2)) main = paste(main, "|", "RISK", "=", val.RISK, ",", "AC", "=", val.AC, ",", "RAT", "=", val.RAT, ",", "RES.1", "=", val.RES.1)
   
 }
 
 #----------
 
 x = seq(d.subset$grad.max[1], d.subset$grad.min[1], length.out = d.subset$nr[1])
 
 unit = d.subset$unit[1]
 
 vec.np = c(d.subset$np.high, d.subset$np.mid, d.subset$np.low) / ( unit^2 * 100 )
 max.y = ceiling( max( vec.np ) ) + margin.up
 
 if(create.plot == TRUE) plot(x, rep(0, length(x)), xlim = c(min(x), max(x)), ylim = c(0, max.y), type = "n", bty = "n", xlab = "Elevation (km)", ylab = ylab, main = main, cex.main = 1.2)
 
 #----------
 
 add.s.batch = function(
   d.subset = d.subset, wch.plot = wch.plot, wch.cat = 1,
   x = x, 
   col = "pink", pch = 15, cex = 0.7, lty = lty.mean, add.mean = TRUE, add.cline = TRUE) {
   
   d.subset = d.subset[ d.subset[, wch.plot] == wch.cat, ]
   
   r = unique(d.subset$batch); nr = length(r)
   
   if(add.mean == TRUE) {
     
     # identifies all data points (y-axis), plots them, and calculates elevation-specific average
     
     get.p = function(d = "r1.rdata") { load(d); y.mamo = x$np / ( unit^2 * 100 ); y.mamo }
     
     dp = rbind()
     
     for(i in 1 : nr) {
       
       p_ = list()
       for(j in 1 : n.sim) {
         p_[[j]] = get.p( d = paste("s", r[i], ".", j, ".Rdata", sep = "") )
         #points(x, p_[[j]], col = col, pch = pch, cex = cex)	
         dp = rbind(dp, p_[[j]])
       }# j
       
     } # i
     
     dp = as.data.frame(dp); nc = dim(dp)[2]; mp = rep(NA, nc)
     for(j in 1 : nc) mp[j] = mean(dp[,j])
     
     # plots a fitting cline or line for elevation-specific average
     if(add.cline == TRUE) wrap.cline(x = x, y = mp, col = col, lty = lty, lwd = 3) else 
       segments(min(x), mean(mp), max(x), mean(mp), col = col, lty = lty, lwd = 3)
     
   } else {
     
     for(i in 1 : nr) {
       for(j in 1 : n.sim) {
         add.s(d = paste("s", r[i], ".", j, ".Rdata", sep = ""), col.p = col, pch = 15, cex = 0.7, col.l = col, lty = lty, add.cline = add.cline)
       }
     }
     
   } # if(add.mean == TRUE)
   
 } # function
 
 #----------
 
 n.cat = length( unique( d.subset[, wch.plot] ) )
 cat.i = sort( unique( d.subset[, wch.plot] ) )
 
 for(i in 1 : n.cat) {
   
   add.s.batch(
     d.subset = d.subset, wch.plot = wch.plot, wch.cat = cat.i[i], x = x, col = col.cat[i], lty = lty.mean, 
     add.mean = TRUE, add.cline = add.cline.cat[i]
   )
   
 }
 
 if(add.leg == TRUE) legend(leg.x, leg.y, leg.text, pch = leg.pch, lty = leg.lty, cex = 0.7, col = leg.col, text.col = "black", bg = "white")
 
 #----------
 
 list( d.subset = d.subset, test.mamo = max(d$test.disp.breed) )
 
} # end of f.plot.composite.2

#---

library(compiler)
f.plot.composite.2 = cmpfun(f.plot.composite.2)

#---

f.plot.single.sp = function(
 output.dir = "C:/Programs/MAMO/RUN/akep.10",
 ylab = "# pairs AKEPA / ha", fig.title = "",
 margin.up = 0){
 
 var.excl = rep("RISK", 3); val.excl = c(1, 3, 5)
 main = NA
 wch.plot = "RISK"; col.cat =  c("blue", "orange")
 add.leg = FALSE
 
 f.plot.composite.2(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
                    ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
                    wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
                    val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
                    add.leg = add.leg, create.plot = TRUE, lty.mean = 3)
 
 val.excl = c(1, 2, 4)
 
 f.plot.composite.2(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
                    ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
                    wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
                    val.RISK = NA, val.AC = 2, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
                    add.leg = add.leg, create.plot = FALSE, lty.mean = 1)
 
 f.plot.composite.2(output.dir = output.dir, var.excl = var.excl, val.excl = val.excl, 
                    ylab = ylab, fig.title = fig.title, main = main, margin.up = margin.up,
                    wch.plot = wch.plot, col.cat = col.cat, add.cline.cat = rep("TRUE", 2),
                    val.RISK = NA, val.AC = 2, val.RAT = 2, val.RES.1 = 1, val.RES.2 = 1,
                    add.leg = add.leg, create.plot = FALSE, lty.mean = 2)
 
}

#---

#layout(matrix(c(1:2), 1, 2)); par(mar = c(2.2,2.2,0,0)+0.5)
layout(matrix(c(1:2), 1, 2)); par(mar = c(3.5,3.5,0,0)+0.5, cex.lab = 1.5)

f.plot.single.sp(
 output.dir = "C:/Programs/MAMO/RUN/akip.10",
 ylab = "# pairs / ha", fig.title = "AKIP",
 margin.up = -0.92)
legend(1, 0.06, c("2015", "2100"), bty = "n", cex = 1.2, lwd = 2, col = c("blue", "orange"), text.col = "black", lty = c(1,1), merge = TRUE)
legend(1, 0.08, c("No action (000)", "Lower Malaria Risk (m00)", "Lower Rat Predation (0r0)"), bty = "n", cex = 1.2, lwd = 2, col = "black", text.col = "black", lty = c(1,3,2), merge = TRUE)

f.plot.single.sp(
 output.dir = "C:/Programs/MAMO/RUN/akep.10",
 ylab = "", fig.title = "AKEP",
 margin.up = 0)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# BLOCK 27 

# A series of add-on / patch code used to create tables and Figs for the community-level paper

#---

# = modified version of BLOCK 20 to create Appendix S2 Table S2 of the community-level paper

# Effect size

f.effect.size.2 = function(){
  
  # f.delta functions based on info provided in http://www.ats.ucla.edu/stat/r/faq/deltamethod.htm and http://www.inside-r.org/packages/cran/msm/docs/deltamethod
  # https://stats.idre.ucla.edu/r/faq/how-can-i-estimate-the-standard-error-of-transformed-regression-parameters-in-r-using-the-delta-method/
  # I checked that vcov(m) works for fixed-effect models (returns matrix of variance - covariance of fixed effects) by comparing the output of summary(m), returning Correlation of Fixed Effects, and vcov(m), calculating correlation as cov(i,j)/ (sqrt(var(i))*sqrt(var(j))
  
  f.delta.2 = function(m = m1.p, b = 2){ # b: risk = 2, hab = 3, rat = 4
    a = 1 # intercept
    cf = c(fixef(m)[a], fixef(m)[b])
    #cf
    k = matrix(NA, 2, 2)
    v = vcov(m)
    k[1,1] = v[a,a]
    k[1,2] = k[2,1] = v[a,b]
    k[2,2] = v[b,b]
    deltamethod(~ (x1+x2) / x1, cf, k)
  }
  
  f.delta.3 = function(m = m1.p, b = 2, p = 3){ # b: risk = 2, hab = 3, rat = 4
    a = 1 # intercept
    cf = c(fixef(m)[a], fixef(m)[b], fixef(m)[p])
    #cf
    k = matrix(NA, 3, 3)
    v = vcov(m)
    k[1,1] = v[a,a]
    k[2,2] = v[b,b]
    k[3,3] = v[p,p]
    k[1,2] = k[2,1] = v[a,b]
    k[1,3] = k[3,1] = v[a,p]
    k[2,3] = k[3,2] = v[b,p]
    deltamethod(~ (x1+x2+x3) / x1, cf, k)
  }
  
  f.delta.4 = function(m = m1.p) deltamethod(~ (x1+x2+x3+x4) / x1, fixef(m), vcov(m))
  
  #---
  
  m = list(m1.p, m2.p, m3.p, m1.f, m2.f, m3.f)
  
  fixef_ = list(); fixef_[[1]] = fixef(m1.p); fixef_[[2]] = fixef(m2.p); fixef_[[3]] = fixef(m3.p)
  fixef_[[4]] = fixef(m1.f); fixef_[[5]] = fixef(m2.f); fixef_[[6]] = fixef(m3.f)
  
  se_ = list(); se_[[1]] = summary(m1.p)$coef[,2]; se_[[2]] = summary(m2.p)$coef[,2]
  se_[[3]] = summary(m3.p)$coef[,2]; se_[[4]] = summary(m1.f)$coef[,2]
  se_[[5]] = summary(m2.f)$coef[,2]; se_[[6]] = summary(m3.f)$coef[,2]
  
  #---
  
  es = as.data.frame(matrix(NA, nr = 6, nc = 30))
  colnames(es) = c(
    "Time", "response", 
    "eff.m00", "eff.0r0", "eff.00h", "eff.mr0", "eff.m0h", "eff.0rh", "eff.mrh", 
    "se.eff.m00", "se.eff.0r0", "se.eff.00h", "se.eff.mr0", "se.eff.m0h", "se.eff.0rh", "se.eff.mrh", 
    "CI1.eff.m00", "CI2.eff.m00",
    "CI1.eff.0r0", "CI2.eff.0r0",
    "CI1.eff.00h", "CI2.eff.00h",
    "CI1.eff.mr0", "CI2.eff.mr0",
    "CI1.eff.m0h", "CI2.eff.m0h",
    "CI1.eff.0rh", "CI2.eff.0rh",
    "CI1.eff.mrh", "CI2.eff.mrh"
  )
  
  es$Time = c(rep("Present", 3), rep("Future", 3))
  
  es$response = rep(c("n.w1", "n.w2", "n.w3"), 2)
  
  for(i in 1: 6) {
    
    es$eff.m00[i] = ( fixef_[[i]][2] + fixef_[[i]][1] ) / fixef_[[i]][1]
    es$eff.0r0[i] = ( fixef_[[i]][4] + fixef_[[i]][1] ) / fixef_[[i]][1]
    es$eff.00h[i] = ( fixef_[[i]][3] + fixef_[[i]][1] ) / fixef_[[i]][1]
    es$eff.mr0[i] = ( fixef_[[i]][4] + fixef_[[i]][2] + fixef_[[i]][1] ) / fixef_[[i]][1]
    es$eff.m0h[i] = ( fixef_[[i]][3] + fixef_[[i]][2] + fixef_[[i]][1] ) / fixef_[[i]][1]	
    es$eff.0rh[i] = ( fixef_[[i]][4] + fixef_[[i]][3] + fixef_[[i]][1] ) / fixef_[[i]][1]		
    es$eff.mrh[i] = ( fixef_[[i]][4] + fixef_[[i]][3] + fixef_[[i]][2] + fixef_[[i]][1] ) / fixef_[[i]][1]			
    
    es$se.eff.m00[i] = f.delta.2(m[[i]], 2)
    es$se.eff.0r0[i] = f.delta.2(m[[i]], 4)
    es$se.eff.00h[i] = f.delta.2(m[[i]], 3)	
    es$se.eff.mr0[i] = f.delta.3(m[[i]], 2, 4)
    es$se.eff.m0h[i] = f.delta.3(m[[i]], 2, 3)
    es$se.eff.0rh[i] = f.delta.3(m[[i]], 3, 4)
    es$se.eff.mrh[i] = f.delta.4(m[[i]])
    
  }
  
  # debugged entirely up to here. Calculation and plot of confidence intervals debugged using the ch function.
  
  for(i in 1: 6) {
    
    es$CI1.eff.m00[i] = es$eff.m00[i] - 2 * es$se.eff.m00[i]; es$CI2.eff.m00[i] = es$eff.m00[i] + 2 * es$se.eff.m00[i]
    es$CI1.eff.0r0[i] = es$eff.0r0[i] - 2 * es$se.eff.0r0[i]; es$CI2.eff.0r0[i] = es$eff.0r0[i] + 2 * es$se.eff.0r0[i]
    es$CI1.eff.00h[i] = es$eff.00h[i] - 2 * es$se.eff.00h[i]; es$CI2.eff.00h[i] = es$eff.00h[i] + 2 * es$se.eff.00h[i]
    es$CI1.eff.mr0[i] = es$eff.mr0[i] - 2 * es$se.eff.mr0[i]; es$CI2.eff.mr0[i] = es$eff.mr0[i] + 2 * es$se.eff.mr0[i]
    es$CI1.eff.m0h[i] = es$eff.m0h[i] - 2 * es$se.eff.m0h[i]; es$CI2.eff.m0h[i] = es$eff.m0h[i] + 2 * es$se.eff.m0h[i]
    es$CI1.eff.0rh[i] = es$eff.0rh[i] - 2 * es$se.eff.0rh[i]; es$CI2.eff.0rh[i] = es$eff.0rh[i] + 2 * es$se.eff.0rh[i]
    es$CI1.eff.mrh[i] = es$eff.mrh[i] - 2 * es$se.eff.mrh[i]; es$CI2.eff.mrh[i] = es$eff.mrh[i] + 2 * es$se.eff.mrh[i]
    
  }
  
  #nd = 3
  #for(i in 3: 30) es[,i] = round(es[,i], nd)
  
  es
  
}

#---

# Appendix S2 Table S2 of the community-level paper
es = f.effect.size.2()
es

#setwd("C:/programs/MAMO/SAUV")
#write.table(es, "effect.size.txt", sep = "\t", dec = ".", row.names = F)

#---

# Plot and compare alternative scenarios

# DEBUG
Time = "Present"
y = "n.w1"
ylim = c(0.5, 3.5)
n.scenarios = 7
labels.scenarios = c("m00", "0r0", "00h", "mr0", "m0h", "0rh", "mrh")
main = "w1 - present"
add.x.labels = TRUE
create.plot = TRUE
col.plot = "black"
d.x = 0

#---

f.plot.effect.size = function(es = es, Time = "Present",
                              y = "n.w1", ylim = c(0.5, 3.5), n.scenarios = 7, labels.scenarios = c("m00", "0r0", "00h", "mr0", "m0h", "0rh", "mrh"), main = "w1 - present", add.x.labels = TRUE, create.plot = TRUE, col.plot = "black", lty.plot = 1, d.x = 0){
  
  # Figure
  
  if(create.plot == TRUE) plot(1:n.scenarios, seq(ylim[1], ylim[2], length.out = n.scenarios), type = "n", ylim = ylim, xaxt = "n", bty = "n", pch = 20, cex = 0.3, col = "black", xlab = "Scenarios", ylab = "Effect Size",
                               main = main, cex.main = 1)
  
  if(add.x.labels == TRUE) {
    axis(1, at = 1 : n.scenarios, labels = labels.scenarios, cex.axis = 1) } else {
      axis(1, at = 1 : n.scenarios, labels = rep("", n.scenarios)) 
    }
  
  abline(h = 1, col = "black", lwd = 2, lty = 2)
  #
  
  if(y == "n.w1" & Time == "Present") i = 1
  if(y == "n.w2" & Time == "Present") i = 2
  if(y == "n.w3" & Time == "Present") i = 3
  if(y == "n.w1" & Time == "Future") i = 4
  if(y == "n.w2" & Time == "Future") i = 5
  if(y == "n.w3" & Time == "Future") i = 6
  
  points(1+d.x, es$eff.m00[i], pch = 19, cex = 1, lwd = 2, col = col.plot)
  points(2+d.x, es$eff.0r0[i], pch = 19, cex = 1, lwd = 2, col = col.plot)
  points(3+d.x, es$eff.00h[i], pch = 19, cex = 1, lwd = 2, col = col.plot)
  points(4+d.x, es$eff.mr0[i], pch = 19, cex = 1, lwd = 2, col = col.plot)
  points(5+d.x, es$eff.m0h[i], pch = 19, cex = 1, lwd = 2, col = col.plot)
  points(6+d.x, es$eff.0rh[i], pch = 19, cex = 1, lwd = 2, col = col.plot)
  points(7+d.x, es$eff.mrh[i], pch = 19, cex = 1, lwd = 2, col = col.plot)
  
  #
  q1 = c(es$CI1.eff.m00[i], es$CI2.eff.m00[i])
  q2 = c(es$CI1.eff.0r0[i], es$CI2.eff.0r0[i])
  q3 = c(es$CI1.eff.00h[i], es$CI2.eff.00h[i])
  q4 = c(es$CI1.eff.mr0[i], es$CI2.eff.mr0[i])
  q5 = c(es$CI1.eff.m0h[i], es$CI2.eff.m0h[i])
  q6 = c(es$CI1.eff.0rh[i], es$CI2.eff.0rh[i])
  q7 = c(es$CI1.eff.mrh[i], es$CI2.eff.mrh[i])
  
  z = 1+d.x
  w.b = 0.05
  
  segments(z, q1[1], z, q1[2], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q1[1], z+w.b, q1[1], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q1[2], z+w.b, q1[2], lwd = 2, col = col.plot, lty = lty.plot)
  z = z+1
  
  segments(z, q2[1], z, q2[2], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q2[1], z+w.b, q2[1], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q2[2], z+w.b, q2[2], lwd = 2, col = col.plot, lty = lty.plot)
  z = z+1
  
  segments(z, q3[1], z, q3[2], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q3[1], z+w.b, q3[1], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q3[2], z+w.b, q3[2], lwd = 2, col = col.plot, lty = lty.plot)
  z = z+1
  
  segments(z, q4[1], z, q4[2], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q4[1], z+w.b, q4[1], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q4[2], z+w.b, q4[2], lwd = 2, col = col.plot, lty = lty.plot)
  z = z+1
  
  segments(z, q5[1], z, q5[2], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q5[1], z+w.b, q5[1], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q5[2], z+w.b, q5[2], lwd = 2, col = col.plot, lty = lty.plot)
  z = z+1
  
  segments(z, q6[1], z, q6[2], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q6[1], z+w.b, q6[1], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q6[2], z+w.b, q6[2], lwd = 2, col = col.plot, lty = lty.plot)
  z = z+1
  
  segments(z, q7[1], z, q7[2], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q7[1], z+w.b, q7[1], lwd = 2, col = col.plot, lty = lty.plot)
  segments(z-w.b, q7[2], z+w.b, q7[2], lwd = 2, col = col.plot, lty = lty.plot)
  
} # function

#---

# Plot not used for publication

layout(matrix(c(1:2), 1, 2)); par(mar = c(2.3,2.3,0,0)+0.5, cex.main = 1)
ylim.1 = 0.95
abs.d.x = 0.08
f.plot.effect.size(es = es, Time = "Present", y = "n.w1", ylim = c(ylim.1, 3), main="PRESENT", add.x.labels = TRUE, create.plot = TRUE, col.plot = "black", d.x = -abs.d.x)
f.plot.effect.size(es = es, Time = "Present", y = "n.w3", ylim = c(ylim.1, 3), main="PRESENT", add.x.labels = FALSE, create.plot = FALSE, col.plot = "grey", d.x = abs.d.x)
legend(1, 2.5, bty = "n", cex = 1, c("Equal weights", "Weight = 1 / N"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "grey")) 
f.plot.effect.size(es = es, Time = "Future", y = "n.w1", ylim = c(ylim.1, 3), main="FUTURE", add.x.labels = TRUE, create.plot = TRUE, col.plot = "black", d.x = -abs.d.x)
f.plot.effect.size(es = es, Time = "Future", y = "n.w3", ylim = c(ylim.1, 3), main="FUTURE", add.x.labels = FALSE, create.plot = FALSE, col.plot = "grey", d.x = abs.d.x)

#---

ff = function(){
  ylim.1 = 0.95
  abs.d.x = 0.15
  f.plot.effect.size(es = es, Time = "Present", y = "n.w1", ylim = c(ylim.1, 3), main=NA, add.x.labels = TRUE, create.plot = TRUE, col.plot = "dodgerblue", lty.plot = 1, d.x = -abs.d.x)
  f.plot.effect.size(es = es, Time = "Present", y = "n.w3", ylim = c(ylim.1, 3), main=NA, add.x.labels = FALSE, create.plot = FALSE, col.plot = "dodgerblue3", lty.plot = 1, d.x = abs.d.x)
  f.plot.effect.size(es = es, Time = "Future", y = "n.w1", ylim = c(ylim.1, 3), main=NA, add.x.labels = FALSE, create.plot = FALSE, col.plot = "darkorange", lty.plot = 1, d.x = -abs.d.x)
  f.plot.effect.size(es = es, Time = "Future", y = "n.w3", ylim = c(ylim.1, 3), main=NA, add.x.labels = FALSE, create.plot = FALSE, col.plot = "darkorange3", lty.plot = 1, d.x = abs.d.x)
  abline(v = 1.5, col = "black", lwd = 1, lty = 1)
  abline(v = 2.5, col = "black", lwd = 1, lty = 1)
  abline(v = 3.5, col = "black", lwd = 1, lty = 1)
  abline(v = 4.5, col = "black", lwd = 1, lty = 1)
  abline(v = 5.5, col = "black", lwd = 1, lty = 1)
  abline(v = 6.5, col = "black", lwd = 1, lty =1)
  abline(v = 1, col = "black", lwd = 1, lty = 3)
  abline(v = 2, col = "black", lwd = 1, lty = 3)
  abline(v = 3, col = "black", lwd = 1, lty = 3)
  abline(v = 4, col = "black", lwd = 1, lty = 3)
  abline(v = 5, col = "black", lwd = 1, lty = 3)
  abline(v = 6, col = "black", lwd = 1, lty =3)
  abline(v = 7, col = "black", lwd = 1, lty =3)
  legend(1, 3, c("2015 - W=1", "2015 - W=1/N", "2100 - W=1", "2100 - W=1/N"), cex = 0.7, lty = rep(1, 4), lwd = rep(2, 4), col = c("dodgerblue", "dodgerblue3", "darkorange", "darkorange3"), bg = "white") 
  
}

par(mfrow = c(1,1))
ff()

#---

# Check of calculation and plot of CI, knowing that es[,1:16] is fully debugged

ch = function(i = 1, cl = "m00", col.ch = "red"){
  
  ff()
  
  if(cl == "m00") { m = "eff.m00" ; s = "se.eff.m00" }
  if(cl == "0r0") { m = "eff.0r0" ; s = "se.eff.0r0" }
  if(cl == "00h") { m = "eff.00h" ; s = "se.eff.00h" }
  if(cl == "mr0") { m = "eff.mr0" ; s = "se.eff.mr0" }
  if(cl == "m0h") { m = "eff.m0h" ; s = "se.eff.m0h" }
  if(cl == "0rh") { m = "eff.0rh" ; s = "se.eff.0rh" }
  if(cl == "mrh") { m = "eff.mrh" ; s = "se.eff.mrh" }
  
  abline(h = es[i, m], lty = 1, col = col.ch)
  abline(h = es[i, m] - 2 * es[i, s], lty = 2, col = col.ch)
  abline(h = es[i, m] + 2 * es[i, s], lty = 2, col = col.ch)
  
}

i = 1 #3, 4, 6
ch(i, "m00")
ch(i, "0r0")
ch(i, "00h")
ch(i, "mr0")
ch(i, "m0h")
ch(i, "0rh")
ch(i, "mrh")

#---

ff()

#---

es.pub = cbind(es[, 1:2], round(es[, 3:9], 2), round(es[, 17:30], 2))

#

x = es.pub[, 1:9]
ci = es.pub[, 10:23]
colnames(x) = c("Time", "Y", "m00",  "0r0", "00h", "mr0", "m0h", "0rh", "mrh")

for(i in 1:6) {
  x[i, "m00"] = paste(x[i, "m00"], " [", ci[i, "CI1.eff.m00"], ",", ci[i, "CI2.eff.m00"], "]")
  x[i, "0r0"] = paste(x[i, "0r0"], " [", ci[i, "CI1.eff.0r0"], ",", ci[i, "CI2.eff.0r0"], "]")
  x[i, "00h"] = paste(x[i, "00h"], " [", ci[i, "CI1.eff.00h"], ",", ci[i, "CI2.eff.00h"], "]")
  x[i, "mr0"] = paste(x[i, "mr0"], " [", ci[i, "CI1.eff.mr0"], ",", ci[i, "CI2.eff.mr0"], "]")
  x[i, "m0h"] = paste(x[i, "m0h"], " [", ci[i, "CI1.eff.m0h"], ",", ci[i, "CI2.eff.m0h"], "]")
  x[i, "0rh"] = paste(x[i, "0rh"], " [", ci[i, "CI1.eff.0rh"], ",", ci[i, "CI2.eff.0rh"], "]")
  x[i, "mrh"] = paste(x[i, "mrh"], " [", ci[i, "CI1.eff.mrh"], ",", ci[i, "CI2.eff.mrh"], "]")
}

#setwd("C:/Alban/Alban HAWAII/HAKALAU/paper")
#write.table(x, "Table B3.txt", sep = "\t", dec = ".", row.names = F)

# debugged visually x versus es.pub 7/20/2016

#---

# Figure 4 of community-level paper

fff = function( col.f = c("dodgerblue", "dodgerblue3", "darkorange", "darkorange3") ){
  ylim.1 = 0.95
  abs.d.x = 0.15
  f.plot.effect.size(es = es, Time = "Present", y = "n.w1", ylim = c(ylim.1, 3), main=NA, add.x.labels = TRUE, create.plot = TRUE, col.plot = col.f[1], lty.plot = 1, d.x = -abs.d.x)
  f.plot.effect.size(es = es, Time = "Present", y = "n.w3", ylim = c(ylim.1, 3), main=NA, add.x.labels = FALSE, create.plot = FALSE, col.plot = col.f[2], lty.plot = 1, d.x = abs.d.x)
  f.plot.effect.size(es = es, Time = "Future", y = "n.w1", ylim = c(ylim.1, 3), main=NA, add.x.labels = FALSE, create.plot = FALSE, col.plot = col.f[3], lty.plot = 1, d.x = -abs.d.x)
  f.plot.effect.size(es = es, Time = "Future", y = "n.w3", ylim = c(ylim.1, 3), main=NA, add.x.labels = FALSE, create.plot = FALSE, col.plot = col.f[4], lty.plot = 1, d.x = abs.d.x)
  abline(v = 1.5, col = "black", lwd = 1, lty = 1)
  abline(v = 2.5, col = "black", lwd = 1, lty = 1)
  abline(v = 3.5, col = "black", lwd = 1, lty = 1)
  abline(v = 4.5, col = "black", lwd = 1, lty = 1)
  abline(v = 5.5, col = "black", lwd = 1, lty = 1)
  abline(v = 6.5, col = "black", lwd = 1, lty =1)
  abline(v = 1, col = "black", lwd = 1, lty = 3)
  abline(v = 2, col = "black", lwd = 1, lty = 3)
  abline(v = 3, col = "black", lwd = 1, lty = 3)
  abline(v = 4, col = "black", lwd = 1, lty = 3)
  abline(v = 5, col = "black", lwd = 1, lty = 3)
  abline(v = 6, col = "black", lwd = 1, lty =3)
  abline(v = 7, col = "black", lwd = 1, lty =3)
  legend(1, 3, c("2015 - W=1", "2015 - W=1/N", "2100 - W=1", "2100 - W=1/N"), cex = 0.7, lty = rep(1, 4), lwd = rep(2, 4), col = col.f, bg = "white") 
  
}

par(mfrow = c(1,1), mar = c(4.3,4.3,0,0)+0.5)
fff( col.f = c("blue", "purple", "orange", "red") )

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
