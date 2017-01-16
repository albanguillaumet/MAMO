#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO CALIBRATE MAMO

f.calibr = function( 

# species
sp = "IIWI",
output.dir = "C:/Alban/Alban HAWAII/HAKALAU/R code/CALIBRATION/iiwi-dum",
y.obs = c(5.39163173, 4.84949943, 5.11686774, 4.70852507, 3.66411526, 1.75865820, 0.24413093, 0.00000000, 0.00887159, 0.00000000), # Fit to expected distribution # y.obs = # pairs / ha
d = read.table("C:/Alban/Alban HAWAII/HAKALAU/Data/Model calibration/Parameters/table.IIWI.1.txt", header = T, sep = "\t", dec = "."), # Parameter range/values

# Spatial structure
nr = 10, nc = 2, grad = c(1900, 1000), unit = 1,

# Time frame
T = 60, Tm = 5, SD.fledg = 0,

# Survival (ad = annual, juv = from fledging to breeding age)
n.s.ad = 3, input.direct.s.ad = FALSE, s.ad.direct = NA,
n.rat.s = 1,

# Reproduction and habitat quality
n.fec = 3, input.direct.fec = FALSE, fec.direct = NA,
n.fec.1 = 3, input.direct.fec.1 = FALSE, fec.1.direct = NA,
paired.s.ad.fec = FALSE,
n.rat.f = 1,
n.K.b = 3,
K.nb.1 = list(K.nb.1.2003, K.nb.1.2004, K.nb.1.avg),
n.K.nb.1 = 3,
K.nb.2 = K.nb.2,
reproduction.malaria = "simple",

# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, alpha.1 = alpha.1,
n.Sm.ac = 1, input.direct.Sm.ac = FALSE, Sm.ac.direct = NA,

# Movements
n.gamma.mov = 2, input.direct.gamma.mov = FALSE, gamma.mov.direct = NA,
calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
n.R.ter = 1,
n.fidelity.ad = 2,
n.m.natal = 2,
n.psi.DD = 2,

# Other options 
add.cline = FALSE,

# Simulations 
n.sim = 1, design = "simple", batch = 1

) {

f.s = function(x = s, n = batch.n) {
name = paste("s", n, sep ="")
save(x, file = paste(name, ".rdata", sep = ""))
}
	
#------------------

setwd(output.dir)

if(dim(d)[1] == 1) rownames(d) = sp else rownames(d) = d$species

# Time frame
t.b = d[sp, "t.b"]
f.nb.1 = d[sp, "f.nb.1"]
min.fledg = d[sp, "min.fledg"]
peak.fledg = d[sp, "peak.fledg"]

# Survival 

if(input.direct.s.ad == FALSE) s.ad_ = seq(d[sp, "s.ad_min"], d[sp, "s.ad_max"], length.out = n.s.ad)
if(input.direct.s.ad == TRUE) s.ad_ = s.ad.direct

rat.s_ = seq(d[sp, "rat.s_min"], d[sp, "rat.s_max"], length.out = n.rat.s)

# Reproduction and habitat quality
if(input.direct.fec == FALSE) fec_ = seq(d[sp, "fec_min"], d[sp, "fec_max"], length.out = n.fec)
if(input.direct.fec == TRUE) fec_ = fec.direct

if(input.direct.fec.1 == FALSE) fec.1_ = seq(d[sp, "fec.1_min"], d[sp, "fec.1_max"], length.out = n.fec.1)
if(input.direct.fec.1 == TRUE) fec.1_ = fec.1.direct

rat.f_ = seq(d[sp, "rat.f_min"], d[sp, "rat.f_max"], length.out = n.rat.f)
K.b_ = seq(d[sp, "K.b_min"], d[sp, "K.b_max"], length.out = n.K.b)
if(sp == "IIWI" | sp == "APAP") { 
n.K.nb.1 = n.K.nb.1 
	if(n.K.nb.1 == 3) K.nb.1_ = K.nb.1
	if(n.K.nb.1 == 2) { 
		K.nb.1_ = list()
			if(d[sp, "K.nb.1_min"] == 2003) K.nb.1_[[1]] = K.nb.1[[1]]
			if(d[sp, "K.nb.1_min"] == 2004) K.nb.1_[[1]] = K.nb.1[[2]]
			if(d[sp, "K.nb.1_min"] == 2003.5) K.nb.1_[[1]] = K.nb.1[[3]]
			if(d[sp, "K.nb.1_max"] == 2003) K.nb.1_[[2]] = K.nb.1[[1]]
			if(d[sp, "K.nb.1_max"] == 2004) K.nb.1_[[2]] = K.nb.1[[2]]
			if(d[sp, "K.nb.1_max"] == 2003.5) K.nb.1_[[2]] = K.nb.1[[3]]
	}
} else { K.nb.1_ = list(NA); n.K.nb.1 = 1 }

# Malaria parameters
if(sp == "IIWI" | sp == "APAP") {  
	alpha.b = alpha.b; alpha.nb.1 = alpha.nb.1; alpha.nb.2 = alpha.nb.2 
}	else {
	alpha.b = alpha.1; alpha.nb.1 = alpha.nb.2 = NA 
}

if(input.direct.Sm.ac == FALSE) Sm.ac_ = seq(d[sp, "Sm.ac_min"], d[sp, "Sm.ac_max"], length.out = n.Sm.ac)
if(input.direct.Sm.ac == TRUE) Sm.ac_ = Sm.ac.direct

# Movements
if(sp == "IIWI" | sp == "APAP") 	{ 
	n.gamma.mov = n.gamma.mov 
		if(input.direct.gamma.mov == FALSE) gamma.mov_ = seq(d[sp, "gamma.mov_min"], d[sp, "gamma.mov_max"], length.out = n.gamma.mov)
		if(input.direct.gamma.mov == TRUE) gamma.mov_ = gamma.mov.direct
} else { 
	n.gamma.mov = 1
	gamma.mov_ = -10
}
R.ter_ = seq(d[sp, "R.ter_min"], d[sp, "R.ter_max"], length.out = n.R.ter)
fidelity.ad_ = seq(d[sp, "fidelity.ad_min"], d[sp, "fidelity.ad_max"], length.out = n.fidelity.ad)
m.natal_ = seq(d[sp, "m.natal_min"], d[sp, "m.natal_max"], length.out = n.m.natal)
psi.DD_ = seq(d[sp, "psi.DD_min"], d[sp, "psi.DD_max"], length.out = n.psi.DD)

#------------------

batch.n = 1

for(i in 1 : n.s.ad) {
	for(j in 1 : n.rat.s) { 
		for(k in 1 : n.fec) {
			for(l in 1 : n.fec.1) {
				for(m in 1 : n.rat.f ) {
					for(n in 1 : n.K.b) {
						for(o in 1 : n.K.nb.1) {
							for(p in 1 : n.Sm.ac) {
								for(q in 1 : n.gamma.mov) {
									for(r in 1 : n.R.ter) {				
										for(t in 1 : n.fidelity.ad) {								
											for(u in 1 : n.m.natal) {			
												for(v in 1 : n.psi.DD)  {		

													if(fec.1_[l] > fec_[k]) next
													if(rat.s_[j] > rat.f_[m]) next
													
													if(paired.s.ad.fec == TRUE) {
														if(i != k) next
													}
													
													s = simul_mamo( 				
														# Spatial structure
														nr = nr, nc = nc, grad = grad, unit = unit,
														# Time frame
														T = T, Tm = Tm, SD.fledg = SD.fledg,
														t.b = t.b, f.nb.1 = f.nb.1, min.fledg = min.fledg, peak.fledg = peak.fledg, 
														# Initial conditions
														init.1 = round( (K.b_[n] / 2), 0 ), init.2 = round( (K.b_[n] / 2), 0 ),
														# Survival (ad = annual, juv = from fledging to breeding age)
														s.ad = s.ad_[i], s.juv = (s.ad_[i] / 2), 
														rat.s = rat.s_[j], 
														# Reproduction and habitat quality
														fec = fec_[k], 
														fec.1 = fec.1_[l],
														rat.f = rat.f_[m], 
														K.b = round( K.b_[n], 0 ), thr.DD = round( (K.b_[n] / 2), 0 ),
														K.nb.1 = K.nb.1_[[o]],
														K.nb.2 = K.nb.2, 
														reproduction.malaria = reproduction.malaria,											
														# Malaria parameters
														alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2,
														Sm.ac = Sm.ac_[p],
														# Movements
														gamma.mov = gamma.mov_[q],
														R.ter = R.ter_[r],
														fidelity.ad = fidelity.ad_[t],
														m.natal = m.natal_[u], SD.natal = m.natal_[u],
														psi.DD = psi.DD_[v],			
														# Fit to expected distribution
														y.obs = y.obs,						
														n.sim = n.sim, batch = batch.n, output.dir = output.dir	
													)

													f.s(s, batch.n)
													 
													batch.n = batch.n +1
							
}}}}}}}}}}}}}#

} # f.calibr

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
f.calibr = cmpfun(f.calibr)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------