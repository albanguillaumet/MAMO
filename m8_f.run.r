#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

f.run = function( 

# species-specific 
sp = "AKIP",
envp = akip.calibr$envp,
output.dir = "C:/Alban/Alban HAWAII/HAKALAU/R code/RUN/AKIP.11",
y.obs =  c(NA, y.obs.AKIP), # pairs / ha
d = read.table("C:/Alban/Alban HAWAII/HAKALAU/Data/Model calibration/Parameters/table.calib.txt", header = T, sep = "\t", dec = "."),

# Spatial structure
nr = 11, nc = 3, grad = c(2000, 1000), unit = 1,

# Time frame
T = 60, Tm = 5, SD.fledg = 0,

# Survival (ad = annual, juv = from fledging to breeding age)
n.RAT = 2,
mg.RAT = list( c(rep(0.3, 4), rep(1, 7)), rep(1, 11) ),

# Reproduction and habitat quality
K.nb.1 = list(K.nb.1.2003_nc.11, K.nb.1.2004_nc.11, K.nb.1.avg_nc.11),
n.K.NB.1 = 3,
mg.K.NB.1 = list( rep(1, 11), c(rep(1.5, 4), rep(1, 7)), c(rep(2, 4), rep(1, 7)) ), # n.K.NB.1 = 1; mg.K.NB.1 = rep(1, 11) # no management
K.nb.2 = K.nb.2_nc.11,
n.K.NB.2 = 3,
mg.K.NB.2 = list( rep(1, 11), c(rep(1.5, 4), rep(1, 7)), c(rep(2, 4), rep(1, 7)) ), # n.K.NB.2 = 1; mg.K.NB.2 = rep(1, 11) # no management
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
n.sim = 1, design = "simple", batch = 1

) {
	
f.s = function(x = s, n = batch.n) {
name = paste("s", n, sep ="")
save(x, file = paste(name, ".rdata", sep = ""))
}

#------------------

rownames(d) = d$species

setwd(output.dir)

n.envp = dim(envp)[1]

if(sp == "IIWI" | sp == "APAP") { 
	# Identify which K.nb.1 has been used in envp
	v = envp$K.nb.1.high
	wch.K.nb.1 = rep(NA, n.envp)
	for(i in 1 : n.envp) {
		if(round(v[i], 5) == 3.23333) wch.K.nb.1[i] = 1
		if(round(v[i], 5) == 9.58000) wch.K.nb.1[i] = 2
		if(round(v[i], 5) == 6.40667) wch.K.nb.1[i] = 3
	}
 } else { 
	wch.K.nb.1 = rep(1, n.envp)
}

#i
if(sp == "IIWI" | sp == "APAP") {  
	alpha.b_ = list( rep(0, nr), (alpha.b * mg.RISK), alpha.b, (alpha.b.2100 * mg.RISK), alpha.b.2100 )
	alpha.nb.1_ = list( rep(0, nr), (alpha.nb.1 * mg.RISK), alpha.nb.1, (alpha.nb.1.2100 * mg.RISK), alpha.nb.1.2100 )
	alpha.nb.2_ = list( rep(0, nr), (alpha.nb.2 * mg.RISK), alpha.nb.2, (alpha.nb.2.2100 * mg.RISK), alpha.nb.2.2100 )
} else {
	alpha.b_ = list( rep(0, nr), (alpha.1 * mg.RISK), alpha.1, (alpha.1.2100 * mg.RISK), alpha.1.2100 )
	alpha.nb.1_ = alpha.nb.2_ = list( rep(NA, nr),  rep(NA, nr), rep(NA, nr), rep(NA, nr), rep(NA, nr) )
}

#j
Sm.ac_ = list(); for(i in 1 : n.AC) Sm.ac_[[i]] =  envp$Sm.ac * evol.AC[i]

#---

batch.n = 1

vec.RISK = vec.AC = vec.RAT = vec.RES.1 = vec.RES.2 = vec.envp = c()

for(i in 1 : n.RISK) {
	for(j in 1 : n.AC) { 
		for(k in 1 : n.RAT) {
			for(l in 1 : n.K.NB.1) {
				for(m in 1 : n.K.NB.2) {
					for(e in 1 : n.envp) {
								
								for(z in 1 : n.sim) {
									vec.RISK = c(vec.RISK, i)
									vec.AC = c(vec.AC, j)
									vec.RAT = c(vec.RAT, k)
									vec.RES.1 = c(vec.RES.1, l)
									vec.RES.2 = c(vec.RES.2, m)
									vec.envp = c(vec.envp, e)
								}
	
								s = simul_mamo( 				
									# Spatial structure
									nr = nr, nc = nc, grad = grad, unit = unit,
									# Time frame
									T = T, Tm = Tm, SD.fledg = SD.fledg,
									t.b = d[sp, "t.b"], f.nb.1 = d[sp, "f.nb.1"], min.fledg = d[sp, "min.fledg"], peak.fledg = d[sp, "peak.fledg"], 
									# Initial conditions
									init.1 = round( (envp$K.b[e] / 2), 0 ), init.2 = round( (envp$K.b[e] / 2), 0 ),
									# Survival (ad = annual, juv = from fledging to breeding age)
									s.ad = envp$s.ad[e], s.juv = (envp$s.ad[e] / 2),
									# Reproduction and habitat quality
									fec = envp$fec[e],
									fec.1 = envp$fec.1[e], 
									K.b = envp$K.b[e], thr.DD = round( (envp$K.b[e] / 2), 0 ),
									reproduction.malaria = reproduction.malaria,		
									# Movements
									gamma.mov = envp$gamma.mov[e],
									calc.gamma.d = calc.gamma.d,
									n.sim.disp = n.sim.disp,
									R.ter = envp$R.ter[e],
									fidelity.ad = envp$fidelity.ad[e],
									m.natal = envp$m.natal[e], SD.natal = envp$m.natal[e],
									psi.DD = envp$psi.DD[e],			
									# i
									alpha.b = alpha.b_[[i]], alpha.nb.1 = alpha.nb.1_[[i]], alpha.nb.2 = alpha.nb.2_[[i]],
									# j
									Sm.ac = min( Sm.ac_[[j]][e], 1 ),
									# k 
									rat.s = mg.RAT[[k]] * envp$rat.s[e], 
									rat.f = mg.RAT[[k]] * envp$rat.f[e], 		
									# l
									K.nb.1 = mg.K.NB.1[[l]] * K.nb.1[[ wch.K.nb.1[e] ]],
									# m
									K.nb.2 = mg.K.NB.2[[m]] * K.nb.2,
									# Other options 
									add.cline = add.cline,									
									# Fit to expected distribution
									y.obs = y.obs,
									# Simulations 
									n.sim = n.sim, design = design, batch = batch.n, output.dir = output.dir	
								)

								f.s(s, batch.n)
								 
								batch.n = batch.n +1
								
}}}}}}

n.comb = n.RISK * n.AC * n.RAT * n.K.NB.1 * n.K.NB.2 * n.envp

factors = as.data.frame(cbind(
RISK = vec.RISK, AC = vec.AC, RAT.S = vec.RAT, RAT.F = vec.RAT, RES.1 = vec.RES.1, RES.2 = vec.RES.2, ENVP = vec.envp, SIM = rep(1:n.sim, n.comb)
))

write.table(factors,  file = "factors.txt", sep = "\t", col.names = TRUE, row.names = FALSE)

factors

} # f.run	


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
f.run = cmpfun(f.run)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
