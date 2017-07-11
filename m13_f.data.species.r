#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

f.data.species = function(
sp = "AKIP", w = w,
directory = c("C:/programs/MAMO/RUN/akip.10", "C:/programs/MAMO/RUN/akip.11"),
vec.hab = c(0, 1), table.risk = table.risk
){

n.data = length(directory)
data = list()

for(i in 1 : n.data) {

	setwd(directory[i])
	
	d.sim = read.table("t.sim.txt", header = T, sep = "\t", dec = ".")
	d.factors = read.table("factors.txt", header = T, sep = "\t", dec = ".")

	d = d.sim
	d = d[!is.na(d[,1]),]
	n = dim(d)[1]

	if(dim(d.factors)[1] != n) stop("d.sim and d.factors have different lengths")

	d.f = cbind(d, d.factors)
	
	#---
	
	mg.hab = rep(vec.hab[i], n)
	
	tr = table.risk
	
	risk = mg.risk = rep(NA, n)
	
		for(j in 1 : n) {
		
			risk[j] = tr$RISK_analysis[ d.f$RISK[j] ] 
			mg.risk[j] = tr$mg.RISK[ d.f$RISK[j] ] 
			
		}
	
	ac =  d.f$AC
	mg.rat = d.f$RAT.S - 1 
	mg.res.1 = d.f$RES.1 - 1
	mg.res.2 = d.f$RES.2 - 1
	envp = d.f$ENVP
	sim = d.f$SIM
	data[[i]] = cbind(mg.hab, risk, mg.risk, ac, mg.rat, mg.res.1, mg.res.2, envp, sim, d.f)
	
}

d = data[[1]]
if(n.data > 1) { for(k in 2 : n.data) d = rbind(d, data[[k]]) }

n = dim(d)[1]
w1 = rep(w[sp, "w1"], n)
w2 = rep(w[sp, "w2"], n)
w3 = rep(w[sp, "w3"], n)

d = cbind(w1, w2, w3, d)

n.min = min(d$np.metapop); n.max = max(d$np.metapop)
n.st = (d$np.metapop - n.min) / (n.max - n.min) 

d = cbind(n.st, d)

d
}

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
f.data.species = cmpfun(f.data.species)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
