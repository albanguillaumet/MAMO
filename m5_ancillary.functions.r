#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

simul_mamo = function(

# Spatial structure
nr = 10, nc = 2, 
grad = c(1900, 1000),
unit = 1,

# Time frame
T = 60, Tm = 5, t.b = 242, f.nb.1 = 0.5,
min.fledg = 100, peak.fledg = (2/3), SD.fledg = 40,

# Initial conditions
init.1 = 275, init.2 = 275, 

# Survival (ad = annual, juv = from fledging to breeding age)
s.ad = 0.729, rat.s = 0,
s.juv = (0.729/2),

# Reproduction and habitat quality
fec = 3, fec.1 = 3, rat.f = 0, 
K.b = 550, thr.DD = 275,
K.nb.1 = K.nb.1.avg, K.nb.2 =  K.nb.2,
reproduction.malaria = "simple",

# Malaria parameters (daily except Sm.ac)
alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, Sm.ac = 0.16, 

# Movements
gamma.mov = 0.159, calc.gamma.d = "fast.risky",
n.sim.disp = 10000,
R.ter = 0.0234,
fidelity.ad = 0.2,
m.natal = 0.3, SD.natal = 0.3, psi.DD = 1,

# Other options 
add.cline = FALSE,

# Fit to expected distribution
# y.obs = # pairs / ha
y.obs = c(4.83508206, 4.63024348, 4.79647277, 4.43224374, 3.44144990, 1.63243871, 0.22050535, 0.00000000, 0.00574044, 0.00000000),

# Simulations 
n.sim = 2, design = "simple", batch = 1, output.dir = "C:/Alban/Alban HAWAII/HAKALAU/R code/CALIBRATION/IIWI"

) {
	
setwd(output.dir)

		
	if(design == "simple") { 
		nr_s = nc_s = grad.max_s = grad.min_s = unit_s = 
		T_s = Tm_s = t.b_s = t.nb.1_s = t.nb.2_s = 
		min.fledg_s = peak.fledg_s = SD.fledg_s  = 
		init.1_s = init.2_s = 
		s.ad_s = s.ad.d_s = rat.s_s = s.juv_s = s.juv.d_s =
		fec_s = fec.1_s = rat.f_s = K.b_s = thr.DD_s = 
		K.nb.1.high_s = K.nb.1.mid_s = K.nb.1.low_s = K.nb.2.high_s = K.nb.2.mid_s = K.nb.2.low_s = 
		reproduction.malaria_s = 
		alpha.b.low_s = alpha.nb.1.low_s = alpha.nb.2.low_s = Sm.ac_s = 
		gamma.mov_s = calc.gamma.d_s = n.sim.disp_s = R.ter_s = fidelity.ad_s = m.natal_s =  SD.natal_s = psi.DD_s =  
		add.cline_s = 
		gamma.d.breed_s = test.disp.breed_s =
		disp.breed.p1_s = disp.natal.p1_s = disp.natal.DD.p1_s = disp.nb.1.p1_s = disp.nb.2.p1_s = 
		np.high_s = np.mid_s = np.low_s = np.metapop_s = gr.high_s = gr.mid_s = gr.low_s =
		r.hm_s =  r.ml_s = 
		cline.5_s = cline.2_s = cline.8_s =
		m.elev_s = m.elev.np_s = m.elev.K.nb.1_s = m.elev.K.nb.2_s = m.elev.gr_s = m.elev.alpha.b_s = m.elev.alpha.nb.1_s = m.elev.alpha.nb.2_s = 
		data.fit_s = 
		rep(NA, n.sim) 
		} 

#------------------

	for(s in 1: n.sim) {

		x  = mamo(
		nr = nr, nc = nc, 
		grad = grad,
		unit = unit,
		T = T, Tm = Tm, t.b = t.b, f.nb.1 = f.nb.1, 
		min.fledg = min.fledg, peak.fledg = peak.fledg, SD.fledg = SD.fledg,
		init.1 = init.1, init.2 = init.2,
		s.ad = s.ad, rat.s = rat.s, s.juv = s.juv,
		fec = fec, fec.1 = fec.1, rat.f = rat.f, K.b = K.b, thr.DD = thr.DD, 
		K.nb.1 = K.nb.1, K.nb.2 = K.nb.2,
		reproduction.malaria = reproduction.malaria,
		alpha.b = alpha.b, alpha.nb.1 = alpha.nb.1, alpha.nb.2 = alpha.nb.2, Sm.ac = Sm.ac, 
		gamma.mov = gamma.mov, calc.gamma.d = calc.gamma.d,
		n.sim.disp = n.sim.disp,
		R.ter = R.ter,
		fidelity.ad = fidelity.ad, 
		m.natal = m.natal, SD.natal = SD.natal, psi.DD = psi.DD,
		add.cline = add.cline
		)
		
		save(x, file = paste("s", batch, ".", s, ".rdata", sep = ""))
		
#------------------
		
	#  Simulation results table
	
	t = read.table("t.sim.txt", sep = "\t", header = TRUE) 
	
	get.empty = function(x) which(is.na(x))[1]
	l = get.empty(t[,1])
	
	t[l, "total.sim"] = l
	t[l, "batch"] = batch
	t[l, "sim"] = s
	t[l, "design"] = design
	
		if(design == "simple") { 
			nr_s[s] = t[l, "nr"] = x$nr
			nc_s[s] = t[l, "nc"] = x$nc
			grad.max_s[s] = t[l, "grad.max"] = max(x$grad)
			grad.min_s[s] = t[l, "grad.min"] = min(x$grad)
			unit_s[s] = t[l, "unit"] = x$unit
			T_s[s] = t[l, "T"] = x$T		
			Tm_s[s] = t[l, "Tm"] = x$Tm
			t.b_s[s] = t[l, "t.b"] = x$t.b[1,1]
			t.nb.1_s[s] = t[l, "t.nb.1"] = x$t.nb.1[1,1]
			t.nb.2_s[s] =  t[l, "t.nb.2"] = x$t.nb.2[1,1]
			min.fledg_s[s] = t[l, "min.fledg"] = x$min.fledg[1,1]
			peak.fledg_s[s] = t[l, "peak.fledg"] = x$peak.fledg[1,1]
			SD.fledg_s[s] = t[l, "SD.fledg"] = x$SD.fledg[1,1]
			init.1_s[s] = t[l, "init.1"] = x$init.1[1,1]
			init.2_s[s] = t[l, "init.2"] = x$init.2[1,1]
			s.ad_s[s] = t[l, "s.ad"] = x$s.ad[1,1]
			s.ad.d_s[s] = t[l, "s.ad.d"] = x$s.ad.d[1,1]
			rat.s_s[s] = t[l, "rat.s"] = x$rat.s[1,1]
			s.juv_s[s] = t[l, "s.juv"] = x$s.juv[1,1]
			s.juv.d_s[s] = t[l, "s.juv.d"] = x$s.juv.d[1,1]			
			fec_s[s] = t[l, "fec"] = x$fec[1,1]			
			fec.1_s[s] = t[l, "fec.1"] = x$fec.1[1,1]		
			rat.f_s[s] = t[l, "rat.f"] = x$rat.f[1,1]			
			K.b_s[s] = t[l, "K.b"] = x$K.b[1,1]
			thr.DD_s[s] = t[l, "thr.DD"] = x$thr.DD[1,1]
			K.nb.1.high_s[s] = t[l, "K.nb.1.high"] = x$K.nb.1.high
			K.nb.1.mid_s[s] = t[l, "K.nb.1.mid"] = x$K.nb.1.mid
			K.nb.1.low_s[s] = t[l, "K.nb.1.low"] = x$K.nb.1.low
			K.nb.2.high_s[s] = t[l, "K.nb.2.high"] = x$K.nb.2.high	
			K.nb.2.mid_s[s] = t[l, "K.nb.2.mid"] = x$K.nb.2.mid
			K.nb.2.low_s[s] = t[l, "K.nb.2.low"] = x$K.nb.2.low
			reproduction.malaria_s[s] = t[l, "reproduction.malaria"] = x$reproduction.malaria
			alpha.b.low_s[s] = t[l, "alpha.b.low"] = x$alpha.b.low
			alpha.nb.1.low_s[s] = t[l, "alpha.nb.1.low"] = x$alpha.nb.1.low
			alpha.nb.2.low_s[s] = t[l, "alpha.nb.2.low"] = x$alpha.nb.2.low
			Sm.ac_s[s] =  t[l, "Sm.ac"] = x$Sm.ac
			gamma.mov_s[s] = t[l, "gamma.mov"] = x$gamma.mov
			calc.gamma.d_s[s] = t[l, "calc.gamma.d"] = x$calc.gamma.d
			n.sim.disp_s[s] = t[l, "n.sim.disp"] = x$n.sim.disp	
			R.ter_s[s] = t[l, "R.ter"] = x$R.ter	
			fidelity.ad_s[s] = t[l, "fidelity.ad"] = x$fidelity.ad
			m.natal_s[s] = t[l, "m.natal"] = x$m.natal
			SD.natal_s[s] = t[l, "SD.natal"] = x$SD.natal
			psi.DD_s[s] = t[l, "psi.DD"] = x$psi.DD	
			add.cline_s[s] = t[l, "add.cline"] = x$add.cline
			gamma.d.breed_s[s] = t[l, "gamma.d.breed"] = x$gamma.d.breed
			test.disp.breed_s[s] = t[l, "test.disp.breed"] = x$test.disp.breed	
			disp.breed.p1_s[s] = t[l, "disp.breed.p1"] = x$disp.breed[[1]][1,1]
			disp.natal.p1_s[s] = t[l, "disp.natal.p1"] = x$disp.natal[[1]][1,1]
			disp.natal.DD.p1_s[s] = t[l, "disp.natal.DD.p1"] = ifelse( is.matrix(x$disp.natal.DD[[1]]), x$disp.natal.DD[[1]][1,1], x$disp.natal.DD )
			disp.nb.1.p1_s[s] = t[l, "disp.nb.1.p1"] = ifelse( is.matrix(x$disp.nb.1[[1]]), x$disp.nb.1[[1]][1,1], x$disp.nb.1 )
			disp.nb.2.p1_s[s] = t[l, "disp.nb.2.p1"] = ifelse( is.matrix(x$disp.nb.2[[1]]), x$disp.nb.2[[1]][1,1], x$disp.nb.2 )
			np.high_s[s] = t[l, "np.high"] = x$np.high
			np.mid_s[s] = t[l, "np.mid"] = x$np.mid
			np.low_s[s] = t[l, "np.low"] = x$np.low
			np.metapop_s[s] = t[l, "np.metapop"] = x$np.metapop
			gr.high_s[s] = t[l, "gr.high"] = x$gr.high
			gr.mid_s[s] = t[l, "gr.mid"] = x$gr.mid
			gr.low_s[s] = t[l, "gr.low"] = x$gr.low
			r.hm_s[s] = t[l, "r.hm"] = x$r.hm 
			r.ml_s[s] = t[l, "r.ml"] = x$r.ml
			cline.5_s[s] = t[l, "cline.5"] = x$cline.5
			cline.2_s[s] = t[l, "cline.2"] = x$cline.2
			cline.8_s[s] = t[l, "cline.8"] = x$cline.8
			m.elev_s[s] = t[l, "m.elev"] = x$m.elev
			m.elev.np_s[s] = t[l, "m.elev.np"] = x$m.elev.np
			m.elev.K.nb.1_s[s] = t[l, "m.elev.K.nb.1"] = x$m.elev.K.nb.1
			m.elev.K.nb.2_s[s] = t[l, "m.elev.K.nb.2"] = x$m.elev.K.nb.2
			m.elev.gr_s[s] = t[l, "m.elev.gr"] = x$m.elev.gr
			m.elev.alpha.b_s[s] = t[l, "m.elev.alpha.b"] = x$m.elev.alpha.b
			m.elev.alpha.nb.1_s[s] = t[l, "m.elev.alpha.nb.1"] = x$m.elev.alpha.nb.1	
			m.elev.alpha.nb.2_s[s] = t[l, "m.elev.alpha.nb.2"] = x$m.elev.alpha.nb.2

			y.sim = x$np / (x$unit^2 *100 ) # = # pairs / ha; e.g., if unit = 1, x$unit^2 = 1 (100 ha/ Km2); if unit = 0.5, x$unit^2 = 0.25 (25 ha/ 0.25 Km2);
			data.fit_s[s] = t[l, "data.fit"] = - sum( (y.sim - y.obs)^2 )
	
		} # if(design == "simple") 

	write.table(t, "t.sim.txt", sep = "\t", dec = ".", row.names = FALSE)
			
	} # for(s in 1: nsim) 

#------------------

	if(design == "simple") { 

		l = list(
		
			n.sim = n.sim, design = design, 
			
			nr_s = nr_s,
			nc_s = nc_s,
			grad.max_s = grad.max_s,
			grad.min_s = grad.min_s,
			unit_s = unit_s,
			T_s = T_s,	
			Tm_s = Tm_s,
			t.b_s = t.b_s,
			t.nb.1_s = t.nb.1_s,
			t.nb.2_s = t.nb.2_s,
			min.fledg_s = min.fledg_s,
			peak.fledg_s = peak.fledg_s,
			SD.fledg_s = SD.fledg_s,
			init.1_s = init.1_s,
			init.2_s = init.2_s,
			s.ad_s = s.ad_s,
			s.ad.d_s = s.ad.d_s,
			rat.s_s = rat.s_s,
			s.juv_s = s.juv_s,
			s.juv.d_s = s.juv.d_s,			
			fec_s = fec_s,			
			fec.1_s = fec.1_s,		
			rat.f_s = rat.f_s,
			K.b_s = K.b_s,
			thr.DD_s = thr.DD_s,
			K.nb.1.high_s = K.nb.1.high_s,
			K.nb.1.mid_s = K.nb.1.mid_s,
			K.nb.1.low_s = K.nb.1.low_s,
			K.nb.2.high_s = K.nb.2.high_s,
			K.nb.2.mid_s = K.nb.2.mid_s,
			K.nb.2.low_s = K.nb.2.low_s,
			reproduction.malaria_s = reproduction.malaria_s,
			alpha.b.low_s = alpha.b.low_s,
			alpha.nb.1.low_s = alpha.nb.1.low_s,
			alpha.nb.2.low_s = alpha.nb.2.low_s,
			Sm.ac_s =  Sm.ac_s,
			gamma.mov_s = gamma.mov_s,
			calc.gamma.d_s = calc.gamma.d_s,
			n.sim.disp_s = n.sim.disp_s,
			R.ter_s = R.ter_s,
			fidelity.ad_s = fidelity.ad_s,
			m.natal_s = m.natal_s,
			SD.natal_s = SD.natal_s,
			psi.DD_s = psi.DD_s,
			add.cline_s = add.cline_s,
			gamma.d.breed_s = gamma.d.breed_s,
			test.disp.breed_s = test.disp.breed_s,
			disp.breed.p1_s = disp.breed.p1_s,
			disp.natal.p1_s = disp.natal.p1_s,
			disp.natal.DD.p1_s = disp.natal.DD.p1_s,
			disp.nb.1.p1_s = disp.nb.1.p1_s,
			disp.nb.2.p1_s = disp.nb.2.p1_s,
			np.high_s = np.high_s,
			np.mid_s = np.mid_s,
			np.low_s = np.low_s,
			np.metapop_s = np.metapop_s,
			gr.high_s = gr.high_s,
			gr.mid_s = gr.mid_s,
			gr.low_s = gr.low_s,
			r.hm_s = r.hm_s,
			r.ml_s = r.ml_s,
			cline.5_s = cline.5_s,
			cline.2_s = cline.2_s,
			cline.8_s = cline.8_s,
			m.elev_s = m.elev_s,
			m.elev.np_s = m.elev.np_s,
			m.elev.K.nb.1_s = m.elev.K.nb.1_s,
			m.elev.K.nb.2_s = m.elev.K.nb.2_s,
			m.elev.gr_s = m.elev.gr_s,
			m.elev.alpha.b_s = m.elev.alpha.b_s,
			m.elev.alpha.nb.1_s = m.elev.alpha.nb.1_s,
			m.elev.alpha.nb.2_s = m.elev.alpha.nb.2_s,
			data.fit_s = data.fit_s
			
		)
		
	} # if(design == "simple") 
	
} # end of simul_mamo

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
simul_mamo = cmpfun(simul_mamo)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# calculates the parameters of a sigmoid cline based on x and y, and plots the cline on a pre-existing plot with parameters col, lty and lwd.

wrap.cline = function(x = x, y = y, col = "black", lty = 1, lwd = 1) {

	a = min(y); b = max(y); c = mean(x)
	
		f.cline = function(y = y, x = x, a = a, b = b, c = c, K = 25){
		nls( y ~ a + (b - a) / (1 + exp( - K * (x - c) ) ), algorithm = "port", control = c(maxiter = 10000, tol = 1e-05, minFactor = 1/2048, warnOnly = TRUE), start = list(a = a, b = b, c = c, K = K))
		}
		
	cline.fit = f.cline(y = y, x = x, a = a, b = b, c = c, K = 25)
	p = coef(cline.fit)

	x.sig = seq(max(x), min(x), length = 10000)
	
		sigmoid = function(a = p[1], b = p[2], c = p[3], K = p[4], x = x.sig) a + (b-a) / (1 + exp( - K * (x - c) ) )
	
	y.sig = sigmoid(a = p[1], b = p[2], c = p[3], K = p[4])
	lines(x.sig, y.sig, col = col, lty = lty, lwd = lwd)
	
	h.5 = mean(y.sig[x.sig > (p[3]-0.005) & x.sig < (p[3]+0.005)])
	y.sig.st = y.sig - min(y.sig) 
	y.sig.st = y.sig.st / max(y.sig.st)
	m.2 = mean(x.sig[y.sig.st > 0.195 & y.sig.st < 0.205])
	h.2 = mean(y.sig[x.sig > (m.2-0.005) & x.sig < (m.2+0.005)])
	m.8 = mean(x.sig[y.sig.st > 0.795 & y.sig.st < 0.805])
	h.8 = mean(y.sig[x.sig > (m.8-0.005) & x.sig < (m.8+0.005)])

	cline.5 = p[3]; cline.2 = m.2; cline.8 = m.8
	
	list(h.5 = h.5, cline.5 = cline.5, cline.2 = cline.2, cline.8 = cline.8)
	
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# plot the simulated data points d on a pre-existing plot with parameters col.p, pch and cex 
# facultatively (add.cline), calculate the corresponding sigmoid cline plotted with parameters col.l, lwd, lty

add.s = function(d = "r1.rdata", col.p = "black", pch = 15, cex = 1.2, add.cline = FALSE, col.l = "black", lwd = 1, lty = 2) {
load(d)
x.mamo = x$grad
y.mamo = x$np / (x$unit^2 *100) # density = # pairs per ha
points(x.mamo, y.mamo, col = col.p, pch = pch, cex = cex)
if(add.cline == TRUE) { ws = wrap.cline(x = x.mamo, y = y.mamo, col = col.l, lwd = lwd, lty = lty) ; ws }
}

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------