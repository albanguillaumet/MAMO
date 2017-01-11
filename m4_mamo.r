#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

mamo = function(

# Spatial structure
nr = 10, nc = 2, 
grad = c(1900, 1000),
unit = 1,

# Time frame
T = 30, Tm = 5, t.b = 242, f.nb.1 = 0.5,
min.fledg = 100, peak.fledg = (2/3), SD.fledg = 40,

# Initial conditions
init.1 = 290, init.2 = 290, 

# Survival (ad = annual, juv = from fledging to breeding age)
s.ad = 0.78, rat.s = 0.05, 
s.juv = (0.78 / 2), 

# Reproduction and habitat quality
fec = 2.444, fec.1 = 2.444, rat.f = 0.10, K.b = 580, thr.DD = 290,
K.nb.1 = rep(c(0,1,2), 3), K.nb.2 = rep(c(0,2,1), 3),
reproduction.malaria = "simple",

# Malaria parameters (daily except Sm.ac)
alpha.b = rep(c(0, 0.0005, 0.05), 3), alpha.nb.1 = rep(c(0.001, 0.004, 0.2), 3), alpha.nb.2 = rep(c(0.001, 0.004, 0.2), 3), Sm.ac = 0.16, 

# Movements
gamma.mov = 10, calc.gamma.d = "fast.risky", # "slow.robust"
n.sim.disp = 10000,
R.ter = 0.0234, # sqrt(1/(pi*580)) : Pi*R^2 (territoire d'un couple) * npairs in a km^2 (Density) = Aire du quadrat (1km^2)
fidelity.ad = 0.95,
m.natal = 0.3, SD.natal = 0.3, psi.DD = 1,

# Other options 
add.cline = FALSE

){

#library(compiler)
#enableJIT(3)

#---------------------------------------------------------------------------------------------------------------- 
# INITIALISATION OF THE MODEL
#---------------------------------------------------------------------------------------------------------------- 

# Libraries and functions

options(digits = 6)

library(truncnorm)

# function to calculate when (what date) an e-bird got infected by malaria

f.when = function(dr = 0.1, n = 242) {
t = 1 : n
p = rep(0, n)
	for(i in 1 : n) p[i] = (1 - dr)^(t[i]-1) * dr
which(rmultinom(1, 1, p) == 1) 	# note: p is internally (within rmultinom) normalized to sum 1
}

# function to calculate the fledging date

if(SD.fledg == 0 | reproduction.malaria == "simple")  {

f.fledg = function(p.b, t.b, min.fledg, peak.fledg, SD.fledg) { 
r = peak.fledg[p.b] * t.b[p.b]
round(r, 0)
}

} else {

f.fledg = function(p.b, t.b, min.fledg, peak.fledg, SD.fledg) { 
r = rtruncnorm(1, a = min.fledg[p.b], b = t.b[p.b] - 14, mean = peak.fledg[p.b] * t.b[p.b], sd = SD.fledg[p.b])
round(r, 0)
}

}

# function to visualize the effect of gamma.d, controlling the effect of distance

f.e = function(g.d, d){ exp(-10^-g.d * d) }

# function to calculate reduction in fecundity due to density-dependent effects

f.DD = function(npair = 1000, thr = 500, K = 800, f = 0.7, sj = 0.5, sa = 0.8) {
slope = ( 1 - (sa + f * sj) ) / (f * sj * (K - thr))
c  = NA
 if(npair < thr) c = 1
 if(npair >= thr)  c = 1 + slope * (npair - thr)
 if(c < 0) c = 0
 c }

fc = function(k, disp = disp.b, n = 1000) { 
x = rep(NA, n); for(i in 1 : n) x[i] = which(rmultinom(1, 1, prob = disp[[k]]) == 1)
u = rep(0, 9); for(i in 1 : 9)  u[i] = length(x[x == i]); u = as.matrix(u); rownames(u) = 1:9
u = order(u, decreasing = TRUE)
list(disp = disp[[k]], most = u)
}
	
# patch-specific dispersion and 'migration' probabilities
f.disp = function(nr, nc, n.patch, unit, R.ter, n, gamma.mov, season, type, method, m, s, matrix.K) {

	# Calculates the fraction of e-birds dispersing out of the patch 
	if(season == "breeding") {

		calcul.disp = function(unit, R.ter, n, type, method, m, s) {

			x = unit * runif(n)
			y  = unit * runif(n)
			plot(x, y, pch = "+", xlim = c(-0.1, unit + 0.1), ylim = c(-0.1, unit + 0.1), col = "red", xlab = "", ylab = "", bty = "n")
			abline(h = 0, lty = 2, col = "red");  abline(h = unit, lty = 2, col = "red");  abline(v = 0, lty = 2, col = "red");  abline(v = unit, lty = 2, col = "red")

			if(type == "natal" & method == "LN") {
				mu = function(m,s) log(m^2 / sqrt(s^2 + m^2)) # location parameter = mu2 = function(m,s) log(m / sqrt( 1 + (s^2 / m^2)))
				sigma = function(m,s) sqrt(log(1+(s^2/m^2))) # scale parameter
				m.l = mu(m = m, s = s); s.l = sigma(m = m, s = s)
			}

			if(type == "breed") {
				thr = R.ter
				lambda = -log(m) / thr
				} else {
				thr = lambda = NA
			}

			# Calculates dispersion out of the patch 
			if(type == "natal" & method == "LN") {
			r.m = rlnorm(n, meanlog = m.l, sdlog = s.l)
			m.EXP = NA
			long.dist.EXP = plnorm(3, meanlog = m.l, sdlog = s.l, lower.tail = FALSE)
			}
			if(type == "breed") {
			r.m = rexp(n, rate = lambda)
			m.EXP = pexp(thr, rate = lambda, lower.tail = FALSE) 
			long.dist.EXP = pexp(3, rate = lambda, lower.tail = FALSE)
			}
			
			t.m  = 2 * pi * runif(n)
			x.m = x + r.m * cos(t.m) ; y.m = y + r.m * sin(t.m)
			segments(x, y, x.m, y.m, col = "orange")

			test.m = ( x.m > unit | x.m < 0 | y.m > unit  | y.m < 0) 
			disp = length( test.m[test.m == TRUE] ) / n
			
			list(thr.breed = thr, lambda.breed = lambda, m.EXP = m.EXP, long.dist.EXP = long.dist.EXP, X = r.m, mean.X = mean(r.m), SD.X = sd(r.m), disp = disp)

		}
		
		cd = calcul.disp(unit = unit, R.ter = R.ter, n = n, type = type, method = method, m = m, s = s) 

		thr.breed = cd$thr.breed
		lambda.breed = cd$lambda.breed
		beta.breed = (1 / lambda.breed)
		m.EXP = cd$m.EXP
		long.dist.EXP = cd$long.dist.EXP
		X = cd$X
		if(type == "breed") m.X = ( length( X[ X > thr.breed ] ) / n ) else m.X = NA
		mean.X = cd$mean.X
		SD.X = cd$SD.X
		fraction.disp = cd$disp

	} else {
		thr.breed = lambda.breed = beta.breed = m.EXP = long.dist.EXP = X = m.X = mean.X = SD.X = fraction.disp = NA
	}

	# section calculates the values of d and gamma.d and use them to calculate the distribution of dispersers and migrants across the landscape

	# calculus of d = distances between center of (square) patches of length (L, l) = unit

	calcul.dist = function(nr, nc, n.patch) { 

		rc = as.data.frame(matrix(NA, n.patch, 2)); colnames(rc) = c("row", "col")

		f = 1
			for(i in 1 : nc) {
				for(j in 1 : nr)  {
				rc$col[f] = i
				rc$row[f] = j
				f = f + 1
				}
			}

		fp = function(a,b) (a^2+b^2)^0.5

		d = matrix(NA, n.patch, n.patch)
			for(i in 1 : n.patch) {
				for(j  in 1 : n.patch) {
				d[i,j] = d[j,i] = fp( abs( (rc$row[i] - rc$row[j]) * unit) , abs( (rc$col[i] - rc$col[j]) * unit) )
				}
			}
		d 
	}

	d = calcul.dist(nr, nc, n.patch) 
	
	#---

	dxy = disp = list()

		for(i in 1 : n.patch) {
		
			# builds the matrix of distance from each patch to focal patch i
			dxy[[i]] = matrix(d[i,], nr, nc)		
			
			# calculates gamma.d
			if(season == "breeding" & calc.gamma.d == "slow.robust") {
			
				g.k = seq(-5, 5, length = 10000)
				y.k = rep(NA, length(g.k))

				for(k in 1 : length(g.k)) {
					den.k = sum( exp( -10^( - g.k[k] ) * as.vector( dxy[[i]] ) ) ) 
					y.k[k] = abs( (1 / den.k) - (1 - fraction.disp) ) 
				}

				plot(g.k, y.k)
				gamma.d = g.k[y.k == min(y.k)]
				abline(v = gamma.d)				
			}  
			
			if(season == "breeding" & calc.gamma.d == "fast.risky") {
			
				# calculates gamma.d
				f = function(g) {
				d.num = dxy[[i]][ dxy[[i]] == 0 ] # 0
				d.den = as.vector(dxy[[i]])
				num = sum( exp( -10^(- g ) * d.num ) ) # 1
				den = sum( exp( -10^(- g ) * d.den ) ) 
				abs( (num / den) - (1 - fraction.disp) ) 
				}
				gamma.d = optimize(f, interval = c(-5, 5))$minimum
			}
			
			if(season != "breeding") gamma.d = gamma.mov
			
			# disp calculates the probability to move in each patch (matrix form) for an e-bird of patch i
			disp[[i]] = exp( -10^(-gamma.d) * dxy[[i]] ) / sum( exp( -10^(-gamma.d) * dxy[[i]] ) )
				
			# integrates landscape heterogeneity in resources
			disp[[i]] = ( disp[[i]] * matrix.K ) / sum( ( disp[[i]] * matrix.K ) )
				
		} # for(i in 1 : n.patch)

	list(
	thr.breed = thr.breed, lambda.breed = lambda.breed, 
	beta.breed = beta.breed, mean.X = mean.X, SD.X = SD.X, 
	m = m, m.EXP = m.EXP,  m.X = m.X, long.dist.EXP = long.dist.EXP,
	X = X, 
	fraction.disp = fraction.disp, gamma.d = gamma.d, 
	disp = disp
	) 
}

# function to calculate the probability (rate) of infection on a time-period t, based on a known (constant) daily risk

f.alpha = function(daily.rate, t) {
lambda = log( 1 / (1 - daily.rate) )
( 1 - exp(- lambda * t) )
}
	
#-------------------------------------------------------------------------------------------------------------------------

# Loading simulation parameters

# Spatial structure
patch = matrix(1 : (nr * nc), nr, nc)
n.patch = nr * nc

# Time frame
t.b.unique = t.b
t.b = matrix(t.b, nr, nc)
min.fledg = matrix(min.fledg, nr, nc)
peak.fledg = matrix(peak.fledg, nr, nc)
SD.fledg = ifelse( reproduction.malaria == "simple", 0, SD.fledg )
SD.fledg = matrix(SD.fledg, nr, nc)
t.nb = 365 - t.b
f.nb.1 = matrix(f.nb.1, nr, nc)
t.nb.1 = round(f.nb.1 * t.nb, 0)
t.nb.2 = t.nb - t.nb.1

# Initial conditions
init.1 = matrix(init.1, nr, nc)
init.2 = matrix(init.2, nr, nc)

# Reproduction parameters
fec = matrix(fec, nr, nc) 
fec.1 = matrix(fec.1, nr, nc) 
rat.f = matrix(rat.f, nr, nc) 

# Habitat quality and density-dependence
K.b = matrix(K.b, nr, nc)
thr.DD = matrix(thr.DD, nr, nc) 
K.nb.1 = matrix(K.nb.1, nr, nc)
K.nb.2 = matrix(K.nb.2, nr, nc)

# Malaria parameters
alpha.b = matrix(alpha.b, nr, nc)
alpha.nb.1 = matrix(alpha.nb.1, nr, nc)
alpha.nb.2 = matrix(alpha.nb.2, nr, nc)

# Survival 
s.ad = matrix(s.ad, nr, nc)
s.ad.d = s.ad^(1 / 365)
s.juv = matrix(s.juv, nr, nc)
length.yr.juv = 365 - ( peak.fledg * t.b ); s.juv.d = s.juv^( 1 / length.yr.juv ) 
rat.s = matrix(rat.s, nr, nc) 

# Dispersal

# patch-specific dispersion and migration probabilities

disp.breed = f.disp(nr, nc, n.patch, unit, R.ter, n = n.sim.disp, gamma.mov = NA, season = "breeding", type = "breed", method = NA, m = (1 - fidelity.ad), s = NA, matrix.K = K.b)
gamma.d.breed = disp.breed$gamma.d # estimate for last patch
test.disp.breed = abs( ( 1 - disp.breed$disp[[1]][1,1] ) - disp.breed$fraction.disp ) # should be ~ 0
disp.breed = disp.breed$disp 
disp.natal = f.disp(nr, nc, n.patch, unit, R.ter = NA, n = n.sim.disp, gamma.mov = NA, season = "breeding", type = "natal", method = "LN", m = m.natal, s = SD.natal, matrix.K = K.b)$disp 
disp.natal.DD = disp.nb.1 = disp.nb.2 = NA
if(psi.DD > 1) disp.natal.DD = f.disp(nr, nc, n.patch, unit, R.ter = NA, n = n.sim.disp, gamma.mov = NA, season = "breeding", type = "natal", method = "LN", m = (m.natal * psi.DD), s = (SD.natal * psi.DD), matrix.K = K.b)$disp 
if( length( unique( as.vector( t.b.unique ) ) ) == 1 & unique( as.vector( t.b.unique ) ) < 365 ) {
disp.nb.1 = f.disp(nr, nc, n.patch, unit, R.ter = NA, n = NA, gamma.mov, season = "non-breeding", type = NA, method = NA, m = NA, s = NA, matrix.K = K.nb.1)$disp 
disp.nb.2 = f.disp(nr, nc, n.patch, unit, R.ter = NA, n = NA, gamma.mov, season = "non-breeding", type = NA, method = NA, m = NA, s = NA, matrix.K = K.nb.2)$disp 
}

#---------------------------------------------------------------------------------------------------------------- 
# RUN THE SIMULATION
#---------------------------------------------------------------------------------------------------------------- 

n.1 = n.ad = n.juv = n.pairs = array(0, c(nr, nc, T))

mal = list(); mal$y = mal$a = list()

for(t in 1 : T) {

	if(t == 1) {
		n.1[,,t] = init.1
		n.ad[,,t] = init.2	
		n.pairs[,,t] = n.1[,,t] + n.ad[,,t]
	}

	mal.y = mal.a = list()
	for(s in 1 : n.patch) mal.y[[s]] = mal.a[[s]] = NA # initiate the list
	
		for(i in 1 : n.patch) { 

			# Identify the breeding patch
			p.b = i
			
			# Estimate malaria infection and survival during the breeding & non-breeding season
				
				if(t <= Tm) { 
				
					# breeding season
					
					# Calculate the fledging date of offspring
					fledg.t = round(peak.fledg[p.b] * t.b[p.b], 0)
					t.b.juv = t.b[p.b] - fledg.t + 1
								
					n.1_t.i = n.1[,,t][p.b]
					n.ad_t.i = n.ad[,,t][p.b]
					n.pairs_t.i =  n.1_t.i + n.ad_t.i
					fec_t.i = (n.1_t.i * fec.1[p.b] + n.ad_t.i * fec[p.b] ) / n.pairs_t.i
					
					E.juv_t.i = fec_t.i * f.DD(npair = n.pairs_t.i, thr = thr.DD[p.b], K = K.b[p.b], f = fec_t.i, sj = s.juv[p.b], sa = s.ad[p.b]) 
					n.juv_t.i = sum( rpois(n.pairs_t.i, E.juv_t.i) )
					
					n.juv_t.i = n.juv[,,t][p.b] = rbinom(1, n.juv_t.i, s.juv.d[p.b]^t.b.juv) 
				
					n.1_t.i = rbinom(1, n.1_t.i, s.ad.d[p.b]^t.b[p.b])    

					n.ad_t.i  = rbinom(1, n.ad_t.i, s.ad.d[p.b]^t.b[p.b])    
																						
					if(t.b[p.b] < 365) {
					# non-breeding season
					
						# first period
						
							#juveniles
							if(n.juv_t.i > 0) {
								p.nb.1_juv = which(rmultinom(n.juv_t.i, 1, prob = disp.nb.1[[p.b]]) == 1) 
								p.nb.1_juv = p.nb.1_juv - seq(0, n.patch*(n.juv_t.i - 1), length.out = n.juv_t.i)

								surv.juv = rbinom(n.juv_t.i, 1, s.juv.d[p.nb.1_juv]^t.nb.1[p.nb.1_juv]) 
								n.juv_t.i = sum(surv.juv)
								p.nb.1_juv = p.nb.1_juv[surv.juv == 1]
							}
							
							#1-yr old
							if(n.1_t.i > 0) { 
								p.nb.1_1 = which(rmultinom(n.1_t.i, 1, prob = disp.nb.1[[p.b]]) == 1) 
								p.nb.1_1 = p.nb.1_1 - seq(0, n.patch*(n.1_t.i - 1), length.out = n.1_t.i)
								
								surv.1 = rbinom(n.1_t.i, 1, s.ad.d[p.nb.1_1]^t.nb.1[p.nb.1_1]) 
								n.1_t.i = sum(surv.1)
								p.nb.1_1 = p.nb.1_1[surv.1 == 1]
							}
														
							#ad
							if(n.ad_t.i > 0) { 
								p.nb.1_ad = which(rmultinom(n.ad_t.i, 1, prob = disp.nb.1[[p.b]]) == 1) 
								p.nb.1_ad = p.nb.1_ad - seq(0, n.patch*(n.ad_t.i - 1), length.out = n.ad_t.i)
								
								surv.ad = rbinom(n.ad_t.i, 1, s.ad.d[p.nb.1_ad]^t.nb.1[p.nb.1_ad]) 
								n.ad_t.i = sum(surv.ad)
								p.nb.1_ad = p.nb.1_ad[surv.ad == 1]
							}
							
						# second period
						
							#juveniles					
							if(n.juv_t.i > 0) {
								p.nb.2_juv = rep(NA, n.juv_t.i)
								for(k in 1 : n.juv_t.i) p.nb.2_juv[k] = which(rmultinom(1, 1, prob = disp.nb.2[[ p.nb.1_juv[k] ]]) == 1) 
								
								surv.juv = rbinom( n.juv_t.i, 1, s.juv.d[p.nb.2_juv]^t.nb.2[p.nb.2_juv] )    
								n.juv_t.i = sum(surv.juv)
							}
							
							#1-yr old
							if(n.1_t.i > 0) { 
								p.nb.2_1 = rep(NA, n.1_t.i)
								for(k in 1 : n.1_t.i) p.nb.2_1[k] = which(rmultinom(1, 1, prob = disp.nb.2[[ p.nb.1_1[k] ]]) == 1) 
								
								surv.1 = rbinom( n.1_t.i, 1, s.ad.d[p.nb.2_1]^t.nb.2[p.nb.2_1] )    
								n.1_t.i = sum(surv.1)
							}
														
							#ad
							if(n.ad_t.i > 0) { 
								p.nb.2_ad = rep(NA, n.ad_t.i)
								for(k in 1 : n.ad_t.i) p.nb.2_ad[k] = which(rmultinom(1, 1, prob = disp.nb.2[[ p.nb.1_ad[k] ]]) == 1) 
								
								surv.ad = rbinom( n.ad_t.i , 1, s.ad.d[p.nb.2_ad]^t.nb.2[p.nb.2_ad] )     
								n.ad_t.i = sum(surv.ad)
							}	
						
						} # if(t.b[p.b] < 365)
							
						# malaria status
						if(n.juv_t.i > 0) mal.juv_t.i = rep(0, n.juv_t.i) else mal.juv_t.i = NA
						if(n.1_t.i > 0) mal.1_t.i  = rep(0, n.1_t.i) else mal.1_t.i  = NA
						if(n.ad_t.i > 0) mal.ad_t.i  = rep(0, n.ad_t.i) else mal.ad_t.i  = NA
						
				} # if(t <= Tm)
	
				if(t > Tm) { 
							
				# breeding season
									
				n.1_t.i = n.1[,,t][p.b]
				n.ad_t.i = n.ad[,,t][p.b]
				n.pairs_t.i =  n.1_t.i + n.ad_t.i
				
				if(n.pairs_t.i > 0) {
				
					n.get.m_d = n.get.m_s = 0
				
					fec_t.i = (n.1_t.i * fec.1[p.b] + n.ad_t.i * fec[p.b] ) / n.pairs_t.i
					
					E.juv_t.i = fec_t.i * f.DD(npair = n.pairs_t.i, thr = thr.DD[p.b], K = K.b[p.b], f = fec_t.i, sj = s.juv[p.b], sa = s.ad[p.b]) 
									
					#---
										
					if(n.1_t.i > 0) {
					
						mal.1_t.i = mal$y[[p.b]]
						n1.0 = length(mal.1_t.i[mal.1_t.i == 0])
						get.m = rbinom( 1, n1.0, f.alpha(alpha.b[p.b], t.b[p.b]) )
						surv.ac.1 = 0; if(get.m > 0) {
							surv.ac.1 = rbinom( 1, get.m, Sm.ac )
							n.get.m_s = n.get.m_s + surv.ac.1
							n.get.m_d = n.get.m_d + (get.m - surv.ac.1)
						}			
						
						mal.1_t.i = c( rep( 1, length(mal.1_t.i[mal.1_t.i == 1]) ), rep(1, surv.ac.1), rep( 0, (n1.0 - get.m) ) )
						n.1_t.i = rbinom( 1, length(mal.1_t.i), s.ad.d[p.b]^t.b[p.b] )
				
					}
					
					if(n.1_t.i > 0) mal.1_t.i = sample(mal.1_t.i, n.1_t.i)	else mal.1_t.i = NA
					
					if(n.ad_t.i > 0) {
					
						mal.ad_t.i = mal$a[[p.b]]
						nad.0 = length(mal.ad_t.i[mal.ad_t.i == 0])
						get.m = rbinom( 1, nad.0, f.alpha(alpha.b[p.b], t.b[p.b]) )
						surv.ac.ad = 0; if(get.m > 0) {
							surv.ac.ad = rbinom( 1, get.m, Sm.ac )
							n.get.m_s = n.get.m_s + surv.ac.ad
							n.get.m_d = n.get.m_d + (get.m - surv.ac.ad)
						}
						
						mal.ad_t.i = c( rep( 1, length(mal.ad_t.i[mal.ad_t.i == 1]) ), rep(1, surv.ac.ad), rep( 0, (nad.0 - get.m) ) )
						n.ad_t.i = rbinom( 1, length(mal.ad_t.i), s.ad.d[p.b]^t.b[p.b] )	

					} 
					
					if(n.ad_t.i > 0) mal.ad_t.i = sample(mal.ad_t.i, n.ad_t.i) else mal.ad_t.i = NA
					
					#---
					
					# Reproduction 
					
					if(reproduction.malaria == "simple") {
					
						n.juv_t.i = 0
						mal.juv_t.i = NA
					
						fledg.t = f.fledg(p.b, t.b, min.fledg, peak.fledg, SD.fledg)
						
						n.pairs_t.i = (n.pairs_t.i - n.get.m_d - n.get.m_s) + n.get.m_d * ( ( t.b[p.b] - (fledg.t + 14) ) / t.b[p.b] ) + n.get.m_s * ( (t.b[p.b] - 40) / t.b[p.b] ) # warning: fledg.t + 14 must b < t.b[p.b] and t.b[p.b] must be > 40
						n.pairs_t.i = round(n.pairs_t.i, 0)
												
						if(n.pairs_t.i > 0) {
						
							n.juv_t.i = sum( rpois(n.pairs_t.i, E.juv_t.i) )
													
							t.b.juv = round(t.b[p.b] - fledg.t + 1, 0)
							
							if(n.juv_t.i > 0) {
						
								get.m = rbinom( 1, n.juv_t.i, f.alpha(alpha.b[p.b], t.b.juv) ) 
								surv.ac.juv = 0; if(get.m > 0) surv.ac.juv = rbinom( 1, get.m, Sm.ac )
								mal.juv_t.i = c( rep(1, surv.ac.juv), rep( 0, (n.juv_t.i - get.m) ) )
								n.juv_t.i = rbinom( 1, length(mal.juv_t.i), s.juv.d[p.b]^t.b.juv )
								
							} else mal.juv_t.i = NA
							
							if(n.juv_t.i > 0)	mal.juv_t.i = sample(mal.juv_t.i, n.juv_t.i) else mal.juv_t.i = NA
						
						} # if(n.pairs_t.i > 0)

					} # if(reproduction.malaria == "simple")
					
					if(reproduction.malaria == "complex") {
						
						n.juv_t.i = 0
						mal.juv_t.i = NA
						fledg.t_juv = c() 
						
						if(n.get.m_d > 0) {
						
							f.when_d = replicate( n.get.m_d, f.when(dr = alpha.b[p.b], n =  t.b[p.b]) ) # computation time-consuming
							fledg.t = replicate( n.get.m_d, f.fledg(p.b, t.b, min.fledg, peak.fledg, SD.fledg) )
							n.after.fledg = ifelse( f.when_d > (fledg.t + 14), 1, 0 )
							fledg.t_breed = fledg.t[n.after.fledg == 1]
							n.breed = length(fledg.t_breed)
							
							if(n.breed > 0) {
								juv_t.i = rpois(n.breed, E.juv_t.i)
									if(sum(juv_t.i) > 0) {
										n.juv_t.i = n.juv_t.i + sum(juv_t.i)
										for ( k in 1 : n.breed ) fledg.t_juv = c( fledg.t_juv, replicate(juv_t.i[k],  fledg.t_breed[k]) )
									} # if(sum(juv_t.i) > 0)
							} # if(n.breed > 0)
												
						} # if(n.get.m_d > 0)
						
						if(n.get.m_s > 0) {
						
							f.when_s = replicate( n.get.m_s, f.when(dr = alpha.b[p.b], n =  t.b[p.b]) )
							fledg.t = replicate( n.get.m_s, f.fledg(p.b, t.b, min.fledg, peak.fledg, SD.fledg) )
													
								for(r in 1 : n.get.m_s) {
								
									T.acute = f.when_s[r] : (f.when_s[r] + 40)
									T.critic = (fledg.t[r] - 35) : (fledg.t[r]  + 14)
									n.intersect = length( Reduce(intersect, list(T.acute, T.critic)) )
									
									if(n.intersect == 0) {
										n.offsp = rpois(1, E.juv_t.i)
										n.juv_t.i = n.juv_t.i + n.offsp 
										fledg.t_juv = c( fledg.t_juv, replicate(n.offsp,  fledg.t[r]) )
									} # if(n.intersect == 0) 
									
								} # for(r in 1 : n.get.m_s)		 
								
						} # if(n.get.m_s > 0)
						
						n.rest = (n.pairs_t.i - n.get.m_d - n.get.m_s) 
						if(n.rest > 0) {
						
							fledg.t = replicate( n.rest, f.fledg(p.b, t.b, min.fledg, peak.fledg, SD.fledg) )
							juv_t.i = rpois(n.rest, E.juv_t.i)
							n.juv_t.i = n.juv_t.i + sum(juv_t.i)
							for ( k in 1 : n.rest ) fledg.t_juv = c( fledg.t_juv, replicate(juv_t.i[k],  fledg.t[k]) )
					
						}
										
						fledg.t_juv = unlist(fledg.t_juv)					
																
						if(n.juv_t.i > 0) {		
							t.b.juv = round(t.b[p.b] - fledg.t_juv + 1, 0)
							get.m = rbinom( n.juv_t.i, 1, f.alpha(alpha.b[p.b], t.b.juv) ) 
							surv.ac.juv = ifelse( get.m == 0, 1, rbinom( n.juv_t.i, 1, Sm.ac ) )
							n.juv_t.i = sum(surv.ac.juv)						
						}

						if( n.juv_t.i > 0 ) {
							t.b.juv = t.b.juv[surv.ac.juv == 1] 	
							mal.juv_t.i = get.m[surv.ac.juv == 1]	
							juv_t.i = rbinom( n.juv_t.i, 1, s.juv.d[p.b]^t.b.juv )
							n.juv_t.i = sum(juv_t.i)
						} else mal.juv_t.i = NA
							
						if( n.juv_t.i > 0 ) mal.juv_t.i = mal.juv_t.i[ juv_t.i == 1 ] else mal.juv_t.i = NA
						
					} # if(reproduction.malaria == "complex")
					
					if(n.juv_t.i > 0) {	
						n.juv_t.i = n.juv[,,t][p.b] = rbinom( 1, n.juv_t.i, ( 1 - rat.f[p.b] ) )
						mal.juv_t.i = sample(mal.juv_t.i, n.juv_t.i)
					} else mal.juv_t.i = NA
					 
					if(n.1_t.i > 0) {	
						n.1_t.i = rbinom( 1, n.1_t.i, ( 1 - rat.s[p.b] ) ) 
						mal.1_t.i = sample(mal.1_t.i, n.1_t.i)
					} else mal.1_t.i = NA
					
					if(n.ad_t.i > 0) {	
						n.ad_t.i = rbinom( 1, n.ad_t.i, ( 1 - rat.s[p.b] ) )
						mal.ad_t.i = sample(mal.ad_t.i, n.ad_t.i)
					} else mal.ad_t.i = NA
					
						if(t.b[p.b] < 365) {
						# non-breeding season
						
							# first period
							
								#juveniles
								if(n.juv_t.i > 0) {
																			
									p.nb.1_juv = which(rmultinom(n.juv_t.i, 1, prob = disp.nb.1[[p.b]]) == 1) 
									p.nb.1_juv = p.nb.1_juv - seq(0, n.patch*(n.juv_t.i - 1), length.out = n.juv_t.i)
									
									get.m = ifelse(mal.juv_t.i == 0, rbinom( n.juv_t.i, 1, f.alpha(alpha.nb.1[p.nb.1_juv], t.nb.1[p.b]) ), 0)						
									surv.ac = ifelse(get.m == 1, rbinom(n.juv_t.i, 1, Sm.ac), 1)
									surv.juv = ifelse(surv.ac == 1, rbinom(n.juv_t.i, 1, s.juv.d[p.nb.1_juv]^t.nb.1[p.b] ), 0)
									mal.juv_t.i = ifelse(get.m == 0, mal.juv_t.i, 1)
									mal.juv_t.i = mal.juv_t.i[surv.juv == 1]
									p.nb.1_juv = p.nb.1_juv[surv.juv == 1]
									n.juv_t.i = length(mal.juv_t.i)
								
								} # if(n.juv_t.i > 0)
								
								#1-yr old
								if(n.1_t.i > 0) {
								
									p.nb.1_1 = which(rmultinom(n.1_t.i, 1, prob = disp.nb.1[[p.b]]) == 1) 
									p.nb.1_1 = p.nb.1_1 - seq(0, n.patch*(n.1_t.i - 1), length.out = n.1_t.i)
									
									get.m = ifelse(mal.1_t.i == 0, rbinom( n.1_t.i, 1, f.alpha(alpha.nb.1[p.nb.1_1], t.nb.1[p.b]) ), 0)
									surv.ac = ifelse(get.m == 1, rbinom(n.1_t.i, 1, Sm.ac), 1)
									surv.1 = ifelse(surv.ac == 1, rbinom(n.1_t.i, 1, s.ad.d[p.nb.1_1]^t.nb.1[p.b] ), 0)
									mal.1_t.i = ifelse(get.m == 0, mal.1_t.i, 1)
									mal.1_t.i = mal.1_t.i[surv.1 == 1]
									p.nb.1_1 = p.nb.1_1[surv.1 == 1]
									n.1_t.i = length(mal.1_t.i)				

								} # if(n.1_t.i > 0)
									
									#ad
									if(n.ad_t.i > 0) {
									
										p.nb.1_ad = which(rmultinom(n.ad_t.i, 1, prob = disp.nb.1[[p.b]]) == 1) 
										p.nb.1_ad = p.nb.1_ad - seq(0, n.patch*(n.ad_t.i - 1), length.out = n.ad_t.i)
										
										get.m = ifelse(mal.ad_t.i == 0, rbinom( n.ad_t.i, 1, f.alpha(alpha.nb.1[p.nb.1_ad], t.nb.1[p.b]) ), 0)
										surv.ac = ifelse(get.m == 1, rbinom(n.ad_t.i, 1, Sm.ac), 1)
										surv.ad = ifelse(surv.ac == 1, rbinom(n.ad_t.i, 1, s.ad.d[p.nb.1_ad]^t.nb.1[p.b] ), 0)
										mal.ad_t.i = ifelse(get.m == 0, mal.ad_t.i, 1)
										mal.ad_t.i = mal.ad_t.i[surv.ad == 1]
										p.nb.1_ad = p.nb.1_ad[surv.ad == 1]
										n.ad_t.i = length(mal.ad_t.i)	
										
									} # if(n.ad_t.i > 0)
							
							# second period
							
								#juveniles
								if(n.juv_t.i > 0) {
								
									p.nb.2_juv = rep(NA, n.juv_t.i)
									for(k in 1 : n.juv_t.i) p.nb.2_juv[k] = which( rmultinom( 1, 1, prob = disp.nb.2[[ p.nb.1_juv[k] ]] ) == 1 ) 
									
									get.m = ifelse(mal.juv_t.i == 0, rbinom( n.juv_t.i, 1, f.alpha(alpha.nb.2[p.nb.2_juv], t.nb.2[p.b]) ), 0)
									surv.ac = ifelse(get.m == 1, rbinom(n.juv_t.i, 1, Sm.ac), 1)
									surv.juv = ifelse(surv.ac == 1, rbinom(n.juv_t.i, 1, s.juv.d[p.nb.2_juv]^t.nb.2[p.b] ), 0)
									mal.juv_t.i = ifelse(get.m == 0, mal.juv_t.i, 1)
									mal.juv_t.i = mal.juv_t.i[surv.juv == 1]
									n.juv_t.i = length(mal.juv_t.i)
								
								}
								
								#1-yr old
								if(n.1_t.i > 0) {
								
									p.nb.2_1 = rep(NA, n.1_t.i)
									for(k in 1 : n.1_t.i) p.nb.2_1[k] = which( rmultinom( 1, 1, prob = disp.nb.2[[ p.nb.1_1[k] ]] ) == 1 ) 
									
									get.m = ifelse(mal.1_t.i == 0, rbinom( n.1_t.i, 1, f.alpha(alpha.nb.2[p.nb.2_1], t.nb.2[p.b]) ), 0)
									surv.ac = ifelse(get.m == 1, rbinom(n.1_t.i, 1, Sm.ac), 1)
									surv.1 = ifelse(surv.ac == 1, rbinom(n.1_t.i, 1, s.ad.d[p.nb.2_1]^t.nb.2[p.b] ), 0)
									mal.1_t.i = ifelse(get.m == 0, mal.1_t.i, 1)
									mal.1_t.i = mal.1_t.i[surv.1 == 1]
									n.1_t.i = length(mal.1_t.i)						
								
								}
									
								#ad
								if(n.ad_t.i > 0) {
									p.nb.2_ad = rep(NA, n.ad_t.i)
									for(k in 1 : n.ad_t.i) p.nb.2_ad[k] = which( rmultinom( 1, 1, prob = disp.nb.2[[ p.nb.1_ad[k] ]] ) == 1 ) 

									get.m = ifelse(mal.ad_t.i == 0, rbinom( n.ad_t.i, 1, f.alpha(alpha.nb.2[p.nb.2_ad], t.nb.2[p.b]) ), 0)
									surv.ac = ifelse(get.m == 1, rbinom(n.ad_t.i, 1, Sm.ac), 1)
									surv.ad = ifelse(surv.ac == 1, rbinom(n.ad_t.i, 1, s.ad.d[p.nb.2_ad]^t.nb.2[p.b] ), 0)
									mal.ad_t.i = ifelse(get.m == 0, mal.ad_t.i, 1)
									mal.ad_t.i = mal.ad_t.i[surv.ad == 1]
									n.ad_t.i = length(mal.ad_t.i)	
								}
								
							} # if(t.b[p.b] < 365) 

					} else {
					
						n.juv_t.i = n.1_t.i = n.ad_t.i = 0
						mal.juv_t.i = mal.1_t.i = mal.ad_t.i = NA
						
					}

				} # if(t > Tm)

				# natal and breeding dispersal, initiation of year t+1
				
				if(t < T & n.juv_t.i > 0) {	
				
					if(n.pairs_t.i < K.b[p.b] | psi.DD == 1) disp.juv = which(rmultinom(n.juv_t.i, 1, prob = disp.natal[[p.b]]) == 1) 
					if(n.pairs_t.i >= K.b[p.b] & psi.DD > 1) disp.juv = which(rmultinom(n.juv_t.i, 1, prob = disp.natal.DD[[p.b]]) == 1)
					disp.juv = disp.juv - seq(0, n.patch*(n.juv_t.i - 1), length.out = n.juv_t.i)	
					
					for(k in 1 : n.juv_t.i) { 
						n.1[,,(t+1)][ disp.juv[k] ] = n.1[,,(t+1)][ disp.juv[k] ] + 1
						mal.y[[ disp.juv[k] ]] = c( mal.y[[ disp.juv[k] ]], mal.juv_t.i[k] )
					}	

				}
				
				if(t < T & n.1_t.i > 0) {	
									
					disp.1 = which(rmultinom(n.1_t.i, 1, prob = disp.breed[[p.b]]) == 1)
					disp.1 = disp.1 - seq(0, n.patch*(n.1_t.i - 1), length.out = n.1_t.i)	
					
					for(k in 1 : n.1_t.i) { 
						n.ad[,,(t+1)][ disp.1[k] ] = n.ad[,,(t+1)][ disp.1[k] ] + 1
						mal.a[[ disp.1[k] ]] = c( mal.a[[ disp.1[k] ]], mal.1_t.i[k] )
					}
				
				}

				if(t < T & n.ad_t.i > 0) {	
				
					disp.ad = which(rmultinom(n.ad_t.i, 1, prob = disp.breed[[p.b]]) == 1)
					disp.ad = disp.ad - seq(0, n.patch*(n.ad_t.i - 1), length.out = n.ad_t.i)		
					
					for(k in 1 : n.ad_t.i) { 
						n.ad[,,(t+1)][ disp.ad[k] ] = n.ad[,,(t+1)][ disp.ad[k] ] + 1
						mal.a[[ disp.ad[k] ]] = c( mal.a[[ disp.ad[k] ]], mal.ad_t.i[k] )
					}
				
				}
									
		} # for(i in 1 : n.patch)
		
		if(t < T) {	
			n.pairs[,,(t+1)] = n.1[,,(t+1)] + n.ad[,,(t+1)]
			for(s in 1 : n.patch) { mal.y[[s]] = mal.y[[s]][-1]; mal.a[[s]] = mal.a[[s]][-1] }
			mal$y = mal.y; mal$a = mal.a
		}
					
	} #	for(t in 1 : T) 

#---------------------------------------------------------------------------------------------------------------- 

# OUTPUT

# Number of pairs at the end of the simulation (averaged over the last three years to 'counter' the effect of demographic stochasticity and over the nc patches)

np  = rep(NA, nr)
for(i in 1 : nr) { 
	np[i] =  sum(n.pairs[i,,T] + n.pairs[i,,T-1] + n.pairs[i,,T-2])
 }

np.metapop = round( (sum(np) / 3), 0)
np = round( (np / (3 * nc)), 0)

# Number of pairs for 3 different elevations, and various other statistics
grad = seq(grad[1], grad[2], length = nr) / 1000

np.high = np.mid = np.low = NA
K.nb.1.high = K.nb.1.mid = K.nb.1.low = K.nb.2.high = K.nb.2.mid = K.nb.2.low = NA
alpha.b.low = alpha.nb.1.low = alpha.nb.2.low = NA
	if(grad[1] >= 1.9 & grad[nr] <=1.1) {
		np.high = np[grad == 1.8]; np.mid = np[grad == 1.5]; np.low = np[grad == 1.2]
		K.nb.1.high = K.nb.1[grad == 1.8, 1]; K.nb.1.mid = K.nb.1[grad == 1.5, 1]; K.nb.1.low = K.nb.1[grad == 1.2, 1]
		K.nb.2.high = K.nb.2[grad == 1.8, 1]; K.nb.2.mid = K.nb.2[grad == 1.5, 1]; K.nb.2.low = K.nb.2[grad == 1.2, 1]
		alpha.b.low = alpha.b[grad == 1.2, 1]; alpha.nb.1.low = alpha.nb.1[grad == 1.2, 1]; alpha.nb.2.low = alpha.nb.2[grad == 1.2, 1]
	}
	
gr.ij  = matrix(NA, nr, nc)

for(i in 1 : nr) { 
	for(j in 1 : nc) { 
		gr.ij[i,j] =  sum( (n.pairs[i,j,T] / n.pairs[i,j,T-1]) + (n.pairs[i,j,T-1] / n.pairs[i,j,T-2]) + (n.pairs[i,j,T-2] / n.pairs[i,j,T-3]) ) / 3
	}
 }
 
gr = rep(NA, nr)
col.gr = c()
for(i in 1 : nr) { 
	gr[i] =  sum(gr.ij[i,])  / nc
		if(is.na(gr[i])) col.gr =  c(col.gr, "red") else {
			if(gr[i] >= 1) col.gr = c(col.gr, "purple")
			if(gr[i] < 1 & gr[i] >= 0.9) col.gr =  c(col.gr, "blue")
			if(gr[i] < 0.9 & gr[i] >= 0.8) col.gr =  c(col.gr, "darkorange")
			if(gr[i] < 0.8) col.gr =  c(col.gr, "red") 
			}

 }
 
gr.high = gr.mid = gr.low = NA
	if(grad[1] >= 1.9 & grad[nr] <=1.1) {
		gr.high = gr[grad == 1.8]; gr.mid = gr[grad == 1.5]; gr.low = gr[grad == 1.2]
	}
	
# mean elevation for np, K.nb.1 and K.nb.2, gr, and alphas

m.elev = mean(grad)
m.elev.np = sum(grad * np) / sum(np)
if(t.b[p.b] < 365) m.elev.K.nb.1 = sum(grad * K.nb.1[,1]) / sum(K.nb.1[,1]) else m.elev.K.nb.1 = NA
if(t.b[p.b] < 365) m.elev.K.nb.2 = sum(grad * K.nb.2[,1]) / sum(K.nb.2[,1]) else m.elev.K.nb.2 = NA
m.elev.gr = sum(grad * gr) / sum(gr)
m.elev.alpha.b = sum(grad * alpha.b[,1]) / sum(alpha.b[,1])
m.elev.alpha.nb.1 = sum(grad * alpha.nb.1[,1]) / sum(alpha.nb.1[,1])
m.elev.alpha.nb.2 = sum(grad * alpha.nb.2[,1]) / sum(alpha.nb.2[,1])
	
 # Plot of np = f(elevation) if no difference in patches carrying capacity
 
 if(length(unique(as.vector(K.b))) == 1) {
 
	 plot(grad, np, type = "n", xlab = "Elevation (km)", ylab = "# pairs", ylim = c(0, (1.1 * K.b[1,1])), pch = 19, cex = 2)
	 points(grad, np, pch = 19, cex = 2, col = col.gr)
	 
	 #lines(grad, np)
	 abline(a = K.b[1,1], b = 0, lty = 2, col = "black")
		 
 }
 
 # fitting a cline to figure
 
 if(add.cline == TRUE) {
 
	 a = min(np); b = max(np); c = mean(grad)
	 
	 f.cline = function(y = np, x = grad, a = a, b = b, c = c, K = 25){
		nls( y ~ a + (b - a) / (1 + exp( - K * (x - c) ) ), algorithm = "port", control = c(maxiter = 100, tol = 1e-05, minFactor = 1/1024), start = list(a = a, b = b, c = c, K = K))
	}

	cline.fit = f.cline(y = np, x = grad, a = a, b = b, c = c, K = 25)
	p = coef(cline.fit)

	x.sig = seq(max(grad), min(grad), length = 10000)
	sigmoid = function(a = p[1], b = p[2], c = p[3], K = p[4], x = x.sig) a + (b-a) / (1 + exp( - K * (x - c) ) )
	y.sig = sigmoid(a = p[1], b = p[2], c = p[3], K = p[4])
	lines(x.sig, y.sig)
	h.5 = mean(y.sig[x.sig > (p[3]-0.005) & x.sig < (p[3]+0.005)])
	
	y.sig.st = y.sig - min(y.sig) 
	y.sig.st = y.sig.st / max(y.sig.st)
	m.2 = mean(x.sig[y.sig.st > 0.195 & y.sig.st < 0.205])
	h.2 = mean(y.sig[x.sig > (m.2-0.005) & x.sig < (m.2+0.005)])
	m.8 = mean(x.sig[y.sig.st > 0.795 & y.sig.st < 0.805])
	h.8 = mean(y.sig[x.sig > (m.8-0.005) & x.sig < (m.8+0.005)])
	
	legend("topleft", c("GR +", "GR < 1", "GR <0.9", "GR <0.8 | NA"), pch = rep(19, 4), cex = 1.2, col = c("purple", "blue", "darkorange", "red"), text.col = "black", bg = "white")

	cline.5 = p[3]; cline.2 = m.2; cline.8 = m.8
	
} else {

cline.5 = cline.2 = cline.8 = NA

}

#---------------------------------------------------------------------------------------------------------------- 

	l = list(
		nr = nr, nc = nc, grad = grad, unit = unit, 
		T = T, Tm = Tm, t.b = t.b, t.nb.1 = t.nb.1, t.nb.2 = t.nb.2, 
		min.fledg = min.fledg, peak.fledg = peak.fledg, SD.fledg = SD.fledg,
		init.1 = init.1, init.2 = init.2, 
		s.ad = s.ad, s.ad.d = s.ad.d, rat.s = rat.s, s.juv = s.juv, s.juv.d = s.juv.d, 
		fec = fec, fec.1 = fec.1, rat.f = rat.f, K.b = K.b, thr.DD = thr.DD,
		K.nb.1 = K.nb.1, K.nb.2 = K.nb.2, K.nb.1.high = K.nb.1.high, K.nb.1.mid = K.nb.1.mid, K.nb.1.low = K.nb.1.low, K.nb.2.high = K.nb.2.high, K.nb.2.mid = K.nb.2.mid, K.nb.2.low = K.nb.2.low,
		reproduction.malaria = reproduction.malaria,
		alpha.b = alpha.b, alpha.b.low = alpha.b.low, alpha.nb.1 = alpha.nb.1, alpha.nb.1.low = alpha.nb.1.low, alpha.nb.2 = alpha.nb.2, alpha.nb.2.low = alpha.nb.2.low, Sm.ac = Sm.ac, 
		gamma.mov = gamma.mov, calc.gamma.d = calc.gamma.d, n.sim.disp = n.sim.disp, R.ter = R.ter, fidelity.ad = fidelity.ad, m.natal = m.natal, SD.natal = SD.natal, psi.DD = psi.DD, 
		add.cline = add.cline,
		gamma.d.breed = gamma.d.breed, test.disp.breed = test.disp.breed, 
		disp.breed = disp.breed, disp.natal = disp.natal, disp.natal.DD = disp.natal.DD, disp.nb.1 = disp.nb.1, disp.nb.2 = disp.nb.2, 
		mal.y = mal.y, mal.a = mal.a, 
		n.1 = n.1, n.ad = n.ad, n.juv = n.juv, n.pairs = n.pairs,
		np = np, np.high = np.high, np.mid = np.mid, np.low = np.low, np.metapop = np.metapop, gr = gr, gr.high = gr.high, gr.mid = gr.mid, gr.low = gr.low,
		r.hm = np.high /  np.mid,  r.ml = np.mid /  np.low,
		cline.5 = cline.5, cline.2 = cline.2, cline.8 = cline.8,
		m.elev = m.elev, m.elev.np = m.elev.np, m.elev.K.nb.1 = m.elev.K.nb.1, m.elev.K.nb.2 = m.elev.K.nb.2, m.elev.gr = m.elev.gr, 
		m.elev.alpha.b = m.elev.alpha.b, m.elev.alpha.nb.1 = m.elev.alpha.nb.1, m.elev.alpha.nb.2 = m.elev.alpha.nb.2
	)
	
	return(l)
		
	} # end of mamo
	
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
mamo = cmpfun(mamo)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
