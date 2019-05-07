#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# New version of f.plot.composite (km instead of Km)

f.plot.composite = function(
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

if(create.plot == TRUE) plot(x, rep(0, length(x)), xlim = c(min(x), max(x)), ylim = c(0, max.y), type = "n", bty = "n", xlab = "Elevation (km)", ylab = ylab, main = main, cex.main = 0.7)

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

} # end of f.plot.composite

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
f.plot.composite = cmpfun(f.plot.composite)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
