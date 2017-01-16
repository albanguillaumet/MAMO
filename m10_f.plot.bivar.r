#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO PLOT SIMULATION RUN (TWO VARIABLES SIMULTANEOUSLY) 

f.plot.bivar = function(
output.dir = "C:/Programs/MAMO/RUN/IIWI.1_1", 
var.excl = NA, val.excl = NA,
take.subset = TRUE, val.RISK = NA, val.AC = NA, val.RAT = 1, val.RES.1 = 1, val.RES.2 = 1,
y = "np.high",  ylab = "# pairs", main = "Elevation = 1800 m",
x = "AC", xlab = "Malaria mortality", 
o = "RISK", title.o = "RISK", legend.o = TRUE, l.pos = NA
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

di = d.f

y.max = max(di[,y])

if(take.subset == TRUE) {
	if(!is.na(val.RISK)) di = di[di$RISK == val.RISK,]
	if(!is.na(val.AC)) di = di[di$AC == val.AC,]
	if(!is.na(val.RAT)) di = di[di$RAT.S == val.RAT,]
	if(!is.na(val.RES.1)) di = di[di$RES.1 == val.RES.1,]
	if(!is.na(val.RES.2)) di = di[di$RES.2 == val.RES.2,]
}

 cn = colnames(di) 
 
#----------

# Figure title

if(is.na(main)) main.o = ""
if(!is.na(main) & (x == "RISK" | o == "RISK") & (x == "AC" | o == "AC")) main.o = paste(main, "|", "RAT", "=", val.RAT, ",", "RES.1", "=", val.RES.1, ",", "RES.2", "=", val.RES.2)
if(!is.na(main) & (x == "RISK" | o == "RISK") & (x == "RAT" | o == "RAT"  | x == "RAT.S" | o == "RAT.S" | x == "RAT.F" | o == "RAT.F")  ) main.o = paste(main, "|", "AC", "=", val.AC, ",", "RES.1", "=", val.RES.1, ",", "RES.2", "=", val.RES.2)
if(!is.na(main) & (x == "RISK" | o == "RISK") & (x == "RES.1" | o == "RES.1")  ) main.o = paste(main, "|", "AC", "=", val.AC, ",", "RAT", "=", val.RAT, ",", "RES.2", "=", val.RES.2)
if(!is.na(main) & (x == "RISK" | o == "RISK") & (x == "RES.2" | o == "RES.2")  ) main.o = paste(main, "|", "AC", "=", val.AC, ",", "RAT", "=", val.RAT, ",", "RES.1", "=", val.RES.1)
if(!is.na(main) & (x == "AC" | o == "AC") & (x == "RAT" | o == "RAT"  | x == "RAT.S" | o == "RAT.S" | x == "RAT.F" | o == "RAT.F")  ) main.o = paste(main, "|", "RISK", "=", val.RISK, ",","RES.1", "=", val.RES.1, ",", "RES.2", "=", val.RES.2)
if(!is.na(main) & (x == "AC" | o == "AC") & (x == "RES.1" | o == "RES.1")  ) main.o = paste(main, "|", "RISK", "=", val.RISK, ",", "RAT", "=", val.RAT, ",", "RES.2", "=", val.RES.2)
if(!is.na(main) & (x == "AC" | o == "AC") & (x == "RES.2" | o == "RES.2")  ) main.o = paste(main, "|", "RISK", "=", val.RISK, ",", "RAT", "=", val.RAT, ",", "RES.1", "=", val.RES.1)
if(!is.na(main) & (x == "RAT" | o == "RAT"  | x == "RAT.S" | o == "RAT.S" | x == "RAT.F" | o == "RAT.F")  & (x == "RES.1" | o == "RES.1")  ) main.o = paste(main, "|", "RISK", "=", val.RISK, ",", "AC", "=", val.AC, ",", "RES.2", "=", val.RES.2)
if(!is.na(main) & (x == "RAT" | o == "RAT"  | x == "RAT.S" | o == "RAT.S" | x == "RAT.F" | o == "RAT.F")  & (x == "RES.2" | o == "RES.2")  ) main.o = paste(main, "|", "RISK", "=", val.RISK, ",", "AC", "=", val.AC, ",", "RES.1", "=", val.RES.1)
if(!is.na(main) & (x == "RES.1" | o == "RES.1")  & (x == "RES.2" | o == "RES.2")  ) main.o = paste(main, "|", "RISK", "=", val.RISK, ",", "AC", "=", val.AC, ",", "RAT", "=", val.RAT)

#------------------
 
 if(is.na(l.pos)) {
 
	raw.means.plot(di, col.offset = which(cn == o), col.x = which(cn == x), col.value = which(cn == y), 
	pt.cex = 1.8, ylim = c(0, y.max), xlab = xlab , ylab = ylab, main = main.o, cex.main = 0.7, legend = legend.o, title = title.o)
	
 } else {

	 raw.means.plot(di, col.offset = which(cn == o), col.x = which(cn == x), col.value = which(cn == y), 
	 pt.cex = 1.8, ylim = c(0, y.max), xlab = xlab , ylab = ylab, main = main.o, cex.main = 0.7, legend = legend.o, title = title.o, l.pos = l.pos)
	 
 }
 
ENVP = di[, "ENVP"]
SIM = di[, "SIM"]

d.subset = di

print( summary( lmer( d.subset[, y] ~ d.subset[, x] * d.subset[, o] + (1 | ENVP) + (1 | SIM) ) ) )
print( summary( lmer( d.subset[, y] ~ d.subset[, x] + d.subset[, o] + (1 | ENVP) + (1 | SIM) ) ) )

list( d.subset = d.subset, test.mamo = max(d$test.disp.breed) )

 }

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

library(compiler)
f.plot.bivar = cmpfun(f.plot.bivar)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
