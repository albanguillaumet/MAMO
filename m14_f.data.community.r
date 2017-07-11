f.data.community = function(
AKIP = AKIP, AKEP = AKEP, HCRE = HCRE, IIWI = IIWI, 
ELEP = ELEP, OMAO = OMAO, APAP = APAP, HAAM = HAAM,
var_ = c("mg.hab", "risk", "mg.risk", "ac", "mg.rat", "mg.res.1", "mg.res.2", "envp", "sim")
){

nr = dim(AKIP)[1]
nv = length(var_)
d = as.data.frame( matrix( NA, nr = nr, nc = (nv + 3 + 8) ) )
for(i in 1 : nv) d[, i] = AKIP[, var_[i]]
colnames(d) = c(var_, "n.w1", "n.w2", "n.w3", "AKIP", "AKEP", "HCRE", "IIWI", "ELEP", "OMAO", "APAP", "HAAM")

for(i in 1 : nr) {

	d$n.w1[i] = ( (AKIP$w1[i] * AKIP$n.st[i]) + (AKEP$w1[i] * AKEP$n.st[i]) + (HCRE$w1[i] * HCRE$n.st[i]) + (IIWI$w1[i] * IIWI$n.st[i]) +
					   (ELEP$w1[i] * ELEP$n.st[i]) + (OMAO$w1[i] * OMAO$n.st[i]) + (APAP$w1[i] * APAP$n.st[i]) + (HAAM$w1[i] * HAAM$n.st[i])  )  /
					   ( AKIP$w1[i] + AKEP$w1[i] + HCRE$w1[i] + IIWI$w1[i] + ELEP$w1[i] + OMAO$w1[i] + APAP$w1[i] + HAAM$w1[i] ) 
					   
	d$n.w2[i] = ( (AKIP$w2[i] * AKIP$n.st[i]) + (AKEP$w2[i] * AKEP$n.st[i]) + (HCRE$w2[i] * HCRE$n.st[i]) + (IIWI$w2[i] * IIWI$n.st[i]) +
					   (ELEP$w2[i] * ELEP$n.st[i]) + (OMAO$w2[i] * OMAO$n.st[i]) + (APAP$w2[i] * APAP$n.st[i]) + (HAAM$w2[i] * HAAM$n.st[i])  )  /
					   ( AKIP$w2[i] + AKEP$w2[i] + HCRE$w2[i] + IIWI$w2[i] + ELEP$w2[i] + OMAO$w2[i] + APAP$w2[i] + HAAM$w2[i] ) 		

	d$n.w3[i] = ( (AKIP$w3[i] * AKIP$n.st[i]) + (AKEP$w3[i] * AKEP$n.st[i]) + (HCRE$w3[i] * HCRE$n.st[i]) + (IIWI$w3[i] * IIWI$n.st[i]) +
					   (ELEP$w3[i] * ELEP$n.st[i]) + (OMAO$w3[i] * OMAO$n.st[i]) + (APAP$w3[i] * APAP$n.st[i]) + (HAAM$w3[i] * HAAM$n.st[i])  )  /
					   ( AKIP$w3[i] + AKEP$w3[i] + HCRE$w3[i] + IIWI$w3[i] + ELEP$w3[i] + OMAO$w3[i] + APAP$w3[i] + HAAM$w3[i] ) 			

}

d$AKIP = AKIP$n.st; d$AKEP = AKEP$n.st; d$HCRE = HCRE$n.st; d$IIWI = IIWI$n.st
d$ELEP = ELEP$n.st; d$OMAO = OMAO$n.st; d$APAP = APAP$n.st; d$HAAM = HAAM$n.st

d
}
