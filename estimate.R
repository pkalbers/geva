# Allele age estimation from filtered pairs

rm.dupl = function(x) {
	k = sprintf("%d %d %d %d %d", x$SampleID0, x$Chr0, x$SampleID1, x$Chr1, x$Shared)
	if (any(duplicated(k))) {
		w = which(!duplicated(k))
		x = x[w,]
	}
	x
}

qc.tmrca = function(x, ne) {
	est = (x$Shape / x$Rate) * ne
	est = round(est/10)* 10
	
	con = which(x$Shared == 1)
	dis = which(x$Shared == 0)
	
	if (length(con) < 1) return(NULL)
	if (length(dis) < 1) return(NULL)
	
	if (max(est[con]) > min(est[dis])) {
		rng = sort(unique(est))
		rng = rng[-length(rng)] + (diff(rng) / 2)
		num = sapply(rng, function(y, c, d) { length(which(c > y)) + length(which(d < y)) }, est[con], est[dis])
		lim = rng[which.min(num)]
		
		wc = which(est[con] < lim)
		wd = which(est[dis] > lim)
		
		if (length(wc) == 0) { wc = which.min(est[con]) }
		if (length(wd) == 0) { wd = which.max(est[dis]) }
		
		con = con[wc]
		dis = dis[wd]
	}
	
	c(con, dis)
}

parse.data = function(d, ne) {
	d = lapply(d, function(x) {
		x = rm.dupl(x)
		x$Pass = 0
		i = qc.tmrca(x, ne)
		if (is.null(i)) { 
			return(NULL)
		}
		x$Pass[i] = 1
		x
	})
	d
}

com.post = function(xc, xd, t) {
	cp = rep(0, length(t))
	for (i in 1:length(t)) {
		cc = sum(pgamma(t[i], shape = xc$Shape, rate = xc$Rate, log.p = T, lower.tail = T) )
		cd = sum(pgamma(t[i], shape = xd$Shape, rate = xd$Rate, log.p = T, lower.tail = F) )
		cp[i] = cc + cd
	}
	cp
}

com.post.est = function(xc, xd, lower = (0.01/ne), upper = (1e8/ne), last = 0, n = 101, lim = 1/ne) {
	if (upper - lower < lim) {
		return(last)
	}
	
	t = exp(seq(log(lower), log(upper), length.out = n))
	
	cp = com.post(xc, xd, t)
	
	wx = which.max(cp)
	lo = max(1, wx-2)
	up = min(n, wx+2)
	
	return(com.post.est(xc, xd, t[lo], t[up], t[wx]))
}

run.age.est = function(d, ne) {
	cp = lapply(d, function(x, ne) {
		con = which(x$Shared == 1 & x$Pass == 1)
		dis = which(x$Shared == 0 & x$Pass == 1)
		
		if (length(con) < 1) return(NULL)
		if (length(dis) < 1) return(NULL)
		
		xc = data.frame(Shape = x$Shape[con], Rate = x$Rate[con])
		xd = data.frame(Shape = x$Shape[dis], Rate = x$Rate[dis])
		
		data.frame(MarkerID = x$MarkerID[1],
							 Clock = x$Clock[1],
							 N_Concordant = length(con),
							 N_Discordant = length(dis),
							 PostMode  = com.post.est(xc, xd) * ne)
	}, ne)
	cp
}



args = commandArgs(T)


if (length(args) == 0) {
	cat("Input arguments not specified\n")
	quit(save = "no")
}

if (length(args) != 2) {
	cat("Unexpected number of arguments; 2 required\n")
	quit(save = "no")
}

if (!file.exists(args[1])) {
	cat("Input file does not exist\n")
	quit(save = "no")
}

if (is.na(as.numeric(args[2]))) {
	cat(sprintf("Unable to interpret scaling parameter (Ne): %s\n", args[2]))
	quit(save = "no")
}


cat(sprintf("Input file (pairs)     : %s\n", args[1]))
cat(sprintf("Scaling parameter (Ne) : %s\n", args[2]))

file = args[1]
ne = as.numeric(args[2]) * 2


cat("Reading pairs ... ")
d = read.table(file, header = T, stringsAsFactors = F)
cat("OK\n")


cat("Parsing data ... ")
d = split(d, list(d$MarkerID, d$Clock))
d = parse.data(d, ne)
cat("OK\n")

cat("Estimating age ... ")
p = run.age.est(d, ne)
cat("OK\n")


outfile = sub("^(.+)\\.(pairs\\.txt)$", "\\1.sites2.txt", file)

cat(sprintf("Writing to output file:  %s\n", outfile))
p = Reduce(rbind, p)
write.table(p, file = outfile, append = F, quote = F, sep = " ", row.names = F, col.names = T)
cat("Done\n")
