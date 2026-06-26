# Copyright (c) 2019 Lars Fritsche
# QQ plot Rscript

options(stringsAsFactors=F)
library("optparse")
library("data.table")

stopf <- function(...) stop(sprintf(...), call. = FALSE)

pToLog10 <- function(P, column = "PVALUE") {
	text <- trimws(as.character(P))
	num <- suppressWarnings(as.numeric(text))
	log10P <- rep(NA_real_, length(text))
	valid <- !is.na(num) & is.finite(num) & num > 0 & num <= 1
	log10P[valid] <- -log10(num[valid])

	underflow <- !is.na(num) & num == 0 &
		grepl("^\\+?([0-9]*\\.?[0-9]+)[Ee]-([0-9]+)$", text)
	if(any(underflow)) {
		mantissa <- as.numeric(sub("^\\+?([0-9]*\\.?[0-9]+)[Ee]-([0-9]+)$", "\\1", text[underflow]))
		exponent <- as.numeric(sub("^\\+?([0-9]*\\.?[0-9]+)[Ee]-([0-9]+)$", "\\2", text[underflow]))
		log10P[underflow] <- exponent - log10(mantissa)
	}

	invalid <- which(is.na(log10P) | !is.finite(log10P))
	if(length(invalid) > 0) {
		stopf("Invalid p-values in column '%s': expected 0 < P <= 1. First invalid row(s): %s",
			column, paste(head(invalid, 10), collapse = ","))
	}
	return(log10P)
}

validateLog10P <- function(values, column = "PVALUE") {
	log10P <- suppressWarnings(as.numeric(values))
	invalid <- which(is.na(log10P) | !is.finite(log10P) | log10P < 0)
	if(length(invalid) > 0) {
		stopf("Invalid -log10(P) values in column '%s': expected finite values >= 0. First invalid row(s): %s",
			column, paste(head(invalid, 10), collapse = ","))
	}
	return(log10P)
}

normalizeAlleleFrequency <- function(values, allele_frequency_type = "auto") {
	allele_frequency_type <- tolower(allele_frequency_type)
	if(!(allele_frequency_type %in% c("auto", "maf", "eaf"))) {
		stopf("--allele_frequency_type must be one of auto, MAF, or EAF")
	}

	freq <- suppressWarnings(as.numeric(values))
	invalid <- which(is.na(freq) | !is.finite(freq) | freq < 0 | freq > 1)
	if(length(invalid) > 0) {
		stopf("Invalid allele frequencies in MAF column: expected finite values in [0, 1]. First invalid row(s): %s",
			paste(head(invalid, 10), collapse = ","))
	}

	over_half <- which(freq > 0.5)
	if(allele_frequency_type == "maf" && length(over_half) > 0) {
		stopf("MAF column contains values > 0.5. Use --allele_frequency_type EAF or auto if this column is allele frequency. First row(s): %s",
			paste(head(over_half, 10), collapse = ","))
	}
	if(allele_frequency_type %in% c("auto", "eaf") && length(over_half) > 0) {
		warning(sprintf("Converting %s allele frequencies > 0.5 to minor allele frequency", length(over_half)),
			call. = FALSE)
		freq[over_half] <- 1 - freq[over_half]
	}
	return(freq)
}

# QQ plot function
qqplotdata <- function(logpvector){
	require(data.table)
    o = sort(logpvector,decreasing=T)
    e = -log10(ppoints(length(o)))       
    qqdata <- data.table(o,e)
	# thinning
    qqdata[,o:=round(o,3)]
    qqdata[,e:=round(e,3)]
    keepU <- which(!duplicated(qqdata))
    qqdata <- qqdata[keepU,]
    
    N <- length(logpvector) ## number of p-values
    ## create the confidence intervals
    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)
      
    qqdata[,c975:=sapply(keepU,function(x) -log10(qbeta(0.975,x,N-x+1)))]
    qqdata[,c025:=sapply(keepU,function(x) -log10(qbeta(0.025,x,N-x+1)))]
    return(qqdata)
}

QQplot <- function(res,maintitle="",top.size=0.125,break.top=15,
	maf="MAF",pvalue="PVALUE",log10p=F,sigthreshold="5E-8",DTthreads=1,
	allele_frequency_type="auto"){
	options(stringsAsFactors=F)
	require("plotrix")
	require("data.table")
	require("RColorBrewer")
	setDTthreads(DTthreads)
	if(!is.data.table(res)) res <- as.data.table(res)

	#check columns
	setnames(res,c(maf,pvalue),c("MAF","PVALUE"),skip_absent=T)
	missing <- setdiff(c("MAF","PVALUE"), names(res))
	if(length(missing) > 0) stopf("Missing required column(s): %s", paste(missing, collapse = ", "))
	
	# horizontal lines and corresponding colors
	yLine <- c(-log10(sort(as.numeric(unlist(strsplit(sigthreshold,","))))))
	colLine <- c("red")

	res[,MAF:=normalizeAlleleFrequency(MAF, allele_frequency_type)]
	if(!log10p){
		res[,LOG10P:=pToLog10(PVALUE)]
	} else {
		res[,LOG10P:=validateLog10P(PVALUE)]
	}
	res <- na.omit(res[,.(MAF,LOG10P)])
	res <- res[!is.infinite(LOG10P),]
	if(nrow(res) == 0) stopf("No valid rows remain after QQ input validation")

	minMAF <- min(res$MAF)

	# Determine frequency bins and create variable for binned QQ plot
	freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
	res[,freqbin:=cut(MAF, freqbins,include.lowest=T)]
	freqtable <- table(res$freqbin)
	freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
	freqtable <- freqtable[freqtable > 0]

	## Generate QQ plot data by frequency bin
	fbin <- character(0)
	fN <- integer(0)
	fx <- numeric(0)
	fy <- numeric(0)
	fcol <- character(0)
	legendcol <- character(0)
	conf <- list()
	allcols <- brewer.pal(4,"Set1")
	for(f in 1:length(freqtable)){
		fbin <- c(fbin,names(freqtable)[f])	
		plotdata <- qqplotdata(res[freqbin == names(freqtable)[f],LOG10P])
		fN <- c(fN,freqtable[f])
		fx <- c(fx,plotdata$e)
		fy <- c(fy,plotdata$o)
		fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
		conf[[f]] <- data.table('x'=c(plotdata$e,rev(plotdata$e)),
								'y'=c(plotdata$c975,rev(plotdata$c025)))
		legendcol <- c(legendcol,allcols[f])
	}
	legendtext <- paste0("MAF=",fbin,"; N SNPs=",format(fN,big.mark=",",scientific=FALSE))

	## QQ plot by binned frequencies
	xlim <- c(0,max(fx,na.rm=T))
	ylim <- c(0,max(fy,na.rm=T))
	maxY <- max(fy,na.rm=T)
	
	par(mar=c(5.1,5.1,4.1,1.1))
	# plot version with two axes
	if(maxY > break.top/(1 - top.size)){
		# create pretty y-axis labels
		lab1 <- pretty(c(0,break.top),n=ceiling(12 * (1-top.size)))
		lab1 <- c(lab1[lab1 < break.top],break.top)
		lab2 <- pretty(c(break.top,maxY),n=max(3,floor(12 * top.size)))
		lab2 <- lab2[lab2 > max(lab1)]

		# resulting range of top scale in bottom scale units
		top.range = break.top/(1 - top.size) - break.top
		top.data = max(lab2)-break.top
	
		# function to rescale the top part
		rescale = function(y) { break.top+(y-break.top)/(top.data/top.range)}
		rescaled.y = rescale(fy[fy>break.top])
		plot(0,0,
			ylim=c(min(fy),break.top*(1+top.size)),xlim=xlim,axes=FALSE,
			xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
			ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
			cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
			main=maintitle,pch=19)
	
		# Plot confidence intervals	
		for(p in 1:length(conf)){
			polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
				col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
				border = NA)
		}

		# add points
		points(fx[fy<break.top],fy[fy<break.top],cex=1,col=fcol[fy<break.top],pch=19)

		# identify line & add axis break
		lines(xlim,xlim,col="black",lty = 2)
		axis(1,cex.axis=1.5,cex.lab=1.5)
		par(las=1)
		axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
		par(las=0)
		box()
		par(las=0)
		points(fx[fy>break.top],rescaled.y,cex=1,col=fcol[fy>break.top],pch=19)
		par(las=1)
		axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
		axis.break(axis=2,breakpos=break.top,style="zigzag",brw=0.02)
		axis.break(axis=4,breakpos=break.top,style="zigzag",brw=0.02)
		lines(range(fx),c(break.top,break.top),col = "grey",lty = 6)
		abline(h=ifelse(yLine<break.top,
			yLine,
			rescale(yLine)),
			col=colLine,lwd=1.5,lty=2)
		legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
	# plot version with single y axes
	} else {
		par(mar=c(5.1,5.1,4.1,1.1),las=1)
		axislim <- ceiling(range(xlim,ylim,yLine))
		plot(0,0,
			ylim=axislim,xlim=xlim,axes=T,
			xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
			ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
			cex=1,cex.lab=1.5,cex.axis=1.5,col="transparent",
			main=maintitle,pch=19)
		# Plot confidence intervals
		for(p in 1:length(conf)){
				polygon(conf[[p]]$'x',conf[[p]]$'y',
					col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
					border = NA)
		}
		points(fx,fy,col=fcol,pch=19)
		# identity line & genome-wide significance line
		lines(axislim,axislim,col = "grey",lwd=1.5,lty=2)
		abline(h=yLine,col=colLine,lwd=1.5,lty=2)
		legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
	}
}

option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited; required columns: 'MAF' and 'PVALUE'"),   
  make_option("--prefix", type="character", default="",
    help="Prefix of output files"),   
  make_option("--top.size", type="numeric", default=0.125,
    help="top size = proportion of total length y axis [default=0.125]"),
  make_option("--break.top", type="numeric", default=15,
    help="set axis break at -log10(P) [default=15]"),
  make_option("--width", type="numeric", default=900,
    help="Width QQ plot in pixel [default=900]"),
  make_option("--height", type="numeric", default=900,
    help="Height QQ plot in pixel [default=900]"),
  make_option("--pointsize", type="numeric", default=16,
    help="Point size of plots [default=16]"),
  make_option("--maf", type="character", default="MAF",
    help="name of column with MAF [default='MAF']"),
  make_option("--pvalue", type="character", default="PVALUE",
    help="name of column with p.value [default='PVALUE']"),
  make_option("--log10p", type="logical", default=F,
    help="Input p.value column with -log10(p.value) [default=F]"),    
  make_option("--allele_frequency_type", type="character", default="auto",
    help="Interpret MAF column as auto, MAF, or EAF; auto/EAF fold values >0.5 [default='auto']"),
  make_option("--maintitle", type="character", default="",
    help="Plot title"),
  make_option("--threads", type="numeric", default=1,
    help="DTthreads") 
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

setDTthreads(opt$threads)

gwas <- fread(opt$input,select=c(opt$maf,opt$pvalue),
	colClasses=setNames("character", opt$pvalue))

png(filename = paste0(opt$prefix,"_QQ.png"), width=opt$width,height=opt$height,pointsize=opt$pointsize)
tryCatch({
	QQplot(res=gwas,top.size=opt$top.size,break.top=opt$break.top,
		maf=opt$maf,pvalue=opt$pvalue,log10p=opt$log10p,sigthreshold="5E-8",
		maintitle=opt$maintitle,DTthreads=opt$threads,
		allele_frequency_type=opt$allele_frequency_type)
}, error=function(e) {
	dev.off()
	stop(e)
})
dev.off()
