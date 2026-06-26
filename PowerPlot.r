options(stringsAsFactors=F)

library("optparse")
library("RColorBrewer")
library("pwr")

stopf <- function(...) stop(sprintf(...), call. = FALSE)

parsePositiveIntegers <- function(value, option_name) {
	num <- suppressWarnings(as.integer(strsplit(value,",")[[1]]))
	if(length(num) == 0 || any(is.na(num)) || any(num <= 0)) {
		stopf("%s must contain comma-separated positive integers", option_name)
	}
	return(num)
}

parsePositiveNumerics <- function(value, option_name) {
	num <- suppressWarnings(as.numeric(strsplit(value,",")[[1]]))
	if(length(num) == 0 || any(is.na(num)) || any(!is.finite(num)) || any(num <= 0)) {
		stopf("%s must contain comma-separated positive numeric values", option_name)
	}
	return(num)
}

# Power Analysis
getPower <- function(ncases,ncontrols,alpha,raf.controls){
	# Model two allele-count proportions with potentially unequal case/control sizes.
	h <- pwr.2p2n.test(h = NULL, n1 = 2 * ncases, n2 = 2 * ncontrols,
		sig.level = alpha, power = 0.8, alternative = "two.sided")$h
	
	if(ncases < ncontrols) {
		h <- -h
	 	raf.cases <- sapply(raf.controls,function(p1) (sin(asin(sqrt(p1)) - h/2))^2)	
	} else {
		raf.cases <- sapply(raf.controls,function(p2) (sin(h/2 + asin(sqrt(p2))))^2)	
	}
	raf <- (raf.controls* ncontrols + raf.cases * ncases) / (ncases + ncontrols)
	OR <- (raf.cases / (1-raf.cases)) / (raf.controls / (1-raf.controls))

	maf <- raf
	maf[which(maf > 0.5)] <- 1 - maf[which(maf > 0.5)]	
	h_check <- ES.h(raf.cases,raf.controls)
	
	# some cohen's d values of certain allele frequencies might be crazy
	# only report 80% power if MAC > 0
	rac.controls <- raf.controls * 2 * ncontrols
	rac.cases <- raf.cases * 2 * ncases

	notok <- which(signif(h_check,3) != signif(median(h_check),3)
		| (floor(rac.controls) == 0 & floor(rac.cases) == 0)
		| (ceiling(rac.controls) == 2*ncontrols & ceiling(rac.cases) == 2 * ncases))

	out <- data.frame(raf.controls,raf.cases,raf,OR,maf)
	out[notok,-1] <- NA
	return(out)
}

getQuantPower <- function(nsamples,alpha,af,covariates,trait_sd){
	v <- nsamples - covariates - 2
	if(v <= 0) {
		stopf("--samples must be greater than --covariates + 2")
	}
	f2 <- pwr.f2.test(u = 1, v = v, f2 = NULL, sig.level = alpha, power = 0.8)$f2
	r2 <- f2 / (1 + f2)
	beta <- sqrt(r2 / (2 * af * (1 - af))) * trait_sd
	maf <- af
	maf[which(maf > 0.5)] <- 1 - maf[which(maf > 0.5)]
	return(data.frame(af,beta,maf))
}

# function to scale x axis
transx <- function(x){
    if (is.na(x)){
    	lx <- NA
    } else if (x > 0.5){
        lx <- (abs(log(1-x))-abs(log(0.5)))+log(0.5)
    } else {
        lx <- log(x)
    }
    return(lx)
}

option_list <- list(
  make_option("--prefix", type="character", default="",
    help="Prefix of output files"),   
  make_option("--trait_type", type="character", default="binary",
    help="Trait type: binary or quantitative [default='binary']"),
  make_option("--cases", type="character", default="100,1000,10000",
    help="Number(s) of cases; comma-separated [default='100,1000,10000']"),        
  make_option("--controls", type="character", default="100,1000,10000",
    help="Number(s) of controls; comma-separated [default='100,1000,10000']"),        
  make_option("--samples", type="character", default="100,1000,10000",
    help="Number(s) of samples for quantitative traits; comma-separated [default='100,1000,10000']"),
  make_option("--covariates", type="integer", default=0,
    help="Number of covariates for quantitative-trait power [default=0]"),
  make_option("--trait_sd", type="numeric", default=1,
    help="Trait standard deviation for quantitative-trait beta scale [default=1]"),
  make_option("--minMAF", type="numeric", default=0.0005,
    help="minimal minor frequency threshold [default=0.0005]"),
  make_option("--alpha", type="numeric", default=5E-8,
    help="alpha [default=5E-8]"),
  make_option("--risksnps", type="character", default="",
    help="File with risk SNPs for binary traits; columns RAF,RAF.CONTROLS and OR, tab-delimited [default='']"),
  make_option("--raf.controls", type="logical", default=T,
    help="Plot RAF in controls instead of overall RAF for binary traits [default=T]"),
  make_option("--width", type="numeric", default=900,
    help="Plot width in pixel [default=900]"),
  make_option("--height", type="numeric", default=900,
    help="Plot height in pixel [default=900]"),
  make_option("--pointsize", type="numeric", default=16,
    help="Point size of plots [default=16]"),
  make_option("--stepsize", type="numeric", default=0.00001,
    help="Stepwise frequency increment / smoothness of plot [default=0.00001]"),
  make_option("--ytix", type="character", default="1,1.2,1.5,2,3,5,10",
    help="Y-axis ticks for binary-trait odds ratios [default='1,1.2,1.5,2,3,5,10']"),
  make_option("--beta_tix", type="character", default="0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10",
    help="Y-axis ticks for quantitative-trait beta values [default='0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10']")
)
parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

trait_type <- tolower(opt$trait_type)
if(!trait_type %in% c("binary","quantitative")) stopf("--trait_type must be binary or quantitative")
if(!is.finite(opt$alpha) || opt$alpha <= 0 || opt$alpha >= 1) stopf("--alpha must be in (0, 1)")
if(!is.finite(opt$minMAF) || opt$minMAF <= 0 || opt$minMAF >= 0.5) stopf("--minMAF must be in (0, 0.5)")
if(!is.finite(opt$stepsize) || opt$stepsize <= 0 || opt$stepsize >= 1) stopf("--stepsize must be in (0, 1)")

if(trait_type == "binary") {
	cases <- parsePositiveIntegers(opt$cases, "--cases")
	controls <- parsePositiveIntegers(opt$controls, "--controls")
	if(length(cases) != length(controls)) stopf("--cases and --controls must contain the same number of values")
	ytix <- parsePositiveNumerics(opt$ytix, "--ytix")
	if(length(cases) > 24) stopf("At most 24 case/control pairs can be plotted")
	samplesizes <- data.frame('cases'=cases, 'controls'=controls)
} else {
	samples <- parsePositiveIntegers(opt$samples, "--samples")
	if(is.na(opt$covariates) || opt$covariates < 0) stopf("--covariates must be a non-negative integer")
	if(!is.finite(opt$trait_sd) || opt$trait_sd <= 0) stopf("--trait_sd must be positive")
	beta_tix <- parsePositiveNumerics(opt$beta_tix, "--beta_tix")
	if(length(samples) > 24) stopf("At most 24 sample sizes can be plotted")
	if(any(samples <= opt$covariates + 2)) stopf("--samples must be greater than --covariates + 2")
	if(opt$risksnps != "") stopf("--risksnps is supported only for --trait_type binary")
	samplesizes <- data.frame('samples'=samples)
}

plotdata <- list()

freqs <- seq(0,1,opt$stepsize)

# Maximal 24 unique line/color combinations. Keep the first palette cycle solid
# to avoid dotted-line artifacts on dense transformed curves.
allcols <- rep(brewer.pal(8,"Dark2"),3)
alllines <- rep(c(1,2,5),each=8)

xrange <- c(opt$minMAF,1-opt$minMAF)

for(i in 1:dim(samplesizes)[1]){
	if(trait_type == "binary") {
		minMAF <- 1 / (2*(samplesizes$cases[i] + samplesizes$controls[i]))
	} else {
		minMAF <- 1 / (2*samplesizes$samples[i])
	}
	plotfreqs <- c(minMAF,freqs[which(freqs>minMAF & freqs < 1-minMAF)],1-minMAF)	
	if(trait_type == "binary") {
		plotdata[[i]] <- getPower(ncases=samplesizes$cases[i],ncontrols=samplesizes$controls[i],
						alpha=opt$alpha,plotfreqs)
	} else {
		plotdata[[i]] <- getQuantPower(nsamples=samplesizes$samples[i],alpha=opt$alpha,
						af=plotfreqs,covariates=opt$covariates,trait_sd=opt$trait_sd)
	}
}

if(trait_type == "binary" & opt$risksnps != "" & !file.exists(opt$risksnps)) {
	stopf("Risk SNP file not found: %s", opt$risksnps)
}
if(trait_type == "binary" & opt$risksnps != ""){
	risksnps <- read.table(opt$risksnps,sep="\t",header=T,comment.char="")
	if( !all(c("OR","RAF","RAF.CONTROLS") %in% colnames(risksnps))) stop("Check columns in risk SNP file")
	risksnps$OR <- suppressWarnings(as.numeric(risksnps$OR))
	risksnps$RAF <- suppressWarnings(as.numeric(risksnps$RAF))
	risksnps$RAF.CONTROLS <- suppressWarnings(as.numeric(risksnps$RAF.CONTROLS))
	if(any(is.na(risksnps$OR) | !is.finite(risksnps$OR) | risksnps$OR <= 0) ||
			any(is.na(risksnps$RAF) | !is.finite(risksnps$RAF) | risksnps$RAF < 0 | risksnps$RAF > 1) ||
			any(is.na(risksnps$RAF.CONTROLS) | !is.finite(risksnps$RAF.CONTROLS) |
				risksnps$RAF.CONTROLS < 0 | risksnps$RAF.CONTROLS > 1)) {
		stopf("Risk SNP columns OR, RAF, and RAF.CONTROLS must be numeric; OR > 0 and frequencies in [0, 1]")
	}
	risksnps$RAF <- ifelse(risksnps$OR < 1,1 - risksnps$RAF,risksnps$RAF)
	risksnps$RAF.CONTROLS <- ifelse(risksnps$OR < 1,1 - risksnps$RAF.CONTROLS,risksnps$RAF.CONTROLS)	
	risksnps$OR <- ifelse(risksnps$OR < 1,1/risksnps$OR,risksnps$OR)	
}

png(paste0(opt$prefix,"_PowerAnalysis.png"),width=opt$width,height=opt$height,pointsize=opt$pointsize)
tryCatch({
    par(las=1,mar=c(5.1,5.1,2.1,1.1))

    if(trait_type == "binary") {
	    yaxislab <- log(ytix)
	    names(yaxislab) <- ytix
	    plot(0,0,type="l",ylim=c(log(min(ytix)),log(max(ytix))),xlim=sapply(xrange,transx),main="",
	        xlab=ifelse(opt$raf.controls,"Risk Allele Frequency (% in Controls)",
	        "Overall Risk Allele Frequency (%)"),ylab="Odds Ratio",
	        lwd=1.2,cex=1,cex.lab=1.2,cex.axis=1.2,xaxt="n",bty="l",yaxt="n",col="transparent")
    } else {
	    yaxislab <- log(beta_tix)
	    names(yaxislab) <- beta_tix
	    plot(0,0,type="l",ylim=c(log(min(beta_tix)),log(max(beta_tix))),xlim=sapply(xrange,transx),main="",
	        xlab="Effect Allele Frequency (%)",
	        ylab=ifelse(abs(opt$trait_sd - 1) < sqrt(.Machine$double.eps),"Absolute Beta (SD Units)","Absolute Beta (Trait Units)"),
	        lwd=1.2,cex=1,cex.lab=1.2,cex.axis=1.2,xaxt="n",bty="l",yaxt="n",col="transparent")
    }
    ybottom <- par("usr")[3]
    ytop <- par("usr")[4]

    rect(xleft=par("usr")[1], ybottom, transx(0.005), ytop,col="grey70",border=NA)
    rect(xleft=transx(0.005), ybottom, transx(0.05), ytop,col="grey90",border=NA)
    rect(xleft=transx(0.995), ybottom, par("usr")[2], ytop,col="grey70",border=NA)
    rect(xleft=transx(0.995), ybottom, transx(0.95), ytop,col="grey90",border=NA)

	if(trait_type == "binary") {
		xplot <- ifelse(opt$raf.controls,"raf.controls","raf")
		for(i in 1:dim(samplesizes)[1]){
			xydata <- na.omit(data.frame('x'=sapply(plotdata[[i]][[xplot]],transx),'y'=log(plotdata[[i]]$OR)))
			lines(xydata$x,xydata$y,type="l",lty=alllines[i],lwd=3,col=allcols[i])
		}
	} else {
		for(i in 1:dim(samplesizes)[1]){
			xydata <- na.omit(data.frame('x'=sapply(plotdata[[i]]$af,transx),'y'=log(plotdata[[i]]$beta)))
			lines(xydata$x,xydata$y,type="l",lty=alllines[i],lwd=3,col=allcols[i])
		}
	}

    axisX <- c(0.0001,0.0005,0.005,0.05,0.5)
    axisX <- sort(unique(c(axisX,1-axisX)))
    axis(side=1,at=sapply(axisX,transx),labels=F,tick=T)
    axis(side=1,at=sapply(axisX,transx),labels=axisX*100,tick=F,cex.axis=1.2)
    axis(side=2,at=yaxislab,labels=names(yaxislab),tick=T,cex.axis=1.2)

	if(exists("risksnps")){
		points(sapply(risksnps[[ifelse(opt$raf.controls,"RAF.CONTROLS","RAF")]],transx),log(risksnps$OR),col="black",pch=4)	
	}

	if(trait_type == "binary") {
		legend("topleft",legend=paste(samplesizes$cases,"cases versus",samplesizes$controls,"controls"),
			col= allcols[1:dim(samplesizes)[1]],
			lty = alllines[1:dim(samplesizes)[1]],
			pt.cex=2,title="80% Power",lwd=3)
	} else {
		legend("topleft",legend=paste(samplesizes$samples,"samples"),
			col= allcols[1:dim(samplesizes)[1]],
			lty = alllines[1:dim(samplesizes)[1]],
			pt.cex=2,title="80% Power",lwd=3)
	}
}, error=function(e) {
	dev.off()
	stop(e)
})
dev.off()

print(paste("Plot can be found here:",paste0(opt$prefix,"_PowerAnalysis.png")))
