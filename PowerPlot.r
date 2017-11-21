options(stringsAsFactors=F)

library("optparse")
library("RColorBrewer")
library("pwr")

# Power Analysis
getPower <- function(ncases,ncontrols,alpha,raf.controls){
	# effective samples size in alleles (x2)
	n <- 2*2/(1/ncases+1/ncontrols)
	h <- pwr.2p.test(h = NULL, n = n, sig.level = alpha, power = 0.8, alternative = "two.sided")$h
	
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

	out <- signif(data.frame(raf.controls,raf.cases,raf,OR,maf),4)
	out[notok,-1] <- NA
	return(out)
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
  make_option("--cases", type="character", default="100,1000,10000",
    help="Number(s) of cases; comma-separated [default='100,1000,10000']"),        
  make_option("--controls", type="character", default="100,1000,10000",
    help="Number(s) of controls; comma-separated [default='100,1000,10000']"),        
  make_option("--minMAF", type="numeric", default=0.0005,
    help="minimal minor frequency threshold [default=0.0005]"),
  make_option("--alpha", type="numeric", default=5E-8,
    help="alpha [default=5E-8]"),
  make_option("--risksnps", type="character", default="",
    help="File with risk SNPs; columns RAF,RAF.CONTROLS and OR, tab-delimited [default='']"),
  make_option("--raf.controls", type="logical", default=T,
    help="Plot RAF in controls instead of overall RAF [default=T]"),
  make_option("--width", type="numeric", default=900,
    help="Plot width in pixel [default=900]"),
  make_option("--height", type="numeric", default=900,
    help="Plot height in pixel [default=900]"),
  make_option("--pointsize", type="numeric", default=16,
    help="Point size of plots [default=16]"),
  make_option("--stepsize", type="numeric", default=0.00001,
    help="Stepwise frequency increment / smoothness of plot [default=0.00001]"),
  make_option("--ytix", type="character", default="1,1.2,1.5,2,3,5,10",
    help="Y-axis ticks [default='1,1.2,1.5,2,3,5,10']")
)
parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

samplesizes <- data.frame('cases'=as.integer(strsplit(opt$cases,",")[[1]]),
		   'controls'=as.integer(strsplit(opt$controls,",")[[1]]))

plotdata <- list()

freqs <- seq(0,1,opt$stepsize)

# maximal 24 unique line/color combinations
allcols <- rep(brewer.pal(8,"Dark2"),3)
alllines <- rep(1:6,4)

xrange <- c(opt$minMAF,1-opt$minMAF)

for(i in 1:dim(samplesizes)[1]){
	minMAF <- 1 / (2*(samplesizes$cases[i] + samplesizes$controls[i]))
	plotfreqs <- c(minMAF,freqs[which(freqs>minMAF & freqs < 1-minMAF)],1-minMAF)	
	plotdata[[i]] <- getPower(ncases=samplesizes$cases[i],ncontrols=samplesizes$controls[i],
					alpha=opt$alpha,plotfreqs)
}

if(opt$risksnps != "" & file.exists(opt$risksnps)){
	risksnps <- read.table(opt$risksnps,sep="\t",header=T,comment.char="")
	if( !all(c("OR","RAF","RAF.CONTROLS") %in% colnames(risksnps))) stop("Check columns in risk SNP file")
	risksnps$RAF <- ifelse(risksnps$OR < 1,1 - risksnps$RAF,risksnps$RAF)
	risksnps$RAF.CONTROLS <- ifelse(risksnps$OR < 1,1 - risksnps$RAF.CONTROLS,risksnps$RAF.CONTROLS)	
	risksnps$OR <- ifelse(risksnps$OR < 1,1/risksnps$OR,risksnps$OR)	
}

png(paste0(opt$prefix,"_PowerAnalysis.png"),width=opt$width,height=opt$height,pointsize=opt$pointsize)
    ytix <- as.numeric(strsplit(opt$ytix,",")[[1]])
    yaxislab <- log(ytix)
    names(yaxislab) <- ytix

    par(las=1,mar=c(5.1,5.1,2.1,1.1))
    plot(0,0,type="l",ylim=c(log(min(ytix)),log(max(ytix))),xlim=sapply(xrange,transx),main="",
        xlab=ifelse(opt$raf.controls,"Risk Allele Frequency (% in Controls)",
        "Overall Risk Allele Frequency (%)"),ylab="Odds Ratio",
        lwd=1.2,cex=1,cex.lab=1.2,cex.axis=1.2,xaxt="n",bty="l",yaxt="n",col="transparent")
    ybottom <- par("usr")[3]
    ytop <- par("usr")[4]

    rect(xleft=par("usr")[1], ybottom, transx(0.005), ytop,col="grey70",border=NA)
    rect(xleft=transx(0.005), ybottom, transx(0.05), ytop,col="grey90",border=NA)
    rect(xleft=transx(0.995), ybottom, par("usr")[2], ytop,col="grey70",border=NA)
    rect(xleft=transx(0.995), ybottom, transx(0.95), ytop,col="grey90",border=NA)

	xplot <- ifelse(opt$raf.controls,"raf.controls","raf")

	for(i in 1:dim(samplesizes)[1]){
		xydata <- na.omit(data.frame('x'=sapply(plotdata[[i]][[xplot]],transx),'y'=log(plotdata[[i]]$OR)))
    	lines(xydata$x,xydata$y,type="l",lty=alllines[i],lwd=3,col=allcols[i])
    }

    axisX <- c(0.0001,0.0005,0.005,0.05,0.5)
    axisX <- sort(unique(c(axisX,1-axisX)))
    axis(side=1,at=sapply(axisX,transx),labels=F,tick=T)
    axis(side=1,at=sapply(axisX,transx),labels=axisX*100,tick=F,cex.axis=1.2)
    axis(side=2,at=yaxislab,labels=names(yaxislab),tick=T,cex.axis=1.2)

	if(exists("risksnps")){
		points(sapply(risksnps[[ifelse(opt$raf.controls,"RAF.CONTROLS","RAF")]],transx),log(risksnps$OR),col="black",pch=4)	
	}

	legend("topleft",legend=paste(samplesizes$cases,"cases versus",samplesizes$controls,"controls"),
    	col= allcols[1:dim(samplesizes)[1]],
    	lty = alllines[1:dim(samplesizes)[1]],
    	pt.cex=2,title="80% Power",lwd=3)
dev.off()

print(paste("Plot can be found here:",paste0(opt$prefix,"_PowerAnalysis.png")))
