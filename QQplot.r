options(stringsAsFactors=F)
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

# QQ plot function
qqplotdata <- function(logpvector){
    o = sort(logpvector,decreasing=T)
    e = -log10(ppoints(length(o)))       
    qqdata <- data.frame(o,e)
    qqdata$o <- round(qqdata$o,3)
    qqdata$e <- round(qqdata$e,3)
    keepU <- which(!duplicated(qqdata))
    qqdata <- qqdata[keepU,]
    
    N <- length(logpvector) ## number of p-values
    ## create the confidence intervals
    qqdata$c95 <- NA
    qqdata$c05 <- NA

    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)

    for(i in 1:length(keepU)){
        j <- keepU[i]
        qqdata$c95[i] <- -log10(qbeta(0.95,j,N-j+1))
        qqdata$c05[i] <- -log10(qbeta(0.05,j,N-j+1))
    }
    return(qqdata)
}


option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited"),   
  make_option("--prefix", type="character", default="",
    help="Prefix of output files"),   
  make_option("--top.size", type="numeric", default=0.125,
    help="top size = proportion of total length y axis [default=0.125]"),
  make_option("--break.top", type="numeric", default=15,
    help="set axis break at -log10(P) [default=15]"),
  make_option("--width", type="numeric", default=900,
    help="Width Manhattan plot in pixel [default=900]"),
  make_option("--height", type="numeric", default=900,
    help="Height Manhattan plot in pixel [default=900]"),
  make_option("--pointsize", type="numeric", default=16,
    help="Point size of plots [default=16]"),
  make_option("--maf", type="character", default="MAF",
    help="name of column with MAF [default='MAF']"),
  make_option("--pvalue", type="character", default="PVALUE",
    help="name of column with p.value [default='PVALUE']"),
  make_option("--log10p", type="logical", default=F,
    help="Input p.value column with -log10(p.value) [default=F]"),    
  make_option("--maintitle", type="character", default="",
    help="Plot title")  
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# horizontal lines and corresponding colors
yLine <- c(-log10(5E-8))
colLine <- c("red")

gwas <- fread(opt$input)

if(!opt$log10p) {
	gwas$log10P <- -log10(gwas[[opt$p]])
	ycol <- "log10P"
} else {
	ycol <- opt$pvalue
}

gwas <- na.omit(data.frame(gwas[,c(opt$maf,ycol),with=F]))

minMAF <- min(gwas[[opt$maf]])

# Determine frequency bins and create variable for binned QQ plot
freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
gwas$freqbin <- cut(gwas[[opt$maf]], freqbins,include.lowest=T)
freqtable <- table(gwas$freqbin)
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
	fsnps <- which(gwas$freqbin ==names(freqtable)[f])
	plotdata <- qqplotdata(gwas[[ycol]][fsnps])
	fN <- c(fN,freqtable[f])
	fx <- c(fx,plotdata$e)
	fy <- c(fy,plotdata$o)
	fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
	conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
							'y'=c(plotdata$c95,rev(plotdata$c05)))
	legendcol <- c(legendcol,allcols[f])
}
legendtext <- paste0("MAF=",fbin,"; N SNPs=",format(fN,big.mark=",",scientific=FALSE))


## QQ plot by binned frequencies
png(filename = paste0(opt$prefix,".QQ_Plot.png"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
	xlim <- c(0,max(fx,na.rm=T))
	ylim <- c(0,max(fy,na.rm=T))
	maxY <- max(fy,na.rm=T)
	par(mar=c(5.1,5.1,4.1,1.1))
	# plot version with two axes
	if(maxY > opt$break.top){
		# create pretty y-axis labels
        lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
        lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
        lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
        lab2 <- lab2[lab2 > max(lab1)]

        # resulting range of top scale in bottom scale units
        top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
        top.data = max(lab2)-opt$break.top
        
        # function to rescale the top part
		rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
		rescaled.y = rescale(fy[fy>opt$break.top])
		plot(0,0,
			ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
			xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
			ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
			cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
			main=opt$maintitle,pch=19)
		
		# Plot confidence intervals	
		for(p in 1:length(conf)){
			polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
				col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
				border = NA)
		}

		# add points
		points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)

		# identify line & add axis break
		lines(xlim,xlim,col="black",lty = 2)
		axis(1,cex.axis=1.5,cex.lab=1.5)
		par(las=1)
		axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
		par(las=0)
		box()
		par(las=0)
		points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
		par(las=1)
		axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
		axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
		axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
		lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
		abline(h=ifelse(yLine<opt$break.top,
			yLine,
			rescale(yLine)),
			col=colLine,lwd=1.5,lty=2)
		legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
	# plot version with single y axes
	} else {
		par(mar=c(5.1,5.1,4.1,1.1),las=1)
		axislim <- ceiling(range(xlim,ylim,-log10(gwthreshold)))
		plot(0,0,
			ylim=axislim,xlim=xlim,axes=T,
			xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
			ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
			cex=1,cex.lab=1.5,cex.axis=1.5,col="transparent",
			main=opt$maintitle,pch=19)
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
dev.off()
