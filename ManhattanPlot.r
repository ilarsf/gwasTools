# Manhattan plot by Lars Fritsche 2019

options(stringsAsFactors=F)
library("optparse")
library("data.table")
 
ManhattanPlot <- function(res,top.size=0.125,break.top=15,hitregion="",
	chr="CHROM",pos="POS",pvalue="PVALUE",
	log10p=F,sigthreshold="5E-8",coltop=F,maintitle="",DTthreads=1) {

	options(stringsAsFactors=F)
	require("plotrix")
	require("data.table")
	require("RColorBrewer")
	setDTthreads(DTthreads)

	if(!is.data.table(res)) res <- as.data.table(res)
	setnames(res,c(chr,pos,pvalue),c("CHROM","POS","PVALUE"),skip_absent=T)

	if (hitregion != ""){
		candidateRegions <- read.table(hitregion,sep="\t",header=T,check.names=F,comment.char="")
	} else {
		candidateRegions <- data.frame(
			'CHROM'=character(0),
			'START'=numeric(0),
			'END'=numeric(0),
			'COL'=character(0),
			'LEGENDTEXT'=character(0)
			)
	}

	# horizontal lines and corresponding colors
	yLine <- c(-log10(sort(as.numeric(unlist(strsplit(sigthreshold,","))))))
	colLine <- c("red")

	if(!log10p) {
		res[,LOG10P:=-log10(PVALUE)]
	} else {
		setnames(res,"PVALUE","LOG10P")
	}

	res <- na.omit(res)
	res <- res[!is.infinite(LOG10P),]
	res <- res[order(as.numeric(gsub("chr","",gsub("X",23,CHROM))),POS),]
	Nmarkers <- nrow(res)

	# Thinning; remove 95% of variants with P > 0.05
	prethin1 <- which(res$LOG10P >= -log10(0.05) & res$LOG10P <= 6)
	prethin2 <- which(res$LOG10P < -log10(0.05))
	prethin2 <- prethin2[order(rnorm(length(prethin2)))[1:(length(prethin2)/20)]]
	prethin <- c(prethin1,prethin2)

	# Additional thinning using unique after rounding position and p-value
	thinned <- prethin[which(!duplicated(data.table(res$CHROM[prethin],
					round(res$POS[prethin]/10000),round(res$LOG10P[prethin],1))))]

	# Plot all SNPs with p < 1E-6
	keepTop <- which(res$LOG10P > 6)
	thinned <- sort(c(keepTop,thinned))

	# Prepare plot data / two-colored chromosomes / with fixed gap
	plotdata <- data.table(res[thinned,.(CHROM,POS,LOG10P)],
		'pch'=20,
		'highlightColor'=as.character(NA))
	
	chrs <- c(1:22,"X",23,"Y",24,"XY",25,"MT",26)
	chrs <- c(chrs,paste0("chr",chrs))
	chrs <- chrs[which(chrs %in% plotdata$CHROM)]
	chrNr <- as.numeric(gsub("chr","",as.character(factor(plotdata$CHROM,levels=chrs,labels=1:length(chrs)))))
	chrColors <- c("grey40","grey60")[(1:length(chrs)%%2)+1]
	names(chrColors) <- chrs

	plotdata[,pcol:=chrColors[CHROM]]
	plotdata[,chrNr:=chrNr]

	endPos <- 0
	plotPos <- numeric(0)
	chrLab <- numeric(0)
	chrGAP <- 1E7
	for(chr in chrs){
		chrTemp <- which(plotdata$CHROM == chr)
		chrPOS <- plotdata$POS[chrTemp]-min(plotdata$POS[chrTemp],na.rm=T)+endPos+1
		chrLab <- c(chrLab,mean(chrPOS,na.rm=T))
		endPos <- max(chrPOS,na.rm=T)+chrGAP
		plotPos <- c(plotPos,chrPOS)
	}
	plotdata[,plotPos:=plotPos]
	
	chrs <- gsub("chr","",chrs)

	# update numeric non-autosomal chromosome names
	fixChr <- c(1:22,"X","Y","XY","MT","X","Y","XY","MT")
	names(fixChr) <- c(1:26,"X","Y","XY","MT")
	chrs <- fixChr[chrs]

	# Highlight candidate regions
	if(nrow(candidateRegions)>0){
		a <- 0
		while(a < nrow(candidateRegions)){
			a <- a + 1 
			if(!coltop){
				overlap <- plotdata[CHROM == candidateRegions$CHROM[a] &
							 POS >= candidateRegions$START[a] &
							 POS <= candidateRegions$END[a] &
							 is.na(highlightColor),.I]
			} else {
				overlap <- plotdata[CHROM == candidateRegions$CHROM[a] &
							 POS >= candidateRegions$START[a] &
							 POS <= candidateRegions$END[a] &
							 is.na(highlightColor) &
							 LOG10P >= min(yLine),.I]
			}
			if(length(overlap)==0) next
			plotdata[overlap,highlightColor:=candidateRegions$COL[a]]
			plotdata[overlap,pcol:=NA]
		}
	}

	# Manhattan plot
	par(mar=c(5.1,5.1,4.1,1.1),las=1)
	x = plotdata$plotPos
	y = plotdata$LOG10P
	maxY <- max(y,na.rm=T)

	# Version with two y axes
	if(maxY > break.top/(1 - top.size)){
		# Manhattan plot with two different y axis scales

		# set axis labels of both scales
		lab1 <- pretty(c(0,break.top),n=ceiling(12 * (1-top.size)))
		lab1 <- c(lab1[lab1 < break.top],break.top)
		lab2 <- pretty(c(break.top,maxY),n=max(3,floor(12 * top.size)))
		lab2 <- lab2[lab2 > max(lab1)]

		# resulting range of top scale in bottom scale units
		top.range = break.top/(1 - top.size) - break.top
		top.data = max(lab2)-break.top
		# function to rescale the top part
		rescale = function(y) {break.top+(y-break.top)/(top.data/top.range)}

		# plot bottom part / rescaled
		plot(x[y<break.top],y[y<break.top],ylim=c(0,break.top+top.range),axes=FALSE,
			pch=plotdata$pch[y<break.top], cex=0.9,cex.lab=1.5,cex.axis=1.5, xaxt="n",
			col=plotdata$pcol[y<break.top], ylab=expression(-log[10]*italic(P)), xlab="",bty="n",
			main=gsub("_"," ",paste0(maintitle,"\n",format(Nmarkers,big.mark=",",scientific=F),
			" variants")))
		# plot top part
		points(x[y>break.top],rescale(y[y>break.top]),pch=plotdata$pch[y>break.top],
			col=plotdata$pcol[y>break.top],cex=0.9)

		# plot highlighted regions
		for(hcol in unique(plotdata$highlightColor)){
			topDot <- plotdata[which(plotdata$highlightColor == hcol & y>break.top),]
			if(length(topDot)>0) {
				points(topDot$plotPos,rescale(topDot$LOG10P),pch=20,col=hcol, cex=0.9)
			}
			bottomDot <- plotdata[which(plotdata$highlightColor == hcol & y<=break.top),]
			if(length(bottomDot)>0) {
				points(bottomDot$plotPos,bottomDot$LOG10P,pch=20,col=hcol, cex=0.9)
			}
		}

		# add axes and axis labels
		axis(1,at=chrLab[seq(1,length(chrLab),by=2)],
			labels=chrs[1:length(chrLab)][seq(1,length(chrLab),by=2)],
			las=1,tick=F,cex.axis=1.5,cex.lab=1.5,line=2)
		axis(1,at=chrLab[seq(2,length(chrLab),by=2)],
			labels=chrs[1:length(chrLab)][seq(2,length(chrLab),by=2)],
			las=1,tick=F,cex.axis=1.5,cex.lab=1.5,line=0)
		axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
		axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)

		# plot axis breaks and indicate line of axis break
		box()
		axis.break(axis=2,breakpos=break.top,style="zigzag",brw=0.02)
		axis.break(axis=4,breakpos=break.top,style="zigzag",brw=0.02)
		abline(h=break.top,lwd=1.5,lty=2,col="grey")
		if(length(yLine)>0) {
			for(rl in 1:length(yLine)){
				if(yLine[rl] <= break.top) {
					abline(h=yLine[rl],lwd=1.5,col=colLine[rl],lty=2)
				} else {
					abline(h=rescale(yLine[rl]),lwd=1.5,col=colLine[rl],lty=2)
				}
			}
		}
	# Version with one y axis / no break in axis
	} else {
		plot(x,y,xaxt="n",pch=plotdata$pch,cex=0.9,cex.lab=1.5,cex.axis=1.5,xaxt="n",
			col=plotdata$pcol,ylab=expression(-log[10]*italic(P)),xlab="",bty="o",
			main=gsub("_"," ",paste0(maintitle,"\n",format(Nmarkers,big.mark=",",scientific=F),
				" variants")),
			ylim=c(0,ceiling(max(maxY+1,yLine))))
		axis(1,at=chrLab[seq(1,length(chrLab),by=2)],
			labels=chrs[1:length(chrLab)][seq(1,length(chrLab),by=2)],
			las=1,tick=F,cex.axis=1.5,cex.lab=1.5,line=2)
		axis(1,at=chrLab[seq(2,length(chrLab),by=2)],
			labels=chrs[1:length(chrLab)][seq(2,length(chrLab),by=2)],
			las=1,tick=F,cex.axis=1.5,cex.lab=1.5,line=0)
		# plot highlighted regions on top
		for(hcol in unique(plotdata$highlightColor)){
			extraDot <- plotdata[which(plotdata$highlightColor == hcol),]
			points(extraDot$plotPos,extraDot$LOG10P,pch=20,col=hcol, cex=0.9)
		}
		# genome-wide significance linn
		if(length(yLine)>0) abline(h=yLine,lwd=1.5,col=colLine,lty=2)
	}

	# Add legend if candidate regions are present 
	ltext <- unique(candidateRegions[,c("COL","LEGENDTEXT")])
	if(dim(ltext)[1]>0) legend("topleft",legend=ltext$LEGENDTEXT,col=ltext$COL,pch=15,bty="n")
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
  make_option("--width", type="numeric", default=1600,
    help="Width Manhattan plot in pixel [default=1600]"),
  make_option("--height", type="numeric", default=900,
    help="Height Manhattan plot in pixel [default=900]"),
  make_option("--pointsize", type="numeric", default=16,
    help="Point size of plots [default=16]"),
  make_option("--hitregion", type="character", default="",
    help="File with candidate regions, CHROM;START;END;COL;LEGENDTEXT [default='']"),
  make_option("--chr", type="character", default="CHR",
    help="name of column with chromosome [default='CHR']"),
  make_option("--pos", type="character", default="POS",
    help="name of column with position [default='POS']"),
  make_option("--pvalue", type="character", default="PVALUE",
    help="name of column with p.value [default='PVALUE']"),
  make_option("--log10p", type="logical", default=F,
    help="Input p.value column with -log10(p.value) [default=F]"),    
  make_option("--sigthreshold", type="character", default="5E-8",
    help="Significance threshold [default=5E-8]"),    	
  make_option("--coltop", type="logical", default=F,
    help="Highlight only markers above the significance threshold [default=F]"),    	
  make_option("--maintitle", type="character", default="",
    help="Plot title"),
  make_option("--threads", type="numeric", default=1,
    help="DTthreads") 
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

gwas <- fread(opt$input,select=c(opt$chr,opt$pos,opt$pvalue),col.names=c("CHROM","POS","PVALUE"))
png(filename = paste0(opt$prefix,"_Manhattan.png"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
	ManhattanPlot(res=gwas,top.size=opt$top.size,break.top=opt$break.top,hitregion=opt$hitregion,
		log10p=opt$log10p,sigthreshold=opt$sigthreshold,coltop=opt$coltop,maintitle=opt$maintitle,DTthreads=opt$threads)
dev.off()
