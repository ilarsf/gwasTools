options(stringsAsFactors=F)
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

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
  make_option("--log10p", type="logical", default=T,
    help="Input p.value column with -log10(p.value) [default=T]"),    
  make_option("--maintitle", type="character", default="",
    help="Plot title")  
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

if (opt$hitregion != ""){
	candidateRegions <- read.table(opt$hitregion,sep="\t",header=T,check.names=F)
} else {
	candidateRegions <- data.frame(
		'CHROM'=character(0),
		'START'=numeric(0),
		'END'=numeric(0),
		'COL'=character(0),
		'LEGENDTEXT'=character(0)
		)
}

chrcol <- opt$chr
poscol <- opt$pos

# horizontal lines and corresponding colors
yLine <- c(-log10(5E-8))
colLine <- c("red")

gwas <- fread(opt$input)

if(!opt$log10p) {
	gwas$log10P <- -log10(gwas[[opt$pvalue]])
	ycol <- "log10P"
} else {
	ycol <- opt$pvalue
}

gwas <- na.omit(data.frame(gwas[,c(chrcol,poscol,ycol),with=F]))
Nmarkers <- dim(gwas)[1]

# Thinning; remove 95% of variants with P > 0.05
prethin1 <- which(gwas[[ycol]] >= -log10(0.05) & gwas[[ycol]] <= 6)
prethin2 <- which(gwas[[ycol]] < -log10(0.05))
prethin2 <- prethin2[order(rnorm(length(prethin2)))[1:(length(prethin2)/20)]]
prethin <- c(prethin1,prethin2)

# Additional thinning using unique after rounding position and p-value
thinned <- prethin[which(!duplicated(data.frame(gwas[[chrcol]][prethin],
				round(gwas[[poscol]][prethin]/10000),round(gwas[[ycol]][prethin],1))))]

# Plot all SNPs with p < 1E-6
keepTop <- which(gwas[[ycol]] > 6)
thinned <- sort(c(keepTop,thinned))

# Prepare plot data / two-colored chromosomes / with fixed gap
CHR <- gwas[[chrcol]][thinned]

POS <- gwas[[poscol]][thinned]
log10P <- gwas[[ycol]][thinned]
chrs <- c(1:22,"X",23,"Y",24,"XY",25,"MT",26)
chrs <- chrs[which(chrs %in% CHR)]
chrNr <- as.numeric(as.character(factor(CHR,levels=chrs,labels=1:length(chrs))))
chrColors <- c("grey40","grey60")[(1:length(chrs)%%2)+1]
names(chrColors) <- chrs
plotdata <- data.frame(
    CHR,
    POS,
    log10P,
    'plotPos'=NA,
    'chrNr'= chrNr,
    pch=20,
    highlightColor=NA,
    pcol=chrColors[CHR],
    check.names=F)
endPos <- 0
plotPos <- numeric(0)
chrLab <- numeric(0)
chrGAP <- 1E7
for(chr in chrs){
    chrTemp <- which(CHR == chr)
    chrPOS <- POS[chrTemp]-min(POS[chrTemp],na.rm=T)+endPos+1
    chrLab <- c(chrLab,mean(chrPOS,na.rm=T))
    endPos <- max(chrPOS,na.rm=T)+chrGAP
    plotPos <- c(plotPos,chrPOS)
}
plotdata$plotPos <- plotPos


# Highlight candidate regions
if(dim(candidateRegions)[1]>0){
	a <- 0
	while(a < dim(candidateRegions)[1]){
		a <- a + 1 
		overlap <- which(plotdata$CHR == candidateRegions$CHROM[a] &
						 plotdata$POS >= candidateRegions$START[a] &
						 plotdata$POS <= candidateRegions$END[a] &
						 is.na(plotdata$highlightColor)
						)
		if(length(overlap)==0) next
		plotdata$highlightColor[overlap] <- candidateRegions$COL[a]
		plotdata$pcol[overlap] <- NA
	}
}

# Manhattan plot
png(filename = paste0(opt$prefix,"_Manhattan.png"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
    par(mar=c(5.1,5.1,4.1,1.1),las=1)
    x = plotdata$plotPos
    y = plotdata$log10P
    maxY <- max(y,na.rm=T)

	# Version with two y axes
    if(maxY > opt$break.top){
        # Manhattan plot with two different y axis scales

        # set axis labels of both scales
        lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
        lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
        lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
        lab2 <- lab2[lab2 > max(lab1)]

        # resulting range of top scale in bottom scale units
        top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
        top.data = max(lab2)-opt$break.top
        # function to rescale the top part
        rescale = function(y) {opt$break.top+(y-opt$break.top)/(top.data/top.range)}

        # plot bottom part / rescaled
        plot(x[y<opt$break.top],y[y<opt$break.top],ylim=c(0,opt$break.top+top.range),axes=FALSE,
            pch=plotdata$pch[y<opt$break.top], cex=0.9,cex.lab=1.5,cex.axis=1.5, xaxt="n",
            col=plotdata$pcol[y<opt$break.top], ylab=expression(-log[10]*italic(P)), xlab="",bty="n",
            main=gsub("_"," ",paste0(opt$maintitle,"\n",format(Nmarkers,big.mark=",",scientific=F),
            " variants")))
        # plot top part
        points(x[y>opt$break.top],rescale(y[y>opt$break.top]),pch=plotdata$pch[y>opt$break.top],
            col=plotdata$pcol[y>opt$break.top],cex=0.9)

        # plot highlighted regions
        for(hcol in unique(plotdata$highlightColor)){
            topDot <- plotdata[which(plotdata$highlightColor == hcol & y>opt$break.top),]
            if(length(topDot)>0) {
                points(topDot$plotPos,rescale(topDot$log10P),pch=20,col=hcol, cex=0.9)
            }
            bottomDot <- plotdata[which(plotdata$highlightColor == hcol & y<=opt$break.top),]
            if(length(bottomDot)>0) {
                points(bottomDot$plotPos,bottomDot$log10P,pch=20,col=hcol, cex=0.9)
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
        axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
        axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
        abline(h=opt$break.top,lwd=1.5,lty=2,col="grey")
        if(length(yLine)>0) {
            for(rl in 1:length(yLine)){
                if(yLine[rl] <= opt$break.top) {
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
            main=gsub("_"," ",paste0(opt$maintitle,"\n",format(Nmarkers,big.mark=",",scientific=F),
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
            points(extraDot$plotPos,extraDot$log10P,pch=20,col=hcol, cex=0.9)
        }
        # genome-wide significance linn
        if(length(yLine)>0) abline(h=yLine,lwd=1.5,col=colLine,lty=2)
    }

	# Add legend if candidate regions are present 
    ltext <- unique(candidateRegions[,c("COL","LEGENDTEXT")])
    if(dim(ltext)[1]>0) legend("topleft",legend=ltext$LEGENDTEXT,col=ltext$COL,pch=15,bty="n")
dev.off()
