# Copyright (c) 2019 Lars Fritsche
# Manhattan plot Rscript

options(stringsAsFactors=F)
library("optparse")
library("data.table")
t_white <- grDevices::rgb(t(grDevices::col2rgb("white")), alpha = 100, maxColorValue = 255)
t_black <- grDevices::rgb(t(grDevices::col2rgb("black")), alpha = 75, maxColorValue = 255)

getNearestGene <- function(input,Marker="Marker",chromosome="chromosome",
	position="position",build="37"){
	require(Map2NCBI)
	require(data.table)

	markers <- input[,c(Marker,chromosome,position),with=F]
	names(markers) <- c("Marker","chromosome","position")

	if(build == "37"){
		file.genelist37 <- path.expand("~/GeneList_Homo_sapiens_BUILD.37.3.txt")
	
		if(!file.exists(file.genelist37)){
			genelist37 <- GetGeneList_v11("Homo sapiens",build="BUILD.37.3",savefiles=F,destfile=tempfile())
			1
			y
			fwrite(genelist37,file.genelist37,sep="\t",quote=T)
		}
		genelist37 <- fread(file.genelist37)
		genelist37 <- genelist37[feature_type == "GENE",]
		MappedGenes <- MapMarkers(genelist37, markers, nAut=22, other = c("X"),savefiles=F)
		MappedGenes <- MappedGenes[,.(Marker,chromosome,position,feature_name,`Inside?`)]
		setnames(MappedGenes,c("Marker","chromosome","position","feature_name","Inside?"),
			c(Marker,chromosome,position,"LABEL","RelativeToGene"))
	}

	if(build == "38"){
		file.genelist38 <- path.expand("~/GeneList_Homo_sapiens_BUILD.38.txt")
		if(!file.exists(file.genelist38)){
			genelist38 <- GetGeneList("Homo sapiens",savefiles=F,destfile=tempfile())
			y
			y
			fwrite(genelist38,file.genelist38,sep="\t",quote=T)
		}
		genelist38 <- fread(file.genelist38)
		genelist38 <- genelist38[feature == "gene" &
			seq_type %in% c("chromosome","mitochondrion") &
			`attributes` != "pseudo",]
		MappedGenes <- MapMarkers(genelist38, markers, nAut=22, other = c("X"),savefiles=F)
		MappedGenes <- MappedGenes[,.(Marker,chromosome,position,symbol,`Inside?`)]
		setnames(MappedGenes,c("Marker","chromosome","position","symbol","Inside?"),
			c(Marker,chromosome,position,"LABEL","RelativeToGene"))
	}
	output <- merge(input,MappedGenes,by=c(Marker,chromosome,position))
	return(output)
}

# main function to generate Manhattan plots
ManhattanPlot <- function(res,top.size=0.125,break.top=15,hitregion=NULL,
	chr="CHROM",pos="POS",pvalue="PVALUE",build="37",
	log10p=F,sigthreshold="5E-8",coltop=F,maintitle="",DTthreads=1) {

	options(stringsAsFactors=F)
	require("plotrix")
	require("data.table")
	require("RColorBrewer")
	# devtools::install_github('JosephCrispell/basicPlotteR')
	require(basicPlotteR)
	setDTthreads(DTthreads)

	if(!is.data.table(res)) res <- as.data.table(res)

	setnames(res,c(chr,pos,pvalue),c("CHROM","POS","PVALUE"),skip_absent=T)

	# horizontal lines and corresponding colors
	yLine <- c(-log10(sort(as.numeric(unlist(strsplit(sigthreshold,","))))))
	colLine <- c("red")

	if(!log10p) {
		res[,LOG10P:=-log10(PVALUE)]
	} else {
		setnames(res,"PVALUE","LOG10P")
	}

	res <- na.omit(res[!is.infinite(LOG10P),.(CHROM,POS,LOG10P)])
	res[,numCHR:=as.numeric(gsub("chr","",gsub("^X$|^XY$","23",CHROM)))]

	res <- res[order(numCHR,POS),]
	Nmarkers <- nrow(res)

	if (!is.null(hitregion)){
		candidateRegions <- fread(hitregion,sep="\t",header=T)
	} else {
		hits <- res[LOG10P >= min(yLine),]
		hits[,`:=`(POS0=POS-1,CHROM=gsub("chr","",CHROM))]
		x <- as.numeric(hits$POS)
		y <- hits$numCHR
		start = c(1, which(diff(y) != 0 | diff(x) <= -500000 | diff(x) >= 500000) + 1)
		end = c(start - 1, length(x))
		candidateRegions <- data.table(
			'CHROM'=hits$CHROM[start],
			'START'=hits$POS[start] - 500000,
			'END'=hits$POS[end] + 500000,
			'COL'="blue",
			'MARKER'=1:length(start),
			'POS'=NA)
		candidateRegions[START < 1,START:=1]
		for(r in 1:nrow(candidateRegions)){
			rhits <- hits[CHROM == candidateRegions$CHROM[r] & 
				POS >= candidateRegions$START[r] & POS <= candidateRegions$END[r],]
			candidateRegions$POS[r] <- rhits$POS[which.max(rhits$LOG10P)]
		}

		# add nearest GENENAME
		candidateRegions <- getNearestGene(input=candidateRegions,Marker="MARKER",chromosome="CHROM",position="POS",build=build)
		candidateRegions <- candidateRegions[,.(CHROM,START,END,COL,POS,LABEL,RelativeToGene)]	
	}

	# Thinning; remove 95% of variants with P > 0.05
	prethin1 <- which(res$LOG10P >= -log10(0.05) & res$LOG10P <= 6)
	prethin2 <- which(res$LOG10P < -log10(0.05))
	prethin2 <- prethin2[order(rnorm(length(prethin2)))[1:(length(prethin2)/20)]]
	prethin <- c(prethin1,prethin2)

	# Additional thinning using unique after rounding position and p-value
	thinned <- prethin[which(!duplicated(data.frame(res$CHROM[prethin],
					round(res$POS[prethin]/10000),round(res$LOG10P[prethin],1))))]

	# Plot all SNPs with p < 1E-6
	keepTop <- which(res$LOG10P > 6)
	thinned <- sort(c(keepTop,thinned))

	# Prepare plot data / two-colored chromosomes / with fixed gap
	plotdata <- res[thinned,.(CHROM,POS,LOG10P)]
	plotdata[,`:=`(pch=20,highlightColor=as.character(NA))]
	
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
	plotdata[,x:=plotPos]
	
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
							 is.na(highlightColor),`:=`
							 	(highlightColor = candidateRegions$COL[a],pcol=NA)]
			} else {
				overlap <- plotdata[CHROM == candidateRegions$CHROM[a] &
							 POS >= candidateRegions$START[a] &
							 POS <= candidateRegions$END[a] &
							 is.na(highlightColor) &
							 LOG10P >= min(yLine),`:=`
							 	(highlightColor = candidateRegions$COL[a],pcol=NA)]
			}
		}
	}


	# Manhattan plot
	par(mar=c(5.1,5.1,4.1,1.1),las=1)
	maxY <- plotdata[,max(LOG10P,na.rm=T)]

	# Version with two y axes
	if(maxY > break.top/(1 - top.size)){
		# Manhattan plot with two different y axis scales
		# set axis labels for both scales
		lab1 <- pretty(c(0,break.top),n=ceiling(12 * (1-top.size)))
		lab1 <- c(lab1[lab1 < break.top],break.top)
		lab2 <- pretty(c(break.top,maxY),n=max(3,floor(12 * top.size)))
		lab2 <- lab2[lab2 > max(lab1)]
		lab <- c(lab1,lab2)
		# resulting range of top scale in bottom scale units
		top.range = break.top/(1 - top.size) - break.top
		top.data = max(lab2)-break.top
		# function to rescale the top part
        rescaleY <- function(y) {
            if (y <= break.top) {
                y
            } else {
                break.top + (y - break.top)/(top.data/top.range)
            }
        }
		ylim=c(0,break.top+top.range)
		addbreak <- T
	} else {
	    break.top <- maxY
        rescaleY <- function(y) y
        addbreak <- F
        ylim <- c(0,ceiling(max(maxY+1,yLine)))
		lab <- pretty(ylim)
		lab <- lab[lab < maxY]
	}
	plotdata[,y:=sapply(LOG10P,rescaleY)]

	regionLabels <- merge(plotdata[,.(CHROM,POS,x,y)],
		candidateRegions,
		by=c("CHROM","POS"))

	# plot non-highlighted positions
	plotdata[,plot(x,y,ylim=ylim,axes=FALSE,
		pch=`pch`, cex=0.9,cex.lab=1.5,cex.axis=1.5, xaxt="n",
		col=`pcol`, ylab=expression(-log[10]*italic(P)), xlab="",bty="n",
		main=gsub("_"," ",paste0(maintitle,"\n",format(Nmarkers,big.mark=",",scientific=F),
		" variants")))]

	# add axes and axis labels
	axis(1,at=chrLab[seq(1,length(chrLab),by=2)],
		labels=chrs[1:length(chrLab)][seq(1,length(chrLab),by=2)],
		las=1,tick=F,cex.axis=1.5,cex.lab=1.5,line=2)
	axis(1,at=chrLab[seq(2,length(chrLab),by=2)],
		labels=chrs[1:length(chrLab)][seq(2,length(chrLab),by=2)],
		las=1,tick=F,cex.axis=1.5,cex.lab=1.5,line=0)
	axis(side=2,at=sapply(lab,rescaleY),labels=lab,cex.axis=1.5,cex.lab=1.5)

	if(addbreak){
		# plot axis breaks and indicate line of axis break
		box()
		axis.break(axis=2,breakpos=break.top,style="zigzag",brw=0.02)
		axis.break(axis=4,breakpos=break.top,style="zigzag",brw=0.02)
		abline(h=break.top,lwd=1.5,lty=2,col="grey")
	} 
	
	if(length(yLine)>0) abline(h=sapply(yLine,rescaleY),lwd=1.5,col=colLine,lty=2)	

	# label top regions
	if(nrow(regionLabels)>0){
		for(hcol in unique(regionLabels$COL)){
			plotdata[highlightColor == hcol,points(x,y,pch=20,col=hcol, cex=0.9)]
		}	
		# non-overlapping labels
		regionLabels[,addTextLabels(xCoords = `x`, yCoords = `y`, labels = `NearestGene`, 
			col.label = "black", col.line = t_black, cex.label = 1, col.background = t_white)]
	}
	
	regionLabels[,`:=`(x=NULL,y=NULL)]
	print(regionLabels)
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
    help="File with candidate regions, CHROM;START;END;COL;LABEL [default='']"),
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
