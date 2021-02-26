# Copyright (c) 2019 Lars Fritsche
# Manhattan plot Rscript

options(stringsAsFactors=F)
library("optparse")
require("plotrix")
require("data.table")
require("RColorBrewer")

t_white <- grDevices::rgb(t(grDevices::col2rgb("white")), alpha = 100, maxColorValue = 255)
t_black <- grDevices::rgb(t(grDevices::col2rgb("black")), alpha = 75, maxColorValue = 255)

# read very small pvalues and covert to -lop10P values
Ptolog10 <- function(P) {
    if (as.numeric(P) > 0) {
        log10P <- -log10(as.numeric(P))
    } else if (is.na(P) | as.numeric(P) > 1 | as.numeric(P) < 0) {
        log10P <- NA
    } else {
        part1 <- as.numeric(gsub("(.+)[Ee]-.+", "\\1", P))
        part1 <- -log10(as.numeric(part1) * 0.1)
        part2 <- as.numeric(gsub(".+[Ee]-(.+)", "\\1", P))
        part2 <- part2 - 1
        log10P <- part1 + part2
    }
    return(as.numeric(log10P))
}

# see latest version: /net/junglebook/home/larsf/Rfunctions/function.getNearestGene.r
getNearestGene <- function(input,Marker="Marker",chromosome="chromosome",
	position="position",build="hg19"){
	require(Map2NCBI)
	require(data.table)

	markers <- input[,c(Marker,chromosome,position),with=F]
	names(markers) <- c("Marker","chromosome","position")
	markers[chromosome == "X",chromosome:= 23]

	file.genelist <- path.expand(paste0("~/ExonList_",build,".txt"))

	if(!file.exists(file.genelist)){
		files.url <- paste0("http://hgdownload.soe.ucsc.edu/goldenPath/",build,
			"/database/ncbiRefSeq",c(".sql",".txt.gz"))
		files.temp <- tempfile(fileext=c(".sql",".gz"))
		
		download.file(files.url[2],destfile=files.temp[2])
		download.file(files.url[1],destfile=files.temp[1])

		header <- readLines(files.temp[1])
		header <- header[(grep("^CREATE TABLE",header)+1):(grep("^  KEY",header)[1]-1)]
		genetable <- fread(files.temp[2],col.names=gsub(".+`(.+)`.+","\\1",header))
		genetable[,chrom:=gsub("^chr","",chrom)]
		genetable <- unique(genetable[chrom %in% c(1:22,"X","Y","M") & !grepl("^X",name),
			.(name2,chrom,exonStarts,exonEnds)])
		
		exontable <- apply(genetable,1,function(x) {
			data.table(
				'FeatureName'=x[["name2"]],
				'chromosome'=x[["chrom"]],
				'start'=as.numeric(strsplit(x[["exonStarts"]],",")[[1]]),
				'end'=as.numeric(strsplit(x[["exonEnds"]],",")[[1]]))
		})
		genetable <- unique(rbindlist(exontable))		
		genetable <- genetable[order(chromosome,`start`,`end`),]
		fwrite(genetable,file.genelist,sep="\t",quote=T)
	} else {
		genetable <- fread(file.genelist,sep="\t")
	}

	MappedGenes <- MapMarkers(genetable, markers, nAut=22, other = c("X","Y"),savefiles=F)
	
	MappedGenes <- MappedGenes[,.(Marker,chromosome,position,FeatureName,`Inside?`)]
	MappedGenes[chromosome %in% c(23,25),chromosome:="X"]
	MappedGenes[chromosome == 24,chromosome:="Y"]
	
	setnames(MappedGenes,c("Marker","chromosome","position","FeatureName","Inside?"),
		c(Marker,chromosome,position,"LABEL","RelativeToGene"))

	output <- merge(input,MappedGenes,by=c(Marker,chromosome,position))
	return(output)
}

# main function to generate Manhattan plots
ManhattanPlot <- function(res,top.size=0.125,break.top=15,hitregion=NULL,
	chr="CHROM",pos="POS",pvalue="PVALUE",build="hg19",labelPeaks=T,regionSize=500000,
	log10p=F,sigthreshold="5E-8",coltop=F,maintitle="",DTthreads=1) {

	# devtools::install_github('JosephCrispell/basicPlotteR')
	require(basicPlotteR)
	setDTthreads(DTthreads)

	if(!is.data.table(res)) res <- as.data.table(res)

	setnames(res,c(chr,pos,pvalue),c("CHROM","POS","PVALUE"),skip_absent=T)

	# horizontal lines and corresponding colors
	yLine <- c(-log10(sort(as.numeric(unlist(strsplit(sigthreshold,","))))))
	colLine <- c("red")

	if(!log10p) {
		res[,LOG10P:=-log10(as.numeric(PVALUE))]
		res[as.numeric(PVALUE) == 0,LOG10P:=sapply(PVALUE,Ptolog10)]
	} else {
		res[,PVALUE:=as.numeric(PVALUE)]
		setnames(res,"PVALUE","LOG10P")
	}

	res <- na.omit(res[!is.infinite(LOG10P),.(CHROM,POS,LOG10P)])
	
	chrs <- c(1:22,"X",23,"Y",24,"XY",25,"MT",26)
	chrs <- c(chrs,paste0("chr",chrs))
	chrs <- chrs[which(chrs %in% unique(res$CHROM))]	
	
	res[,`:=`(CHROM=as.character(CHROM),
		POS=as.numeric(POS),
		numCHR=as.numeric(gsub("chr","",as.character(
			factor(CHROM,levels=chrs,labels=1:length(chrs))))))]

	res <- res[order(numCHR,POS),]
	Nmarkers <- nrow(res)

	hits <- res[LOG10P >= min(yLine),]
	
	if (!is.null(hitregion)){
		candidateRegions <- fread(hitregion,sep="\t",header=T)
		candidateRegions[,CHROM:=as.character(CHROM)]
		for(r in 1:nrow(candidateRegions)){
			rhits <- res[CHROM == candidateRegions$CHROM[r] & 
				POS >= candidateRegions$START[r] & POS <= candidateRegions$END[r],]
			candidateRegions$POS[r] <- rhits$POS[which.max(rhits$LOG10P)]
		}		
	} else if(nrow(hits) > 0){
		hits[,`:=`(POS0=POS-1,CHROM=gsub("chr","",CHROM))]
		x <- as.numeric(hits$POS)
		y <- hits$numCHR
		start = c(1, which(diff(y) != 0 | diff(x) <= -regionSize | diff(x) >= regionSize) + 1)
		end = c(start - 1, length(x))
		candidateRegions <- data.table(
			'CHROM'=as.character(hits$CHROM[start]),
			'START'=hits$POS[start] - regionSize,
			'END'=hits$POS[end] + regionSize,
			'COL'="blue",
			'MARKER'=1:length(start),
			'POS'=NA)
		candidateRegions[START < 1,START:=1]
		for(r in 1:nrow(candidateRegions)){
			rhits <- hits[CHROM == candidateRegions$CHROM[r] & 
				POS >= candidateRegions$START[r] & POS <= candidateRegions$END[r],]
			candidateRegions$POS[r] <- rhits$POS[which.max(rhits$LOG10P)]
		}

		if(labelPeaks){
			# add nearest GENENAME
			candidateRegions <- getNearestGene(input=candidateRegions,Marker="MARKER",chromosome="CHROM",position="POS",build=build)
			candidateRegions <- candidateRegions[,.(CHROM,START,END,COL,POS,LABEL,RelativeToGene)]
		} else {
			candidateRegions$LABEL <- NA
		}
		
	} else {
	# empty table if there are no hits
		candidateRegions <- data.table(
			'CHROM'=character(0),
			'START'=numeric(0),
			'END'=numeric(0),
			'COL'=character(0),
			'POS'=numeric(0),
			'LABEL'=numeric(0))
	}

	# Thinning; remove 95% of variants with P > 0.05	
	thinned <- res[
		(LOG10P >= -log10(0.05) & LOG10P <= 6) |
		(LOG10P < -log10(0.05) & sample(c(T,rep(F,19)),.N,replace=T)),
			.(CHROM,POS,LOG10P,numCHR)]

	# Additional thinning using unique after rounding position and p-value
	# keep all SNPs with p < 1E-6
	plotdata <- rbind(
		res[LOG10P > 6,],
		thinned[!duplicated(paste(CHROM,round(POS/1000),round(LOG10P,1))),])
	plotdata <- plotdata[order(numCHR,POS),.(CHROM,POS,LOG10P,numCHR)]

	# Prepare plot data / two-colored chromosomes
	plotdata[,`:=`(
		pch=20,
		highlightColor=as.character(NA),
		pcol = ifelse(numCHR %%2 == 0, "grey40","grey60"))]

	# covert pos on chromosomes to continuous plot positions
	chrGAP <- 1E7
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
	
	# update numeric non-autosomal chromosome names
	chrs <- gsub("chr","",chrs)
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
		by=c("CHROM","POS"),sort=F)

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
		if(labelPeaks)
			regionLabels[,addTextLabels(xCoords = `x`, yCoords = `y`, labels = `LABEL`, 
				col.label = "black", col.line = t_black, cex.label = 1, col.background = t_white)]
	}	
	return(candidateRegions)
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
  make_option("--hitregion", type="character", default=NULL,
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
  make_option("--build", type="character", default="hg19",
    help="Genome build [default='hg19']"),
  make_option("--threads", type="numeric", default=1,
    help="DTthreads") 
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

gwas <- fread(opt$input,select=c(opt$chr,opt$pos,opt$pvalue),col.names=c("CHROM","POS","PVALUE"))
png(filename = paste0(opt$prefix,"_Manhattan.png"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
	candidateRegions <- ManhattanPlot(res=gwas,top.size=opt$top.size,break.top=opt$break.top,hitregion=opt$hitregion,
		log10p=opt$log10p,sigthreshold=opt$sigthreshold,coltop=opt$coltop,maintitle=opt$maintitle,
		build=opt$build,DTthreads=opt$threads)
dev.off()

fwrite(candidateRegions,paste0(opt$prefix,"_Manhattan.txt"),sep="\t",quote=F)