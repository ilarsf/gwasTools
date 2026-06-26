# Copyright (c) 2019 Lars Fritsche
# Manhattan plot Rscript

options(stringsAsFactors=F)
library("optparse")
require("plotrix")
require("data.table")
require("RColorBrewer")

t_white <- grDevices::rgb(t(grDevices::col2rgb("white")), alpha = 100, maxColorValue = 255)
t_black <- grDevices::rgb(t(grDevices::col2rgb("black")), alpha = 75, maxColorValue = 255)

stopf <- function(...) stop(sprintf(...), call. = FALSE)

defaultGeneDataDir <- function() {
	file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
	if(length(file_arg) == 0) return(file.path("data", "gene_positions"))
	script_path <- path.expand(sub("^--file=", "", file_arg[1]))
	return(file.path(dirname(script_path), "data", "gene_positions"))
}

DEFAULT_GENE_DATA_DIR <- defaultGeneDataDir()

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

validatePositiveNumeric <- function(values, column) {
	num <- suppressWarnings(as.numeric(values))
	invalid <- which(is.na(num) | !is.finite(num) | num <= 0)
	if(length(invalid) > 0) {
		stopf("Invalid %s values: expected finite values > 0. First invalid row(s): %s",
			column, paste(head(invalid, 10), collapse = ","))
	}
	return(num)
}

parseSigThresholds <- function(sigthreshold) {
	thresholds <- suppressWarnings(as.numeric(unlist(strsplit(sigthreshold,","))))
	if(length(thresholds) == 0 || any(is.na(thresholds)) || any(!is.finite(thresholds)) ||
			any(thresholds <= 0 | thresholds > 1)) {
		stopf("--sigthreshold must contain comma-separated p-values in (0, 1]")
	}
	return(thresholds)
}

normalizeChromosomes <- function(chrom) {
	raw <- trimws(as.character(chrom))
	canon <- toupper(sub("^chr", "", raw, ignore.case = TRUE))
	canon[canon == "M"] <- "MT"
	canon[canon == "23"] <- "X"
	canon[canon == "24"] <- "Y"
	canon[canon == "25"] <- "XY"
	canon[canon == "26"] <- "MT"
	valid <- c(as.character(1:22), "X", "Y", "XY", "MT")
	invalid <- which(is.na(canon) | canon == "" | !(canon %in% valid))
	if(length(invalid) > 0) {
		bad <- unique(raw[invalid])
		stopf("Unsupported chromosome value(s): %s", paste(head(bad, 10), collapse = ", "))
	}
	return(list(CHROM = canon, numCHR = match(canon, valid),
		levels = valid[valid %in% unique(canon)]))
}

# see latest version: /net/junglebook/home/larsf/Rfunctions/function.getNearestGene.r
getNearestGene <- function(input,Marker="Marker",chromosome="chromosome",
	position="position",build="hg19",gene_table="",gene_data_dir=DEFAULT_GENE_DATA_DIR,
	allow_download=F){
	if(!requireNamespace("Map2NCBI", quietly = TRUE)) {
		stopf("Nearest-gene labels require the optional R package 'Map2NCBI'")
	}
	file.genelist <- if(!is.null(gene_table) && nzchar(gene_table)) {
		path.expand(gene_table)
	} else {
		file.path(path.expand(gene_data_dir), paste0("ExonList_", build, ".txt"))
	}

	markers <- input[,c(Marker,chromosome,position),with=F]
	names(markers) <- c("Marker","chromosome","position")
	markers[chromosome == "X",chromosome:= 23]

	if(!file.exists(file.genelist)){
		if(!allow_download) {
			stopf("Nearest-gene labels require a local gene table: %s. Provide --gene_table, or pass --allow_gene_download T to create it under --gene_data_dir.",
				file.genelist)
		}
		if(!dir.exists(dirname(file.genelist))) {
			if(!dir.create(dirname(file.genelist), recursive = TRUE, showWarnings = FALSE)) {
				stopf("Could not create gene table directory: %s", dirname(file.genelist))
			}
		}
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

	MappedGenes <- Map2NCBI::MapMarkers(genetable, markers, nAut=22, other = c("X","Y"),savefiles=F)
	
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
	chr="CHROM",pos="POS",pvalue="PVALUE",build="hg19",labelPeaks=F,regionSize=500000,
	log10p=F,sigthreshold="5E-8",coltop=F,maintitle="",DTthreads=1,
	thinning_mode="adaptive",thinning_aggressiveness=1.0,thinning_sig_threshold=6,
	custom_thinning_rates=NULL,disable_thinning=F,gene_table="",
	gene_data_dir=DEFAULT_GENE_DATA_DIR,allow_gene_download=F,random_seed=1) {

	# devtools::install_github('JosephCrispell/basicPlotteR')
	require(basicPlotteR)
	setDTthreads(DTthreads)

	if(!is.data.table(res)) res <- as.data.table(res)

	setnames(res,c(chr,pos,pvalue),c("CHROM","POS","PVALUE"),skip_absent=T)
	missing <- setdiff(c("CHROM","POS","PVALUE"), names(res))
	if(length(missing) > 0) stopf("Missing required column(s): %s", paste(missing, collapse = ", "))

	thinning_mode <- tolower(thinning_mode)
	if(!(thinning_mode %in% c("adaptive","random","none"))) {
		stopf("--thinning_mode must be one of adaptive, random, or none")
	}
	if(thinning_mode == "none") disable_thinning <- T
	if(!is.finite(thinning_aggressiveness) || thinning_aggressiveness <= 0) {
		stopf("--thinning_aggressiveness must be a finite value > 0")
	}
	if(!is.finite(thinning_sig_threshold) || thinning_sig_threshold < 0) {
		stopf("--thinning_sig_threshold must be a finite value >= 0")
	}

	# horizontal lines and corresponding colors
	yLine <- c(-log10(sort(parseSigThresholds(sigthreshold))))
	colLine <- c("red")

	if(!log10p) {
		res[,LOG10P:=pToLog10(PVALUE)]
	} else {
		res[,LOG10P:=validateLog10P(PVALUE)]
	}

	chroms <- normalizeChromosomes(res$CHROM)
	res[,`:=`(CHROM=chroms$CHROM,
		POS=validatePositiveNumeric(POS, "POS"),
		numCHR=chroms$numCHR)]
	res <- res[,.(CHROM,POS,LOG10P,numCHR)]
	if(nrow(res) == 0) stopf("No valid rows remain after Manhattan input validation")
	
	chrs <- chroms$levels

	res <- res[order(numCHR,POS),]
	Nmarkers <- nrow(res)

	hits <- res[LOG10P >= min(yLine),]
	
	if (!is.null(hitregion)){
		candidateRegions <- fread(hitregion,sep="\t",header=T)
		missing_regions <- setdiff(c("CHROM","START","END","COL"), names(candidateRegions))
		if(length(missing_regions) > 0) {
			stopf("Candidate-region file is missing required column(s): %s",
				paste(missing_regions, collapse = ", "))
		}
		region_chroms <- normalizeChromosomes(candidateRegions$CHROM)
		candidateRegions[,`:=`(CHROM=region_chroms$CHROM,
			START=validatePositiveNumeric(START, "START"),
			END=validatePositiveNumeric(END, "END"))]
		if(any(candidateRegions$START > candidateRegions$END)) {
			stopf("Candidate-region START must be <= END")
		}
		if(!("LABEL" %in% names(candidateRegions))) candidateRegions[,LABEL:=NA_character_]
		for(r in 1:nrow(candidateRegions)){
			rhits <- res[CHROM == candidateRegions$CHROM[r] & 
				POS >= candidateRegions$START[r] & POS <= candidateRegions$END[r],]
			candidateRegions$POS[r] <- if(nrow(rhits) > 0) rhits$POS[which.max(rhits$LOG10P)] else NA_real_
		}		
	} else if(nrow(hits) > 0){
		hits[,POS0:=POS-1]
		x <- as.numeric(hits$POS)
		y <- hits$numCHR
		start = c(1, which(diff(y) != 0 | diff(x) <= -regionSize | diff(x) >= regionSize) + 1)
		end = c(start - 1, length(x))
		candidateRegions <- data.table(
			'CHROM'=as.character(hits$CHROM[start]),
			'START'=hits$POS[start] - regionSize/2,
			'END'=hits$POS[end] + regionSize/2,
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
			candidateRegions <- getNearestGene(input=candidateRegions,Marker="MARKER",chromosome="CHROM",
				position="POS",build=build,gene_table=gene_table,gene_data_dir=gene_data_dir,
				allow_download=allow_gene_download)
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

	# Improved adaptive thinning algorithm
	# Calculate plot parameters for visual-aware thinning - moved after plotPos calculation
	
	# Prepare plot data / two-colored chromosomes
	res[,`:=`(
		pch=20,
		highlightColor=as.character(NA),
		pcol = ifelse(numCHR %%2 == 0, "grey40","grey60"))]

	# covert pos on chromosomes to continuous plot positions
	chrGAP <- 1E7
	endPos <- 0
	plotPos <- numeric(0)
	chrLab <- numeric(0)
	for(chr in chrs){
		chrTemp <- which(res$CHROM == chr)
		chrPOS <- res$POS[chrTemp]-min(res$POS[chrTemp],na.rm=T)+endPos+1
		chrLab <- c(chrLab,mean(chrPOS,na.rm=T))
		endPos <- max(chrPOS,na.rm=T)+chrGAP
		plotPos <- c(plotPos,chrPOS)
	}
	res[,x:=plotPos]

	# Now calculate plot parameters for visual-aware thinning
	plot_width_bp <- max(1, max(plotPos) - min(plotPos))  # Total genomic span in plot units
	plot_height_log10p <- max(1, max(res$LOG10P) - min(res$LOG10P))  # Total p-value span
	
	# Estimate visual resolution (approximate pixels per point)
	# Assuming typical plot dimensions and point size
	visual_resolution_x <- plot_width_bp / 1600  # bp per pixel horizontally
	visual_resolution_y <- plot_height_log10p / 900  # log10p per pixel vertically
	
	# Point radius in plot units (considering cex=0.9, default point size ~3-4 pixels)
	point_radius_x <- visual_resolution_x * 4
	point_radius_y <- visual_resolution_y * 4
	
	# Adaptive thinning with distance-based approach
	adaptiveThinning <- function(data, sig_threshold = 6, aggressiveness = 1.0, custom_rates = NULL) {
		setorder(data, -LOG10P)  # Start with most significant
		
		# Track statistics
		original_count <- nrow(data)
		stats <- list()
		
		# Always keep highly significant variants (p < 1e-6)
		keep_indices <- which(data$LOG10P >= sig_threshold)
		if(length(keep_indices) > 0) {
			thinned <- data[keep_indices, .(CHROM, POS, LOG10P, numCHR, x)]
			stats[["genome_wide_significant"]] <- list(original = length(keep_indices), kept = length(keep_indices))
		} else {
			thinned <- data.table(CHROM=character(0), POS=numeric(0), 
								  LOG10P=numeric(0), numCHR=numeric(0), x=numeric(0))
			stats[["genome_wide_significant"]] <- list(original = 0, kept = 0)
		}
		
		if(length(keep_indices) > 0) remaining <- data[-keep_indices] else remaining <- data
		
		if(nrow(remaining) == 0) {
			stats[["total"]] <- list(original = original_count, kept = nrow(thinned), 
									proportion_kept = nrow(thinned) / original_count)
			return(list(data = thinned, stats = stats))
		}
		
		# Define thinning strata based on p-value significance
		remaining[, thin_group := cut(LOG10P, 
			breaks = c(-Inf, 1, 2, 3, 4, 5, sig_threshold, Inf),
			labels = c("very_ns", "ns", "weak", "moderate", "strong", "very_strong", "genome_wide"),
			include.lowest = TRUE)]
		
		# Configurable thinning rates based on aggressiveness
		if(!is.null(custom_rates)) {
			# Use custom rates if provided
			thin_rates <- custom_rates
		} else {
			# Base rates adjusted by aggressiveness multiplier
			base_rates <- c("very_ns" = 0.05, "ns" = 0.15, "weak" = 0.35, "moderate" = 0.55, 
							"strong" = 0.75, "very_strong" = 0.9, "genome_wide" = 1.0)
			
			# Adjust rates based on aggressiveness (higher = keep more points)
			# aggressiveness: 0.1 = very aggressive (keep fewer), 2.0 = less aggressive (keep more)
			if(aggressiveness <= 0.5) {
				# More aggressive thinning
				multiplier <- aggressiveness * 0.4  # Range: 0.04 to 0.2
				thin_rates <- pmax(base_rates * multiplier, c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0))
			} else if(aggressiveness >= 1.5) {
				# Less aggressive thinning  
				multiplier <- 1 + (aggressiveness - 1) * 0.8  # Range: 1.0 to 1.8
				thin_rates <- pmin(base_rates * multiplier, c(0.15, 0.35, 0.65, 0.85, 0.95, 1.0, 1.0))
			} else {
				# Normal range adjustment
				multiplier <- 0.5 + aggressiveness * 0.5  # Range: 1.0 to 1.25
				thin_rates <- base_rates * multiplier
			}
			
			names(thin_rates) <- names(base_rates)
		}
		
		for(group in names(thin_rates)) {
			group_data <- remaining[thin_group == group]
			if(nrow(group_data) == 0) {
				stats[[group]] <- list(original = 0, kept = 0)
				next
			}
			
			original_group_count <- nrow(group_data)
			
			if(thin_rates[group] >= 1.0) {
				# Keep all variants in this group
				group_result <- group_data[, .(CHROM, POS, LOG10P, numCHR, x)]
				thinned <- rbind(thinned, group_result)
				stats[[group]] <- list(original = original_group_count, kept = nrow(group_result))
			} else {
				# Distance-based thinning within group
				setorder(group_data, numCHR, POS)
				
				# Adjust grid resolution based on aggressiveness
				grid_multiplier <- 1.0 + (2.0 - aggressiveness) * 0.5  # More aggressive = smaller grids
				min_dist_x <- point_radius_x * 1.5 * grid_multiplier
				min_dist_y <- point_radius_y * 1.5 * grid_multiplier
				
				# Grid-based thinning approach
				group_data[, grid_x := round(x / min_dist_x)]
				group_data[, grid_y := round(LOG10P / min_dist_y)]
				group_data[, grid_key := paste(numCHR, grid_x, grid_y, sep = "_")]
				
				# Within each grid cell, keep more variants
				grid_thinned <- group_data[, {
					if(.N == 1) {
						.SD[, .(CHROM, POS, LOG10P, numCHR, x)]
					} else {
						# Always keep the most significant in this grid cell
						most_sig <- .SD[which.max(LOG10P)]
						others <- .SD[-which.max(LOG10P)]
						
						# More generous sampling based on thinning rate
						n_keep <- max(1, round(nrow(others) * thin_rates[group]))
						if(nrow(others) > 0 && n_keep > 0) {
							sampled <- others[sample(.N, min(n_keep, .N))]
							result <- rbind(most_sig, sampled)
						} else {
							result <- most_sig
						}
						result[, .(CHROM, POS, LOG10P, numCHR, x)]
					}
				}, by = grid_key]
				
				# Select only the needed columns
				grid_result <- grid_thinned[, .(CHROM, POS, LOG10P, numCHR, x)]
				thinned <- rbind(thinned, grid_result)
				stats[[group]] <- list(original = original_group_count, kept = nrow(grid_result))
			}
		}
		
		# Calculate total statistics
		stats[["total"]] <- list(original = original_count, kept = nrow(thinned), 
								proportion_kept = nrow(thinned) / original_count)
		
		return(list(data = thinned, stats = stats))
	}
	
	# Apply thinning based on mode
	if(disable_thinning) {
		cat("\n=== THINNING DISABLED ===\n")
		cat(sprintf("Using all %s variants without thinning\n", format(nrow(res), big.mark=",")))
		cat("=========================\n\n")
		plotdata <- res
		plotdata[,`:=`(
			pch=20,
			highlightColor=as.character(NA),
			pcol = ifelse(numCHR %%2 == 0, "grey40","grey60"))]
	} else if(thinning_mode == "random") {
		# Simple random thinning
		if(!is.null(random_seed) && !is.na(random_seed)) set.seed(random_seed)
		retention_threshold <- min(thinning_sig_threshold, min(yLine))
		must_keep <- which(res$LOG10P >= retention_threshold)
		n_keep <- max(1, round(nrow(res) * thinning_aggressiveness * 0.1))  # aggressiveness as percentage
		candidate_indices <- setdiff(seq_len(nrow(res)), must_keep)
		n_sample <- max(0, min(n_keep - length(must_keep), length(candidate_indices)))
		sampled <- if(n_sample > 0) sample(candidate_indices, n_sample) else integer(0)
		keep_indices <- sort(unique(c(must_keep, sampled)))
		plotdata <- res[keep_indices]
		cat(sprintf("\n=== RANDOM THINNING ===\nKept %s of %s variants (%.1f%%)\n========================\n\n",
			format(nrow(plotdata), big.mark=","), format(nrow(res), big.mark=","),
			(nrow(plotdata)/nrow(res)) * 100))
		plotdata[,`:=`(
			pch=20,
			highlightColor=as.character(NA),
			pcol = ifelse(numCHR %%2 == 0, "grey40","grey60"))]
	} else {
		# Adaptive thinning (default)
		retention_threshold <- min(thinning_sig_threshold, min(yLine))
		thinning_result <- adaptiveThinning(res, sig_threshold = retention_threshold,
											aggressiveness = thinning_aggressiveness, 
											custom_rates = custom_thinning_rates)
		plotdata <- thinning_result$data
		thinning_stats <- thinning_result$stats
		
		# Print thinning statistics
		cat(sprintf("\n=== ADAPTIVE THINNING (Aggressiveness: %.1f) ===\n", thinning_aggressiveness))
		cat(sprintf("Total variants: %s -> %s (%.1f%% kept)\n", 
				format(thinning_stats$total$original, big.mark=","),
				format(thinning_stats$total$kept, big.mark=","),
				thinning_stats$total$proportion_kept * 100))
		cat("\nBy significance group:\n")
		
		significance_labels <- c(
			"very_ns" = "Very non-significant (-log10P < 1)",
			"ns" = "Non-significant (-log10P 1-2)", 
			"weak" = "Weak (-log10P 2-3)",
			"moderate" = "Moderate (-log10P 3-4)",
			"strong" = "Strong (-log10P 4-5)",
			"very_strong" = "Very strong (-log10P 5-6)",
			"genome_wide_significant" = "Genome-wide significant (-log10P ≥ 6)"
		)
		
		for(group in names(significance_labels)) {
			if(group %in% names(thinning_stats) && thinning_stats[[group]]$original > 0) {
				cat(sprintf("  %s: %s -> %s (%.1f%% kept)\n",
					significance_labels[group],
					format(thinning_stats[[group]]$original, big.mark=","),
					format(thinning_stats[[group]]$kept, big.mark=","),
					(thinning_stats[[group]]$kept / thinning_stats[[group]]$original) * 100))
			}
		}
		cat("=============================================\n\n")
		
		plotdata[,`:=`(
			pch=20,
			highlightColor=as.character(NA),
			pcol = ifelse(numCHR %%2 == 0, "grey40","grey60"))]
	}

	if(nrow(candidateRegions) > 0 && "POS" %in% names(candidateRegions)) {
		peakKeys <- unique(candidateRegions[!is.na(POS),.(CHROM,POS)])
		if(nrow(peakKeys) > 0) {
			peakRows <- merge(res, peakKeys, by=c("CHROM","POS"))
			if(nrow(peakRows) > 0) {
				plotKeys <- paste(plotdata$CHROM, plotdata$POS, sep=":")
				peakRows <- peakRows[!(paste(CHROM, POS, sep=":") %in% plotKeys)]
				if(nrow(peakRows) > 0) plotdata <- rbind(plotdata, peakRows, fill=T)
			}
		}
	}
	if(nrow(plotdata) == 0) stopf("No markers remain to plot after thinning")
	
	plotdata <- plotdata[order(numCHR, POS)]

	# Recalculate plot positions for thinned data
	chrGAP <- 1E7
	endPos <- 0
	plotPos <- numeric(0)
	chrLab <- numeric(0)
	for(chr in chrs){
		chrTemp <- which(plotdata$CHROM == chr)
		if(length(chrTemp) > 0) {
			chrPOS <- plotdata$POS[chrTemp]-min(plotdata$POS[chrTemp],na.rm=T)+endPos+1
			chrLab <- c(chrLab,mean(chrPOS,na.rm=T))
			endPos <- max(chrPOS,na.rm=T)+chrGAP
			plotPos <- c(plotPos,chrPOS)
		} else {
			chrLab <- c(chrLab, endPos + 5E6)  # Add placeholder for empty chromosomes
			endPos <- endPos + chrGAP
		}
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
	oddChr <- seq(1,length(chrLab),by=2)
	axis(1,at=chrLab[oddChr],
		labels=chrs[1:length(chrLab)][oddChr],
		las=1,tick=F,cex.axis=1.5,cex.lab=1.5,line=2)
	if(length(chrLab) >= 2) {
		evenChr <- seq(2,length(chrLab),by=2)
		axis(1,at=chrLab[evenChr],
			labels=chrs[1:length(chrLab)][evenChr],
			las=1,tick=F,cex.axis=1.5,cex.lab=1.5,line=0)
	}
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
    help="Tab-delimited candidate regions with CHROM, START, END, COL, optional LABEL [default='']"),
  make_option("--flankingsize", type="integer", default=500000,
    help="Flanking region to determine candidate regions, top hit +/- X bp [default=500000]"),
  make_option("--chr", type="character", default="CHROM",
    help="name of column with chromosome [default='CHROM']"),
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
  make_option("--labelPeaks", type="logical", default=F,
    help="Label significant Manhattan peaks using a local gene table; may require optional Map2NCBI [default=F]"),
  make_option("--gene_table", type="character", default="",
    help="Local exon/gene table for --labelPeaks; overrides --gene_data_dir when provided [default='']"),
  make_option("--gene_data_dir", type="character", default=DEFAULT_GENE_DATA_DIR,
    help="Directory for default ExonList_<build>.txt cache used by --labelPeaks [default='%default']"),
  make_option("--allow_gene_download", type="logical", default=F,
    help="Allow UCSC download to create a missing gene table when labeling peaks [default=F]"),
  make_option("--threads", type="numeric", default=1,
    help="DTthreads"),
  make_option("--thinning_mode", type="character", default="adaptive",
    help="Thinning mode: adaptive, random, or none [default='adaptive']"),
  make_option("--thinning_aggressiveness", type="numeric", default=1.0,
    help="Thinning retention control: 0.25=preview/fewer background points, 1.0=default, 2.0=dense/more background points [default=1.0]"),
  make_option("--thinning_sig_threshold", type="numeric", default=6,
    help="Significance threshold for thinning (-log10P) [default=6]"),
  make_option("--disable_thinning", type="logical", default=F,
    help="Disable all thinning and plot all variants [default=F]"),
  make_option("--custom_thinning_rates", type="character", default=NULL,
    help="Custom thinning rates as comma-separated values: very_ns,ns,weak,moderate,strong,very_strong [default=NULL]"),
  make_option("--random_seed", type="integer", default=1,
    help="Random seed used by --thinning_mode random [default=1]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# Parse custom thinning rates if provided
custom_rates <- NULL
if(!is.null(opt$custom_thinning_rates)) {
	rates_vector <- as.numeric(strsplit(opt$custom_thinning_rates, ",")[[1]])
	if(length(rates_vector) == 6) {
		if(any(is.na(rates_vector)) || any(rates_vector < 0 | rates_vector > 1)) {
			stop("Custom thinning rates must be numeric values between 0 and 1", call. = FALSE)
		}
		custom_rates <- c("very_ns" = rates_vector[1], "ns" = rates_vector[2], 
						  "weak" = rates_vector[3], "moderate" = rates_vector[4],
						  "strong" = rates_vector[5], "very_strong" = rates_vector[6],
						  "genome_wide" = 1.0)
	} else {
		warning("Custom thinning rates must have exactly 6 values, using default rates")
	}
}

input_header <- names(fread(opt$input,nrows=0))
if(!(opt$chr %in% input_header) && opt$chr == "CHROM" && "CHR" %in% input_header) {
	warning("Using legacy chromosome column 'CHR'; prefer --chr CHROM or a CHROM column", call. = FALSE)
	opt$chr <- "CHR"
}
missing_input <- setdiff(c(opt$chr,opt$pos,opt$pvalue), input_header)
if(length(missing_input) > 0) {
	stop(sprintf("Input file is missing required column(s): %s", paste(missing_input, collapse = ", ")),
		call. = FALSE)
}

gwas <- fread(opt$input,select=c(opt$chr,opt$pos,opt$pvalue),
	colClasses=setNames("character", opt$pvalue))
png(filename = paste0(opt$prefix,"_Manhattan.png"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
tryCatch({
	candidateRegions <- ManhattanPlot(res=gwas,top.size=opt$top.size,break.top=opt$break.top,hitregion=opt$hitregion,
		chr=opt$chr,pos=opt$pos,pvalue=opt$pvalue,regionSize=opt$flankingsize,
		log10p=opt$log10p,sigthreshold=opt$sigthreshold,coltop=opt$coltop,maintitle=opt$maintitle,
		build=opt$build,labelPeaks=opt$labelPeaks,gene_table=opt$gene_table,
		gene_data_dir=opt$gene_data_dir,allow_gene_download=opt$allow_gene_download,
		DTthreads=opt$threads,thinning_mode=opt$thinning_mode,
		thinning_aggressiveness=opt$thinning_aggressiveness,thinning_sig_threshold=opt$thinning_sig_threshold,
		custom_thinning_rates=custom_rates,disable_thinning=opt$disable_thinning,
		random_seed=opt$random_seed)
}, error=function(e) {
	dev.off()
	stop(e)
})
dev.off()

if(nrow(candidateRegions) > 0) fwrite(candidateRegions,paste0(opt$prefix,"_ManhattanPlot_",
	round(opt$flankingsize/1000),"kb_CandidateRegions.txt"),sep="\t",quote=F)
