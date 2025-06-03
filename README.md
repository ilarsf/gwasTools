# gwasTools

A collection of R scripts that might be useful to plot GWAS results.

The following R pakages need to be installed for running these Rscripts:

*optparse, data.table, RColorBrewer, plotrix, pwr, Map2NCBI, basicPlotteR, devtools*
*devtools::install_github('JosephCrispell/basicPlotteR')*

## QQ Plot

Frequency binned and thinned QQ plot.

Please check out the required/available Rscript parameters by using the following command

    Rscript QQplot.r --help

<img src="/Example_QQ.png" height="480">

## Manhattan Plot

Fast Manhattan plot script with **adaptive data thinning** for handling large GWAS datasets efficiently. The script intelligently reduces data points while preserving significant associations.

**New Thinning Features:**
- `--thinning_mode`: Choose between "adaptive" (default), "random", or disable with `--disable_thinning T`
- `--thinning_aggressiveness`: Control thinning intensity (0.1=aggressive, 1.0=normal, 2.0=less aggressive)
- `--thinning_sig_threshold`: Significance threshold for keeping variants (-log10P, default=6)

Please check out the required/available Rscript parameters by using the following command

    Rscript ManhattanPlot.r --help
    
<img src="/Example_Manhattan.png" height="480">

## 80% Power Plot

80% Power plots based on Cohen's effect size calculation for proportions as implemented in the 'pwr' R package

    Rscript PowerPlot.r --help

<img src="/Example_PowerAnalysis.png" height="480">

## Minimal/example Input format for QQ and Manhattan plots

|CHROM	|POS	|MAF	|PVALUE	|
|---	|---	|---	|---	|
|1  	|1  	|0.05	|0.99	|
|2  	|2   	|0.15	|0.1	|
|3  	|3  	|0.5	|0.25	|


### Examples (run from command line):

    Rscript QQplot.r \
    --input ExampleGWAS.txt \
    --prefix Example \
    --maf MAF \
    --pvalue PVALUE \
    --maintitle 'An Example QQ plot'


    # GWAS results from http://csg.sph.umich.edu/abecasis/public/amd2015/Fritsche_2015_AdvancedAMD.txt.gz
    Rscript ManhattanPlot.r \
    --input Fritsche_2015_AdvancedAMD.txt.gz \
    --prefix Example \
    --chr Chrom \
    --pos Pos \
    --pvalue GC.Pvalue \
    --coltop T \
    --maintitle 'Age-related macular degeneration (Fritsche et al. 2016)' \
    --threads 8

    # Example with custom thinning for large datasets
    Rscript ManhattanPlot.r \
    --input large_gwas_results.txt \
    --prefix LargeGWAS \
    --thinning_aggressiveness 0.5 \
    --maintitle 'Large GWAS with aggressive thinning'
 
 
    Rscript PowerPlot.r \
    --prefix Example \
    --cases 500,1000,5000 \
    --controls 1000,2000,10000 \
    --minMAF 0.001 \
    --alpha 5E-8
