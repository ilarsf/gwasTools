# gwasTools

A collection of R scripts that might be useful to plot GWAS results.

The following R pakages need to be installed for running these Rscripts:

*optparse, data.table, RColorBrewer, plotrix*


## QQ plot:

Frequency binned and thinned QQ plot.

Please check out the required/available Rscript parameters by using the following command

    Rscript QQplot.r --help

![QQ plot](/Example_QQ.png)

## Manhattan plot:

Fast Manhattan plot script. Please check out the required/available Rscript parameters by using the following command

    Rscript ManhattanPlot.r --help

![Manhattan plot](/Example_Manhattan.png)

## 80% Power plot(s):

80% Power plots based on Cohen's effect size calculation for proportions as implemented in the 'pwr' R package

    Rscript --vanilla PowerPlot.r

![Manhattan plot](/Example_PowerAnalysis.png)

## Minimal/example Input format for QQ and Manhattan plots

|CHROM	|POS	|MAF	|PVALUE	|
|---	|---	|---	|---	|
|1  	|1  	|0.05	|0.99	|
|2  	|2   	|0.15	|0.1	|
|3  	|3  	|0.5	|0.25	|


### Examples

    Rscript --vanilla QQplot.r \
    --input ExampleGWAS.txt \
    --prefix Example \
    --maf MAF \
    --pvalue PVALUE \
    --maintitle 'An Example QQ plot'


    Rscript --vanilla ManhattanPlot.r \
    --input ExampleGWAS.txt \
    --prefix Example \
    --chr CHROM \
    --pos POS \
    --pvalue PVALUE \
    --maintitle 'An Example Manhattan plot'
 
 
    Rscript --vanilla PowerPlot.r \
    --prefix Example \
    --cases 500,1000,5000 \
    --controls 1000,2000,10000 \
    --minMAF 0.001 \
    --alpha 5E-8 
