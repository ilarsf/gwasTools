# gwasTools

A collecion of R scripts that might be useful to plot GWAS results.

The following R pakages need to be installed for running these Rscripts:

*optparse, data.table, RColorBrewer, plotrix*


## Creating a QQ plot:

Check the required/available Rscript parameters by using the following command

    Rscript QQplot.r --help


## Creating a Manhattan plot:

Check the required/available Rscript parameters by using the following command

    Rscript ManhattanPlot.r --help

## Minimal/example Input format (uncomressed / tab delimited)

|CHROM	|POS	|MAF	|PVALUE	|
|---	|---	|---	|---	|
|1  	|1  	|0.05	|0.99	|
|2  	|2   	|0.15	|0.1	|
|3  	|3  	|0.5	|0.25	|


### Examples

    Rscript --vanilla QQplot.r \
    --input ExampleGWAS.txt \
    --prefix ExampleGWAS \
    --maf MAF \
    --p PVALUE \
    --maintitle 'An Example QQ plot'


    Rscript --vanilla ManhattanPlot.r \
    --input ExampleGWAS.txt \
    --prefix ExampleGWAS \
    --chr CHROM \
    --pos POS \
    --p PVALUE \
    --maintitle 'An Example Manhattan plot'
 
