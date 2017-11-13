# gwasTools

A collecion of R scripts that might be useful to plot GWAS results.

The following R pakages need to be installed for running these Rscripts:

*optparse, data.table, RColorBrewer, plotrix*


## Creating a QQ plot:

Check the required options by using the following command

    Rscript QQplot.r --help


## Creating a Manhattan plot:

Check the required options by using the following command

    Rscript ManhattanPlot.r --help
        
### Examples

    Rscript --vanilla QQplot.r \
    --input ExampleGWAS.txt \
    --prefix ExampleGWAS \
    --maf MAF \
    --p LOG10P_GC \
    --log10p T \
    --maintitle 'An Example QQ plot'


    Rscript --vanilla ManhattanPlot.r \
    --input ExampleGWAS.txt \
    --prefix ExampleGWAS \
    --chr CHROM \
    --pos POS \
    --p LOG10P_GC \
    --log10p T \
    --maintitle 'An Example Manhattan plot'
 
