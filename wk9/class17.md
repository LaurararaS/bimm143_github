# Genome Informatics
Laura Sun (PID: A17923552)

# Section 1: Proportion of G/G in a population

Download a csv from Ensemble \<
https://useast.ensembl.org/Homo_sapiens/Variation/Sample?db=core;r=17:39876388-40046389;v=rs8067378;vdb=variation;vf=959672880#373531_tablePanel
\>

Here we read csv file for MXL:

``` r
mxl <- read.csv("wk9/373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
head(mxl)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  NA19648 (F)                       A|A ALL, AMR, MXL      -
    2                  NA19649 (M)                       G|G ALL, AMR, MXL      -
    3                  NA19651 (F)                       A|A ALL, AMR, MXL      -
    4                  NA19652 (M)                       G|G ALL, AMR, MXL      -
    5                  NA19654 (F)                       G|G ALL, AMR, MXL      -
    6                  NA19655 (M)                       A|G ALL, AMR, MXL      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
table(mxl$Genotype..forward.strand.)
```


    A|A A|G G|A G|G 
     22  21  12   9 

``` r
table(mxl$Genotype..forward.strand.)/nrow(mxl)*100
```


        A|A     A|G     G|A     G|G 
    34.3750 32.8125 18.7500 14.0625 

Let’s look at a different population:

``` r
gbr <- read.csv("wk9/373522-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```

``` r
table(gbr$Genotype..forward.strand.)
```


    A|A A|G G|A G|G 
     23  17  24  27 

``` r
round(table(gbr$Genotype..forward.strand.)/nrow(gbr) * 100,2)
```


      A|A   A|G   G|A   G|G 
    25.27 18.68 26.37 29.67 

This variant that is associated with child asthma is more frequent in
the GBR population than MKL population.

Let’s dig further:

## Section 4: Population Scale Analysis \[HOMEWORK\]

One sample is obviously not enough to know what is happening in a
population. You are interested in assessing genetic differences on a
population scale.

How many samples do we have?

``` r
expr <- read.table("wk9/rs8067378_ENSG00000172057.6.txt")
head(expr)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
nrow(expr)
```

    [1] 462

Sample size of each genotype:

``` r
table(expr$geno)
```


    A/A A/G G/G 
    108 233 121 

Median expression levels of each genotype:

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
expr %>% group_by(geno) %>% summarize(median=median(exp))
```

    # A tibble: 3 × 2
      geno  median
      <chr>  <dbl>
    1 A/A     31.2
    2 A/G     25.1
    3 G/G     20.1

``` r
library(ggplot2)
```

Let’s make a boxplot:

``` r
ggplot(expr) + aes(geno, exp, fill=geno) + geom_boxplot(notch=TRUE)
```

![](class17_files/figure-commonmark/unnamed-chunk-12-1.png)

From this plot, we see that the relative expression levels are A\|A \>
A\|G \> G\|G, with a relatively big difference. SNP does effect the
expression of ORMDL3.
