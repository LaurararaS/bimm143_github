# Class11: Structural Bioinformatics
Laura Sun (PID: A16923552)

## AlphaFold Data Base (AFBD)

The EBI maintains the largest databaes of alphafold prediction models
at: https://alphafold.ebi.ac.uk

From last class we saw that the PDB had 244,290 (October 2025)

The total number of protein sequences in UniProtKB is 199,579,901

> Key point: this is a tiny fraction o f sequence space that has
> structural coverage (0.12%)

``` r
244290/199579901 * 100
```

    [1] 0.1224021

AFDB is attempting to address this gap…

There are two “quality scores” from AlphaFold. One for residues
(i.e. each amino acid) called **pLDDT** score, the other called **PAE**
Predicted Aligned Error measures the confidence in the relative position
of two residues (i.e. a score for every pair of residues).

## Generating your own structure predicitons

https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb

![](HIVPR_DIMER_superimposed.png) of all 5 superimposed

And the top ranked model colored by chain.

![](HIVPR_DIMER_m1.png)

pLDDT score for model 1:

![](HIVPR_pLDDT_m1.png)

and model 5

![](HIVPR_pLDDT_m5.png)

## Custom Analysis of Resulting Models in R

Read key result files into R. The first thing I need to know is what my
results directory/folder is called (i.e. its name is different for every
AlphaFold run/job)

``` r
results_dir <- "HIVPR_dimer_23119/"

pdb_files <- list.files(path=results_dir,
                        pattern="*.pdb",
                        full.names = TRUE)

# Print our PDB file names
basename(pdb_files)
```

    [1] "HIVPR_dimer_23119_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb"
    [2] "HIVPR_dimer_23119_unrelaxed_rank_002_alphafold2_multimer_v3_model_4_seed_000.pdb"
    [3] "HIVPR_dimer_23119_unrelaxed_rank_003_alphafold2_multimer_v3_model_1_seed_000.pdb"
    [4] "HIVPR_dimer_23119_unrelaxed_rank_004_alphafold2_multimer_v3_model_5_seed_000.pdb"
    [5] "HIVPR_dimer_23119_unrelaxed_rank_005_alphafold2_multimer_v3_model_3_seed_000.pdb"

``` r
library(bio3d)

m1 <- read.pdb(pdb_files[1])
m1
```


     Call:  read.pdb(file = pdb_files[1])

       Total Models#: 1
         Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)

         Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
         Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)

         Non-protein/nucleic Atoms#: 0  (residues: 0)
         Non-protein/nucleic resid values: [ none ]

       Protein sequence:
          PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
          QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
          ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
          VNIIGRNLLTQIGCTLNF

    + attr: atom, xyz, calpha, call

``` r
head(m1$atom)
```

      type eleno elety  alt resid chain resno insert       x     y     z o     b
    1 ATOM     1     N <NA>   PRO     A     1   <NA> -17.688 3.373 4.414 1 90.56
    2 ATOM     2    CA <NA>   PRO     A     1   <NA> -17.859 2.959 3.020 1 90.56
    3 ATOM     3     C <NA>   PRO     A     1   <NA> -17.000 1.748 2.658 1 90.56
    4 ATOM     4    CB <NA>   PRO     A     1   <NA> -17.422 4.195 2.230 1 90.56
    5 ATOM     5     O <NA>   PRO     A     1   <NA> -16.078 1.391 3.404 1 90.56
    6 ATOM     6    CG <NA>   PRO     A     1   <NA> -16.641 5.012 3.207 1 90.56
      segid elesy charge
    1  <NA>     N   <NA>
    2  <NA>     C   <NA>
    3  <NA>     C   <NA>
    4  <NA>     C   <NA>
    5  <NA>     O   <NA>
    6  <NA>     C   <NA>

``` r
plot(m1$atom$b[m1$calpha], typ="l", ylim=c(0,100))
```

![](class11_files/figure-commonmark/unnamed-chunk-5-1.png)

## Residue conservation from alignment file

Find the large AlphaFold alignment file

``` r
aln_file <- list.files(path=results_dir,
                       pattern=".a3m$",
                        full.names = TRUE)
aln_file
```

    [1] "HIVPR_dimer_23119//HIVPR_dimer_23119.a3m"

Read this into R

``` r
aln <- read.fasta(aln_file[1], to.upper = TRUE)
```

    [1] " ** Duplicated sequence id's: 101 **"
    [2] " ** Duplicated sequence id's: 101 **"

How many sequences are in this alignment

``` r
dim(aln$ali)
```

    [1] 5397  132

We can score residue conservation in the alignment with the conserv()
function.

``` r
sim <- conserv(aln)
```

``` r
plotb3(sim[1:99], ylab="Conservation Score")
```

![](class11_files/figure-commonmark/unnamed-chunk-10-1.png)

``` r
con <- consensus(aln, cutoff = 0.9)
con$seq
```

      [1] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [19] "-" "-" "-" "-" "-" "-" "D" "T" "G" "A" "-" "-" "-" "-" "-" "-" "-" "-"
     [37] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [55] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [73] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
     [91] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [109] "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-" "-"
    [127] "-" "-" "-" "-" "-" "-"
