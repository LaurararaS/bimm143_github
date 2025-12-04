# Class 06: R Functions
Laura Sun (PID: A17923552)

All functions in R have at least 3 things:

- A **name**, we pick this and use it to call the function.
- Input **arguments**, there can by multiple comma separated inputs to
  the function.
- The **body**, lines of R code that do the work of the function.

Our first wee function:

``` r
add <- function(x,y=1){
  x + y
}
```

Let’s test our function

``` r
add(c(1,2,3), y=10)
```

    [1] 11 12 13

``` r
add(10,100)
```

    [1] 110

## A second function

Let’s try something more interesting. Make a sequence generation tool.

The `sample()` function could be useful here.

``` r
sample(1:10, size=3)
```

    [1] 7 8 9

Change this to work with the nucleotides ATCG and return 3 of them.

``` r
sample(c("A", "T", "C", "G"), size=3)
```

    [1] "G" "T" "A"

Or I can name them as n and sample(n, size=3).

With the above code, I can’t generate more than 4 (sample larger than
the population). Try this instead, replace (sampling be with
replacement). Default is FALSE, we set to TRUE.

``` r
n <- c("A", "T", "C", "G")
sample(n, size=10, replace = TRUE)
```

     [1] "G" "A" "A" "C" "A" "C" "A" "A" "C" "T"

Turn this snipet into a function that returns a user specified length
DNA sequence. Let’s call it `generate_DNA`.

``` r
generate_DNA <- function(len=10){
  n <- c("A", "T", "C", "G")
  sample(n, size=len, replace = TRUE)
}
```

or

``` r
generate_DNA <- function(len=10, fasta=FALSE){
  n <- c("A", "T", "C", "G")
  v <- sample(n, size=len, replace = TRUE)
  # make a single element vector
  s <- paste(v, collapse="")
  cat("Well done!\n")
  if(fasta){
    return(s)
  }
  else{
    return(v)
  }
}
```

``` r
s <- generate_DNA(15)
```

    Well done!

``` r
s
```

     [1] "G" "C" "C" "C" "A" "T" "G" "T" "T" "T" "T" "G" "T" "T" "T"

I want the option to return a single element character vector with my
sequence all together like this: “GGAGTAC”.

``` r
s
```

     [1] "G" "C" "C" "C" "A" "T" "G" "T" "T" "T" "T" "G" "T" "T" "T"

``` r
paste(s, collapse="")
```

    [1] "GCCCATGTTTTGTTT"

``` r
generate_DNA(10, fasta=TRUE)
```

    Well done!

    [1] "CTAAAGTTTT"

## A more advanced example

Make a third function that generates protein sequence of a user
specified length and format.

``` r
proseq <- function(len=10, fasta=FALSE){
  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", 
          "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  seq <- sample(aa, size = len, replace = TRUE)
  if (fasta) {
    return(paste(seq, collapse = ""))
  }
  else {
    return(seq)
  }
}
```

``` r
proseq(15, TRUE)
```

    [1] "SSCAKMWPDRTQIFV"

> Generate random protein sequences between lengths 5 and 12 amino
> acids.

Our approach is to do this by brute force calling our function for each
length 5 to 12.

Another approach is to write a `for()` loop to iterate over the input
valued 5 to 12.

A very useful third R specific approach is to use the `saaply()`
function.

``` r
seq_length <- 6:12
for (i in seq_length) {
  cat(">", i, "\n")
  cat(proseq(i))
  cat("\n")
}
```

    > 6 
    D C W L Y R
    > 7 
    S I C P K V P
    > 8 
    E K I S V W L Q
    > 9 
    R S V W A G R D A
    > 10 
    F P T H V E F G Q L
    > 11 
    F M G K V H H L H L V
    > 12 
    E F S R S L S M N N M T

``` r
sapply(6:12, proseq)
```

    [[1]]
    [1] "K" "S" "K" "Q" "Q" "N"

    [[2]]
    [1] "D" "W" "Q" "T" "A" "E" "M"

    [[3]]
    [1] "F" "G" "M" "G" "E" "A" "C" "P"

    [[4]]
    [1] "T" "C" "S" "S" "T" "K" "A" "G" "T"

    [[5]]
     [1] "T" "I" "P" "Q" "R" "C" "R" "D" "T" "I"

    [[6]]
     [1] "H" "G" "V" "I" "A" "I" "D" "C" "W" "V" "P"

    [[7]]
     [1] "Y" "C" "C" "M" "C" "I" "D" "R" "I" "E" "F" "V"

> **Key-point**: writing functions in R is doable but not the easiest
> thing in the world. Starting with a working snippet of come and then
> using LLM tools to improve and generalize your function code is a
> productive approach.
