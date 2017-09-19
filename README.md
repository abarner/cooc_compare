# Supplementary code for analyses of Barner et al.

Barner et al. (in review) compares the estimation of pairwise species interactions from observational co-occurrence data against species interactions estimated experimentally.

Here we provide code to:

1. Replicate the only co-occurrence method not previously implemented in R ('frequency distribution' or the COOC algorithm of Sfenthourakis et al. 2006 [doi: 10.1111/j.1466-822X.2005.00192.x](doi.org/10.1111/j.1466-822X.2005.00192.x)). All other pairwise co-occurrence analyses use methods implemented in existing R packages, detailed in the main text and Appendix S2.   

    **COOC_code_example.md** illustrates an example implementation of this algorithm   
    **COOC_algorithm.R** is the raw R code needed to execute the algorithm

2. Replicate the aggregation of spatially gridded community data, described in Appendix S1.   

    **spatial_aggregate.R** aggregates gridded presence/absence data to increasingly large grain size, using functions written by [Rich FitzJohn](https://gist.github.com/richfitz/11018949)

3. Replicate the analysis that accounts for habitat filtering in species pairwise co-occurrences, described in Appendix S2.   

    **habitat_filter.md** illustrates an example implementation of several environmental-correcton methods 

