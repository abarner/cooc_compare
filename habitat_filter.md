# Example implementation of the environmental correction for co-occurrence



## Correcting pairwise co-occurrence for environmental/habitat filtering

Species pairs may or may not co-occur for a variety of reasons. Although pairwise co-occurrence is often interpreted directly as a signal of pairwise species interactions, co-occurrences themselves may be a result of habitat filtering (e.g., Blois et al. 2014, Morueta-Holme et al. 2016). Several methods have been proposed to first account for the influence of environmental factors in driving species co-occurrences and then examine the remaining "significant" associations. For several methods, accounting for the environment is implemented alongside the co-occurrence method (see the R package `BayesComm` for the JSDM residual covariance method, the R Bioconductor package `ccrepe` for the correlation methods). For the remaining methods (all constraint-based methods and partial correlation), we use the framework proposed by Blois et al. (2014; see also an implementation in Li & Waller 2016). In brief, the framework operates post-hoc on significant pairwise associations to evaluate whether negative co-occurrence is due to differences in habitat or whether positive co-occurrence is due to similarity in habitat (see Figures 1 and 2 in Blois et al. 2014 for a conceptual diagram). Note that the Blois et al. (2014) also considers how dispersal could generate patterns of positive or negative co-occurrence, though this is not considered in the present study nor implemented in this example.

## Example implementation

First, we need a site x species community matrix of the presence/absence of each species at each site (always need fewer species than sites).


```r
n_species <- 100
n_sites <- 250
data <- data.frame(matrix(ifelse(runif(n_species * n_sites) < 0.5, 1, 0), 
                          nrow = n_sites))
# helpful to have data as a data frame for later analysis
colnames(data) <- paste0("sp_", 1:n_species)
rownames(data) <- paste0("site", 1:n_sites)
```

Also need a list of the "significant" positive and negative co-occurrences for each species pair (implement this method _after_ running a co-occurrence analysis). Here we'll implement a simple co-occurrence method (combinatorics) available in R from the `cooccur` package.


```r
library(tidyverse); library(cooccur); library(vegan)
# run the combinatorics method
# note: must transpose our "site x species" data to "species x site" data for function to work
cooccur(t(data), type = "spp_site", spp_names = TRUE)$results %>%
  # select only our variables of interest: species pairs, and their
  # probabilities of being positively associated ("p_gt" < 0.05) or
  # negatively associated ("p_lt" < 0.05)
  # for more information see ?cooccur::prob.table
  select(sp1_name, sp2_name, positive = p_gt, negative = p_lt) %>%
  gather(positive, negative, key = sign, value = probability) %>%
  filter(sign == "positive" & probability < 0.05 |
           sign == "negative" & probability < 0.05) %>%
  # assign "association" based on the sign
  mutate(association = ifelse(sign == "positive", 1, -1)) %>%
  select(-sign, -probability) -> data_pairs
# place positive & negative associations in different data frames
data_pairs_neg <- filter(data_pairs, association < 0)
data_pairs_pos <- filter(data_pairs, association > 0)
```

Finally, need a list of environmental variables relevant to species occurrences in the system of interest.


```r
n_env <- 5 # number of environmental variables
env <- data.frame(matrix(runif(n_env * n_sites), 
                         nrow = n_sites))
# helpful to have data as a data frame for later analysis
colnames(env) <- paste0("env_var_", 1:n_env)
rownames(env) <- paste0("site", 1:n_sites)
```

The goal is to loop through each pair of negatively associated species, find the sites where the pair never co-occurs and determine whether the environmental characteristics are different at those sites. To test whether the environmental characteristics are different, Blois et al. (2014) recommend using ANOVA, but in many cases one might have multiple environmental characteristics of interest and we thus illustrate the use of multivariate Anova (PERMANOVA; Anderson 2001). Finally, if the sites are significantly different, remove that pair of species from the list of negatively associated species.


```r
num_neg <- nrow(data_pairs_neg)
data_pairs_neg_env <- numeric()

for (i in 1:num_neg){
  both_sp <- cbind(data[,data_pairs_neg[i,"sp1_name"]],
                 data[,data_pairs_neg[i,"sp2_name"]])
  sp1_only <- which(both_sp[,1] == 1 & both_sp[,2] == 0)
  sp2_only <- which(both_sp[,1] == 0 & both_sp[,2] == 1)
  env_tmp_1 <- env[sp1_only,]
  env_tmp_2 <- env[sp2_only,]
  if (nrow(env_tmp_1) == 0 | nrow(env_tmp_2) == 0) {
    data_pairs_neg_env[i] <- NA
  } else {
    env_tmp <- bind_rows(list(one = env_tmp_1, two = env_tmp_2), .id = "fac")
    # 'adonis' is the PERMANOVA function
    ad_tmp <- adonis(select(env_tmp, -fac) ~ env_tmp$fac, permutations = 99)
    data_pairs_neg_env[i] <- ad_tmp$aov.tab$`Pr(>F)`[1]
  }
}

# which species associations remain after removing those caused by
# environmental conditions that are significantly different?
data_pairs_neg[which(data_pairs_neg_env > 0.05),]
```

Then, for each pair of positively associated species, test whether the environment at the site of positive co-occurrence is different from the sites where both species never occur and remove those associations from subsequent analysis.


```r
num_pos <- nrow(data_pairs_pos)
data_pairs_pos_env <- numeric()

for (i in 1:num_pos){
  both_sp<-cbind(data[,data_pairs_pos[i,"sp1_name"]],
                 data[,data_pairs_pos[i,"sp2_name"]])
  both_occ <- which(both_sp[,1] == 1 & both_sp[,2] == 1)
  neith_occ <- which(both_sp[,1] == 0 & both_sp[,2] == 0)
  env_tmp_1 <- env[both_occ,]
  env_tmp_2 <- env[neith_occ,]
  if (nrow(env_tmp_1) == 0 | nrow(env_tmp_2) == 0) {
    data_pairs_pos_env[i] <- NA
  } else {
    env_tmp <- bind_rows(list(one = env_tmp_1, two = env_tmp_2), .id = "fac")
    ad_tmp <- adonis(select(env_tmp, -fac) ~ env_tmp$fac, permutations = 99)
    data_pairs_pos_env[i] <- ad_tmp$aov.tab$`Pr(>F)`[1]
  }
}

# which species associations remain after removing those caused by
# environmental conditions that are significantly different?
data_pairs_pos[which(data_pairs_pos_env > 0.05),]
```


## References

Anderson, M. J. 2001. A new method for non-parametric multivariate analysis of variance. Austral Ecology 26:32–46. doi: [10.1111/j.1442-9993.2001.01070.pp.x](https://doi.org/10.1111/j.1442-9993.2001.01070.pp.x)

Blois, J. L., N. J. Gotelli, A. K. Behrensmeyer, J. T. Faith, S. K. Lyons, J. W. Williams, K. L. Amatangelo, A. Bercovici, A. Du, J. T. Eronen, G. R. Graves, N. Jud, C. Labandeira, C. V. Looy, B. McGill, D. Patterson, R. Potts, B. Riddle, R. Terry, A. Tóth, A. Villaseñor, and S. Wing. 2014. A framework for evaluating the influence of climate, dispersal limitation, and biotic interactions using fossil pollen associations across the late Quaternary. Ecography 37:1095–1108. doi: [10.1111/ecog.00779](https://doi.org/10.1111/ecog.00779)

Li, D., and D. Waller. 2016. Long-term shifts in the patterns and underlying processes of plant associations in Wisconsin forests. Global Ecology and Biogeography 25:516-526. doi: [10.1111/geb.12432](https://doi.org/10.1111/geb.12432)

Morueta-Holme, N., B. Blonder, B. Sandel, B. J. McGill, R. K. Peet, J. E. Ott, C. Violle, B. J. Enquist, P. M. Jørgensen, and J.-C. Svenning. 2016. A network approach for inferring species associations from co-occurrence data. Ecography 39:1–12. doi: [10.1111/ecog.01892](https://doi.org/10.1111/ecog.01892)
