
library(vegan)

# this generates example data, or data can be read in as a presence/absence site x species matrix
sites <- 100
species <- 5
p_presence <- 0.9
DATA <- matrix(nrow = sites,
               ncol = species, rbinom(sites * species, 1, p_presence))
colnames(DATA) <- letters[1:species]

# get all possible pairs of species in the data
combpairs <- function(x) {apply(combn(x, 2), 2, paste, sep = "", collapse = "")}
allpairs <- combpairs(colnames(DATA))

# create empty data frame
NR <- length(allpairs)
coocbyrow <- matrix(nrow = NR, ncol = nrow(DATA))
row.names(coocbyrow) <- allpairs
# fill with observed co-occurrences
for (i in 1:nrow(DATA)) {
  # which combinations are present in this site (row)
  temp <- combpairs(colnames(DATA)[which(DATA[i, ] == 1)])
  coocbyrow[, i] <- ifelse(row.names(coocbyrow) %in% temp, 1, 0)
}

# tally pairwise co-occurrences by site
obs.cooc <- apply(X = coocbyrow, MARGIN = 1, FUN = sum)

# compare with a null model of pairwise co-occurrences
nperm <- 500
nul <-permatswap(DATA, fixedmar = "both", mtype = "prab",
                 # use sequential swap algorithm
                 method = "swap", 
                 times = nperm)

# for each one, count the null co-occurrences
expc.cooc <- matrix(nrow = NR, ncol = nperm)
# note: below nested for-loop is slow
for (i in 1:nperm) {
  tmp <- nul$perm[[i]]
  coocbyrow <- matrix(nrow = NR, ncol = nrow(DATA))
  row.names(coocbyrow) <- allpairs
  for (j in 1:nrow(tmp)) {
    temp <- combpairs(colnames(tmp) [which(tmp[j, ] == 1)])
    coocbyrow[, j] <- ifelse(row.names(coocbyrow) %in% temp, 1, 0)
  }
  expc.cooc[, i] <- apply(X = coocbyrow, MARGIN = 1, FUN = sum)
}
row.names(expc.cooc) <- allpairs

# calculate 'significance' of the difference between observed & expected pairwise co-occurrences
