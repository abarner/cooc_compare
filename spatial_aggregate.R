# Code to aggregate the spatially gridded data described in Barner et al. (in review) Appendix S1.
# In brief, in a 25x25 cm2 quadrat, presence data is recorded for every species in each of 25 5x5 cm2
# squares. The goal is to aggregate the data into all possible square combinations of: 10x10 cm2,
# 15x15 cm2, 20x20 cm2, and finally to 25x25 cm2, over many replicated quadrat surveys.

## Functions & aggregation procedure developed by R. FitzJohn: https://gist.github.com/richfitz/11018949

## Create the 25x25 cm2 spatial grid:
nr <- 5
grid <- matrix(seq_len(nr^2), nr, nr, byrow=TRUE)
# Row & column indices for the spatial grid
rc <- index_to_rc(seq_len(nr^2), nr)


## Simulate example community data for one quadrat, 3 species present/absent in each 5x5 grid
data <- data.frame(quad_num = rep(1:2, each = nr^2),
                   grid = rep(seq_len(nr^2), times = 2),
                   sp1 = runif(2 * (nr^2)) < .7,
                   sp2 = runif(2 * (nr^2)) < .8,
                   sp3 = runif(2 * (nr^2)) < .9,
                   stringsAsFactors = FALSE)
# Index the columns of data that record species occurrences
sp_col <- which(substr(colnames(data), start = 1, stop = 2) == "sp")
# Above data is T/F for presence/absence, but we can convert to more traditional 0/1 format
data_pa <- cbind(data[,-c(sp_col)], ifelse(data[,c(sp_col)] == TRUE, 1, 0))


## Function to convert from index within the vector to row/column index and vice-versa:
index_to_rc <- function(i, nr) {
  cbind((i - 1) %/% nr + 1, (i - 1) %% nr + 1)
}
rc_to_index <- function(rc, nr) {
  (rc[,1] - 1) * nr + rc[,2]
}

## Function to index a square sub-array with n rows and columns out of a matrix with nr rows
sub_array <- function(i, n, nr) {
  rc <- index_to_rc(i, nr)
  delta <- cbind(rep(seq_len(n), each=n), seq_len(n)) - 1L
  apply(rc, 1, function(x)
    rc_to_index(cbind(x[1] + delta[,1], x[2] + delta[,2]), nr))
}
f <- function(x) colSums(x) > 0


## Aggregate to 10x10 cm2 squares (or a "2x2" array, if 5x5 cm2 is a 1x1 array)

# Get the indices of the 2x2 array
i2 <- rc_to_index(rc[rowSums(rc <= 4) == 2,], nr)

# Create a matrix where each column is a set of indices to the grid.
sub2 <- sub_array(i2, 2, nr)
# (4 numbers, because the 2x2 array [10x10 cm2] contains four 5x5 cm2 squares,
# 16 columns because there are 16 ways to aggregate 5x5 into 10x10 within a 
# 25x25 cm square)

# Loop through all the quadrats in the data, get the presence/absence of each species in each 
res2 <- vector(mode = "list", length = max(data$quad_num))
# Loop through each quadrat, end up with 9 sub arrays of 3x3 squares
for (i in 1:max(data$quad_num)) {
  temp <- subset(data, quad_num==i) 
  res2[[i]] <- t(apply(sub2, 2, function(x) f(temp[temp$grid[x],-c(1,2)])))
}
# Name each element of the list by its quadrat
names(res2) <- 1:max(data$quad_num)

data2 <- plyr::ldply(res2, data.frame)
colnames(data2)[1] <- "quad_num"
# Index subarrays (each corresponds to a different subarray)
data2$subarray <- c(rep(1:ncol(sub2), times = max(data$quad_num)))


## Repeat: aggregate to 15x15 cm2 squares (a "3x3" array)

i3 <- rc_to_index(rc[rowSums(rc <= 3) == 2,], nr)
sub3 <- sub_array(i3, 3, nr)
res3 <- vector(mode = "list", length = max(data$quad_num))
for (i in 1:max(data$quad_num)) {
  temp <- subset(data, quad_num==i) 
  res3[[i]]<-t(apply(sub3, 2, function(x) f(temp[temp$grid[x],-c(1,2)])))
}
names(res3) <- 1:max(data$quad_num)
data3 <- plyr::ldply(res3, data.frame)
colnames(data3)[1] <- "quad_num"
data3$subarray <- c(rep(1:ncol(sub3), times = max(data$quad_num)))


## Repeat: aggregate to 20x20 cm2 squares (a "4x4" array)

i4 <- rc_to_index(rc[rowSums(rc <= 2) == 2,], nr)
sub4 <- sub_array(i4, 4, nr)
res4 <- vector(mode="list", length = max(data$quad_num))
for (i in 1:max(data$quad_num)) {
  temp <- subset(data, quad_num==i) 
  res4[[i]]<-t(apply(sub4, 2, function(x) f(temp[temp$grid[x],-c(1,2)])))
}
names(res4) <- 1:max(data$quad_num)
data4 <- plyr::ldply(res4, data.frame)
colnames(data4)[1] <- "quad_num"
data4$subarray <- c(rep(1:ncol(sub4), times = max(data$quad_num)))


## Aggregate to 25x25 cm2 square (the full quadrat)

# This is the maximum size that data can be spatially aggregated - 
# can be aggregated more easily with the 'aggregate' function,
# using the presence/absence data

data5 <- aggregate(data_pa[,c(sp_col)],
                  # Aggregate by quadrat identity
                  by = list(quad_num = data_pa$quad_num),
                  FUN = max)

