library(arm)
library(plyr)
library(parallel)
library(rjags)
library(coda)
library(rstan)
load('../data/gradient_setup.rdata')

n <- 4
sight <- sight[1:n]

sizes <- laply(sight, dim)
time.bin <- as.numeric(sizes[1, 2])
prov.size <- sizes[, 1]

taxon <- as.numeric(as.factor(Reduce(c, llply(sight, rownames))))
ntaxa <- length(unique(taxon))

total.names <- sort(unique(Reduce(c, llply(sight, rownames))))

total.occur <- array(0, c(ntaxa, time.bin, length(sight)))
for(ii in seq(length(sight))) {
  total.occur[match(rownames(sight[[ii]]), total.names), , ii] <- sight[[ii]]
}

data <- list(nindiv = nrow(total.occur), 
             nyear = ncol(total.occur), 
             nprov = length(sight), 
             y = total.occur)
save(data, file = '../data/data_dump/occurrence_data.rdata')
with(data, {dump(c('nindiv', 'nyear', 'nprov', 'y'), 
                 file = '../data/data_dump/occurrence_dump.R')})
p_norm <- array(0.5, dim = c(data$nyear, data$nprov))
gamma_norm <- phi_norm <- array(0.5, dim = c(data$nyear - 1, data$nprov))
z <- data$y
psi <- rep(0.5, data$nprov)
dump(c('gamma_norm', 'phi_norm', 'p_norm', 'z', 'psi'), 
     file = '../data/data_dump/hmm_inits.R')


# 80/20 split
train <- sort(sample(nrow(total.occur), (nrow(total.occur) / 100) * 80))
train.occur <- total.occur[train, , ]
test.occur <- total.occur[-(train), , ]

train.data <- list(nindiv = nrow(train.occur), 
                   nyear = ncol(train.occur), 
                   nprov = length(sight), 
                   y = train.occur)

test.data <- list(nindiv = nrow(test.occur), 
                  nyear = ncol(test.occur), 
                  nprov = length(sight), 
                  y = test.occur)
save(train.data, test.data, file = '../data/data_dump/occurrence_holdout.rdata')
with(train.data, {dump(c('nindiv', 'nyear', 'nprov', 'y'), 
                       file = '../data/data_dump/occurrence_train.R')})
p_norm <- array(0.5, dim = c(train.data$nyear, train.data$nprov))
gamma_norm <- phi_norm <- array(0.5, dim = c(train.data$nyear - 1, 
                                             train.data$nprov))
z <- train.data$y
psi <- rep(0.5, train.data$nprov)
dump(c('gamma_norm', 'phi_norm', 'p_norm', 'z', 'psi'), 
     file = '../data/data_dump/hmm_inits_train.R')
