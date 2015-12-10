library(rstan)
library(arm)
library(parallel)
source('../R/gts.r')
source('../R/mung.r')

data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil

lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))

shape <- readShapeSpatial('../data/ne_10m_coastline.shp')  # from natural earth

sight <- space.time(bibr, 
                    bins = 'StageNewOrdSplitNoriRhae20Nov2013', 
                    gts = rev(as.character(lump[, 2])),
                    cuts = 'Chang',
                    bot = 'Trem',
                    shape = shape)

sizes <- laply(sight, dim)
time.bin <- sizes[1, 2]
prov.size <- sizes[, 1]

taxon <- as.numeric(as.factor(Reduce(c, llply(sight, rownames))))
ntaxa <- length(unique(taxon))

# R rows
# C columns
# T taxa
# P provinces
#
# matrix[R, C] sight;
# vector[R] taxon;
# vector[P] prov;
data <- list(R = sum(sizes[, 1]), C = time.bin, 
             T = ntaxa, P = length(prov.size),
             sight = Reduce(rbind, sight),
             taxon = taxon,
             prov = prov.size)

with(data, {stan_rdump(list = c('R', 'C', 'T', 'P', 'sight', 'taxon', 'prov'),
                       file = '../data/data_dump/sight_info.data.R')})
