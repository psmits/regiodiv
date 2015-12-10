library(rstan)
library(arm)
library(parallel)
source('../R/gts.r')
source('../R/mung.r')

get.post <- function(occ.filt, flat = TRUE) {
  out <- list()
  for(ii in seq(length(occ.filt))) {
    trys <- occ.filt[[ii]]
    a.tax.post <- b.tax.post <- c()
    a.bck.post <- b.bck.post <- c()
    for(jj in seq(nrow(trys))) {
      if(flat) {
        a.tax.prior <- b.tax.prior <- 1
        a.bck.prior <- b.bck.prior <- 1
      } else {
        if(jj == 1) {
          a.tax.prior <- b.tax.prior <- 1
          a.bck.prior <- b.bck.prior <- 1
        } else {
          a.tax.prior <- a.tax.post[jj - 1]
          b.tax.prior <- b.tax.post[jj - 1]
          a.bck.prior <- a.bck.post[jj - 1]
          b.bck.prior <- b.bck.post[jj - 1]
        }
      }

      a.tax.post[jj] <- trys[jj, 2] + a.tax.prior
      b.tax.post[jj] <- trys[jj, 3] + b.tax.prior
      a.bck.post[jj] <- trys[jj, 4] + a.bck.prior
      b.bck.post[jj] <- trys[jj, 5] + b.bck.prior
    }
    out[[ii]] <- data.frame(a.tax.post, b.tax.post,
                            a.bck.post, b.bck.post)
  }
  out
}


data.file <- list.files('../data', pattern = 'Occs')
fossil <- read.csv(paste0('../data/', data.file))
bibr <- fossil

lump.file <- list.files('../data', pattern = 'lump')
lump <- read.csv(paste0('../data/', lump.file))

occ.filt <- environ.occ(bibr, gts = rev(as.character(lump[, 2])))
occ.post <- get.post(occ.filt)  # flat = FALSE for fully updating
