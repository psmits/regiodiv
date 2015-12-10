library(plyr)
library(rstan)
# processing coda into stan-like posterior draws
process.coda <- function(post, data) {
  post <- Reduce(rbind, post)
  vars <- unique(llply(str_split(colnames(post), pattern = '\\['), 
                       function(x) x[1]))

  # now break it apart
  new.post <- list()
  for(ii in seq(length(vars))) {
    temp <- post[, str_detect(colnames(post), 
                              paste0(vars[[ii]], '\\['))]
    comma.detect <- unique(str_count(colnames(temp), pattern = '\\,'))
    col.detect <- ncol(temp)
    
    if(ncol(temp) == 0) {
      temp <- post[, str_detect(colnames(post), vars[[ii]])]
    }

    if(comma.detect == 0 && col.detect == data$nprov) { 
      # goes by province ONLY
      colnames(temp) <- NULL
      new.post[[ii]] <- temp
    } else if(comma.detect == 1  && col.detect > data$nprov) { 
      # goes by time AND prov
      hold <- array(dim = c(nrow(temp), col.detect / data$nprov, data$nprov))
      
      colnames(temp) <- NULL
      for(jj in seq(data$nprov)) {
        if(jj == 1) {
          grab <- seq(from = 1, to = col.detect / data$nprov)
        } else {
          grab <- seq(from = ((col.detect / data$nprov) * (jj - 1)) + 1, 
                      to = (col.detect / data$nprov) * jj)
        }
        hold[, , jj] <- temp[, grab]
      }
      new.post[[ii]] <- hold
    } else if(comma.detect == 0 && 
              (col.detect == data$nyear | 
               col.detect == (data$nyear - 1))) { 
      # goes by time ONLY
      colnames(temp) <- NULL
      new.post[[ii]] <- temp
    } else if(is.null(dim(temp))) {
      names(temp) <- NULL
      new.post[[ii]] <- temp
    }
  }
  names(new.post) <- unlist(vars)
  new.post
}


# data for posterior predictive checks given the observed sample sizes
# post is the extraction from stan object
# data is the data used to fit the initial model
#   P: number of provinces 
#   prov: vector of province "start points"
#   C: number of time points
#   sight: sighting record (just need initial, i think)
posterior.turnover <- function(post, data) {
  # simulate for each province
  region <- list()
  for(jj in seq(data$nprov)) {
    rand <- sample(nrow(post[[1]]), 1)
    ## diversification process
    time.step <- data$nyear

    # simulate for each taxon
    life.time <- matrix(nrow = data$nindiv, ncol = data$nyear)
    life.time[, 1] <- data$y[, 1, jj]
    for(nn in seq(data$nindiv)) {
      for(ii in 2:time.step) {
        if(life.time[nn, ii - 1] == 1) {
          life.time[nn, ii] <- rbinom(1, size = 1, 
                                      prob = post$phi[rand, ii - 1, jj])
        } else {
          life.time[nn, ii] <- rbinom(1, size = 1, 
                                      prob = post$gamma[rand, ii - 1, jj])
        }
      }
    }

    sample.time <- life.time
    # observation process
    for(nn in seq(data$nindiv)) {
      for(ii in seq(from = 2, to = time.step)) {
        sample.time[nn, ii] <- life.time[nn, ii] * 
        rbinom(1, size = 1, prob = post$p[rand, ii, jj])
      }
    }

    region[[jj]] <- list(true = life.time, obs = sample.time)
  }
  region 
}


# count how many taxa actually seen from the simulation
# out is just a single provincial simulation
# set can be 'true' or 'observed'
check.count <- function(out, set = 'true') {
  chck <- c()
  if(set == 'true') {
    for(ii in seq(nrow(out$true))) {
      chck[ii] <- any(out$true[ii, ] != 0)
    }
  } else if (set == 'observed') {
    for(ii in seq(nrow(out$true))) {
      chck[ii] <- any(out$observed[ii, ] != 0)
    }
  }
  chck
}

# extract certain value types from the posterior
macro.prob <- function(data, post, ww = 'gamma') {
  samp <- sample(4000, 1)
  regions <- list()
  for(jj in seq(data$nprov)) {
    hold <- c()
    for(cc in seq(data$nyear - 1)) {
      hold[cc] <- post[[ww]][samp, cc, jj]
    }
    regions[[jj]] <- hold
  }
  regions
}

# calculate, recursively, occpuancy probability
occ.prob <- function(post, data) {
  samp <- sample(4000, 1)
  regions <- list()
  for(jj in seq(data$nprov)) {
    hold <- c()
    for(cc in seq(data$nyear)) {
      if(cc == 1) {
        psi <- post$psi[samp, jj]
        hold[cc] <- psi
      } else {
        psi <- hold[cc - 1]
        phi <- post$phi[samp, cc - 1, jj]
        gam <- post$gamma[samp, cc - 1, jj]
        hold[cc] <- psi * phi + (1 - psi) * gam
      }
    }
    regions[[jj]] <- hold
  }
  regions
}

growth.rate <- function(occupancy, data = data) {
  est <- llply(seq(data$nprov), function(y) 
               Reduce(rbind, llply(occupancy, function(x) x[[y]])))
  regions <- list()
  for(jj in seq(data$nprov)) {
    out <- matrix(ncol = ncol(est[[1]]) - 1, nrow = nrow(est[[1]]))
    for(ii in seq(nrow(est[[1]]))) {
      hold <- c()
      for(cc in seq(from = 1, to = data$nyear - 1)) {
        hold[cc] <- est[[jj]][ii, cc + 1] / est[[jj]][ii, cc]
      }
      out[ii, ] <- hold
    }
    regions[[jj]] <- out
  }
  regions
}
