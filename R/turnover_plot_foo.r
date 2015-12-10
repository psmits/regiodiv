macro.plot <- function(posterior, 
                       time = lump,
                       label = 'Pr(z(i,t - 1) = 0 | z(i, t) = 1)',
                       filename = '../doc/gradient/figure/turnover.png',
                       subtract = FALSE,
                       process = TRUE, 
                       log = FALSE,
                       rate = FALSE) {

  if(process) {
    est.turn <- llply(seq(data$nprov), function(y) 
                      Reduce(rbind, llply(posterior, function(x) x[[y]])))
  } else {
    est.turn <- posterior
  }
  est.mean <- melt(laply(est.turn, colMeans))
  names(est.mean) <- c('prov', 'year', 'div')
  if(subtract) {
    est.mean$div <- 1 - est.mean$div
  }

  est.turn <- Map(function(x) {
                  rownames(x) <- seq(nrow(x))
                  x}, est.turn)
  est.turn <- melt(est.turn)
  names(est.turn) <- c('sim', 'year', 'div', 'prov')
  est.turn$prov <- factor(est.turn$prov)
  if(subtract) {
    est.turn$div <- 1 - est.turn$div
  }

  # province names
  est.turn$prov <- mapvalues(est.turn$prov, unique(est.turn$prov), 
                             c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  est.turn$prov <- factor(est.turn$prov, levels = 
                          c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  est.mean$prov <- mapvalues(est.mean$prov, unique(est.mean$prov), 
                             c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  est.mean$prov <- factor(est.mean$prov, levels = 
                          c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))

  # make in years
  time.slice <- time[seq(from = 6, to = data$nyear + 4), ]
  est.turn$year <- mapvalues(est.turn$year, 
                             unique(est.turn$year), time.slice[, 3])
  est.mean$year <- mapvalues(est.mean$year, 
                             unique(est.mean$year), time.slice[, 3])

  if(rate) {
    width <- rev(diff(rev(lump[seq(from = 5, to = data$nyear + 4), 3])))   
    split.turn <- split(est.turn, est.turn$year)   
    for(ii in seq(length(split.turn))) {   
      ss <- split.turn[[ii]]$div   
      est.turn$div <- -log(ss) / width[ii]    
    }
  }

  turn.est <- ggplot(est.turn, aes(x = year, y = div, group = sim))
  turn.est <- turn.est + geom_line(alpha = 0.01) 
  if(!rate) {
    turn.est <- turn.est + geom_line(data = est.mean,
                                     mapping = aes(x = year, 
                                                   y = div, 
                                                   group = NULL), 
                                     colour = 'blue')
  }
  turn.est <- turn.est + scale_x_reverse() + facet_grid(prov ~ .)
  if(log) {
    turn.est <- turn.est + scale_y_continuous(trans = log10_trans())
  }
  turn.est <- turn.est + labs(x = 'time', y = label)
  ggsave(plot = turn.est, filename = filename,
         width = 10, height = 5)
}


# diversity graphs
diversity.plot <- function(mat, lump = lump, 
                           filename = '../doc/gradient/figure/est_birth.png',
                           ylab = 'log estimated births',
                           log = TRUE) {
  birth.mat <- Map(function(x) {
                   colnames(x) <- seq(ncol(x))
                   x}, mat)

  birth.mean <- llply(birth.mat, function(x) cbind(seq(ncol(x)), colMeans(x)))
  birth.mean <- Map(function(x, y) cbind(x, prov = y), x = birth.mean, y = 1:4)
  birth.mean <- Reduce(rbind, birth.mean)
  colnames(birth.mean) <- c('year', 'div', 'prov')
  birth.mean <- data.frame(birth.mean)

  birth.melt <- melt(birth.mat)
  names(birth.melt) <- c('sim', 'year', 'div', 'prov')
  birth.melt$prov <- factor(birth.melt$prov)


  # province names
  birth.melt$prov <- mapvalues(birth.melt$prov, unique(birth.melt$prov), 
                               c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  birth.melt$prov <- factor(birth.melt$prov, levels = 
                            c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  birth.mean$prov <- mapvalues(birth.mean$prov, unique(birth.mean$prov), 
                               c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))
  birth.mean$prov <- factor(birth.mean$prov, levels = 
                            c('N. Temp', 'N. Trop', 'S. Trop', 'S. Temp'))

  # make in years
  time.slice <- lump[seq(from = 6, to = data$nyear + 4), ]
  birth.melt$year <- mapvalues(birth.melt$year, 
                               unique(birth.melt$year), time.slice[, 3])
  birth.mean$year <- mapvalues(birth.mean$year, 
                               unique(birth.mean$year), time.slice[, 3])

  prov.bir <- ggplot(birth.melt, aes(x = year, y = div, group = sim))
  prov.bir <- prov.bir + geom_hline(yintercept = 0, colour = 'grey')
  prov.bir <- prov.bir + geom_line(alpha = 0.01)
  prov.bir <- prov.bir + geom_line(data = birth.mean,
                                   mapping = aes(x = year, 
                                                 y = div,
                                                 group = NULL),
                                   colour = 'blue')
  prov.bir <- prov.bir + facet_grid(prov ~ .)
  prov.bir <- prov.bir + scale_x_reverse()
  if(log) {
    prov.bir <- prov.bir + scale_y_continuous(trans = log10_trans())
  }
  prov.bir <- prov.bir + labs(x = 'time', y = ylab)
  ggsave(plot = prov.bir, filename = filename,
         width = 10, height = 5)
}
