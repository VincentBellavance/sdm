#
#
#
#
#

filter_abs <- function(obs, perc) {

  abs <- obs[obs$occurrence == 0,]

  abs <- abs[sample(nrow(abs), nrow(abs)*perc),]

  obs <- rbind(obs[obs$occurrence == 1,], abs)

  return(obs)

}