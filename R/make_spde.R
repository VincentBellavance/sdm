make_spde <- function(mesh) {
  spde <- inla.spde2.pcmatern(mesh=mesh,
                              alpha=2,
                              prior.range=c(1000, 0.95),
                              prior.sigma=c(15, 0.05))
  return(spde)
}

