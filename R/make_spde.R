make_spde <- function(mesh) {
  spde <- inla.spde2.pcmatern(mesh=mesh,
                              alpha=2,
                              prior.range=c(2500, 0.95),
                              prior.sigma=c(20, 0.05))
  return(spde)
}

