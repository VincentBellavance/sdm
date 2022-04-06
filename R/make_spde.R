make_spde <- function(mesh) {
  spde <- inla.spde2.pcmatern(mesh=mesh,
                              alpha=2,
                              prior.range=c(100, 0.01),
                              prior.sigma=c(0.5, 0.01))
  return(spde)
}

