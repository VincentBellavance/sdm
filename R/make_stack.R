make_stack <- function(mesh, obs) {

        # Step 2 - Define the stochastic partial differential equation object
        spde <- inla.spde2.pcmatern(mesh=mesh,
                                    alpha=2,
                                    prior.range=c(100, 0.01),
                                    prior.sigma=c(0.5, 0.01))


        # Step 3 - Index matrix
        field <- inla.spde.make.index("field", n.spde=mesh$n)


        # Step 4 - A matrix
        xy <- as.matrix(cbind(obs$lon, obs$lat))
        colnames(xy) <- c("x", "y")
        ## For estimation
        AEst <- inla.spde.make.A(mesh, loc=xy)
        ## For prediction
        APred <- inla.spde.make.A(mesh)


        # Step 5 - Organise the A matrix into a list
        ## For estimation
        AEstlist <- list(AEst)
        ## For prediction
        APredlist <- list(APred)


        # Step 6 - Organise the effects (field variables)
        ## Estimation
        effectEst <- list(i = 1:spde$n.spde)
        ## Prediction
        effectPred <- list(i = 1:spde$n.spde)


        # Step 7 - Build stack
        ## Stack for estimation
        StackEst <- inla.stack(data=list(occurrence = as.data.frame(obs)[,"occurrence"]),
                               A = AEstlist,
                               effects = effectEst,
                               tag="est")
        ## Stack for prediction
        StackPred <- inla.stack(data=list(occurrence = NA),
                                A=APredlist,
                                effects=effectPred,
                                tag="pred")
        ## Organise StackEst and StackPred into a single stack
        Stack <- inla.stack(StackEst,StackPred)

	return(Stack)
}
