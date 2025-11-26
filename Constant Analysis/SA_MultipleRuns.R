# This is where we run the simulation multiple times to do the Sobol SA

multiRun_CulexOnly <- function(lhs_sample) {
  return(mapply(singleSim_OnlyCulex,
                lhs_sample[,1], lhs_sample[,2],
                lhs_sample[,3], lhs_sample[,4],
                lhs_sample[,5], lhs_sample[,6]))
}



multiRun_AedesOnly <- function(lhs_sample) {
  return(mapply(singleSim_OnlyAedes,
                lhs_sample[,1], lhs_sample[,2],
                lhs_sample[,3], lhs_sample[,4],
                lhs_sample[,5], lhs_sample[,6]))
}