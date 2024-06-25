calculate_uHe <- function(ns, Hes){
  uHes <- c()
  for(i in 1:length(ns)){
    uHes[i] <- (2*ns[i]/(2*ns[i]-1))*Hes[i]
  }
  uHe <- mean(uHes,  na.rm=TRUE) 
  return(uHe)
}