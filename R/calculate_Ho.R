calculate_Ho <- function(gt){
  heteros <- sum(gt==1, na.rm=TRUE)
  total <- sum(!is.na(gt))
  Ho <- heteros/total
  return(Ho)
}