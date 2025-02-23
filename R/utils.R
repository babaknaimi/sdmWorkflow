

# return an Evaluation object given the sdmModels as input
# get the ensemble prediction for all records, and cacluate evaluations against species data:
.getEval <- function(m,n='species',...) {
  
  en <- as.data.frame(m@data,n)
  
  obs <- en[,n]
  
  en <- ensemble(m, en,...)[,1]
  evaluates(obs,en)
}
#---------


# x can be a list with multiple groups of variables, or a vector with the names of variables:

# .getComb <- function(x,n=3) {
#   if (is.list(x)) {
#     if (length(x) < n) stop('The length of list should be equal or greater than n!')
#     if (length(x) > n) {
#       
#     }
#   }
# }