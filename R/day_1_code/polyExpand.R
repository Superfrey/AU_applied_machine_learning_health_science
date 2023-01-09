polyExpand <-
function(X,k){
  #POLYEXPAND Perform polynomial expansion of input array.
  #
  # X is an input array, and k is the order of the polynomial expansion.
  #
  # E.g. if the function call is Xp = polyExpand(X,3) the output will be
  # Xp = [X X.^2 X.^3];
  # 13Dec2023 PMR
  Xp = data.frame(matrix(ncol = 0, nrow = nrow(X)))
  for (n in 1:k){
    dum = X^n
    names(dum) = paste0(names(X), paste0('_pol',n))
    Xp = cbind(Xp,dum)
  }
  return(Xp)
}
