MatrixFusion <- function(true_lab,pre_lab){
  # true_lab is the true community label
  # pre_lab is the predicted community label
  r = max(true_lab)
  s = max(pre_lab)
  A = matrix(0,r,s)
  for (i in 1:r){
    for(j in 1:s){
      A[i,j] = length(intersect(which(true_lab==i),which(pre_lab==j)))
    }
  }
  return(A)
}
