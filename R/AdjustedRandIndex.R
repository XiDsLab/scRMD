
CalculateARI <- function(true_lab,pre_lab){
  # true_lab is the true community label
  # pre_lab is the predicted community label
  # return the adjusted rand index
  A = MatrixFusion(true_lab,pre_lab)
  a = apply(A,1,sum)
  b = apply(A,2,sum)
  n = length(true_lab)
  c1 = sum(a*(a-1))/2
  c2 = sum(b*(b-1))/2
  c0 = sum(A*(A-1))/2
  ari = (n*(n-1)*c0/2-c1*c2)/(n*(n-1)*(c1+c2)/4-c1*c2)
  return(ari)
}