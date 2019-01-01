svt <- function(A, thresh, econ = 0) {
  # econ is a flag that determines whether to use fast.svd in the "corpcor" package
  # or to use svds in the "RSpectra" pacakge
  if (econ == 0){
    res <- fast.svd(A)
    A.svt <- res$u %*% diag(pmax(res$d - thresh, 0)) %*% as.matrix(t(res$v))
  }
  else if (econ == 1) {
    A.d <- svd(A, 0, 0)$d
    k <- sum(A.d > thresh)
    res <- svds(A, k)
    if (k > 1) {
      d <- diag(pmax(res$d - thresh, 0))
    } else {
      d <- res$d
    }
    A.svt <- res$u %*% d %*% t(res$v)
  }
  return(list(A.svt=A.svt, r=k))
}