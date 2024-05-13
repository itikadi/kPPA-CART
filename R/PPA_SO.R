#' Internal function to run kurtosis PPA
#' @noRd
#' @keywords internal
PPA_SO <- function(Xorig, ndim = 2, nguess = 4, orthflag = 1) {

  nsamp <- nrow(Xorig)
  nvars <- ncol(Xorig)

  set.seed(0)

  PPOUT <- list(T = array(0, dim = c(nsamp, nguess, ndim)),
                V = array(0, dim = c(nvars, nguess, ndim)),
                K = matrix(0, nrow = nguess, ncol = ndim),
                CFlag = matrix(FALSE, nrow = nguess, ncol = ndim),
                imin = numeric(ndim))

  Tstor <- array(0, dim = c(nsamp, nguess, ndim))
  Pstor <- array(0, dim = c(nvars, nguess, ndim))

  Xmu <- colMeans(Xorig)
  Xmc <- Xorig - matrix(Xmu, nrow = nsamp, ncol = nvars, byrow = TRUE)
  rk <- matrix(rankMatrix(Xmc), nrow = 1)
  svd_result <- svd(Xmc)
  U0 <- svd_result$u[,1:rk]
  S0 <- diag(svd_result$d[1:rk])
  V0 <- svd_result$v[,1:rk]

  X0 <- U0 %*% S0
  Xnow <- X0

  for (i in 1:ndim) {
    ppaOne <- ppaone(Xnow, nguess)

    PPOUT$K[,i] <- ppaOne$PPAinfo$k
    PPOUT$CFlag[,i] =ppaOne$PPAinfo$conflag
    PPOUT$imin[i]=ppaOne$imin;
    Pstor[,,i]=ppaOne$PPAinfo$P;
    Tstor[,,i]=ppaOne$PPAinfo$T;

    Tk <- c()
    Pk <- c()
    if(i > 1){
      for (k in 1:(i - 1)) {
        Tk <- cbind(Tk, Tstor[, PPOUT$imin[k], k])
        Pk <- cbind(Pk, Pstor[, PPOUT$imin[k], k])
      }
    }

    Q <- V0 %*% solve(S0) %*% t(U0)

    for (j in 1:nguess) {
      T <- cbind(Tk, Tstor[, j, i])
      V <- Q %*% T
      T <- T + as.numeric(Xmu %*% V)
      tem <- sqrt(colSums(V^2))
      V <- V / (matrix(1, nrow = nrow(V), ncol = ncol(V)) * tem)
      T <- T / (matrix(1, nrow = nrow(T), ncol = ncol(T)) * tem)
      PPOUT$T[, j, i] <- T[, i]
      PPOUT$V[, j, i] <- V[, i]
    }

    if (i < ndim) {
      T <- cbind(Tk, Tstor[, PPOUT$imin[i], i])
      P <- cbind(Pk, Pstor[, PPOUT$imin[i], i])
      Xnow <- Xnow - T[, i] %*% t(P[, i])

      if (orthflag != 0 && i > 1) {
        Ttemp <- apply(T, 1, prod)
        w <- pinv(Xnow) %*% Ttemp
        w <- w / sqrt(sum(w^2))
        t <- Xnow %*% w
        p <- (t(Xnow) %*% t) / as.numeric((t(t) %*% t))
        Xnow <- Xnow - t %*% t(p)
      }
    }
  }

  T <- matrix(0, nrow = nsamp, ncol = ndim)
  V <- matrix(0, nrow = nvars, ncol = ndim)
  kurt <- numeric(ndim)


  for (i in 1:ndim) {
    T[, i] <- PPOUT$T[, PPOUT$imin[i], i]
    V[, i] <- PPOUT$V[, PPOUT$imin[i], i]
    kurt[i] <- PPOUT$K[PPOUT$imin[i], i]
  }


  return(list(T = T, V = V, kurt = kurt, PPOUT = PPOUT))
}
