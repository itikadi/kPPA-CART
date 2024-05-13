#' Internal function to run kurtosis PPA
#' @noRd
ppaone <- function(X, guess) {
  r <- nrow(X)
  c <- ncol(X)
  maxcount <- 10000
  convFlag <- rep(FALSE, guess)
  cc <- rankMatrix(X)
  convlimit <- (1e-10) * cc
  svd_result <- svd(X)
  U <- svd_result$u[,1:cc]
  S <- diag(svd_result$d[1:cc])
  Vj <- svd_result$v[,1:cc]
  X <- X %*% Vj
  wall <- matrix(0, nrow = cc, ncol = guess)
  Mat2 <- diag(S)^2
  VM <- matrix(0, nrow = cc * cc, ncol = r)

  for (i in 1:r) {
    tem <- X[i,] %*% t(X[i,])
    VM[, i] <- as.vector(tem)
  }

  for (k in 1:guess) {
    w <- rnorm(cc)
    w <- w / sqrt(sum(w^2))
    oldw <- w
    count <- 0

    while (TRUE) {
      count <- count + 1
      x <- X %*% w
      Mat1 <- rowSums(VM %*%(x*x))
      Mat1 <- matrix(Mat1, nrow = cc, ncol = cc)
      w <- solve(Mat1, (Mat2 * w))
      w <- w / sqrt(sum(w^2))
      L1 <- (t(w) %*% oldw)^2

      if ((1 - L1) < convlimit) {
        convFlag[k] <- TRUE
        break
      } else if (count > maxcount) {
        convFlag[k] <- FALSE
        break
      }

      w <- w + 0.5 * oldw
      w <- w / sqrt(sum(w^2))
      oldw <- w
    }

    wall[, k] <- w
  }

  PPAInfoT <- X%*%wall
  Pj <- ((t(X) %*% PPAInfoT))/matrix(rep(1,dim(X)[2])) %*% colSums(PPAInfoT**2)
  PPAinfo <- list(T = X %*% wall, V = Vj %*% (t(X) %*% PPAInfoT / sum(PPAInfoT^2)),
                  k = kurtosis(PPAInfoT, na.rm = TRUE), conflag = convFlag,
                  P = Vj %*% Pj)

  kmin <- min(PPAinfo$k)
  indxmin <- which(PPAinfo$k == kmin)
  Tmin <- PPAInfoT[, indxmin, drop = FALSE]

  return(list(Tmin = Tmin, kmin = kmin, imin = indxmin, PPAinfo = PPAinfo))
}
