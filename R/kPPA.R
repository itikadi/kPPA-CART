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

#' Function to Execute KPPACart
#'
#' This function executes iterative KPPA with classification and regression trees as outlined in the paper "...".
#'
#' @param X Data matrix of samples by features
#' @param n_iterations The number of iterations to run (n_iterations * n_features > n_total_features)
#' @param n_features The number of features to choose in each iteration (n_iterations * n_features > n_total_features)
#' @param k_dim Dimensions used in MDS
#' @param exp_clusters Number of expected clusters
#' @param n_cores Number of cores for multi-core execution
#'
#' @return list of scores, important features/data, and cluster assignments for each sample
#' @export
KPPACart <- function(X,n_features=100,
                        n_iterations=3000,
                        k_dim=10,
                        exp_clusters=4,
                        kppa_dim = 2,
                        n_cores = 4){

  # create array of importances and features
  importances <- c()
  features <- c()
  kurtosis <- c()

  # set up pararlleization
  n.cores <- n_cores

  # create cluster
  compute.clust <- makePSOCKcluster(n.cores)

  # register it to be used by %dopar%
  registerDoSNOW(compute.clust)

  # check if it is registered (optional)
  print("")
  print(paste("Cluster detected by doParallel: ", foreach::getDoParRegistered()))

  # for loop
  res <- foreach(it=1:n_iterations, .packages = c("randomForest", "Matrix", "moments")) %dopar% {

    # set seed
    set.seed(it)

    # randomly sample features
    samp <- sample(1:dim(X)[1], n_features)

    # select random features
    samp <- X[samp,]

    # do an initial kPPA run
    orig_mds <- cmdscale(dist(t(samp)), k = k_dim)
    orig_ppa <- PPA_SO(orig_mds, kppa_dim)

    # extract kurtosis
    orig_kurt <- orig_ppa$kurt
    # kurtosis <- append(kurtosis, sum(orig_kurt))

    # kluster based on 4 groups in this case
    klust <- kmeans(orig_ppa$T, exp_clusters, nstart=exp_clusters)

    # samp contains actual data
    rf <- randomForest(x = t(samp), y = factor(klust$cluster))

    # get feature importance
    # features <- append(features,rownames(rf$importance))
    # importances <- append(importances,as.vector(rf$importance))


    # return everything
    return(list(sum(orig_kurt),rownames(rf$importance),as.vector(rf$importance)))

  }

  print("DONE!")

  # rework data
  for(i in 1:n_iterations){
    features <- append(features, res[[i]][2])
    importances <- append(importances, res[[i]][3])
    kurtosis <- append(kurtosis, res[[i]][1])
  }

  # unlist
  features <- unlist(features)
  importances <- unlist(importances)
  kurtosis <- unlist(kurtosis)

  # create dataframe from importances and features
  imp_df <- data.frame(importances,features)

  # group by features and sort descending
  imp_df <- imp_df %>%
    group_by(features) %>%
    dplyr::summarize(total_importance = mean(importances), sd_importance = sd(importances)) %>%
    #arrange(total_importance)
    arrange(desc(total_importance))

  # check dim df
  if(dim(X)[1] != dim(imp_df)[1]){
    print("Please run more iterations or use more features per iteration.")
    print(dim(imp_df)[1])
    print(dim(X)[1])
  }

  # exctract n_top features
  n_top <- n_features
  top <- imp_df$features[1:n_top]

  # get data
  top_data <- X[rownames(X) %in% top,]

  # keep track of solutions
  solutions = c()

  # run PPA multiple times and find the solution that is most similar
  for(i in 1:10){
    # do MDS
    mds.top <- cmdscale(dist(t(top_data)), k = k_dim)
    # apply kPPA to MDS
    kppa.top <- PPA_SO(mds.top, ndim=kppa_dim)
    # add solution to array
    solutions <- append(solutions, list(kppa.top$T))
  }

  # create empty matrxi of required size
  similarity <- matrix(rep(0,100),10,10)

  # calculate similarity between each solution and it
  for(i in 1:10){
    for(j in 1:10){
      # calculate mds of both solutions of interest
      mds.1 <- cmdscale(dist(solutions[i][[1]]), k = 2)
      mds.2 <- cmdscale(dist(solutions[j][[1]]), k = 2)
      similarity[i,j] = sum(abs(mds.1-mds.2))
    }
  }

  # calculate rowSums
  row.sum <- rowSums(similarity)

  # select index of lowest
  low.idx <- which.min(row.sum)

  # eslect solution
  best.solution <- solutions[low.idx][[1]]

  # create clusters so we can see if they overlap
  klust <- kmeans(best.solution, exp_clusters, nstart=exp_clusters)

  # return best solution adn best data
  return(
    list(
      T = best.solution,
      BestData = top_data,
      assignedClusters = klust$cluster,
      allData = imp_df
    )
  )

}

