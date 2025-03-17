#' Function to Execute KPPACart
#'
#' This function executes iterative KPPA with classification and regression trees as outlined in the paper "...".
#'
#' @param X Data matrix of features by samples
#' @param n_iterations The number of iterations to run (n_iterations * n_features > n_total_features)
#' @param n_features The number of features to choose in each iteration (n_iterations * n_features > n_total_features)
#' @param k_dim Dimensions used in MDS
#' @param exp_clusters Number of expected clusters
#' @param n_cores Number of cores for multi-core execution
#' @param kppa_dim Dimensions used in each KPPA iteration. The higher, the slower the algorithm runs.
#'
#' @return list of scores, important features/data, and cluster assignments for each sample
#' @export
KPPACart <- function(X,n_features=100,
                        n_iterations=NULL,
                        k_dim=10,
                        exp_clusters=4,
                        kppa_dim = 2,
                        n_cores = 4){

  # if n_iterations is NULL, set it to the number of features
  if(is.null(n_iterations)){
    n <- dim(X)[1]/n_features
    res <- sum(1/(1:n))  # Harmonic series summation
    result <-ceiling(res * n)  # Ceiling function (round up to the nearest integer)
    n_iterations <- result * n_features
    cat("Number of iterations set to: ", n_iterations, "\n")
  }

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
  if(foreach::getDoParRegistered() == FALSE){
    stop("Cluster not registered. \n Please make sure your clusters are accessible to R.")
  }else{
    cat("Clusters registered.\n")
  }

  # start message
  cat("Executing KPPACart.\n")

  # Create a progress bar
  pb <- progress_bar$new(
    format = "  Progress [:bar] :percent eta: :eta",
    total = n_iterations,
    width = 60
  )

  # Create a wrapper function to update the progress bar
  progress_wrapper <- function(...) {
    pb$tick()
    return(TRUE)
  }

  # Register the progress bar to be used with %dopar%
  opts <- list(progress = progress_wrapper)


  # Initialize a list to store the results
  res <- list()

  # for loop
  res <- foreach(it=1:n_iterations, .options.snow = opts, .packages = c("randomForest", "Matrix", "moments", "KPPACart")) %dopar% {

    # set seed
    set.seed(it)

    # randomly sample features
    samp <- sample(1:dim(X)[1], n_features)

    # select random features
    samp <- X[samp,]

    # do an initial kPPA run
    orig_mds <- cmdscale(dist(t(samp)), k = k_dim)
    orig_ppa <- KPPACart:::PPA_SO(orig_mds, kppa_dim)

    # extract kurtosis
    orig_kurt <- orig_ppa$kurt

    # kluster based on 4 groups in this case
    klust <- kmeans(orig_ppa$T, exp_clusters, nstart=exp_clusters)

    # samp contains actual data
    rf <- randomForest(x = t(samp), y = factor(klust$cluster))

    # return everything
    return(list(sum(orig_kurt),rownames(rf$importance),as.vector(rf$importance)))

  }

  # close cluster
  stopCluster(compute.clust)

  # message that is done
  cat("Done executing KPPACart. \nProcessing results.\n")

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
    dplyr::group_by(features) %>%
    dplyr::summarize(total_importance = mean(importances), sd_importance = sd(importances)) %>%
    #arrange(total_importance)
    dplyr::arrange(desc(total_importance))

  # check dim df
  if(dim(X)[1] != dim(imp_df)[1]){
    warning("There are features that have not been sampled due to the number of iterations and features per iteration. \n",
            "Consider increasing the number of iterations or features per iteration. \n",
            "This may affect the final results. \n",
            "Sampled: ", dim(imp_df)[1], " Total: ", dim(X)[1])
  }

  # exctract n_top features
  n_top <- n_features
  top <- imp_df$features[1:n_top]

  # get data
  top_data <- X[top,]

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

  # select solution
  best.solution <- solutions[low.idx][[1]]

  # create clusters so we can see if they overlap
  klust <- kmeans(best.solution, exp_clusters, nstart=exp_clusters)

  # run random forest with bestData
  rf <- randomForest(x = t(top_data), y = factor(klust$cluster))

  # return best solution adn best data
  return(
    list(
      T = best.solution,
      bestData = top_data,
      assignedClusters = klust$cluster,
      allData = imp_df
    )
  )

}

