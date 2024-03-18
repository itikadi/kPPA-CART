#' Function that loads sample data
#'
#' This function loads sample data for mouse fraelty.
#' X is RNA Seq data and SexID and AgeID are numerical factors for
#' different ages and sexes.
#'
#' @return A list of X, Age, and Sex
#' @export
kPPA_sample_data <- function(){
  # read in sex and age
  SexID <- read.csv("sample/SexID.csv")
  AgeID <- read.csv("sample/AgeID.csv")
  # read in RNA seq data
  X <- read.csv("sample/X.csv")
  # return data
  return(list(X = t(X), Age = factor(AgeID$X2), Sex = factor(SexID$X1)))
}
