#' @title Function generating the dynamic geometric grid
#' @param t Number of observations
#' @return A numeric vector containing the generated grid values
#' @examples
#' get_grid(10)
#' @export
get_grid <- function(t) {
  # if(t==1){
  #   return(c(1))
  # }else{
  #   return(1:(t-1))}
  upper1 <- floor(log((t - 1) / 3, base = 2)) + 1
  upper2 <- floor(log(t - 1, base = 2)) - 1
  len <- 1 + max(0, upper1) + max(0, upper2)
  Gt <- rep(NA, len)
  Gt[1] <- 1
  counter <- 2
  if (upper1 >= 1) {
    for (j in 1:upper1) {
      gL <- 2^j + (t - 1) %% 2^(j - 1)
      Gt[counter] <- gL
      counter <- counter + 1
      if (j <= upper2) {
        Gt[counter] <- gL + 2^(j - 1)
        counter <- counter + 1
      }
    }
  }

  return(Gt)
}

#' @title Function generating penalty and threshold values for
#' the mean change detector
#' @param t An integer representing the sample size
#' @param p An integer representing the data dimension
#' @return A list containing the following elements:
#' \describe{
#'   \item{as}{A numeric vector of threshold values}
#'   \item{nu_as}{Mean centering terms}
#'   \item{r_values}{A numeric vector of penalty values}
#'   \item{ss}{A numeric vector of sparsity values}
#' }
#' @examples
#' t <- 100
#' p <- 10
#' result <- mean_penalty_and_threshold(t, p)
#' print(result)
#'
#' @export
mean_penalty_and_threshold <- function(t, p) {
  max_s <- min(sqrt(p * log(t)), p)
  log2ss <- 0:floor(log(max_s, base = 2))
  ss <- 2^(log2ss)
  ss <- c(p, rev(ss))
  as <- ss[]
  as[2:length(ss)] <- sqrt(2 * log(exp(1) * p * log(t) / ss[2:length(ss)]^2))
  as[1] <- 0
  nu_as <- 1 + as * exp(dnorm(as, log = TRUE) -
    pnorm(as, lower.tail = FALSE, log.p = TRUE))



  r_values <- nu_as[]


  r_values[2:length(ss)] <- ss[2:length(ss)] *
    log(1 + sqrt(p * log(t)) / ss[2:length(ss)]) + log(t)

  r_values[1] <- (sqrt(p * log(t)) + log(t))



  return(list(as, nu_as, r_values, ss))
}
