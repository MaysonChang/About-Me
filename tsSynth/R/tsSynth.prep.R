#' @title Construct a list of matrices from panel dataset to be loaded into other functions needed.
#' @description  This function extracts relevant data objects from the given dataframe and produce a list of matrices necessary for running appropriate synthetic control method.
#' @details Users should import a dataframe ("data") with required wide format, identify the columns associated with treatment and control units respectively, time variable, the pre-treatment and post-treatment time periods.
#' @param data The dataframe with required wide format.
#' @param treatment.identifier A scalar identifying the column number of treatment unit in the dataframe or a character string giving the column name of treatment unit in the dataframe.
#' @param control.identifier A scalar identifying the column numbers of control units in the dataframe or a vector of character strings giving the column names of control units in the dataframe.
#' @param time.variable A scalar identifying column number or column-name character string associated with period (time) data. The time variable has to be numeric.
#' @param time.prior A numeric vector identifying the row numbers corresponding to pre-treatment periods.
#' @param time.post A numeric vector identifying the row numbers corresponding to post-treatment periods.
#' @export
#' @return A list with a series of components prepared for running \code{synth.meth}.
#' \item{Y1}{A matrix of treatment unit data for pre-treatment periods.}
#' \item{X1}{A matrix of control units data for pre-treatment periods.}
#' \item{Y2}{A matrix of treatment unit data for post-treatment periods.}
#' \item{X2}{A matrix of control units data for post-treatment periods.}
#' \item{time}{A matrix of period(time) data.}
#' \item{T}{A scalar identifying the column number of period(time) data.}
#' \item{T1}{A scalar identifying the number of pre-treatment periods.}
#' \item{T2}{A scalar identifying the number of post-treatment periods.}
#' \item{n}{A scalar identifying the column number of control units data.}
#' @references
#'   Li K, Shankar V (2020) Estimating the causal effect of a digitally native retailer opening a new store: a new two-step synthetic control method. Working Paper, University of Texas at Austin, Austin, TX.
#'
#'   Li K, Shankar V (2020) Statistical inference for average treatment effects estimated by synthetic control methods. \emph{Journal of American Statistical Association}.
#' @examples
#'   ##First example: wide-format toy dataset sourced from Package "Synth".
#'
#'   #load data from the package.
#'   data(synth.data)
#'
#'   #extract relevant components necessary for running synth.meth()
#'   #from wide-format panel data.
#'   tsSynth.prep.out =
#'     tsSynth.prep(
#'      synth.data,
#'      treatment.identifier = 9,
#'      controls.identifier = c(2:8),
#'      time.variable = 1,
#'      time.prior = c(1984:1989),
#'      time.post = c(1990:1998)
#'      )
#'
#'   ## Second example: The economic impact of terrorism in the
#'   ## Basque country using data from Abadie and Gardeazabal (2003)
#'
#'   #load data from the package.
#'   data(basque)
#'
#'   #extract relevant components necessary for running synth.meth()
#'   #from wide-format panel data.
#'   tsSynth.prep.out =
#'     tsSynth.prep(
#'      data = basque,
#'      treatment.identifier = 5,
#'      controls.identifier = c(2:4,6:19),
#'      time.variable = "year",
#'      time.predictors.prior = c(1964:1969),
#'      time.post = c(1970:1985)
#'      )
#'


tsSynth.prep = function(data = NULL,treatment.identifier = NULL,
                      controls.identifier = NULL, time.variable = NULL,
                      time.prior = NULL,time.post = NULL)
{
  if (is.data.frame(data) == FALSE) {
    stop("\n No data.frame supplied in data.\n")
  }
  if (mode(treatment.identifier) == "character") {
    treatment.identifier <- which(names(data) == treatment.identifier)
  }
  if (length(treatment.identifier) != 1) {
    stop("\n please specify a single treated unit\n")
  }
  if (is.null(treatment.identifier) == TRUE || mode(data[, treatment.identifier]) !=
      "numeric") {
    stop("\n treatment.identifier not found as numeric variable in data.\n")
  }
  if (length(controls.identifier) < 2) {
    stop("\n please specify at least two control units\n")
  }
  if (mode(controls.identifier) == "character") {
    control.no = c()
    for (i in 1:length(controls.identifier)) {
      control.no = c(control.no, which(names(data) == controls.identifier[i]))
    }
    controls.identifier = control.no
  }
  if (is.null(controls.identifier) == FALSE) {
    if (sum(duplicated(controls.identifier)) > 0) {
      stop("\n duplicate control units in controls.identifier \n")
    }
    for (i in controls.identifier) {
      if (mode(data[, i]) != "numeric") {
        stop(paste("\n control unit", names(data)[i], "not found as numeric variable in data \n"))
      }
    }
  }
  else {
    for (i in controls.identifier) {
      if (is.null(controls.identifier[i]) == TRUE) {
        stop(paste("\n control unit", names(data)[i], "not found as numeric variable in data \n"))
      }
    }
  }
  if (names(data)[treatment.identifier] %in% names(data)[controls.identifier] == TRUE) {
    stop("\n treated unit among controls\n")
  }
  if (mode(time.variable) == "character") {
    time.variable = which(names(data) == time.variable)
  }
  if (is.null(time.variable) == TRUE || mode(data[, time.variable]) !=
      "numeric") {
    stop("\n time.variable not found as numeric variable in data.\n")
  }
  if (length(time.variable) != 1) {
    stop(" Please specify only one time.variable\n")
  }
  if (sum(is.null(time.prior)) > 0) {
    stop("time.prior missing")
  }
  if (mode(time.prior) != "numeric") {
    stop(paste("\n time.prior not found as numeric variable\n"))
  }
  if (sum(duplicated(time.prior)) > 0) {
    stop(paste("\n duplicates in time.prior \n"))
  }
  if (length(time.prior) < 1) {
    stop(paste("\n specificy at least one period in time.prior \n"))
  }
  for (p in time.prior) {
    if (p %in% unique(data[, time.variable]) == FALSE)
      stop(paste("\n time period ", p, " from time.prior
                   not found in time.variable\n",
                 sep = ""))
  }
  if (sum(is.null(time.post)) > 0) {
    stop("time.post missing")
  }
  if (mode(time.post) != "numeric") {
    stop(paste("\n time.post not found as numeric variable\n"))
  }
  if (sum(duplicated(time.post)) > 0) {
    stop(paste("\n duplicates in time.post \n"))
  }
  if (length(time.post) < 1) {
    stop(paste("\n specificy at least one period in time.post \n"))
  }
  for (p in time.post) {
    if (p %in% unique(data[, time.variable]) == FALSE)
      stop(paste("\n time period ", p, " from time.post
                   not found in time.variable\n",
                 sep = ""))
  }

  Y = as.matrix(data[, treatment.identifier])
  X = as.matrix(data[, controls.identifier])
  Y1 = as.matrix(Y[c(which(data[, time.variable] == time.prior[1]):which(data[, time.variable] == time.prior[length(time.prior)])), ])
  if (sum(is.na(Y1)) == nrow(Y1)) {
    stop("\n Treated unit is missing for all pre-treatment time periods \n")
  }
  for (i in 1:nrow(Y1)) {
    if (is.na(Y1[i])) {
      stop(paste("\n Treated unit is missing for period",
                 rownames(Y1)[i], "in pre-treatment interval"))
    }
  }
  X1 = as.matrix(X[c(which(data[, time.variable] == time.prior[1]):which(data[, time.variable] == time.prior[length(time.prior)])), ])
  for (i in 1:ncol(X1)) {
    if (sum(is.na(X1[, i])) == nrow(X1)) {
      stop(paste("\n Control unit", names(X1)[i],
                 "has missing data for all pre-treatment time periods \n"))
    }
    for (j in 1:nrow(X1)) {
      if (is.na(X1[j, i]) == TRUE) {
        cat(paste("\n Missing data- control unit(pre-treatment):",
                  names(X1)[i], "; for period:", rownames(X1)[j], "\n"))
      }
    }
  }
  Y2 = as.matrix(Y[c(which(data[, time.variable] == time.post[1]):which(data[, time.variable] == time.post[length(time.post)])), ])
  if (sum(is.na(Y2)) == nrow(Y2)) {
    stop("\n Treated unit is missing for all post-treatment time periods \n")
  }
  for (i in 1:nrow(Y2)) {
    if (is.na(Y2[i])) {
      stop(paste("\n Treated unit is missing for period",
                 rownames(Y2)[i], "in post-treatment interval"))
    }
  }
  X2 = as.matrix(X[c(which(data[, time.variable] == time.post[1]):which(data[, time.variable] == time.prior[length(time.post)])), ])
  for (i in 1:ncol(X2)) {
    if (sum(is.na(X2[, i])) == nrow(X2)) {
      stop(paste("\n Control unit", names(X1)[i],
                 "has missing data for all post-treatment time periods \n"))
    }
    for (j in 1:nrow(X2)) {
      if (is.na(X2[j, i]) == TRUE) {
        cat(paste("\n Missing data- control unit(post-treatment):",
                  names(X2)[i], "; for period:", rownames(X2)[j], "\n"))
      }
    }
  }
  time = as.matrix(data[, time.variable])
  T1 = as.numeric(length(time.prior))
  T2 = as.numeric(length(time.post))
  T = T1 + T2
  n = as.numeric(ncol(X))

  output <- list(Y1 = Y1, X1 = X1, Y2 = Y2, X2 = X2, time.overall = time, time.prior = time.prior,
                 time.post = time.post, T = T, T1 = T1, T2 = T2, n = n)
  return(invisible(output))
}
