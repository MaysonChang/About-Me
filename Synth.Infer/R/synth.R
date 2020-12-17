#' @title Construct a list of matrices from panel dataset to be loaded into synth.meth().
#' @description  This function extracts relevant data objects from the given dataframe and produce a list of matrices necessary for running synth.meth().
#' @details User should import a dataframe ("data") with required wide format, identify the columns associated with treatment and control units respectively, time variable and the pre-treatment time period.
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
#'   Li K (2019) Statistical inference for average treatment effects estimated by synthetic control methods. \emph{Journal of American Statistical Association}.
#' @examples
#'   ##First example: wide-format toy dataset sourced from Package "Synth".
#'
#'   #load data from the package.
#'   data(synth.data)
#'
#'   #extract relevant components necessary for running synth.meth()
#'   #from wide-format panel data.
#'   synth.prep.out =
#'     synth.prep(
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
#'   synth.prep.out =
#'     synth.prep(
#'      data = basque,
#'      treatment.identifier = 5,
#'      controls.identifier = c(2:4,6:19),
#'      time.variable = "year",
#'      time.prior = c(1964:1969),
#'      time.post = c(1970:1985)
#'      )
#'


synth.prep = function(data = NULL,treatment.identifier = NULL,
                      controls.identifier = NULL, time.variable = NULL,
                      time.prior = NULL, time.post = NULL)
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
  Y1 = as.matrix(Y[c(which(data[, time.variable] == time.predictors.prior[1]):which(data[, time.variable] == time.predictors.prior[length(time.predictors.prior)])), ])
  if (sum(is.na(Y1)) == nrow(Y1)) {
    stop("\n Treated unit is missing for all pre-treatment time periods \n")
  }
  for (i in 1:nrow(Y1)) {
    if (is.na(Y1[i])) {
      stop(paste("\n Treated unit is missing for period",
                 rownames(Y1)[i], "in pre-treatment interval"))
    }
  }
  X1 = as.matrix(X[c(which(data[, time.variable] == time.predictors.prior[1]):which(data[, time.variable] == time.predictors.prior[length(time.predictors.prior)])), ])
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
  X2 = as.matrix(X[c(which(data[, time.variable] == time.post[1]):which(data[, time.variable] == time.post[length(time.post)])), ])
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

#' @title Implement synthetic control method to estimate the average treatment effect.
#' @description This function can routinely search for the set of weights that generate the best fitting combination of the control units chosen to approximate counterfactual outcomes of treatment unit and thus estimate SC ATE with long panels under general conditions. This function require the users to input four matrices as its main arguments. they are named Y1, X1, Y2, and X2 accordingly. Y1 and Y2 contain values of the treatment unit over the pre-intervention and post-intervention periods respectively. X1 and X2 contain the the values for the control units over the pre-intervention and post-intervention periods respectively.
#' @details Creating relevant matrices to be loaded into \code{synth.meth} can be tedious. Users are recommended to apply a preparatory function called "synth.prep" that allows the user to easily create all inputs required for synth.meth(). For more information about the synthetic control method, please refer to the paper Statistical Inference for Average Effects Estimated by Synthetic Control Methods (Kathleen T. Li, 2019).
#' @param synth.prep.object The object that comes from running synth.prep(). This objects contains all information about four necessary matrices. Therefore, once users import this object, there is no need to specify any of four matrices manually.
#' @param meth A character string indicating the type of applied method. Possible values are "ordinary" (the default, SC method) and "modified" (MSC method).
#' @param Y1 A matrix of treatment unit data for pre-treatment periods.
#' @param X1 A matrix of control units data for pre-treatment periods.
#' @param Y2 A matrix of treatment unit data for post-treatment periods.
#' @param X2 A matrix of control units data for post-treatment periods.
#' @param time.overall A numeric vector identifying the row numbers corresponding to overall periods
#' @export
#' @importFrom pracma fmincon
#' @return A list that contains the weights on control units and other components necessary for running \code{synth.infer}.
#' \item{b_SC}{A vector of weights across the controls.}
#' \item{ATT_S}{The value of estimated average treatment effect.}
#' \item{sigma_SC}{MSPE from optimization over weights on control units.}
#' \item{sigma2_v}{The variance of estimated average treatment effect.}
#' \item{r2_SC}{The value of R-squared.}
#' \item{r2_bar_SC}{The value of adjusted R-squared.}
#' \item{dataset}{A list that contains a series of matrices and scalars for future use.}
#' @references
#'   Li K (2019) Statistical inference for average treatment effects estimated by synthetic control methods. \emph{Journal of American Statistical Association}.
#' @examples
#'   ##In order to run synth.meth() properly, we recommended users to run
#'   ##synth.prep for component-extraction first.
#'
#'   #load data from the package.
#'   data(synth.data)
#'
#'   #extract relevant components necessary for running synth.meth()
#'   #from wide-format panel data.
#'   synth.prep.out =
#'     synth.prep(
#'      data = synth.data,
#'      treatment.identifier = 9,
#'      controls.identifier = c(2:8),
#'      time.variable = 1,
#'      time.prior = c(1984:1989),
#'      time.post = c(1990:1998)
#'      )
#'   #now run the synth.meth command to identify the optimized weights
#'   synth.meth.out = synth.meth(synth.prep.out)
#'
#'   #or users can choose MSC as estimation method.
#'   synth.meth.out2 = synth.meth(synth.prep.out, meth = "modified")

synth.meth = function(synth.prep.obj = NULL, meth = "ordinary",
                      Y1 = NULL, X1 = NULL, Y2 = NULL, X2 = NULL,
                      time.variable = NULL, time.prior = NULL,
                      time.post = NULL)
{
  if (is.null(synth.prep.obj) == FALSE) {
    cat("Y1, X1, Y2, X2 and other parameters all come directly from synth.prep object.\n\n")
    Y1 = synth.prep.obj[[1]]
    X1 = synth.prep.obj[[2]]
    Y2 = synth.prep.obj[[3]]
    X2 = synth.prep.obj[[4]]
  }
  else {
    cat("Y1, X1, Y2, X2 were individually inputted (not synth.prep object.)\n\n")
  }
  store = list(Y1 = Y1, X1 = X1, Y2 = Y2, X2 = X2)
  for (i in 1:4) {
    if (is.null(store[[i]])) {
      stop(paste("\n", names(store)[i], "is missing \n"))
    }
    if (sum(is.na(store[[i]])) > 0) {
      stop(paste("\n NAs in", names(store)[i], "\n"))
    }
    if (is.matrix(store[[i]]) == FALSE) {
      stop(paste("\n", names(store)[i], "is not a matrix object\n"))
    }
  }
  if (ncol(Y1) != 1) {
    stop("\n Please specify only one treated unit: Y1 has to have ncol= 1")
  }
  if (ncol(X1) < 2) {
    stop("\n Please specify at least two control units: X1 has to have ncol >= 2")
  }
  if (ncol(Y2) != 1) {
    stop("\n Please specify only one treated unit: Y2 has to have ncol= 1 ")
  }
  if (ncol(X2) < 2) {
    stop("\n Please specify at least two control units: X2 has to have ncol >= 2")
  }
  if (nrow(Y1) != nrow(X1)) {
    stop("\n Different number of periods for treated and controls: nrow(Y1) unequal nrow(X1)")
  }
  if (nrow(Y2) != nrow(X2)) {
    stop("\n Different number of periods for treated and controls: nrow(Y2) unequal nrow(X2)")
  }
  if (ncol(X1) != ncol(X2)) {
    stop("\n Different number of units for two-period controls: ncol(X1) unequal ncol(X2)")
  }
  if (nrow(X1) == 0 || nrow(X2) == 0) {
    stop("No periods specified for X1 or X2 (or both). Please specify at least one period for both")
  }
  if (nrow(Y1) == 0 || nrow(Y2) == 0) {
    stop("No periods specified for Y1 or Y2 (or both). Please specify at least one period for both")
  }

  if (is.null(synth.prep.obj) == FALSE) {
    time.overall = synth.prep.obj[[5]]
    time.prior = synth.prep.obj[[6]]
    time.post = synth.prep.obj[[7]]
    T = synth.prep.obj[[8]]
    T1 = synth.prep.obj[[9]]
    T2 = synth.prep.obj[[10]]
    n = synth.prep.obj[[11]]
  }
  else {
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
      if (p %in% unique(data[ ,time.variable]) == FALSE)
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
      if (p %in% unique(data[ ,time.variable]) == FALSE)
        stop(paste("\n time period ", p, " from time.post
                   not found in time.variable\n",
                   sep = ""))
    }
  }
  time = as.matrix(data[, time.variable])
  T1 = as.numeric(length(time.prior))
  T2 = as.numeric(length(time.post))
  T = T1 + T2
  n = as.numeric(ncol(X))

  requireNamespace("pracma", quietly = TRUE)
  if (meth == "ordinary") {
    lb0 = matrix(0, n, 1)
    Aeq0 = matrix(1, 1, n)
    beq0 = 1
    b0 = matrix(1, n, 1)/n
    fn = function(b)
    {
      t(Y1-X1%*%b)%*%(Y1-X1%*%b)
    }
    res = fmincon(b0, fn, Aeq = Aeq0, beq = beq0, lb = lb0)
    b_SC = as.matrix(res$par)
    Y2_SC = X2%*%b_SC
    colnames(Y2_SC) = "Counterfactual Outcome"
    Y1_SC = X1%*%b_SC
    colnames(Y1_SC) = "Weighted Average of Controls"

    ATT_SC = mean(Y2-Y2_SC)
    sigma_SC = sqrt( mean( (Y1-Y1_SC)^2 ) )
    sigma2_v =  mean((Y2-Y2_SC - mean(Y2-Y2_SC))^2 )

    u1_SC = Y1 - Y1_SC
    r2_SC = 1 - sum(u1_SC^2)/sum((Y1 - mean(Y1))^2)
    r2_bar_SC = 1 - (1 - r2_SC)*T1/(T1-n)

    cat("\n****************", "\n****************",
        "\n\n The SC estimated ATT:", ATT_SC, "\n\n Optimized unit weights:\n",
        round(as.numeric(b_SC), 10), "\n\n MSE in pretreatment time period:",
        sigma_SC, "\n\n Variance of the SC estimated ATT:", sigma2_v, "\n\n The R-squared:", r2_SC, "\n\n The adjusted R-squared:",
        r2_bar_SC, "\n\n")

    dataset = list(Y1 = Y1, Y2 = Y2, Y1_SC = Y1_SC, Y2_SC = Y2_SC, X1 = X1, X2 = X2, T = T,
                   T1 = T1, T2 = T2, time.overall = time.overall, time.prior = time.prior,
                   n = n)
    optimize.out = list(b_SC = b_SC, ATT_SC = ATT_SC, sigma_SC = sigma_SC,
                        sigma2_v = sigma2_v, r2_SC = r2_SC, r2_bar_SC = r2_bar_SC, dataset = dataset)
  }
  else if (meth == "modified") {
    CCC = matrix(1,T,1)
    X1 = as.matrix(cbind(CCC[c(1:nrow(X1)), ],X1))
    X2 = as.matrix(cbind(CCC[c(1:nrow(X2)), ],X2))
    n = ncol(X1)

    lb0 = matrix(0, n, 1)
    lb0[1] = -Inf
    b0 = matrix(1, n, 1)/n
    fn = function(b)
    {
      t(Y1-X1%*%b)%*%(Y1-X1%*%b)
    }
    res = fmincon(b0, fn, lb = lb0)
    b_MSC = as.matrix(res$par)
    Y2_MSC = X2%*%b_MSC
    colnames(Y2_MSC) = "Counterfactual Outcome"
    Y1_MSC = X1%*%b_MSC
    colnames(Y1_MSC) = "Weighted Average of Controls"

    ATT_MSC = mean(Y2-Y2_MSC)
    sigma_MSC = sqrt( mean( (Y1-Y1_MSC)^2 ) )
    sigma2_v =  mean((Y2-Y2_MSC - mean(Y2-Y2_MSC))^2 )

    u1_MSC = Y1 - Y1_MSC
    r2_MSC = 1 - sum(u1_MSC^2)/sum((Y1 - mean(Y1))^2)
    r2_bar_MSC = 1 - (1 - r2_MSC)*(T1-1)/(T1-n)

    cat("\n****************", "\n****************",
        "\n\n The MSC estimated ATT:", ATT_MSC, "\n\n Optimized unit weights:\n",
        round(as.numeric(b_MSC), 10), "\n\n MSE in pretreatment time period:",
        sigma_MSC, "\n\n Variance of the MSC estimated ATT:", sigma2_v, "\n\n The R-squared:", r2_MSC, "\n\n The adjusted R-squared:",
        r2_bar_MSC, "\n\n")

    dataset = list(Y1 = Y1, Y2 = Y2, Y1_MSC = Y1_MSC, Y2_MSC = Y2_MSC, X1 = X1, X2 = X2, T = T,
                   T1 = T1, T2 = T2, time.overall = time.overall, time.prior = time.prior, n = n)
    optimize.out = list(b_MSC = b_MSC, ATT_MSC = ATT_MSC, sigma_MSC = sigma_MSC,
                        sigma2_v = sigma2_v, r2_MSC = r2_MSC, r2_bar_MSC = r2_bar_MSC, dataset = dataset)
  }
  else {
    stop("\n No Corresponding estimation method available. \n")
  }
  return(invisible(optimize.out))
}

