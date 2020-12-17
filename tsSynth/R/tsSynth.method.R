#' @title Implement appropriate synthetic control method to estimate the average treatment effect.
#' @description This function can routinely search for the set of weights that generate the best fitting combination of the control units chosen to approximate counterfactual outcomes of treatment unit and thus estimate ATT with long panels under general conditions. This function require the users to input four matrices as its main arguments. they are named Y1, X1, Y2, and X2 accordingly. Y1 and Y2 contain values of the treatment unit over the pre-intervention and post-intervention periods respectively. X1 and X2 contain the the values for the control units over the pre-intervention and post-intervention periods respectively.
#' @details Creating relevant matrices to be loaded into \code{tsSynth.meth} can be tedious. Users are recommended to apply a preparatory function called \code{tsSynth.prep} that allows the user to easily create all inputs required for \code{tsSynth.meth}. After extracting the required matrices, it is advisable to run \code{tsSynth.test} to determine which model is the most appropriate for the corresponding treatment and controls.
#' @param synth.prep.object The object that comes from running \code{tsSynth.prep}. This objects contains all information about four necessary matrices. Therefore, once users import this object, there is no need to specify any of four matrices manually.
#' @param meth A character string indicating the type of applied method. Possible values are "SC" (the default, original SC method), "MSCa" (MSC method with weights sum to one), "MSCb" (MSC method with zero intercept) and "MSC" (MSC method with non-negativity).
#' @param Y1 A matrix of treatment unit data for pre-treatment periods.
#' @param X1 A matrix of control units data for pre-treatment periods.
#' @param Y2 A matrix of treatment unit data for post-treatment periods.
#' @param X2 A matrix of control units data for post-treatment periods.
#' @param time.variable A scalar identifying column number or column-name character string associated with period (time) data.
#' @param time.prior  A numeric vector identifying the row numbers corresponding to pre-treatment periods.
#' @param time.post A numeric vector identifying the row numbers corresponding to post-treatment periods.
#' @export
#' @importFrom pracma fmincon
#' @return A list that contains the weights on control units and other components.
#' \item{b_SC}{A vector of weights across the controls.}
#' \item{ATT_SC}{The value of estimated average treatment effect on treated unit.}
#' \item{sigma_SC}{MSPE from optimization over weights on control units.}
#' \item{sigma2_v}{The variance of estimated average treatment effect.}
#' \item{r2_SC}{The value of R-squared.}
#' \item{r2_bar_SC}{The value of adjusted R-squared.}
#' \item{dataset}{A list that contains a series of matrices and scalars for future use.}
#' @references
#'   Li K, Shankar V (2020) Estimating the causal effect of a digitally native retailer opening a new store: a new two-step synthetic control method. Working Paper, University of Texas at Austin, Austin, TX.
#'
#'   Li K (2019) Statistical inference for average treatment effects estimated by synthetic control methods. \emph{Journal of American Statistical Association}.
#' @examples
#'   ##In order to run tsSynth.meth() properly, we recommend users to run
#'   ##tsSynth.prep for components extraction first.
#'
#'   #load data from the package.
#'   data(synth.data)
#'
#'   #extract relevant components necessary for running tsSynth.meth()
#'   #from wide-format panel data.
#'   tsSynth.prep.out =
#'     tsSynth.prep(
#'      data = synth.data,
#'      treatment.identifier = 9,
#'      controls.identifier = c(2:8),
#'      time.variable = 1,
#'      time.prior = c(1984:1989),
#'      time.post = c(1990:1998)
#'      )
#'
#'   #then run tsSynth.test() to determine the appropriate method based
#'   #on the testing result.
#'   tsSynth.test.out =
#'     tsSynth.test(tsSynth.prep.out, m = 4)
#'
#'   #the testing result indicates that original SC model is the most
#'   #appropriate.
#'   #now run tsSynth.meth function to apply original SC model
#'   #to identify the optimized weights
#'   tsSynth.meth.out =
#'     tsSynth.meth(
#'      tsSynth.prep.out,
#'      meth = "SC"
#'      )
#'
#'   #or users can choose MSCa model (with weights sum to one).
#'   tsSynth.meth.out2 =
#'     tsSynth.meth(
#'      tsSynth.prep.out,
#'      meth = "MSCa"
#'      )
#'
#'   #or apply MSCb model (with zero intercept).
#'   tsSynth.meth.out3 =
#'     tsSynth.meth(
#'      tsSynth.prep.out,
#'      meth = "MSCb"
#'      )
#'
#'   #or apply MSC model (with non-negativity).
#'   tsSynth.meth.out4 =
#'     tsSynth.meth(
#'      tsSynth.prep.out,
#'      meth = "MSC",
#'      )

tsSynth.meth = function(tsSynth.prep.obj = NULL, meth = "ordinary",
                        Y1 = NULL, X1 = NULL, Y2 = NULL,X2 = NULL,
                        time.variable = NULL, time.prior = NULL,
                        time.post = NULL)
{
  if (is.null(tsSynth.prep.obj) == FALSE) {
    cat("Y1, X1, Y2, X2 and other parameters all come directly from tsSynth.prep object.\n\n")
    Y1 = tsSynth.prep.obj[[1]]
    X1 = tsSynth.prep.obj[[2]]
    Y2 = tsSynth.prep.obj[[3]]
    X2 = tsSynth.prep.obj[[4]]
  }
  else {
    cat("Y1, X1, Y2, X2 were individually inputted (not tsSynth.prep object.)\n\n")
  }
  store <- list(Y1 = Y1, X1 = X1, Y2 = Y2, X2 = X2)
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

  if (is.null(tsSynth.prep.obj) == FALSE) {
    time.overall = tsSynth.prep.obj[[5]]
    time.prior = tsSynth.prep.obj[[6]]
    time.post = tsSynth.prep.obj[[7]]
    T = tsSynth.prep.obj[[8]]
    T1 = tsSynth.prep.obj[[9]]
    T2 = tsSynth.prep.obj[[10]]
    n = tsSynth.prep.obj[[11]]
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
    T1 = as.numeric(length(time.prior))
    T2 = as.numeric(length(time.post))
    T = T1 + T2
    n = as.numeric(ncol(X1))
  }

  requireNamespace("pracma", quietly = TRUE)
  if (meth == "SC") {
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
    sigma2_v =  mean((Y2-Y2_SC - mean(Y2-Y2_SC))^2 )#var of ATT_SC

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
                   time.post = time.post, n = n)
    optimize.out = list(b_SC = b_SC, ATT_SC = ATT_SC, sigma_SC = sigma_SC,
                        sigma2_v = sigma2_v, r2_SC = r2_SC, r2_bar_SC = r2_bar_SC, dataset = dataset)
    return(invisible(optimize.out))
    }
  else {
    if (meth == "MSCa") {
      CCC = matrix(1,T,1)
      X1 = as.matrix(cbind(CCC[c(1:nrow(X1)), ],X1))
      X2 = as.matrix(cbind(CCC[c(1:nrow(X2)), ],X2))
      n = ncol(X1)
      lb0 = matrix(0, n, 1)
      lb0[1] = -Inf
      Aeq0 = matrix(1, 1, n)
      Aeq0[1] = 0
      beq0 = 1
      b0 = matrix(1, n, 1)/n
      fn = function(b)
      {
        t(Y1-X1%*%b)%*%(Y1-X1%*%b)
      }
      res = fmincon(b0, fn, Aeq = Aeq0, beq = beq0, lb = lb0)

    }
    else if (meth == "MSCb") {
      lb0 = matrix(0, n, 1)
      b0 = matrix(1, n, 1)/n
      fn = function(b)
      {
        t(Y1-X1%*%b)%*%(Y1-X1%*%b)
      }
      res = fmincon(b0, fn, lb = lb0)
    }
    else if (meth == "MSC") {
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
    }
    else {
      stop("\n No Corresponding estimation method available. \n")
    }

    b_MSC = as.matrix(res$par)
    Y2_MSC = X2%*%b_MSC
    colnames(Y2_MSC) = "Counterfactual Outcome"
    Y1_MSC = X1%*%b_MSC
    colnames(Y1_MSC) = "Weighted Average of Controls"

    ATT_MSC = mean(Y2-Y2_MSC)
    sigma_MSC = sqrt( mean( (Y1-Y1_MSC)^2 ) )
    sigma2_v =  mean((Y2-Y2_MSC - mean(Y2-Y2_MSC))^2 )#var of ATT_MSC

    u1_MSC = Y1 - Y1_MSC
    r2_MSC = 1 - sum(u1_MSC^2)/sum((Y1 - mean(Y1))^2)
    r2_bar_MSC = 1 - (1 - r2_MSC)*(T1-1)/(T1-n)

    cat("\n****************", "\n****************",
        "\n\n The MSC estimated ATT:", ATT_MSC, "\n\n Optimized unit weights:\n",
        round(as.numeric(b_MSC), 10), "\n\n MSE in pretreatment time period:",
        sigma_MSC, "\n\n Variance of the MSC estimated ATT:", sigma2_v, "\n\n The R-squared:", r2_MSC, "\n\n The adjusted R-squared:",
        r2_bar_MSC, "\n\n")

    dataset = list(Y1 = Y1, Y2 = Y2, Y1_MSC = Y1_MSC, Y2_MSC = Y2_MSC, X1 = X1, X2 = X2, T = T,
                   T1 = T1, T2 = T2, time.overall = time.overall, time.prior = time.prior,
                   time.post = time.post, n = n)
    optimize.out = list(b_MSC = b_MSC, ATT_MSC = ATT_MSC, sigma_MSC = sigma_MSC,
                        sigma2_v = sigma2_v, r2_MSC = r2_MSC, r2_bar_MSC = r2_bar_MSC, dataset = dataset)
    return(invisible(optimize.out))
  }
}

