#' @title Conduct subsampling method to produce statistical inference for average treatment effects estimated by synthetic method.
#' @description This function provides users with formal inference methods that can be used to derive confidence intervals for the SC or MSC ATT estimators and conduct inference whereas the standard bootstrap cannot.
#' @details To run this function, we need users to first supply data object obtained by running the \code{synth.meth}. In this function, we provide two methods, including subsampling and parametric bootstrap, for statistical inference. Users can use the argument "type" to specify the selected approach.
#' @param synth.meth.obj The object that comes from running \code{synth.meth}. This objects contains all matrices and parameters necessary for running the function \code{synth.infer}.
#' @param m The value of subsample size. Noted that the value of m should not be larger than the number of pre-treatment periods.
#' @param type A character string indicating the type of applied inference method. Possible values are "subsam" (the default, subsampling method) and "par.boot" (parametric bootstrap method).
#' @param conf A scalar or vector containing the confidence level(s) of the required interval(s)
#' @param nb The number of subsampling-bootstrap replication.
#' @export
#' @importFrom pracma fmincon randn
#' @return A list which contains a series of statistics derived from the inference method and other components.
#' \item{ci_set}{A matrix of intervals across all the confidence level(s) required.}
#' \item{cr_mi}{The value of the least lower bound of intervals.}
#' \item{cr_max}{The value of the greatest upper bound of intervals.}
#' \item{ATT_order}{A matrix that contains a set of subsampling ATT's sorted by ascending order.}
#' \item{ATT_SC}{The value of estimated average treatment effect on treated.}
#' \item{nb}{The number of subsampling-bootstrap replication.}
#' @references
#'   Li K (2019) Statistical inference for average treatment effects estimated by synthetic control methods. \emph{Journal of American Statistical Association}.
#' @examples
#'   ##In order to conduct inference properly, we recommended users to directly use the list object
#'   ##derived from running synth.meth() and choose appropriate subsample size m.
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
#'   #then run the synth.meth command to identify the optimized weights estimated by SC method.
#'   synth.meth.out = synth.meth(synth.prep.out)
#'
#'   #directly use this object to conduct formal inference process and then
#'   #obtain confidence intervals in different levels.
#'   synth.infer.out = synth.infer(synth.meth.out, m = 4, conf = c(0.99, 0.95, 0.90))


synth.infer = function(synth.meth.obj = NULL, m = NULL,type = "subsam",
                       conf = 0.95, nb = 1000)
{
  if (is.null(synth.meth.obj) == FALSE) {
    stop("\n No list object supplied in synth.meth.obj.\n")
    b_SC = synth.meth.obj[[1]]
    Y1 = synth.meth.obj[[7]][[1]]
    X1 = synth.meth.obj[[7]][[5]]
    X2 = synth.meth.obj[[7]][[6]]
    T = synth.meth.obj[[7]][[7]]
    T1 = synth.meth.obj[[7]][[8]]
    T2 = synth.meth.obj[[7]][[9]]
    n = synth.meth.obj[[7]][[12]]
    ATT_SC = synth.meth.obj[[2]]
    sigma2_v = synth.meth.obj[[4]]
  }
  if (is.null(m) == TRUE) {
    stop("\n m missing. \n")
  }
  if (mode(m) != "numeric") {
    stop("\n m not found as a numeric value. \n")
  }
  if (length(m) != 1) {
    stop("\n please specify a value for m. \n")
  }
  if (m > T1) {
    stop("\n m not found smaller or equal(bootstrapping) to T1")
  }
  if (mode(nb) != "numeric") {
    stop("\n nb not found as a numeric value. \n")
  }
  if (length(nb) != 1) {
    stop("\n please specify a value for nb. \n")
  }

  A_star_set = matrix(nrow = nb, ncol = 1)
  colnames(A_star_set) = c("A_star")
  set.seed(2020)
  requireNamespace("pracma", quietly = TRUE)
  e1_star = sqrt(sigma2_v)*randn(T2,nb)
  for(i in 1:nb)
  {
    if (type == "subsam") {
      ZZ1 = cbind(X1,Y1)
      Zm = ZZ1[sample(nrow(ZZ1), m, replace = TRUE), ]
      Xm = Zm[, 1:n]
      Ym = Zm[, n+1]
    }
    else if (type == "par.boot") {
      Xm = X1[T1-m+1:T1,]
      Ym = Xm%*%b_SC$par + sigma_SC*randn(m,1)
    }
    else {
      stop("\n No Corresponding inference method is available. \n")
    }
    lb0 = matrix(0, n, 1)
    Aeq0 = matrix(1, 1, n)
    beq0 = 1
    b0 = matrix(1, n, 1)/n
    fn = function(b)
    {
      t(Ym-Xm%*%b)%*%(Ym-Xm%*%b)
    }
    if (labels(synth.meth.obj)[1] == "b_SC") {
      res_m = fmincon(b0, fn, Aeq = Aeq0, beq = beq0, lb = lb0)
    }
    else if (labels(synth.meth.obj)[1] == "b_MSC") {
      lb0[1] = -Inf
      res = fmincon(b0, fn, lb = lb0)
    }
    bm_SC = as.matrix(res_m$par)

    A1_star = - mean( X2%*%(bm_SC - b_SC) )*sqrt((T2*m)/T1)
    A2_star = sqrt(T2)*mean( e1_star[ , i] )
    A_star = A1_star + A2_star
    A_star_set[i, 1] = A_star
  }
  ATT_order = sort(A_star_set[,1], decreasing = FALSE)/sqrt(T2)

  if (length(conf) < 1) {
    stop("\n please specify at least a value for conf. \n")
  }
  if (mode(conf) != "numeric") {
    stop("\n conf not found as numeric. \n")
  }
  for (i in 1:length(conf)) {
    if (conf[i] > 1 || conf[i] < 0) {
      stop("\n Please specify proper confidence level within [0,1].")
    }
  }
  cr_min = ATT_SC - ATT_order[nb]
  cr_max = ATT_SC - ATT_order[1]
  ci_set = matrix(0,length(conf),2)
  colnames(ci_set) = c("Lower Bound of Confidence Interval","Upper Bound of Confidence Interval")
  cat("\n**************************", "\n************************** \n")
  for (i in 1:length(conf)) {
    ci_lb = ATT_SC - ATT_order[round(nb*(1-(1-conf[i]))/2)]
    ci_ub = ATT_SC - ATT_order[round(nb*(1-conf[i])/2)]
    ci_set[i,1] = ci_lb
    ci_set[i,2] = ci_ub
    cat(conf[i]*100,"% Confidence Interval of ATT is",
        "[", ci_lb, ",", ci_ub, "] \n\n")
  }

  synth.infer.out = list(ci_set = ci_set, cr_min = cr_min, cr_max = cr_max,
                    ATT_order = ATT_order, ATT_SC = ATT_SC, nb = nb)
  return(invisible(synth.infer.out))
}

