#' @title Apply a two-step hypothesis test to determine the appropriate method based on the testing result.
#' @description  This function provides a formal testing procedure for the parallel trends assumption in the first step, and the application of an appropriate synthetic control method used for estimating accurate ATT of treated unit in the second step.
#' @details Users are recommended to supply data object obtained by running the \code{synth.prep}, or import data of both treated and control units over the pretreatment periods manually. Moreover, a subsampling method serves to yield consistent estimators of the test statistic, and thus users need to choose appropriate subsample size and times of subsampling repetition.
#' @param synth.prep.obj The object that comes from running \code{synth.prep}. This objects contains matrices necessary for running this function. Therefore, once users import this object, there is no need to specify Y1, X1 manually.
#' @param m The value of subsample size. Noted that the value of m should not be larger than the number of pre-treatment periods.
#' @param nb The number of subsampling-bootstrap replication.
#' @param Y1 A matrix that contains values of the treatment unit over the pre-intervention periods.
#' @param X1 A matrix that contains values of the control units over the pre-intervention periods.
#' @export
#' @importFrom pracma fmincon
#' @return The decision on which method is the most appropriate to apply.
#' @references
#'   Li K, Shankar V (2020) Estimating the causal effect of a digitally native retailer opening a new store: a new two-step synthetic control method. Working Paper, University of Texas at Austin, Austin, TX.
#'
#'   Li K, Shankar V (2020) Statistical inference for average treatment effects estimated by synthetic control methods. \emph{Journal of American Statistical Association}.
#' @examples
#'   ##In order to run synth.test() properly, we recommended users
#'   ##to run synth.prep for component-extraction first.
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
#'      time.predictors.prior = c(1984:1989)
#'      )
#'
#'   ##apply the two-step hypothesis test to determine the most appropriate
#'   ##method based on the results of the test.
#'   synth.test.out =
#'     synth.test(synth.prep.out,m = 4)
#'


tsSynth.test = function(tsSynth.prep.obj = NULL, m = NULL, nb = 1000,
                      Y1 = NULL, X1 = NULL) {
  if (is.null(tsSynth.prep.obj) == FALSE) {
    cat("Y1 and X1 come directly from tsSynth.prep object.\n\n")
    Y1 = tsSynth.prep.obj[[1]]
    X1 = tsSynth.prep.obj[[2]]
    T1 = tsSynth.prep.obj[[9]]
    n = tsSynth.prep.obj[[11]]
  }
  else {
    cat("Y1 and X1 were individually inputted (not tsSynth.prep object.)\n\n")
    store <- list(Y1 = Y1, X1 = X1)
    for (i in 1:2) {
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
    if (nrow(Y1) != nrow(X1)) {
      stop("\n Different number of periods for treated and controls: nrow(Y1) unequal nrow(X1)")
    }
    T1 = as.numeric(nrow(Y1))
    n = as.numeric(ncol(X1))
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

  requireNamespace("pracma", quietly = TRUE)
  n = n+1
  CCC = matrix(1, nrow = T1, ncol = 1)
  X1 = cbind(CCC,X1)

  Rt = cbind(
    matrix(c(0,1), nrow = 2, ncol = 1),
    rbind(matrix(1,nrow = 1,ncol = n-1),matrix(0,nrow = 1,ncol = n-1))
  )
  qt = matrix(c(1,0),nrow = 2,ncol = 1)

  lb0 = matrix(0, n, 1)
  lb0[1] = -Inf
  b0 = matrix(1, n, 1)/n
  fn = function(b)
  {
    t(Y1-X1%*%b)%*%(Y1-X1%*%b)
  }
  res = fmincon(b0, fn, lb = lb0)
  b_t1 = as.matrix(res$par)

  dt_set = list() #set up an empty list to store matrices
  var_hat = 0

  for(i in 1:nb)
  {
    ZZ1 = cbind(X1,Y1)
    Zm = ZZ1[sample(nrow(ZZ1), m, replace = TRUE), ]
    Xm = Zm[, 1:n]
    Ym = Zm[, n+1]

    lb0 = matrix(0, n, 1)
    lb0[1] = -Inf
    b0 = matrix(1, n, 1)/n
    fn = function(b)
    {
      t(Ym-Xm%*%b)%*%(Ym-Xm%*%b)
    }
    res = fmincon(b0, fn, lb = lb0)
    b_m = as.matrix(res$par)

    dgap = b_m-b_t1
    var_hat = var_hat+dgap%*%t(dgap)

    dt_s = sqrt(m)*(Rt%*%(b_m-b_t1))
    dt_set[[i]] = dt_s

  }
  var_hat = m/nb*var_hat
  v = solve(Rt%*%var_hat%*%t(Rt))
  s_set = matrix(nrow = nb, ncol = 1)
  for(i in 1:nb) {
    s_m = t(dt_set[[i]])%*%v%*%(dt_set[[i]])
    s_set[i,1] = s_m
  }
  s_order = sort(s_set[,1], decreasing = FALSE)
  s_t1 = as.numeric(T1*(t(Rt%*%b_t1-qt)%*%v%*%(Rt%*%b_t1-qt)))
  pJ = mean(s_order>s_t1)

  if (pJ >0.05) {
    cat("\n****************", "\n****************",
        "\n\n Fail to reject the joint hypothesis, p-value:",
        pJ,"\n\n Recommend to apply the original SC method.\n\n")
  }
  else {
    Ra = Rt[1,]
    Rb = Rt[2,]
    qa = 1
    qb = 0
    da = Ra%*%b_t1-qa
    db = Rb%*%b_t1-qb
    s_t1_a = as.numeric(T1*(da^2))
    s_t1_b = as.numeric(T1*(db^2))
    s_set_a = matrix(nrow = nb, ncol = 1)
    s_set_b = matrix(nrow = nb, ncol = 1)

    for(i in 1:nb)
    {
      ZZ1 = cbind(X1,Y1)
      Zm = ZZ1[sample(nrow(ZZ1), m, replace = TRUE), ]
      Xm = Zm[, 1:n]
      Ym = Zm[, n+1]

      lb0 = matrix(0, n, 1)
      lb0[1] = -Inf
      b0 = matrix(1, n, 1)/n
      fn = function(b)
      {
        t(Ym-Xm%*%b)%*%(Ym-Xm%*%b)
      }
      res = fmincon(b0, fn, lb = lb0)
      b_m = as.matrix(res$par)

      s_m_a = m*(Ra%*%(b_m-b_t1))^2
      s_set_a[i,1] = s_m_a
      s_m_b = m*(Rb%*%(b_m-b_t1))^2
      s_set_b[i,1] = s_m_b
    }
    s_order_a = sort(s_set_a[,1], decreasing = FALSE)
    s_order_b = sort(s_set_b[,1], decreasing = FALSE)
    pa = mean(s_order_a>s_t1_a)
    pb = mean(s_order_b>s_t1_b)

    if (pa < 0.05 & pb > 0.05) {
      cat("\n****************", "\n****************",
          "\n\n Reject 'weights sum to one' hypothesis while fail to reject",
          "\n 'zero intercept' hypothesis.\n",
          "\n p-value for 'weights sum to one' hypothesis:",pa,
          "\n p-value for 'without intercept' hypothesis:",pb,
          "\n\n Recommend to apply the MSC method with zero intercept.\n\n")
    }
    else if (pa > 0.05 & pb < 0.05) {
      cat("\n****************", "\n****************",
          "\n\n Fail to reject 'weights sum to one' hypothesis while reject",
          "\n 'zero intercept' hypothesis.\n",
          "\n p-value for 'weights sum to one' hypothesis:",pa,
          "\n p-value for 'without intercept' hypothesis:",pb,
          "\n\n Recommend to apply the MSC method with weights sum to one.\n\n")
    }
    else if (pa < 0.05 & pb < 0.05) {
      cat("\n****************", "\n****************",
          "\n\n Reject both 'weights sum to one' hypothesis and 'zero intercept'",
          "\n hypothesis.\n",
          "\n p-value for 'weights sum to one' hypothesis:",pa,
          "\n p-value for 'without intercept' hypothesis:",pb,
          "\n\n Recommend to apply the MSC method.\n\n")
    }
  }
}
