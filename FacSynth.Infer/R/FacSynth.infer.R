#' @title Conduct statistical inference for factor model approach to estimate average treatment effect on treated.
#' @description This function provides users with a formal inference theory for using the factor model approach to estimate ATT, along with a new criterion proposed to more accurately select the number of factors.
#' @details To run this function, we need users to supply data object obtained by running the \code{FacSynth.prep}, or users should manually input relevant matrices and numeric values to fill in the corresponding arguments.
#' @param Facsynth.prep.obj The object that comes from running \code{FacSynth.prep}. This objects contains all matrices and parameters necessary for running the function \code{FacSynth.infer}. Therefore, once users import this object, there is no need to specify any of matrices and parameters manually.
#' @param y1 A matrix that contains values of the treatment unit over the pre-intervention periods.
#' @param x1 A matrix that contains values of the control units (plus an intercept) over the pre-intervention periods.
#' @param y2 A matrix that contains values of the treatment unit over the post-intervention periods.
#' @param x2 A matrix that contains values of the control units (plus an intercept) over the post-intervention periods.
#' @param time.variable A scalar identifying column number or column-name character string associated with period (time) data. The time variable has to be numeric.
#' @param time.prior The number of pre-treatment periods.
#' @param time.post The number of post-treatment periods.
#' @export
#' @return A list which contains ATT derived from factor model approach, a series of confidence intervals, R-square and adjusted R-square.
#' \item{ATT Estimated by FM}{The value of ATT estimated by factor model.}
#' \item{Confidence Interval}{A matrix of intervals across all the confidence level(s) required.}
#' \item{R-squared}{The value of R-squared.}
#' \item{Adjusted R-squared}{The value of adjusted R-squared.}
#' @references
#'   Li K, Sonnier G (2020) Statistical inference for factor model approach to estimate causal effects in quasi-experimental settings. Working Paper, University of Texas at Austin, Austin, TX.
#' @examples
#'   ##In order to conduct inference properly, we recommended users to directly use the list object
#'   ##derived from running FacSynth.prep() and determine the stationarity of the outcome.
#'
#'   #load data from the package.
#'   data(simdata_wide)
#'
#'   #extract relevant components necessary for running FacSynth.infer()
#'   #from wide-format panel data.
#'   FacSynth.prep.out =
#'     FacSynth.prep(
#'      simdata_wide,
#'      treatment.identifier = 2,
#'      controls.identifier = c(7:50),
#'      time.variable = 1,
#'      time.prior = c(1:20),
#'      time.post = c(21:30)
#'     )
#'
#'   #then directly use this object to conduct formal inference process and then
#'   #obtain estimated ATT and its confidence intervals in different levels.
#'   FacsSynth.infer.out =
#'     FacSynth.infer(
#'      FacSynth.prep.out,
#'      STNRY = FALSE,
#'      conf = c(0.99, 0.95, 0.90)
#'     )

FacSynth.infer = function(FacSynth.prep.obj = NULL, STNRY = FALSE,
                          conf = 0.95, y1 = NULL, x1 = NULL, y2 = NULL,
                          x2 = NULL, time.variable = NULL, time.prior = NULL,
                          time.post = NULL) {
  if (is.null(FacSynth.prep.obj) == FALSE) {
    cat("Y1, X1, Y2, X2 and other parameters all come directly from synth.prep object.\n\n")
    y1 = FacSynth.prep.obj[[1]]
    x1 = FacSynth.prep.obj[[2]]
    y2 = FacSynth.prep.obj[[3]]
    x2 = FacSynth.prep.obj[[4]]
  }
  else {
    cat("Y1, X1, Y2, X2 were individually inputted (not FacSynth.prep object.)\n\n")
  }
  store <- list(y1 = y1, x1 = x1, y2 = y2, x2 = x2)
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
  if (ncol(y1) != 1) {
    stop("\n Please specify only one treated unit: Y1 has to have ncol= 1")
  }
  if (ncol(x1) < 3) {
    stop("\n Please specify at least two control units")
  }
  if (ncol(y2) != 1) {
    stop("\n Please specify only one treated unit: Y2 has to have ncol= 1 ")
  }
  if (ncol(x2) < 3) {
    stop("\n Please specify at least two control units")
  }
  if (nrow(y1) != nrow(x1)) {
    stop("\n Different number of periods for treated and controls: nrow(Y1) unequal nrow(X1)")
  }
  if (nrow(y2) != nrow(x2)) {
    stop("\n Different number of periods for treated and controls: nrow(Y2) unequal nrow(X2)")
  }
  if (ncol(x1) != ncol(x2)) {
    stop("\n Different number of units for two-period controls: ncol(X1) unequal ncol(X2)")
  }
  if (nrow(x1) == 0 || nrow(x2) == 0) {
    stop("No periods specified for X1 or X2 (or both). Please specify at least one period for both")
  }
  if (nrow(y1) == 0 || nrow(y2) == 0) {
    stop("No periods specified for Y1 or Y2 (or both). Please specify at least one period for both")
  }

  if (is.null(FacSynth.prep.obj) == FALSE) {
    time.overall = FacSynth.prep.obj[[5]]
    time.prior = FacSynth.prep.obj[[6]]
    time.post = FacSynth.prep.obj[[7]]
    T = FacSynth.prep.obj[[8]]
    T1 = FacSynth.prep.obj[[9]]
    T2 = FacSynth.prep.obj[[10]]
    n = FacSynth.prep.obj[[11]]
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
  x = as.matrix(rbind(x1,x2))
  if (STNRY == FALSE) {
    criti = 11
  }
  else if (STNRY == TRUE) {
    criti = 10 # Treat data as I(1), using Bai's (2004) criterion.
  }
  else {
    stop("\n 'STNRY' be a logical flag.\n")
  }
  rmax = 10 # Maximum number of factors allowed for
  DEMEAN = 1 # 1 means demean the data
  m_N = 80
  m_T = 80
  t1_one = matrix(1,T1,1)
  t2_one = matrix(1,T2,1)
  nbpiid = NBPIID(x = x, kmax = rmax, jj = criti, DEMEAN = DEMEAN,
                  m_N = m_N, m_T = m_T) # Select the number of factors by MBN
  nfactor = nbpiid$ic1
  X = demean(x)
  XX = X%*%t(X)
  XX.svd = svd(XX) # Conduct the single value decomposition
  eigval = diag(XX.svd$d)
  Fhat0 = XX.svd$u
  Fhat1 = XX.svd$v
  F_hat = Fhat0
  X1 = cbind(t1_one, F_hat[c(1:T1),c(1:nfactor)]) # Factors for pre-treatment periods
  X2 = cbind(t2_one, F_hat[c((T1+1):T),c(1:nfactor)]) # Factors for post-treatment periods
  beta_hat = solve(t(X1)%*%X1)%*%(t(X1)%*%y1) # Factor loadings using pre-treatment data
  u1_hat = y1 - X1%*%beta_hat
  y2_factor = X2%*%beta_hat # Estimated counterfactual outcome
  y1hat = X1%*%beta_hat # In-sample-fit

  ATT_FM = mean(y2 - y2_factor) # ATT estimated by factor model
  t_ratio = T2/T1
  sigma_e_hat = mean((y1-y1hat)^2)
  eta_hat = as.matrix(apply(X2, 2, mean))
  psi_hat = t(X1)%*%X1/T1

  Omega_1_hat = sigma_e_hat*t(eta_hat)%*%solve(psi_hat)%*%eta_hat
  v1 = y2 - y2_factor
  Omega_2_hat = mean(u1_hat^2)
  Omega_hat = (T2/T1)*Omega_1_hat + Omega_2_hat

  ATT_std = sqrt(Omega_hat)/sqrt(T2)
  ATT_std0 = sqrt(T2)*ATT_FM/sqrt(Omega_hat)

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
  CI_set = matrix(0,length(conf),2)
  colnames(CI_set) = c("Lower Bound of Confidence Interval","Upper Bound of Confidence Interval")
  cat("\n**************************", "\n************************** \n",
      "\n ATT estimated by FM is", ATT_FM, "\n")
  for (i in 1:length(conf)) {
    CI_lb = ATT_FM - qnorm(1-(1-conf[i])/2)*ATT_std
    CI_ub = ATT_FM + qnorm(1-(1-conf[i])/2)*ATT_std
    CI_set[i,1] = CI_lb
    CI_set[i,2] = CI_ub
    cat("\n", conf[i]*100,"% Confidence Interval of ATT is",
        "[", CI_lb, ",", CI_ub, "] \n")
  }

  r2_FM = 1 - sum(u1_hat^2)/sum((y1 - mean(y1))^2)
  r2_bar_FM = 1 - (1 - r2_FM)*(T1 - 1)/(T1 - (nfactor + 1))
  cat(" R-square is", r2_FM, "\n\n", "Adjusted R-square is", r2_bar_FM, "\n\n")

  FacSynth.infer.output = list("ATT Estimated by FM" = ATT_FM,
                               "Confidence Interval" = CI_set,
                               "R-squared" = r2_FM,
                               "Adjusted R-squared" = r2_bar_FM)
  return(invisible(FacSynth.infer.output))
}

NBPIID = function(x = NULL, kmax = NULL, jj = NULL, DEMEAN = NULL,
                  m_N = NULL, m_T = NULL) {
  T = nrow(x)
  N = ncol(x)
  NT = N*T
  NT1 = N+T
  CT = matrix(0,1,kmax)
  temp = 1
  ii = c()
  while (temp <= kmax) {
    ii = c(ii,temp)
    temp = temp+1
  }
  ii = t(as.matrix(ii))
  if (jj == 1) {
    CT[1,] = log(NT/NT1)*ii*NT1/NT
  }
  else if (jj == 2) {
    CT[1,] = (NT1/NT)*log(min(N,T))*ii
  }
  else if (jj == 3) {
    CT[1,] = ii*log(min(N,T))/min(N,T)
  }
  else if (jj == 4) {
    CT[1,] = 2*ii/T
  }
  else if (jj == 5) {
    CT[1,] = log(T)*ii/T
  }
  else if (jj == 6) {
    CT[1,] = 2*ii*NT1/NT
  }
  else if (jj == 7) {
    CT[1,] = log(NT)*ii*NT1/NT
  }
  else if (jj == 10) {
    CT[1,] = ((N+m_N)*(T+m_T)/NT)*log(NT/NT1)*ii*NT1/NT
  }
  else if (jj == 11) {
    CT[1,] = (N*T/NT)*(T/(4*log(log(T))))*log(NT/NT1)*ii*NT1/NT
  }
  if (DEMEAN == 1) {
    X = demean(x)
  }
  else if (DEMEAN == 0) {
    X = x
  }
  IC1 = matrix(0, nrow(CT), kmax+1)
  Sigma = matrix(0,1,kmax+1)
  XX = X%*%t(X)
  XX.svd = svd(XX)
  eigval = diag(XX.svd$d)
  Fhat0 = XX.svd$u
  Fhat1 = XX.svd$v

  for (i in kmax:1) {
    Fhat = Fhat0[,c(1:i)]
    lambda = t(Fhat)%*%X
    chat = Fhat%*%lambda
    ehat = X - chat
    Sigma[1,i] = mean(apply(ehat*ehat/T,2,sum))
    IC1[,i] = Sigma[1,i] + CT[1,i]*Sigma[1,kmax]
  }
  Sigma[1,kmax+1] = mean(apply(X*X/T,2,sum))
  IC1[ ,kmax+1] = Sigma[1,kmax+1]
  ic1 = t(minindc(t(IC1)))
  ic1 = as.numeric(ic1*(ic1 <= kmax))
  Fhat = as.matrix(Fhat0[ ,c(1:ic1)])
  lamda = t(Fhat)%*%X
  chat = Fhat%*%lamda
  NBPIID.output = list(ic1 = ic1, chat = chat, Fhat = Fhat)
  return(NBPIID.output)
}

minindc = function(x = NULL) {
  ncols = ncol(x)
  nrows = nrow(x)
  pos = matrix(0,ncols,1)
  seq = seqa(1,1,nrows)

  for (i in 1:ncols) {
    dum = min(x[ ,i])
    dum1 = seq*(x[ ,i]-dum == 0)
    pos[i,1] = sum(dum1)
  }
  return(invisible(pos))
}

seqa = function(a = NULL, b = NULL, c = NULL) {
  requireNamespace("pracma", quietly = TRUE)
  x = t(linspace(a,(a+(c-1)*b),c))
  return(invisible(x))
}

demean = function(X = NULL) {
  requireNamespace("pracma", quietly = TRUE)
  T = nrow(X)
  N = ncol(X)
  m = t(as.matrix(apply(X, 2, mean)))
  x = X - repmat(m,T,1)
  return(invisible(x))
}



