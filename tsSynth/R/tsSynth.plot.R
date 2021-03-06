#' @title Plot graphs relevant to the outcome trajectories for treatment and constructed control.
#' @description  This function allows users to plot the outcome trajectories between the treated units and the synthetic control unit, and gaps in outcome between these two groups. The users can specify some attributes of the plots, like labels for axis, main title and legend text.
#' @details User should import a dataframe ("data") with required wide format, identify the columns associated with treatment and control units respectively, time variable and the pre-treatment time periods.
#' @param tsSynth.meth.obj Output list created by \code{tsSynth.meth}.
#' @param type A character string indicating the type of plot. Possible values are "path" (the default, outcome trajectories) and "gaps" (gaps in outcome trajectories).
#' @param tr.intake Optional scalar to indicate the time period when the intervention occurs.
#' @param Ylab Optional label for Y axis.
#' @param Xlab Optional label for X axis.
#' @param Legend Optional legend text. Only applicable when plotting outcome trajectories.
#' @param Main Optional main title.
#' @export
#' @import ggplot2
#' @return The plot of two groups' trajectories or their gaps.
#' @references
#'   Li K, Shankar V (2020) Estimating the causal effect of a digitally native retailer opening a new store: a new two-step synthetic control method. Working Paper, University of Texas at Austin, Austin, TX.
#'
#'   Li K (2019) Statistical inference for average treatment effects estimated by synthetic control methods. \emph{Journal of American Statistical Association}.
#' @examples
#'   ##In order to run tsSynth.plot() properly, we recommend users to run
#'   ##tsSynth.prep for components extraction first.
#'
#'   #load data from the package.
#'   data(synth.data)
#'
#'   #extract relevant components necessary for running tsSynth.meth()
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
#'   #then run the tsSynth.meth command to identify the optimized weights
#'   tsSynth.meth.out =
#'     tsSynth.meth(
#'      tsSynth.prep.out,
#'      meth = "SC"
#'      )
#'
#'   #now plot the outcome trajectories for treatment and constructed control.
#'   tsSynth.path.plot =
#'     tsSynth.plot(
#'      tsSynth.meth.out,
#'      type = "path",
#'      tr.intake = 1990,
#'      Ylab = "Production",
#'      Xlab = "Year",
#'      Legend = c("Treatment","Control")
#'      )
#'
#'    #or plot the gaps in outcome trajectories for two groups.
#'    tsSynth.gaps.plot =
#'     tsSynth.plot(
#'      tsSynth.meth.out,
#'      type = "gaps",
#'      tr.intake = 1990,
#'      Ylab = "Production Gaps",
#'      Xlab = "Year",
#'      Main = "Gaps in Outcome between Treatment and Constructed Control"
#'      )
#'


tsSynth.plot = function(tsSynth.meth.obj = NULL, type = "path", tr.intake = NULL, Ylab = "Outcome",
                      Xlab = "Time", Legend = c("Treatment","Control"),
                      Main = "Trajectories for Treatment and Constructed Control")
{
  if (is.list(tsSynth.meth.obj) == FALSE) {
    stop("\n Output list created by tsSynth.meth() is required.\n")
  }
  if (is.character(Legend) == FALSE) {
    stop(stop("\n Legend not found as character vector.\n"))
  }
  if (is.character(Xlab) == FALSE) {
    stop(stop("\n Xlab not found as character scalar.\n"))
  }
  if (is.character(Ylab) == FALSE) {
    stop(stop("\n Ylab not found as character scalar.\n"))
  }
  if (is.character(Main) == FALSE) {
    stop(stop("\n Main not found as character scalar.\n"))
  }
  Y1 = tsSynth.meth.obj[[7]][[1]]
  Y2 = tsSynth.meth.obj[[7]][[2]]
  Y1_SC = tsSynth.meth.obj[[7]][[3]]
  Y2_SC = tsSynth.meth.obj[[7]][[4]]
  T1 = tsSynth.meth.obj[[7]][[8]]
  T2 = tsSynth.meth.obj[[7]][[9]]
  time.overall = tsSynth.meth.obj[[7]][[10]]
  time.prior = tsSynth.meth.obj[[7]][[11]]
  time.post = tsSynth.meth.obj[[7]][[12]]

  requireNamespace("ggplot2", quietly = TRUE)
  Y_SC = rbind(Y1_SC,Y2_SC)
  Y = rbind(Y1,Y2)
  if (type == "path") {
    outcome = cbind(Y,Y_SC)
    if (length(Legend) != 2) {
      stop("\n Length of legend text limited to 2.\n")
    }
    colnames(outcome) = c("Treatment","Control")
    date = time.overall[c(which(time.overall == time.prior[1]):which(time.overall == time.post[length(time.post)])),]
    plot.obj = data.frame(cbind(date, outcome))
    if (is.null(tr.intake) == FALSE) {
      plot_SC = ggplot(plot.obj, aes(x=date)) +
        geom_line(aes(y = Treatment, color="firebrick"), size=1) +
        geom_line(aes(y = Control, color="steelblue"), size=1) +
        geom_vline(xintercept = tr.intake, size = 1, color = "tan") +
        labs(title = Main, y = Ylab, x = Xlab) +
        scale_color_identity(
          name = "Groups",
          breaks = c("firebrick", "steelblue"),
          labels = c(Legend[1], Legend[2]),
          guide = "legend")
    }
    else {
      plot_SC = ggplot(plot.obj, aes(x=date)) +
        geom_line(aes(y = Treatment, color="firebrick"), size=1) +
        geom_line(aes(y = Control, color="steelblue"), size=1) +
        geom_vline(xintercept = date[which(date == date[length(time.prior)+1])], size = 1, color = "tan") +
        labs(title = Main, y = Ylab, x = Xlab) +
        scale_color_identity(
          name = "Groups",
          breaks = c("firebrick", "steelblue"),
          labels = c(Legend[1], Legend[2]),
          guide = "legend")
    }
  }
  else if (type == "gaps") {
    outcome = Y-Y_SC
    date = time.overall[c(which(time.overall == time.prior[1]):which(time.overall == time.post[length(time.post)])),]
    plot.obj = data.frame(cbind(date,outcome))
    if (is.null(tr.intake) == FALSE) {
      plot_SC = ggplot(plot.obj, aes(x=date)) +
        geom_line(aes(y=outcome), size=1, color="steelblue") +
        geom_vline(xintercept = tr.intake, size = 1, color = "tan") +
        geom_hline(yintercept=0, size=1, linetype="dotted", color="grey0") +
        labs(title = Main, y = Ylab, x = Xlab)
    }
    else {
      plot_SC = ggplot(plot.obj, aes(x=date)) +
        geom_line(aes(y=outcome), size=1, color="steelblue") +
        geom_vline(xintercept=date[which(date == date[length(time.prior)+1])], size=1, color="tan") +
        geom_hline(yintercept=0, size=1, linetype="dotted", color="grey0") +
        labs(title = Main, y=Ylab, x=Xlab)
    }
  }
  else {
    stop("\n No matching type of plot available.\n")
  }
}

