% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samIPW.R
\name{PH.samIPW}
\alias{PH.samIPW}
\title{IPW estimates for Cox's PH model using Samuelsen (1997)‘s weight}
\usage{
PH.samIPW(
  formula,
  data,
  id,
  case,
  control,
  m0,
  t0,
  yes.matching = FALSE,
  control.matching = list(Mdat = NULL, aM = NULL),
  control.accuracy = list(type = "FPR", u0 = 0.05)
)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a \code{~} operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.}

\item{data}{a \code{data.frame} in which to interpret the variables named in the formula.}

\item{id}{a vector which identifies the subjects. The length of 'id' should be the same as the number of rows of 'data'.}

\item{case}{a vector which indicates whether each subject is a sub-cohort case. The value is 1 if the subject is a sub-cohort case, and it is 0 if the subject is not a sub-cohort case. The length of 'case' should be the same as the number of rows of 'data'.}

\item{control}{a vector which indicates whether each subject is a sub-cohort control. The value is 1 if the subject is a sub-cohort control, and it is 0 if the subject is not a sub-cohort control. The length of 'control' should be the same as the number of rows of 'data'.}

\item{m0}{a numerical value which is the number of controls selected for each case.}

\item{t0}{a numerical value, which is the prediction time period, used for creating the binary outcome status: =1 if the event time is smaller than \code{t0}, = 0 otherwise.}

\item{yes.matching}{logical indicating whether the controls are selected with or without matching; the default value is \code{FALSE}.}

\item{control.matching}{a list which includes two elements: \code{Mdat} and \code{aM}. \code{Mdat} is a \code{data.frame} which includes the matching variables, and its default value is \code{NULL} if the selection is without matching. \code{aM} is a numerical vector which includes the matching constraints for each matching variable, and its default value is \code{NULL} if the selection is without matching. The length of 'aM' should be the same as the number of columns of \code{Mdat}.}

\item{control.accuracy}{a list which includes the arguments used in estimating the accuracy measures. It includes two elements: \code{type} and \code{u0}. \code{type} is a character which could be "FPR" or "TPR", and the default value is "FPR". \code{u0} is a numerical value for the accuracy measure specified in \code{type}, and the default value is 0.05.}
}
\value{
a list of two elements: the first element \code{coef} includes IPW estimates of the model parameters, and the second element \code{acc} includes IPW estimates of the accuracy parameters.
}
\description{
IPW estimates of model parameters and accuracy parameters for Cox's proportional hazards (PH) model using Samuelsen (1997)‘s weight.
}
\examples{
\code{data("myexample")}
\code{PH.samIPW(formula=Surv(time,status)~marker1+marker2,
      data=myexample$data,
      id=myexample$id,
      case=myexample$case,
      control=myexample$control,
      m0=myexample$m0,t0=1)}

}
