#' IPW estimates for Cox's PH model using Samuelsen (1997)‘s weight
#'
#' @description IPW estimates of model parameters and accuracy parameters for Cox's proportional hazards (PH) model using Samuelsen (1997)‘s weight.
#'
#' @param formula a formula object, with the response on the left of a \code{~} operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.
#' @param data a \code{data.frame} in which to interpret the variables named in the formula.
#' @param id a vector which identifies the subjects. The length of 'id' should be the same as the number of rows of 'data'.
#' @param case a vector which indicates whether each subject is a sub-cohort case. The value is 1 if the subject is a sub-cohort case, and it is 0 if the subject is not a sub-cohort case. The length of 'case' should be the same as the number of rows of 'data'.
#' @param control a vector which indicates whether each subject is a sub-cohort control. The value is 1 if the subject is a sub-cohort control, and it is 0 if the subject is not a sub-cohort control. The length of 'control' should be the same as the number of rows of 'data'.
#' @param m0 a numerical value which is the number of controls selected for each case.
#' @param t0 a numerical value, which is the prediction time period, used for creating the binary outcome status: =1 if the event time is smaller than \code{t0}, = 0 otherwise.
#' @param yes.matching logical indicating whether the controls are selected with or without matching; the default value is \code{FALSE}.
#' @param control.matching a list which includes two elements: \code{Mdat} and \code{aM}. \code{Mdat} is a \code{data.frame} which includes the matching variables, and its default value is \code{NULL} if the selection is without matching. \code{aM} is a numerical vector which includes the matching constraints for each matching variable, and its default value is \code{NULL} if the selection is without matching. The length of 'aM' should be the same as the number of columns of \code{Mdat}.
#' @param control.accuracy a list which includes the arguments used in estimating the accuracy measures. It includes two elements: \code{type} and \code{u0}. \code{type} is a character which could be "FPR" or "TPR", and the default value is "FPR". \code{u0} is a numerical value for the accuracy measure specified in \code{type}, and the default value is 0.05.
#'
#' @export
#'
#' @return a list of two elements: the first element \code{coef} includes IPW estimates of the model parameters, and the second element \code{acc} includes IPW estimates of the accuracy parameters.
#'
#' @examples
#' \code{data("myexample")}
#' \code{PH.samIPW(formula=Surv(time,status)~marker1+marker2,
#'       data=myexample$data,
#'       id=myexample$id,
#'       case=myexample$case,
#'       control=myexample$control,
#'       m0=myexample$m0,t0=1)}
#'
PH.samIPW <- function(formula,data,id,case,control,m0,t0,
                      yes.matching=FALSE, control.matching=list(Mdat=NULL,aM=NULL),
                      control.accuracy=list(type="FPR",u0=0.05)){

  varnm.all = all.vars(formula)
  exist.yes = is.element(varnm.all,colnames(data))

  ## warnings
  if (!all(exist.yes)) stop("\"data\" does not have variables:",paste(varnm.all[!exist.yes],collapse=", "))
  if (!is.null(id) & length(id)!=nrow(data)) stop("The length of \"id\" should be the same as the number of rows in \"data\".")
  if (length(case)!=nrow(data)) stop("The length of \"case\" should be the same as the number of rows in \"data\".")
  if (length(control)!=nrow(data)) stop("The length of \"control\" should be the same as the number of rows in \"data\".")
  if (length(unique(case))!=2 | !all(is.element(c(0,1),unique(case)))) stop("The values of \"case\" shoud be 0 or 1.")
  if (length(unique(case))!=2 | !all(is.element(c(0,1),unique(case)))) stop("The values of \"case\" should be 0 or 1.")

  if (is.null(id)) id=as.character(seq(nrow(data)))
  N = nrow(data)
  if (yes.matching){
    M.dat = control.matching$Mdat
    a.M = control.matching$aM
    if (length(a.M)!=ncol(M.dat)) stop("The length of \"aM\" should be the same as the number of columns in \"Mdat\".")
  } else {M.dat = NULL; a.M=NULL}
  wgt.ncc = samWGT(xi=data[,varnm.all[1]],di=data[,varnm.all[2]],id=id,case=case,control=control,m0=m0,M.dat=M.dat,a.M=a.M)
  wgt.cen = WGT.CEN(Ti=data[,varnm.all[1]],Di=data[,varnm.all[2]],t0=t0,Vii.ptb=NULL)
  workdat = data; workdat$w.ncc = wgt.ncc; workdat$w.cen = wgt.cen
  workdat = workdat[workdat$w.ncc>0,]
  fit.m = coxph(formula, data=workdat,robust=T,weights=w.ncc)

  risk.score = predict(fit.m,type="lp")
  acc = ACC.FUN(1*(workdat[,varnm.all[1]]<t0),risk.score,wgtk=workdat$w.ncc*workdat$w.cen,
                type=control.accuracy$type, u0=control.accuracy$u0)
  return(list(coef=fit.m$coefficients,acc=acc))
}

#' Double IPW estimates for a time-dependent GLM using Samuelsen (1997)‘s weight
#'
#' @description Double IPW estimates of model and accuracy parameters for a time-dependent generalized linear model (GLM) using Samuelsen (1997)‘s weight.
#'
#' @param formula a formula object, with the response on the left of a \code{~} operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.
#' @param data a \code{data.frame} in which to interpret the variables named in the formula.
#' @param family a description of the error distribution and link function to be used in the generalized linear model, which is the argument used in the \code{\link{glm}} function. The default is \code{binominal(link=logit)}.
#' @param id a vector which identifies the subjects. The length of 'id' should be the same as the number of rows of 'data'.
#' @param case a vector which indicates whether each subject is a sub-cohort case. The value is 1 if the subject is a sub-cohort case, and it is 0 if the subject is not a sub-cohort case. The length of 'case' should be the same as the number of rows of 'data'.
#' @param control a vector which indicates whether each subject is a sub-cohort control. The value is 1 if the subject is a sub-cohort control, and it is 0 if the subject is not a sub-cohort control. The length of 'control' should be the same as the number of rows of 'data'.
##'@param m0 a numerical value which is the number of controls selected for each case.
#' @param t0 a numerical value, which is the prediction time period, used for creating the binary outcome status: =1 if the event time is smaller than \code{t0}, = 0 otherwise.
#' @param yes.matching logical indicating whether the controls are selected with or without matching; the default value is \code{FALSE}.
#' @param control.matching a list which includes two elements: \code{Mdat} and \code{aM}. \code{Mdat} is a \code{data.frame} which includes the matching variables, and its default value is \code{NULL} if the selection is without matching. \code{aM} is a numerical vector which includes the matching constraints for each matching variable, and its default value is \code{NULL} if the selection is without matching. The length of 'aM' should be the same as the number of columns of \code{Mdat}.
#' @param control.accuracy a list which includes the arguments used in estimating the accuracy measures. It includes two elements: \code{type} and \code{u0}. \code{type} is a character which could be "FPR" or "TPR", and the default value is "FPR". \code{u0} is a numerical value for the accuracy measure specified in \code{type}, and the default value is 0.05.
#'
#' @export
#'
#' @return a list of two elements: the first element \code{coef} includes the IPW estimates of the model parameters, and the second element \code{acc} includes the IPW estimates of the accuracy parameters.
#'
#' @examples
#' \code{data("myexample")}
#' \code{GLM.samIPW(formula=Surv(time,status)~marker1+marker2,
#'       data=myexample$data, family=binomial(link=logit),
#'       id=myexample$id,
#'       case=myexample$case,
#'       control=myexample$control,
#'       m0=myexample$m0,t0=1)}

GLM.samIPW <- function(formula,data,family=binomial(link=logit),id,case,control,m0,t0,
                      yes.matching=FALSE, control.matching=list(Mdat=NULL,aM=NULL),
                      control.accuracy=list(type="FPR",u0=0.05)){

  varnm.all = all.vars(formula)
  exist.yes = is.element(varnm.all,colnames(data))

  ## warnings
  if (!all(exist.yes)) stop("\"data\" does not have variables:",paste(varnm.all[!exist.yes],collapse=", "))
  if (!is.null(id) & length(id)!=nrow(data)) stop("The length of \"id\" should be the same as the number of rows in \"data\".")
  if (length(case)!=nrow(data)) stop("The length of \"case\" should be the same as the number of rows in \"data\".")
  if (length(control)!=nrow(data)) stop("The length of \"control\" should be the same as the number of rows in \"data\".")
  if (length(unique(case))!=2 | !all(is.element(c(0,1),unique(case)))) stop("The values of \"case\" shoud be 0 or 1.")
  if (length(unique(case))!=2 | !all(is.element(c(0,1),unique(case)))) stop("The values of \"case\" should be 0 or 1.")

  if (is.null(id)) id=as.character(seq(nrow(data)))
  N = nrow(data)
  if (yes.matching){
    M.dat = control.matching$Mdat
    a.M = control.matching$aM
    if (length(a.M)!=ncol(M.dat)) stop("The length of \"aM\" should be the same as the number of columns in \"Mdat\".")
  } else {M.dat = NULL; a.M=NULL}

  wgt.ncc = samWGT(xi=data[,varnm.all[1]],di=data[,varnm.all[2]],id=id,case=case,control=control,m0,M.dat=M.dat,a.M=a.M)
  wgt.cen = WGT.CEN(Ti=data[,varnm.all[1]],Di=data[,varnm.all[2]],t0=t0,Vii.ptb=NULL)
  workdat = data; workdat$wgt = wgt.ncc*wgt.cen
  workdat$resp = 1*(workdat[,varnm.all[1]]<t0);
  workdat = workdat[,-match(varnm.all[c(1,2)],colnames(workdat))]
  workdat = workdat[workdat$wgt>0,]
  formula = update(formula,resp~.)
  fit.m = glm(formula,data=workdat,family=family,weights=wgt)
  risk.score = predict(fit.m,type="response")
  acc = ACC.FUN(workdat$resp,risk.score,wgtk=workdat$wgt,
                type=control.accuracy$type, u0=control.accuracy$u0)
  return(list(coef=fit.m$coefficients,acc=acc))
}
