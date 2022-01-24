#' IPW estimates for Cox PH model using the new weight
#'
#' @description IPW estimates of model parameters and accuracy parameters for Cox's proportional hazards (PH) model with perturbed counterparts.
#'
#' @param formula a formula object, with the response on the left of a \code{~} operator, and the terms on the right. The response must be a survival object as returned by the \code{Surv} function.
#' @param data a \code{data.frame} in which to interpret the variables named in the formula.
#' @param id character, which is the variable name in the \code{data} corresponding to the column that identifies the subjects.
#' @param case character, which is the variable name in the \code{data} corresponding to the column that indicates whether each subject is a sub-cohort case. The value is 1 if the subject is a sub-cohort case, and it is 0 if the subject is not a sub-cohort case.
#' @param control character, which is the variable name in the \code{data} corresponding to the column that indicates whether each subject is a sub-cohort control. The value is 1 if the subject is a sub-cohort control, and it is 0 if the subject is not a sub-cohort control. #' @param m0 a numerical value which is the number of controls selected for each case.
#' @param t0 a numerical value, which is the prediction time period, used for creating the binary outcome status: =1 if the event time is smaller than \code{t0}, = 0 otherwise.
#' @param weight.type character which indicate which sampling weight to be used and whose value could be
#' \describe{
#' \item{"sam"}{the modified Samuelsen's weight by assigning the inverse selection probability to event cases;}
#' \item{"new"}{the modified Samuelsen's weight by assigning 0 to event cases;}
#' \item{"HT"}{the Horvitz-Thompson type of weight, which is the default value.}
#' }
#'
#' @param yes.match logical indicating whether the controls are selected with or without matching; the default value is \code{FALSE}.
#' @param control.matching a list which includes two elements: \code{Mdat} and \code{aM}. \code{Mdat} is a \code{data.frame} which includes the matching variables, and its default value is \code{NULL} if the selection is without matching. \code{aM} is a numerical vector which includes the matching constraints for each matching variable, and its default value is \code{NULL} if the selection is without matching. The length of 'aM' should be the same as the number of columns of \code{Mdat}.
#' @param control.accuracy a list which includes the arguments used in estimating the accuracy measures. It includes two elements: \code{type} and \code{u0}. \code{type} is a character which could be "FPR" or "TPR", and the default value is "FPR". \code{u0} is a numerical value for the accuracy measure specified in \code{type}, and the default value is 0.05.
#' @param yes.ptb logical indicating if obtaining the perturbed counterparts of the IPW estimates, and the default value is \code{FALSE}.
#' @param control.ptb a list of elements that control the perturbation procedure. \code{seed} is a numerical value used in the function \code{\link{set.seed}}, and its default value is 123. \code{n.ptb} is a numerical value indicating the number of perturbed counterparts generated, and its default value is 100. \code{CaseID} is character, which is the variable name in the \code{data} corresponding to the column that listed the case that each subject that is selected in the sub-cohort for.
#'
#' @return a list of two elements: \code{estimates} and \code{ptb}. \code{estimates} is a list of IPW estimates of the model parameters (\code{coef}) and accuracy parameters (\code{acc}). If \code{yes.ptb=TRUE}, \code{ptb} is a list of the perturbed counterparts for the IPW estimates of the model parameters (\code{coef}) and accuracy parameters (\code{acc}); otherwise, \code{ptb = NULL}.
#'
#' @import survival dplyr
#'
#' @examples
#' \code{data("myexample")}
#' \code{PH.IPW(formula=Surv(time,status)~marker1+marker2,
#'       data=myexample$data,
#'       id="id",
#'       case="case",
#'       control="control",
#'       m0=3,t0=1,
#'       yes.match=T,control.matching=list(Mdat=myexample$Mdat,aM=myexample$aM),
#'       yes.ptb=TRUE,control.ptb=list(n.ptb=10,CaseID="CaseID"))}

PH.IPW <- function(formula,data,id,case,control,m0,t0,weight.type="HT",yes.match=FALSE,control.matching=list(Mdat=NULL,aM=NULL),control.accuracy=list(type="FPR",u0=0.05),yes.ptb=FALSE,control.ptb=list(seed=123,n.ptb=500,CaseID=NULL)){

  varnm.all = all.vars(formula)
  exist.yes = is.element(varnm.all,colnames(data))

  ## warnings
  if (!all(exist.yes)) stop("\"data\" does not have variables:",paste(varnm.all[!exist.yes],collapse=", "))
  if (is.null(id)) stop("\"id\" should be provided") else {
    if (!is.element(id,colnames(data))) stop("data should include the column with the name provided in \"id\"") else id = as.character(data[,id])
  }
  if (is.null(case)) stop("\"case\" should be provided") else {
    if (!is.element(case,colnames(data))) stop("data should include the column with the name provided in \"case\"") else {
      case = data[,case]
      if (length(unique(case))!=2 | !all(is.element(c(0,1),unique(case)))) stop("The values of \"case\" shoud be 0 or 1.")
    }
  }

  if (is.null(control)) stop("\"control\" should be provided") else{
    if (!is.element(control,colnames(data))) stop("data should include the column with the name provided in \"control\"") else {
      control = data[,control]
      if (length(unique(control))!=2 | !all(is.element(c(0,1),unique(control)))) stop("The values of \"control\" shoud be 0 or 1.")
    }
  }


  N = nrow(data)
  if (yes.match){
    M.dat = control.matching$Mdat
    a.M = control.matching$aM
    if (length(a.M)!=ncol(M.dat)) stop("The length of \"aM\" should be the same as the number of columns in \"Mdat\".")
  } else {M.dat = NULL; a.M=NULL}
  if (yes.ptb){
    CaseID = control.ptb$CaseID
    if (is.null(CaseID)) stop("\"control\" should be provided") else{
      if (!is.element(CaseID,colnames(data))) stop("data should include the column with the name provided in \"control\"") else CaseID = as.character(data$CaseID)
    }
    n.ptb = control.ptb$n.ptb
    set.seed(control.ptb$seed); Vii.ptb = array(rexp(N*n.ptb),dim=c(N,n.ptb))
    if (weight.type=="sam") stop("Perturbation resampling is not provided for \"sam\" weight") else {
      if (weight.type=="new")  wgt.ncc = newWGT(xi=data[,varnm.all[1]],di=data[,varnm.all[2]],id=id,case=case,control=control,m0,M.dat=M.dat,a.M=a.M,yes.ptb=yes.ptb,control.list=list(CaseID=CaseID, Vii.ptb = Vii.ptb,yes.parallel=FALSE)) else wgt.ncc = htWGT(xi=data[,varnm.all[1]],di=data[,varnm.all[2]],id=id,case=case,control=control,m0,M.dat=M.dat,a.M=a.M,yes.ptb=yes.ptb,control.list=list(CaseID=CaseID, Vii.ptb = Vii.ptb,yes.parallel=FALSE))
    }
    wgt.cen = WGT.CEN(Ti=data[,varnm.all[1]],Di=data[,varnm.all[2]],t0=t0,Vii.ptb=Vii.ptb)
    fit.all = sapply(seq(n.ptb+1),function(k){
      workdat = data; workdat$w.ncc = wgt.ncc[,k]; workdat$w.cen = wgt.cen[,k]
      workdat = workdat[workdat$w.ncc>0,]
      fit.m = coxph(formula,data=workdat,robust=T,weights=w.ncc)
      risk.score = predict(fit.m,type="lp")
      acc = ACC.FUN(1*(workdat[,varnm.all[1]]<t0),risk.score,wgtk=workdat$w.ncc*workdat$w.cen,
                    type=control.accuracy$type, u0=control.accuracy$u0)
      return(list(coef=fit.m$coefficients,acc=acc))
    },simplify = FALSE)
    coef.all = do.call(cbind,lapply(fit.all,function(fit)fit$coef)); acc.all = do.call(cbind,lapply(fit.all,function(fit)fit$acc))
    return(list(estimates = list(coef=coef.all[,1],acc=acc.all[,1]),
                ptb = list(coef=coef.all[,-1],acc=acc.all[,-1])
    ))
  } else {
    switch(weight.type,
           "sam"={
             wgt.ncc = samWGT(xi=data[,varnm.all[1]],di=data[,varnm.all[2]],id=id,case=case,control=control,m0,M.dat=M.dat,a.M=a.M)
           },
           "new"={
             wgt.ncc = newWGT(xi=data[,varnm.all[1]],di=data[,varnm.all[2]],id=id,case=case,control=control, m0, M.dat=M.dat, a.M=a.M, yes.ptb=FALSE)
           },
           "HT"={
             wgt.ncc = htWGT(xi=data[,varnm.all[1]],di=data[,varnm.all[2]],id=id,case=case,control=control, m0, M.dat=M.dat, a.M=a.M, yes.ptb=FALSE)
           }
           )

    wgt.cen = WGT.CEN(Ti=data[,varnm.all[1]],Di=data[,varnm.all[2]],t0=t0,Vii.ptb=NULL)
    workdat = data; workdat$w.ncc = wgt.ncc; workdat$w.cen = wgt.cen
    workdat = workdat[workdat$w.ncc>0,]
    fit.m = coxph(formula, data=workdat,robust=T,weights=w.ncc)
    risk.score = predict(fit.m,type="lp")
    acc = ACC.FUN(1*(workdat[,varnm.all[1]]<t0),risk.score,wgtk=workdat$w.ncc*workdat$w.cen,
                  type=control.accuracy$type, u0=control.accuracy$u0)
    return(list(estimates=list(coef=fit.m$coefficients,acc=acc),
                ptb=NULL))
  }
}

