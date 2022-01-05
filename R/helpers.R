#' Calculate the censoring weights
#'
#' @description Calculate the weights accounting for missing binary outcomes due to censoring.
#'
#' @param Ti a vector of censored event times.
#' @param Di a vector of binary values (1 or 0) indicating whether the true event time is observed.
#' @param t0 a numerical value which is the prediction time period, used for creating the binary outcome status. The value is 1 if the event time is smaller than t0; otherwise, it is 0.
#' @param Vii.ptb a vector of perturbation variables. The default is \code{NULL} if the perturbation procedure is not implemented for estimating the variance of the IPW estimators.
#'
#' @export
#'
#' @return If \code{Vii.ptb=NULL}, it returns a single-column matrix which includes the censoring weights; otherwise, it returns a matrix with the first column being the censoring weights and the remaining columns being their perturbed counterparts.
#'
#'
WGT.CEN <- function(Ti,Di,t0,Vii.ptb=NULL){
  ## ================================================================ ##
  ## ============== KM Estimator of Censoring Survival ============== ##
  ## ================================================================ ##
  Ghat.FUN <- function(tt, Ti, Di,type='fl',wgt=NULL){
    surv.fit = survfit(Surv(Ti,1-Di)~1,se.fit=F,type=type,weights=wgt)
    surv.ti = surv.fit$time; surv.fit.surv = surv.fit$surv
    surv.til = sort(surv.fit$time); surv.fit.surv = surv.fit.surv[order(surv.ti)]
    # tmpind1 = sum.I(tt,">=",surv.til) + 1
    tmpind = sapply(tt,function(u){sum(surv.til<=u)}) + 1
    c(1,surv.fit.surv)[tmpind]
  }
  N = length(Ti)
  tt <- c(t0,Ti[Di==1 & Ti<=t0])
  Ghat.all = Ghat.FUN(tt,Ti,Di)
  Ghat.t0 = Ghat.all[1]; Ghat.Ti1 = Ghat.all[-1]
  Wi <- rep(0,length(Ti))
  if (sum(Ti<=t0 & Di==1)>0) Wi[Ti <= t0 & Di==1] <- 1/Ghat.Ti1
  if (sum(Ti>t0)>0) Wi[Ti >  t0] <- 1/Ghat.t0
  if (is.null(Vii.ptb)) return(matrix(Wi,byrow=F,ncol=1,nrow=length(Wi))) else {
    Ghat.all.ptb = apply(Vii.ptb,2,function(x){Ghat.FUN(tt,Ti,Di,wgt=x)})
    Ghat.t0.ptb = Ghat.all.ptb[1,]; Ghat.Ti.ptb = Ghat.all.ptb[-1,]
    Wi.ptb <- matrix(0,nrow=length(Ti),ncol=ncol(Vii.ptb))
    if (sum(Ti<=t0 & Di==1)>0) Wi.ptb[Ti<=t0 & Di==1,] <- 1/Ghat.Ti.ptb
    if (sum(Ti>t0)>0) Wi.ptb[Ti>t0,] <- 1/Ghat.t0.ptb
    return(cbind(Wi,Wi.ptb))
  }
}

#' Calculate the new weight by reformularizing the NCC sampling as a two-stage stratified sampling
#'
#' @description This proposed weight is the Horvitz-Thompson type of weight. 
#'
#' @param xi a vector of censored event times.
#' @param di a vector of binary values (1 or 0) indicating whether the true event time is observed.
#' @param id a vector which identifies the subjects.
#' @param case a vector which indicates whether each subject is a sub-cohort case. The value is 1 if the subject is a sub-cohort case, and it is 0 if the subject is not a sub-cohort case.
#' @param control a vector which indicates whether each subject is a sub-cohort control. The value is 1 if the subject is a sub-cohort control, and it is 0 if the subject is not a sub-cohort control.
##'@param m0 a numerical value which is the number of controls selected for each case.
#' @param M.dat a \code{data.frame} which includes the matching variables. The default value is \code{NULL} if the selection is without matching.
#' @param a.M a numerical vector which includes the matching constraint for each matching variable. The default value is \code{NULL} if the selection is without matching.
#' @param yes.ptb logical indicating if obtaining the perturbed counterparts of the IPW estimates, and the default value is \code{FALSE}.
#' @param control.list a list including two elements: \code{CaseID} and \code{Vii.ptb} that control the perturbation procedure. \code{CaseID} is a vector in which listed the case that each subject that is selected in the sub-cohort for. If \code{id} is the same as \code{CaseID}, it indicates this subject is one of the sub-cohort cases. \code{Vii.ptb} is a vector of perturbation variables, and the default value is \code{NULL} if the perturbation procedure is not implemented for estimating the variance of the IPW estimators.
#'
#' @export
#'
#' @return If \code{yes.ptb=FALSE}, it returns a single-column matrix which includes the NCC weights; otherwise, it returns a matrix with the first column being the NCC weights and the remaining columns being their perturbed counterparts.
#'

htWGT <- function(xi,di,id,case,control, m0, M.dat=NULL, a.M=NULL, yes.ptb=FALSE, control.list=list(CaseID=NULL, Vii.ptb=NULL,yes.parallel=FALSE)){
  N = length(id)
  VTM <- function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
  }
  if (!is.null(M.dat)) {p.M = ncol(M.dat)}
  
  ti1 = xi[case==1]; n1 = length(ti1); pi1=n1/sum(di==1); case.id = id[case==1]
  Ii1j = VTM(xi,n1) >= ti1 # n1xN matrix I(xi>=tj)
  if (!is.null(M.dat)){
    M1 = M.dat[case==1,,drop=F]
    for (ll in 1:p.M){Ii1j = Ii1j*(abs(VTM(M.dat[,ll],n1)-M1[,ll])<=a.M[ll])}
  } # N1 x N matrix
  n.Ri1 = apply(Ii1j,1,sum); names(n.Ri1)=case.id# size of risk set for each case;
  junk = Ii1j*log(1-m0/pmax(n.Ri1-1,m0)) # n1 x n0 matrix: sum over all the selected cases
  P0hati = 1- exp(apply(junk,2,sum,na.rm=T))
  subcohort = case + (1-case)*control
  wgt = rep(0,N); wgt[di==1 & subcohort==1] = 1/(pi1+(1-pi1)*P0hati[di==1 & subcohort==1]); wgt[di==0 & subcohort==1] = 1/P0hati[di==0 & subcohort==1]
  
  if (!yes.ptb) return(wgt) else {
    CaseID = control.list$CaseID; Vii.ptb = control.list$Vii.ptb
    if (pi1<1) pi1.ptb = apply(Vii.ptb[case==1,],2,sum)/apply(Vii.ptb[di==1,],2,sum) else pi1.ptb = rep(1,ncol(Vii.ptb))
    subcohort.ptb = subcohort*Vii.ptb
    V0ij = data.frame(case=CaseID,control=id) %>% 
      filter(!is.na(case)) %>%
      filter(case!=control) %>%
      arrange(case)
    V0ij.ptb = matrix(rexp(nrow(V0ij)*ncol(Vii.ptb)),nrow=nrow(V0ij),ncol=ncol(Vii.ptb))
    control.id = id[di==0 & control==1]
    control.ptb = lapply(control.id,function(u){
      1-apply(1-V0ij.ptb[V0ij[,2]==u,,drop=F],2,prod)
    }) %>% do.call(rbind,.)
    if (control.list$yes.parallel){
      p0hati.ptb = foreach(u=control.id,.combine = rbind) %dopar%{
        id.cases = case.id[Ii1j[,id==u]==1]
        junk = lapply(id.cases,function(v){
          1-apply(V0ij.ptb[V0ij[,1]==v,,drop=F],2,sum)/(n.Ri1[names(n.Ri1)==v]-1)})
        junk = do.call(rbind,junk)
        1-apply(junk,2,prod)
      }
    } else{
      p0hati.ptb = lapply(control.id,function(u){
        id.cases = case.id[Ii1j[,id==u]==1]
        junk = lapply(id.cases,function(v){
        1-apply(V0ij.ptb[V0ij[,1]==v,,drop=F],2,sum)/(n.Ri1[names(n.Ri1)==v]-1)})
        junk = do.call(rbind,junk)
        1-apply(junk,2,prod)
      }) %>% do.call(rbind,.)
    }
    wgt.ptb = matrix(0,nrow=N,ncol=ncol(Vii.ptb))
    wgt.ptb[di==1 & subcohort==1] = subcohort.ptb[di==1 & subcohort==1,]/VTM(pi1.ptb,sum(di==1 & subcohort==1))
    wgt.ptb[di==0 & control==1] = control.ptb/p0hati.ptb
    return(cbind(wgt,wgt.ptb))
  }
  
  return(wgt)
}

#' Calculate the proposed weight which weights the event controls with the inverse selection probability
#'
#' @description This proposed weight modifies the method in Samuelsen (1997) by weighting the event controls with the inverse selection probability. 
#'
#' @param xi a vector of censored event times.
#' @param di a vector of binary values (1 or 0) indicating whether the true event time is observed.
#' @param id a vector which identifies the subjects.
#' @param case a vector which indicates whether each subject is a sub-cohort case. The value is 1 if the subject is a sub-cohort case, and it is 0 if the subject is not a sub-cohort case.
#' @param control a vector which indicates whether each subject is a sub-cohort control. The value is 1 if the subject is a sub-cohort control, and it is 0 if the subject is not a sub-cohort control.
##'@param m0 a numerical value which is the number of controls selected for each case.
#' @param M.dat a \code{data.frame} which includes the matching variables. The default value is \code{NULL} if the selection is without matching.
#' @param a.M a numerical vector which includes the matching constraint for each matching variable. The default value is \code{NULL} if the selection is without matching.
#'
#' @export
#'
#' @return a vector of NCC weights.
#'
#' @references
#' Samuelsen, S. (1997). A psudolikelihood approach to analysis of nested case-control studies. Biometrika, 84(2):379–394.
#'

samWGT <- function(xi,di,id,case,control,m0,M.dat=NULL,a.M=NULL){
  N = length(id)
  VTM <- function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
  }
  if (!is.null(M.dat)) {p.M = ncol(M.dat)}
  ti1 = xi[case==1]; n1 = length(ti1); case.id = id[case==1]
  Ii1j = VTM(xi,n1) >= ti1 # n1xN matrix I(xi>=tj)
  if (!is.null(M.dat)){
    M1 = M.dat[case==1,,drop=F]
    for (ll in 1:p.M){Ii1j = Ii1j*(abs(VTM(M.dat[,ll],n1)-M1[,ll])<=a.M[ll])}
  } # N1 x N matrix
  n.Ri1 = apply(Ii1j,1,sum); names(n.Ri1)=case.id# size of risk set for each case;
  junk = Ii1j*log(1-m0/pmax(n.Ri1-1,m0)) # n1 x n0 matrix: sum over all the selected cases
  P0hati = 1- exp(apply(junk,2,sum,na.rm=T))
  wgtI = rep(0,N); wgtI[case==1] = 1; wgtI[case==0 & control==1] = 1/P0hati[case==0 & control==1]
  return(wgtI)
}

#' Calculate the new weight which weights the event controls with 0
#'
#' @description This proposed weight modifies the method in Samuelsen (1997) by weighting the event controls with 0.
#'
#' @param xi a vector of censored event times.
#' @param di a vector of binary values (1 or 0) indicating whether the true event time is observed.
#' @param id a vector which identifies the subjects.
#' @param case a vector which indicates whether each subject is a sub-cohort case. The value is 1 if the subject is a sub-cohort case, and it is 0 if the subject is not a sub-cohort case.
#' @param control a vector which indicates whether each subject is a sub-cohort control. The value is 1 if the subject is a sub-cohort control, and it is 0 if the subject is not a sub-cohort control.
##'@param m0 a numerical value which is the number of controls selected for each case.
#' @param M.dat a \code{data.frame} which includes the matching variables. The default value is \code{NULL} if the selection is without matching.
#' @param a.M a numerical vector which includes the matching constraint for each matching variable. The default value is \code{NULL} if the selection is without matching.
#' @param yes.ptb logical indicating if obtaining the perturbed counterparts of the IPW estimates, and the default value is \code{FALSE}.
#' @param control.list a list including two elements: \code{CaseID} and \code{Vii.ptb} that control the perturbation procedure. \code{CaseID} is a vector in which listed the case that each subject that is selected in the sub-cohort for. If \code{id} is the same as \code{CaseID}, it indicates this subject is one of the sub-cohort cases. \code{Vii.ptb} is a vector of perturbation variables, and the default value is \code{NULL} if the perturbation procedure is not implemented for estimating the variance of the IPW estimators.
#'
#' @export
#'
#' @return If \code{yes.ptb=FALSE}, it returns a single-column matrix which includes the NCC weights; otherwise, it returns a matrix with the first column being the NCC weights and the remaining columns being their perturbed counterparts.
#'
newWGT <- function(xi,di,id,case,control, m0, M.dat=NULL, a.M=NULL, yes.ptb=FALSE, control.list=list(CaseID=NULL, Vii.ptb=NULL,yes.parallel=FALSE)){
  N = length(id)
  VTM <- function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
  }
  if (!is.null(M.dat)) {p.M = ncol(M.dat)}

  ti1 = xi[case==1]; n1 = length(ti1); pi1=n1/sum(di==1); case.id = id[case==1]
  Ii1j = VTM(xi,n1) >= ti1 # n1xN matrix I(xi>=tj)
  if (!is.null(M.dat)){
    M1 = M.dat[case==1,,drop=F]
    for (ll in 1:p.M){Ii1j = Ii1j*(abs(VTM(M.dat[,ll],n1)-M1[,ll])<=a.M[ll])}
  } # N1 x N matrix
  n.Ri1 = apply(Ii1j,1,sum); names(n.Ri1)=case.id# size of risk set for each case;
  junk = Ii1j*log(1-m0/pmax(n.Ri1-1,m0)) # n1 x n0 matrix: sum over all the selected cases
  P0hati = 1- exp(apply(junk,2,sum,na.rm=T))
  wgt = rep(0,N); wgt[case==1] = 1/pi1; wgt[di==0 & control==1] = 1/P0hati[di==0 & control==1]

  if (!yes.ptb) return(wgt) else {
    CaseID = control.list$CaseID; Vii.ptb = control.list$Vii.ptb
    if (pi1<1) pi1.ptb = apply(Vii.ptb[case==1,],2,sum)/apply(Vii.ptb[di==1,],2,sum) else pi1.ptb = rep(1,ncol(Vii.ptb))
    case.ptb = case*Vii.ptb
    V0ij = data.frame(case=CaseID,control=id) %>% 
      filter(!is.na(case)) %>%
      filter(case!=control) %>%
      arrange(case)
    
    V0ij.ptb = matrix(rexp(nrow(V0ij)*ncol(Vii.ptb)),nrow=nrow(V0ij),ncol=ncol(Vii.ptb))
    control.id = id[di==0 & control==1]
    control.ptb = lapply(control.id,function(u){
      1-apply(1-V0ij.ptb[V0ij[,2]==u,,drop=F],2,prod)
    }) %>% do.call(rbind,.)
    if (control.list$yes.parallel){
      p0hati.ptb = foreach(u=control.id,.combine = rbind) %dopar%{
        id.cases = case.id[Ii1j[,id==u]==1]
        junk = lapply(id.cases,function(v){
          1-apply(V0ij.ptb[V0ij[,1]==v,,drop=F],2,sum)/(n.Ri1[names(n.Ri1)==v]-1)})
        junk = do.call(rbind,junk)
        1-apply(junk,2,prod)
      }
    } else{
      p0hati.ptb = lapply(control.id,function(u){
        id.cases = case.id[Ii1j[,id==u]==1]
        junk = lapply(id.cases,function(v){
          1-apply(V0ij.ptb[V0ij[,1]==v,,drop=F],2,sum)/(n.Ri1[names(n.Ri1)==v]-1)})
        junk = do.call(rbind,junk)
        1-apply(junk,2,prod)
      }) %>% do.call(rbind,.)
    }

    wgt.ptb = matrix(0,nrow=N,ncol=ncol(Vii.ptb))
    wgt.ptb[case==1] = case.ptb[case==1,]/VTM(pi1.ptb,sum(case==1))
    wgt.ptb[di==0 & control==1] = control.ptb/p0hati.ptb
    return(cbind(wgt,wgt.ptb))
  }
}

#' Calculate the accuracy measures
#'
#' @description Calculate the (i) AUC, and (ii) TPR, PPV, and NPV at the cut-off value that leads to FPR being the given value, or FPR, PPV, and NPNV at the cut-off value that leads to TPR being the given value.
#'
#' @param yk a vector of binary outcomes.
#' @param ck a vector of risk scores.
#' @param wgtk a vector of weights.
#' @param type a character which could be "FPR" or "TPR", and the default value is "FPR".
#' @param u0 a numerical value for the accuracy measure specified in \code{type}, and the default value is 0.05.
#'
#' @export
#'
#' @return a vector of IPW estimates of the accuracy measures.

ACC.FUN <- function(yk,ck,wgtk,type="FPR",u0=0.05){

  sum.I <- function(yy,FUN,Yi,Vi=NULL){
    if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
    # for each distinct ordered failure time t[j], number of Xi < t[j]
    pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
    if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
    if (!is.null(Vi)) {
      ## if FUN contains '=', tmpind is the order of decending
      if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
      ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
      Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
      return(rbind(0,Vi)[pos+1,])
    } else return(pos)
  }

  scl = sort(ck); nv = length(ck)
  TPR.cl = sum.I(scl,"<",ck,wgtk*(yk==1))/sum(wgtk*(yk==1))
  FPR.cl = sum.I(scl,"<",ck,wgtk*(yk==0))/sum(wgtk*(yk==0))
  PPV.cl = sum.I(scl,"<",ck,wgtk*(yk==1))/sum.I(scl,"<",ck,wgtk)
  NPV.cl = sum.I(scl,">=",ck,wgtk*(yk==0))/sum.I(scl,">=",ck,wgtk)
  AUC = sum(sum.I(ck,"<",ck,wgtk*(yk==1))*wgtk*(yk==0))/(sum(wgtk*(yk==1))*sum(wgtk*(yk==0)))

  acc.cl = data.frame("cutoff"=scl,"FPR"=FPR.cl,"TPR"=TPR.cl,"NPV"=NPV.cl, "PPV"=PPV.cl)

  ## transform from c by TPR/FPR
  acc.ul = acc.cl; ind0 = match(type,names(acc.cl));
  ul = acc.ul[,ind0]; acc.ul = acc.ul[order(ul),]; ul = sort(ul); # sorted by ul
  if(ind0==1){ indu = sum.I(u0,">=",ul) } else { indu = sum.I(u0,">",ul) }
  c.u0 = acc.ul[indu,1];
  return(c("AUC"=AUC, unlist(acc.ul[indu,-c(1,ind0)])))
}
