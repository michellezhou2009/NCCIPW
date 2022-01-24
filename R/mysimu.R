#' Function to simulation NCC data
#'
#' @param N sample size.
#' @param yes.match logical indicating whether the controls are selected with or without matching. The default value is \code{FALSE}.
#' @param m0 a numerical value which is the number of controls selected for each case.
#' @param pi1 the percentage of events that are selected as cases
#'
#' @export
#'
#' @return data a \code{data.frame} with the following variables: \code{time} (censored event time), \code{status} (censoring indicator), \code{marker1} and \code{marker2} (two markers), \code{id} (subjects' ID), \code{case} (indicator for sub-cohort case), \code{control} (indicator for sub-cohort control), \code{CaseID} (indicator which sub-cohort case that each subject is selected in the sub-cohort for).
#'
#' If \code{yes.match=TRUE}, the list includes another two elements: \code{aM} (a vector of the matching constraint for each matching variable, i.e., \code{a0}) and \code{Mdat} (a \code{data.frame} including two matching variables).
#'

mysimu <- function(N, yes.match=FALSE, m0, pi1){
  VTM <- function(vc, dm){matrix(vc, ncol=length(vc), nrow=dm, byrow=T)}
  gam.true=c(0.5,0.5)
  marker1 = rnorm(N); yi = cbind(marker1,marker2=marker1+rnorm(N))
  ti = -c(yi%*%gam.true)+log(-log(runif(N))); ti = exp(ti/2+1.5);
  ai = 2
  ci = pmin(rgamma(N,2,2)+0.1,ai)
  xi = pmin(ti,ci); di = 1*(ti<=ci)
  id=paste("S",1:N,sep="")

  N1 = sum(di); n1 = round(N1*pi1)
  if(yes.match){
    Zi = 1*(pnorm(yi[,1]+rnorm(N))>0.5);
    Zi = cbind(Zi,pnorm(yi[,2]+rnorm(N))*5); Zi = round(Zi)
    a0 = c(0,1)
    colnames(Zi)=c("Mvar1","Mvar2")
  } else {Zi = NULL}
  index.case = sort(sample((1:N)[di==1],size=n1,replace=F));
  vi.case = vi.control= rep(0,N); vi.case[index.case]=1
  case.id = rep(NA,N); case.id[vi.case==1] = id[vi.case==1]

  for(l in 1:length(index.case)){
    tmp.index = index.case[l]; riskset.index = xi>xi[tmp.index];
    ## =========================================================================== ##
    ## if matching, additional constraint of |Zi - Zl| <= a0 needs to be satisfied ##
    ## =========================================================================== ##
    if(yes.match){riskset.index = riskset.index & (apply(abs(Zi-VTM(Zi[tmp.index,],N))<=VTM(a0,N),1,prod)==1)}
    riskset.index = which(riskset.index==T); nl = length(riskset.index);
    ## =========================================================================== ##
    ## if riskset is empty, no control would be selected                           ##
    ## if riskset size < m.match, only select # of available for Finite Population ##
    ## =========================================================================== ##
    if (nl!=0){
      if (nl==1)  control.index = riskset.index else control.index = sample(riskset.index,min(m0,nl));
      vi.control[control.index]=1
      case.id[control.index] = id[tmp.index]
    }
  }

  if (yes.match) return(list(data=data.frame(time=xi,status=di,yi,id=id,case=vi.case,control=vi.control,CaseID=case.id),Mdat=Zi,aM=a0)) else return(list(data=data.frame(time=xi,status=di,yi,id=id,case=vi.case,control=vi.control,CaseID=case.id)))
}

