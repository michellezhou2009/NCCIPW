#' Function to simulation NCC data
#'
#' @param N sample size.
#' @param yes.matching logical indicating whether the controls are selected with or without matching. The default value is \code{FALSE}.
#' @param a0 a numerical vector which includes the matching constraint for each matching variable. The default value is \code{NULL} if \code{yes.matching=FALSE}.
#' @param m0 a numerical value which is the number of controls selected for each case.
#' @param pi1 the percentage of events that are selected as cases
#'
#' @export
#'
#' @return a list of elements:
#' \describe{
#' \item{data}{a \code{data.frame} with the following variables: \code{time} (censored event time), \code{status} (censoring indicator), \code{marker1} and \code{marker2} (two markers).}
#' \item{id}{a vector which identifies the subjects.}
#' \item{case}{a vector which indicates whether each subject is a sub-cohort case.}
#' \item{control}{a vector which indicates whether each subject is a sub-cohort control.}
#' \item{m0}{number of controls selected for each case.}
#' \item{matchid}{a matrix in which the number of rows is the same as the number of selected sub-cohort cases, and the number of columns is \code{m0}+1. In each row, the first element is the ID of a case and the remaining elements are the IDs of the controls selected for the given case.}
#' }
#' If \code{yes.matching=TRUE}, the list includes another two elements: \code{aM} (a vector of the matching constraint for each matching variable, i.e., \code{a0}) and \code{Mdat} (a \code{data.frame} including two matching variables).
#'

mysimu <- function(N, yes.matching=FALSE, a0=NULL, m0, pi1){
  VTM <- function(vc, dm){matrix(vc, ncol=length(vc), nrow=dm, byrow=T)}
  gam.true=c(0.5,0.5)
  marker1 = rnorm(N); yi = cbind(marker1,marker2=marker1+rnorm(N))
  ti = -c(yi%*%gam.true)+log(-log(runif(N))); ti = exp(ti/2+1.5);
  ci = pmin(rgamma(N,2,2)+0.1,runif(N,0.5,2))
  xi = pmin(ti,ci); di = 1*(ti<=ci)
  id=paste("S",1:N,sep="")

  N1 = sum(di); n1 = round(N1*pi1)
  if(yes.matching){
    Zi = 1*(pnorm(yi[,1]+rnorm(N))>0.5);
    Zi = cbind(Zi,pnorm(yi[,2]+rnorm(N))*5); Zi = round(Zi)
    if(is.null(a0)){a0 = rep(0,ncol(Zi))}
    colnames(Zi)=c("Mvar1","Mvar2")
  } else {Zi = NULL}
  index.case = sort(sample((1:N)[di==1],size=n1,replace=F));
  vi.case = vi.control= rep(0,N); vi.case[index.case]=1
  matchid = matrix(NA,nrow=n1,ncol=m0+1);
  for(l in 1:length(index.case)){
    tmp.index = index.case[l]; riskset.index = xi>xi[tmp.index];
    matchid[l,1] = id[tmp.index]
    ## =========================================================================== ##
    ## if matching, additional constraint of |Zi - Zl| <= a0 needs to be satisfied ##
    ## =========================================================================== ##
    if(yes.matching){riskset.index = riskset.index & (apply(abs(Zi-VTM(Zi[tmp.index,],N))<=VTM(a0,N),1,prod)==1)}
    riskset.index = which(riskset.index==T); nl = length(riskset.index);
    # cat("case ID:",as.character(id[tmp.index]),"has",nl,"controls\n")
    ## =========================================================================== ##
    ## if riskset is empty, no control would be selected                           ##
    ## if riskset size < m.match, only select # of available for Finite Population ##
    ## =========================================================================== ##
    if (nl!=0){
      if (nl==1)  control.index = riskset.index else control.index = sample(riskset.index,min(m0,nl));
      matchid[l,(1:length(control.index))+1] = id[control.index]
      # print(matchid[control.index])
      vi.control[control.index]=1
    }
  }
  if (yes.matching) return(list(data=data.frame(time=xi,status=di,yi),id=id,case=vi.case,control=vi.control,m0=m0,matchid=matchid,Mdat=Zi,aM=a0)) else return(list(data=data.frame(time=xi,status=di,yi),id=id,case=vi.case,control=vi.control,m0=m0,matchid=matchid))
}

