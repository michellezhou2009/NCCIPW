# NCCIPW

This R package is to implement the proposed method in the manuscript titled ["A new weighting method when not all the events are selected as cases in a nested case-control study"](https://arxiv.org/abs/2104.02665).

## Instrall `NCCIPW`

To install this R package, you need to first install the `devtool` package via
```{r}
install.packages("devtools")
```
To install the `NCCIPW` package,
```{r}
library(devtools)
install_github("michellezhou2009/NCCIPW")
library(NCCIPW)
```

## Example

An example data `myexample` is included in the package, and it is a simulated data from the simulation setting included in the manuscript. The following gives the R code to fit a time-dependent generalized linear model, and the output includes the new IPW estimates of the model parameters (i.e., regression coefficients) and accuracy parameters (including AUC and others) as well as their perturbed counterparts. 

```{r}
data("myexample")
GLM.IPW(formula=Surv(time,status)~marker1+marker2,
      data=myexample$data, family=binomial(link=logit),
      id=myexample$id,
      case=myexample$case,
      control=myexample$control,
      m0=myexample$m0,t0=1,
      yes.matching=T,control.matching=list(Mdat=myexample$Mdat,aM=myexample$aM),
      yes.ptb=TRUE,control.ptb=list(n.ptb=100,matchid=myexample$matchid))
```

The following gives the R code to fit a Cox proportional hazards model, and the output includes the new IPW estimates of the model parameters (i.e., regression coefficients) and accuracy parameters (including AUC and others) as well as their perturbed counterparts. 
```{r}
PH.IPW(formula=Surv(time,status)~marker1+marker2,
      data=myexample$data,
      id=myexample$id,
      case=myexample$case,
      control=myexample$control,
      m0=myexample$m0,t0=1,
      yes.matching=T,control.matching=list(Mdat=myexample$Mdat,aM=myexample$aM),
      yes.ptb=TRUE,control.ptb=list(n.ptb=100,matchid=myexample$matchid))
```
