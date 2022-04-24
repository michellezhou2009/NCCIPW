# NCCIPW

This R package is to implement the proposed method in the manuscript titled "New weighting methods when cases are only a subset of
events in a nested case-control study."

## Instrall `NCCIPW`

To install this R package, you need to first install the `devtool` package via
```{r}
install.packages("devtools")
```
To install the `NCCIPW` package,
```{r}
devtools::install_github("michellezhou2009/NCCIPW")
library(NCCIPW)
```

## Example

An example data `myexample` is included in the package, and it is a simulated data from the simulation setting in the manuscript. The following gives the R code to fit a time-dependent generalized linear model, and the output includes the IPW estimates using the Horvitz-Thompson's weight for the model parameters (i.e., regression coefficients) and accuracy parameters (including AUC and others) as well as their perturbed counterparts. 

```{r}
data("myexample")
GLM.IPW(formula=Surv(time,status)~marker1+marker2,
      data=myexample$data,
      id="id",
      case="case",
      control="control",
      m0=3,t0=1, 
      weight.type = "HT",      
      yes.match=T,control.matching=list(Mdat=myexample$Mdat,aM=myexample$aM),
      yes.ptb=TRUE,control.ptb=list(n.ptb=10,CaseID="CaseID"))
```

The following gives the R code to fit a Cox proportional hazards model, and the output includes the IPW estimates using the Horvitz-Thompson's weight for the model parameters (i.e., regression coefficients) and accuracy parameters (including AUC and others) as well as their perturbed counterparts. 
```{r}
PH.IPW(formula=Surv(time,status)~marker1+marker2,
      data=myexample$data,
      id="id",
      case="case",
      control="control",
      m0=3,t0=1, 
      weight.type = "HT",      
      yes.match=T,control.matching=list(Mdat=myexample$Mdat,aM=myexample$aM),
      yes.ptb=TRUE,control.ptb=list(n.ptb=10,CaseID="CaseID"))
```
