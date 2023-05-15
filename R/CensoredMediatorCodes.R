#' Title CensoredMediator
#' Description Analyze a mediation model with a censored mediator
#' @param data Input dataset file in the .csv format.
#' @param outcome.type The type of the outcome in the mediation model. outcome.type="continuous" or "binary".
#' @param casecontrol If the study is a case-control study.
#' #If outcome.type="continuous": casecontrol=FALSE
#' #If outcome.type="binary": casecontrol=FALSE if it is not a case-control study, casecontrol=TRUE if it is a case-control study
#' @param prev  Prevalance of the case in the general population.
#' #Set it as 0 when it is not a case-control study.
#' #Set it as a positive number when it is a case-control study.
#' #If outcome.type="continuous": casecontrol=FALSE; prev=0
#' #If outcome.type="binary": casecontrol=FALSE; prev=0 if it is not a case-control study, casecontrol=TRUE; prev!=0 if it is a case-control study
#' @param printout Print the estimated coefficients for different paths.
#' @import MASS
#' @import dplyr
#' @import survival
#' @import parameters
#'
#' @return list(alpha.hat,beta.opt,IE.results)
#' @export
CensoredMediator <- function(data,
                             outcome.type,
                             casecontrol,
                             prev,
                             printout)
{

  ## If it is a case-control study, the prevalence of the case in general population needs to be provided (prev).
  ## The inverse disease prevalence weight is calculated.
  if (outcome.type=="binary" & isTRUE(casecontrol) & prev!=0){
    n0<-sum(data$Y==0,na.rm=T);n1<-sum(data$Y==1,na.rm=T);fd<-prev;
    w0<-n1*(1-fd)/(n0*fd)
    data$sweight<-rep(1,dim(data)[1]);data$sweight[data$Y==0]<-w0
  }


  ## Estimation of the path from exposure (X) to mediator (M)
  newdata<-data
  my.surv<-Surv(newdata$M,1-newdata$delta)
  my.fit<-survfit(my.surv~1)
  time<-c(0,summary(my.fit)$time)
  discen<-c(1,summary(my.fit)$surv)

  data.surv<-subset(newdata, newdata$delta > 0)
  n.surv<-dim(data.surv)[1]
  kmweight<-rep(0,n.surv)
  for(kk in 1:n.surv)
  {kmweight[kk]<-min(discen[time<=data.surv$M[kk]])}

  X<-as.matrix(cbind(rep(1,n.surv),data.surv$X,data.surv[,-which(names(data.surv) %in% c('Y','X','M','delta','sweight'))]))
  ## case control study
  if (isTRUE(casecontrol)){nkmweight<-diag(data.surv$sweight*1/kmweight)}
  ## not case control study
  if (!isTRUE(casecontrol)){nkmweight<-diag(1/kmweight)}
  alpha.hat<-solve(t(X)%*%nkmweight%*%X)%*%t(X)%*%nkmweight%*%data.surv$M


  ## Estimation of the path from mediator (M) to outcome (Y)
  rr<-myBeta(alpha.hat,data,outcome.type,casecontrol)
  beta.opt<-rr[[1]];tmp<-rr[[2]]

  ## calculate the IE and PM
  IE.results<-myIE(alpha.hat,beta.opt,tmp,data,outcome.type)

  if(isTRUE(printout)){
    cat('Coefficients of M~X:\n\n')
    alpha.hat<-t(alpha.hat)
    colnames(alpha.hat)<-setdiff(names(beta.opt),'M')
    print(round(alpha.hat[1,],3))
    cat('\n\n')
    cat('Coefficients of Y~X+M:\n\n')
    print(round(beta.opt,3))
    cat('\n\n')
  }

  return(list(alpha.hat,beta.opt,IE.results))

}





#' Title myBeta
#' Description:        Calculate the coefficients for Y~X+M. Called by "CensoredMediator".
#' @param alpha Coefficients estimated for M~X.
#' @param data Input dataset.
#' @param outcome.type The type of the outcome in the mediation model. outcome.type="continuous" or "binary".
#' @param casecontrol If the study is a case-control study.
#' If outcome.type="continuous": casecontrol=FALSE
#' If outcome.type="binary":  casecontrol=FALSE if it is not a case-control study  casecontrol=TRUE if it is a case-control study
#' @return list(beta.opt,tmp)
#' @export
myBeta <- function(alpha,data,outcome.type,casecontrol)
{
  # Estimate the complete case coefficients (initial values for optimization)
  data.cc<-data[data$delta==1,]
  if(outcome.type=='binary'& isTRUE(casecontrol)){
    cc.model<-glm(Y~.-delta-sweight,family="binomial",data=data.cc)
    initial<-cc.model$coefficients}
  if(outcome.type=='binary'& !isTRUE(casecontrol)){
    cc.model<-glm(Y~.-delta,family="binomial",data=data.cc)
    initial<-cc.model$coefficients}
  if(outcome.type=="continuous"){
    cc.model<-lm(Y~.-delta,data=data.cc)
    initial<-cc.model$coefficients
    initial<-c(initial,sd(data$Y))}

  # Estimate the distribution of residuals
  if(outcome.type=='binary'){
    data$R<-data$M-as.matrix(cbind(rep(1,dim(data)[1]),data[,-which(names(data) %in% c('Y','sweight','delta','M'))]))%*%alpha
  }
  if(outcome.type=='continuous'){
    data$R<-data$M-as.matrix(cbind(rep(1,dim(data)[1]),data[,-which(names(data) %in% c('Y','delta','M'))]))%*%alpha
  }
  R.model<-survfit(Surv(data$R,data$delta)~1)
  upper<-c(1,summary(R.model)$surv[1:(length(summary(R.model)$time)-1)])
  lower<-summary(R.model)$surv
  d.t<-upper-lower
  t<-summary(R.model)$time
  tmp<-data.frame(t,d.t)

  opt.model<-optim(initial, loglikelihood, data=data, alpha=alpha, tmp=tmp, outcome.type=outcome.type, casecontrol=casecontrol)
  if(outcome.type=='binary'){beta.opt<-opt.model$par}
  if(outcome.type=='continuous'){beta.opt<-opt.model$par[1:(length(opt.model$par)-1)]}
  return(list(beta.opt,tmp))

}




#' Title loglikelihood
#' Description:        Calculate the coefficients for Y~X+M. Called by "myBeta".
#' @param x  Initial values for parameters to be optimized over.
#' @param data Input dataset.
#' @param alpha Coefficients estimated for M~X.
#' @param tmp AFT model residuals by the Kaplan-Meier estimator.
#' @param outcome.type The type of the outcome in the mediation model.  outcome.type="continuous" or "binary".
#' @param casecontrol If the study is a case-control study.
#'                                    If outcome.type="continuous": casecontrol=FALSE
#'                                    If outcome.type="binary":
#'                                              casecontrol=FALSE if it is not a case-control study
#'                                              casecontrol=TRUE if it is a case-control study
#'
#' @return negative loglikelihood
#' @export
loglikelihood <- function(x,data,alpha,tmp,outcome.type,casecontrol)
{
  n_x<-length(x)


  # The first part of the likelihood: the observed individuals
  data.delta1<-data[data$delta==1,]
  if (outcome.type=="continuous"){
    v.exp<-exp(-(data.delta1$Y-(as.matrix(cbind(rep(1,dim(data.delta1)[1]),data.delta1[,-which(names(data.delta1) %in% c('Y','delta','R'))]))%*%x[1:(n_x-1)]))^2/(2*x[n_x]^2))
    LL1<-sum(log(v.exp/sqrt(2*pi*x[n_x]^2)))
  }
  if (outcome.type=="binary"){
    v.exp<-exp(as.matrix(cbind(rep(1,dim(data.delta1)[1]),data.delta1[,-which(names(data.delta1) %in% c('Y','delta','sweight','R'))]))%*%x)
  }
  ## binary outcome and case control study
  if (outcome.type=='binary' & isTRUE(casecontrol)){
    LL.y1<-data.delta1$sweight*log(v.exp/(1+v.exp));LL.y0<-data.delta1$sweight*log(1/(1+v.exp))
    LL1<-sum(LL.y1*data.delta1$Y)+sum(LL.y0*(1-data.delta1$Y))}
  ## binary outcome but not case control study
  if (outcome.type=='binary' & !isTRUE(casecontrol)){
    LL.y1<-log(v.exp/(1+v.exp));LL.y0<-log(1/(1+v.exp))
    LL1<-sum(LL.y1*data.delta1$Y)+sum(LL.y0*(1-data.delta1$Y))}


  # The second part of the likelihood: the censored individuals
  data.delta0<-data[data$delta==0,];
  nn<-dim(data.delta0)[1];
  LL<-rep(0,nn);
  for (i in 1:nn)
  {
    if(outcome.type=="continuous"){
      index<-tmp$t>=(data.delta0$M[i]-(as.matrix(cbind(1,data.delta0[i,-which(names(data.delta0) %in% c('Y','M','delta','R'))]))%*%alpha))[,1]}
    if(outcome.type=="binary"){
      index<-tmp$t>=(data.delta0$M[i]-(as.matrix(cbind(1,data.delta0[i,-which(names(data.delta0) %in% c('Y','M','delta','sweight','R'))]))%*%alpha))[,1]}
    if(sum(index,na.rm=TRUE)!=0)
    {
      tmp.new<-tmp[index,];
      if(outcome.type=="binary"){
        M.tmp<-tmp.new$t+(as.matrix(cbind(1,data.delta0[i,-which(names(data.delta0) %in% c('Y','M','delta','sweight','R'))]))%*%alpha)[,1]
        data.tmp<-data.delta0[i,-which(names(data.delta0) %in% c('Y','delta','sweight','R'))]
        input.tmp<-cbind(rep(1,dim(tmp.new)[1]),data.tmp[rep(seq_len(nrow(data.tmp)),each=dim(tmp.new)[1]),])
        input.tmp$M<-M.tmp
        v.exp<-exp((as.matrix(input.tmp)%*%x)[,1])

        if (data.delta0$Y[i]==1){LL.y<-v.exp/(1+v.exp)}
        if (data.delta0$Y[i]==0){LL.y<-1/(1+v.exp)}
        if(isTRUE(casecontrol)){LL[i]<-data.delta0$sweight[i]*log(sum(LL.y*tmp.new$d.t))}
        if(!isTRUE(casecontrol)){LL[i]<-log(sum(LL.y*tmp.new$d.t))}
      }
      if(outcome.type=="continuous"){
        M.tmp<-tmp.new$t+(as.matrix(cbind(1,data.delta0[i,-which(names(data.delta0) %in% c('Y','M','delta','R'))]))%*%alpha)[,1]
        data.tmp<-data.delta0[i,-which(names(data.delta0) %in% c('Y','delta','R'))]
        input.tmp<-cbind(rep(1,dim(tmp.new)[1]),data.tmp[rep(seq_len(nrow(data.tmp)),each=dim(tmp.new)[1]),])
        input.tmp$M<-M.tmp
        v.exp<-exp(-(data.delta0$Y[i]-((as.matrix(input.tmp)%*%x[1:(n_x-1)])[,1]))^2/(2*x[n_x]^2))
        LL.y<-v.exp/(sqrt(2*pi*x[n_x]^2))
        LL[i]<-log(sum(LL.y*tmp.new$d.t))
      }

    }
  }
  LL0<-sum(LL)
  #LL0<-0

  -(LL0+LL1)
}




#' Title myIE
#' ## Description:        Calculate the indirect effect (IE) and percentage mediated. Called by "CensoredMediator".
#' @param alpha Coefficients estimated for M~X.
#' @param beta Coefficients estimated for Y~X+M.
#' @param tmp AFT model residuals by the Kaplan-Meier estimator.
#' @param data Input dataset.
#' @param outcome.type The type of the outcome in the mediation model. outcome.type="continuous" or "binary".
#'
#' @return data.frame(IE,DE,TE,PM)
#' @export
myIE <- function(alpha,beta,tmp,data,outcome.type) {

  mm<-dim(tmp)[1];
  P00<-P11<-P10<-rep(NA,dim(data)[1])

  for (i in 1:dim(data)[1])
  {
    if (outcome.type=="binary"){
      # X=0 for observed value and X=0 for f(m|x)
      x<-0;xx<-0
      m.value<-(as.matrix(cbind(1,xx,data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))]))%*%alpha)[,1]+tmp$t
      v.exp<-exp(as.matrix(cbind(rep(1,mm),rep(x,mm),m.value,
                                 bind_rows(replicate(1, data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))], simplify = FALSE))
      ))%*%beta)
      P00[i]<-sum((v.exp/(1+v.exp))*tmp$d.t)
      # X=1 for observed value and X=0 for f(m|x)
      x<-1;xx<-0
      m.value<-(as.matrix(cbind(1,xx,data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))]))%*%alpha)[,1]+tmp$t
      v.exp<-exp(as.matrix(cbind(rep(1,mm),rep(x,mm),m.value,
                                 bind_rows(replicate(1, data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))], simplify = FALSE))
      ))%*%beta)
      P10[i]<-sum((v.exp/(1+v.exp))*tmp$d.t)
      # X=1 for observed value and X=1 for f(m|x)
      x<-1;xx<-1
      m.value<-(as.matrix(cbind(1,xx,data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))]))%*%alpha)[,1]+tmp$t
      v.exp<-exp(as.matrix(cbind(rep(1,mm),rep(x,mm),m.value,
                                 bind_rows(replicate(1, data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))], simplify = FALSE))
      ))%*%beta)
      P11[i]<-sum((v.exp/(1+v.exp))*tmp$d.t)

    }

    if (outcome.type=="continuous"){

      # X=0 for observed value and X=0 for f(m|x)
      x<-0;xx<-0
      m.value<-(as.matrix(cbind(1,xx,data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))]))%*%alpha)[,1]+tmp$t
      v<-as.matrix(cbind(rep(1,mm),rep(x,mm),m.value,
                         bind_rows(replicate(1, data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))], simplify = FALSE))
      ))%*%beta
      P00[i]<-sum(v*tmp$d.t)
      # X=1 for observed value and X=0 for f(m|x)
      x<-1;xx<-0
      m.value<-(as.matrix(cbind(1,xx,data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))]))%*%alpha)[,1]+tmp$t
      v<-as.matrix(cbind(rep(1,mm),rep(x,mm),m.value,
                         bind_rows(replicate(1, data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))], simplify = FALSE))
      ))%*%beta
      P10[i]<-sum(v*tmp$d.t)
      # X=1 for observed value and X=1 for f(m|x)
      x<-1;xx<-1
      m.value<-(as.matrix(cbind(1,xx,data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))]))%*%alpha)[,1]+tmp$t
      v<-as.matrix(cbind(rep(1,mm),rep(x,mm),m.value,
                         bind_rows(replicate(1, data[i,-which(names(data) %in% c('Y','M','X','delta','sweight','R'))], simplify = FALSE))
      ))%*%beta
      P11[i]<-sum(v*tmp$d.t)
    }
  }

  DE<-mean(P10-P00)
  IE<-mean(P11-P10)
  TE<-mean(P11-P00)
  PM<-IE/TE

  results<-data.frame(IE,DE,TE,PM)

  return(results)

}



#' Title mybootstrap
#' Calculate the 95% confidence interval for indirect effect (IE) and percentage mediated using bootstrapping.
#'
#' @param data_hold Input dataset.
#' @param IE.star Estimated IE and PM, obtained from "CensoredMediator".
#' @param boot Number of bootstrap samples.
#' @param outcome.type The type of the outcome in the mediation model. outcome.type="continuous" or "binary".
#' @param casecontrol If the study is a case-control study.
#'                                    If outcome.type="continuous": casecontrol=FALSE
#'                                    If outcome.type="binary":
#'                                             casecontrol=FALSE if it is not a case-control study
#'                                              casecontrol=TRUE if it is a case-control study
#' @param prev Prevalance of the case in the general population.
#'                                   Set it as 0 when it is not a case-control study.
#'                                   Set it as a positive number when it is a case-control study.
#'                                   If outcome.type="continuous": casecontrol=FALSE; prev=0
#'                                   If outcome.type="binary":
#'                                             casecontrol=FALSE; prev=0 if it is not a case-control study
#'                                             casecontrol=TRUE; prev!=0 if it is a case-control study
#' @param printout Print the estimated coefficients for different paths. printout=FALSE
#'
#' @return list(CI_IE=CI_IE,CI_PM=CI_PM)
#' @export
mybootstrap <- function(data_hold,IE.star,boot,outcome.type,casecontrol,prev,printout)
{
  alpha.hat<-data.frame(matrix(0,boot,(length(names(data_hold))-2)))
  beta.hat<-data.frame(matrix(0,boot,(length(names(data_hold))-1)))
  IE.results<-data.frame(matrix(0,boot,4))
  nn<-dim(data_hold)[1]

  for (k in 1:boot)
  {
    ## generate bootstrapping data
    index<-sample(1:nn,replace=TRUE)
    data<-data_hold[index,]
    ## obtain the bootstrapping results
    results.tmp<-CensoredMediator(data,outcome.type,casecontrol,prev,printout)
    alpha.hat[k,]<-results.tmp[[1]]
    beta.hat[k,]<-results.tmp[[2]]
    IE.results[k,]<-results.tmp[[3]]
  }
  names(alpha.hat)<-row.names(results.tmp[[1]])
  names(beta.hat)<-names(results.tmp[[2]])
  names(IE.results)<-names(results.tmp[[3]])
  ## calculate the confidence intervals
  ind_star<-IE.star$IE
  ind_boot<-IE.results$IE
  ind_boot_sort<-sort(ind_boot)
  pm_star<-IE.star$PM
  pm_boot<-IE.results$PM
  pm_boot_sort<-sort(pm_boot)

  ### for indirect effect
  p<-t.test(ind_boot,mu=ind_star)$p.value
  ## calculate z0
  CDF_ind_star<-sum(1*(ind_boot<=ind_star))/boot
  if (p<=0.05) {z0<-qnorm(CDF_ind_star);a_dot<-skewness(ind_boot)[1]/6}
  if (p>0.05) {z0<-0;a_dot<-0}
  ## calculate the confidence interval
  z_low<-z0+(z0+qnorm(0.025))/(1-a_dot*(z0+qnorm(0.025)))
  z_high<-z0+(z0+qnorm(0.975))/(1-a_dot*(z0+qnorm(0.975)))
  p_low<-pnorm(as.numeric(z_low))
  p_high<-pnorm(as.numeric(z_high))
  ci_low<-ind_boot_sort[round(boot*p_low)]
  if (!is.na(round(boot*p_low))&&round(boot*p_low)==0){ci_low<-NA}
  if (is.na(round(boot*p_low))){ci_low<-NA}
  ci_high<-ind_boot_sort[round(boot*p_high)]
  if (round(boot*p_high)==0||is.na(round(boot*p_high))){ci_high<-NA}
  ## save the results
  CI_IE<-c(ci_low,ci_high)

  ### for percentage mediated
  z0<-0;a_dot<-0
  ## calculate the confidence interval
  z_low<-z0+(z0+qnorm(0.025))/(1-a_dot*(z0+qnorm(0.025)))
  z_high<-z0+(z0+qnorm(0.975))/(1-a_dot*(z0+qnorm(0.975)))
  p_low<-pnorm(z_low)
  p_high<-pnorm(z_high)
  ci_low<-pm_boot_sort[round(boot*p_low)]
  if (!is.na(round(boot*p_low))&&round(boot*p_low)==0){ci_low<-NA}
  if (is.na(round(boot*p_low))){ci_low<-NA}
  ci_high<-pm_boot_sort[round(boot*p_high)]
  if (round(boot*p_high)==0||is.na(round(boot*p_high))){ci_high<-NA}
  ## save the results
  CI_PM<-c(ci_low,ci_high)

  cat('Indirect effect (IE) and 95% CI:\n\n')
  print(paste(round(ind_star,3),' (',round(CI_IE[1],3),', ',round(CI_IE[2],3),')',sep = ""))
  cat('\n\n')
  cat('Percentage mediated (PM) and 95% CI: \n\n')
  print(paste(round(pm_star,3),' (',round(CI_PM[1],3),', ',round(CI_PM[2],3),')',sep = ""))
  cat('\n\n')

  return(list(CI_IE=CI_IE,CI_PM=CI_PM))
}


