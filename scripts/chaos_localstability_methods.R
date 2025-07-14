# This script contains several functions that are used to:
# - select the best E, tau, and theta
# - get global Lyaponov exponents (LEs), annual LEs, and local stability metrics using the Jacobian method
# - get power spectra for a time series
# - fit best smap model with all hyperparameters free or some fixed
# - fit loglinear and ricker models
# - compute memory statistics
# Adapted from Rogers et al. (2022) Nature Ecology & Evolution

#R2 function
getR2=function(resultsdf, obse="obs", prede="pred") {
  #d=na.omit(select(resultsdf, obse, prede))
  d=na.omit(select(resultsdf, all_of(obse), all_of(prede)))
  R2=1-sum((d[,obse]-d[,prede])^2)/sum((d[,obse]-mean(d[,obse]))^2)
  return(R2)
}

#Calls smap_model with several pre-specified model forms
#Note that models 1&3 and 2&5 produce identical (backtransformed) output
smap_model_options=function(data, hypars=NULL, Efix=NULL, taufix=NULL, thetafix=NULL, y, model, seaspred=F, taumax=12, Emax=6, R2trans=F, h=NULL) {
  if(model==1 | model=="ut-ut") { #"ut-ut", response and predictors both untransformed
    modelresults=smap_model(data, hypars=hypars, y=y, ylog=F, seaspred=seaspred, Efix=Efix, taufix=taufix, thetafix=thetafix, taumax=taumax, Emax=Emax, R2trans=R2trans, h=h)
  }
  if(model==2 | model=="log-log") { #"log-log", response and predictors both log transformed
    modelresults=smap_model(data, hypars=hypars, y=y, ylog=T, seaspred=seaspred, Efix=Efix, taufix=taufix, thetafix=thetafix, taumax=taumax, Emax=Emax, R2trans=R2trans, h=h)
  }
  if(model==3 | model=="fd-ut") { #"fd-ut", response = first difference, predictors untransformed
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="fd", ylog=F, seaspred=seaspred, Efix=Efix, taufix=taufix, thetafix=thetafix, taumax=taumax, Emax=Emax, R2trans=R2trans, h=h)
  }
  if(model==4 | model=="gr-ut") { #"gr-ut", response = growth rate, predictors untransformed
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="gr", ylog=F, seaspred=seaspred, Efix=Efix, taufix=taufix, thetafix=thetafix, taumax=taumax, Emax=Emax, R2trans=R2trans, h=h)
  }
  if(model==5 | model=="gr-log") { #"gr-log", response = growth rate, predictors log transformed
    modelresults=smap_model(data, hypars=hypars, y=y, pgr="gr", ylog=T, seaspred=seaspred, Efix=Efix, taufix=taufix, thetafix=thetafix, taumax=taumax, Emax=Emax, R2trans=R2trans, h=h)
  }
  return(modelresults)
}

#Get best hyperparameters (E, tau, theta)
besthyper=function(
  data=NULL, #a data frame containing the time series
  Efix=NULL, #option to fix E to a particular value
  taufix=NULL, #option to fix tau to a particular value
  thetafix=NULL, #option to fix theta to a particular value
  y, #either a time series vector, or the name of the time series variable in 'data'
  ylog, #whether or not to log the data (T/F)
  pgr="none", #whether response variable should be untransformed "none", first difference "fd", or growth rate "gr"
  seaspred=F, #include seasonal predictors
  returntable=F, #return full model selection table rather than just best model
  taumax=12, #max value of tau to try
  Emax=6, #max value of E to try
  R2trans=F, #do model selection based on transformed y rather than backtransformed y
  h=NULL #alternative steps-ahead model (defaults to h=tau)
  ) {

  #if using first difference as response, ylog must be set to F
  if(ylog==T & pgr=="fd") stop("must use ylog=F for pgr='fd'")

  #if providing data frame containing y, give y as character (column name)
  if(!is.null(data)) {
    ser_or=data[,y]
  } else {
    ser_or=y
  }

  #log transform
  if(ylog==T) ser=log(ser_or) else ser=ser_or
  if(pgr=="gr") ser_log=log(ser_or)

  #get effective ts length (longest string of non-NA values)
  len=rle(!is.na(ser))
  serlen=max(len$lengths[len$values==TRUE])

  #range of E, tau, theta to try
  #max is currently hard coded, but could be changed/added as an argument
  if(!is.null(taufix)) {
    tautry=taufix #set tau to fixed value
  } else {
    tautry=1:taumax
  }
  if(!is.null(Efix)) {
    Etry=Efix #set E to fixed value
  } else {
    Etry=1:Emax
  }

  if(!is.null(thetafix)) {
    thetatest = thetafix #set theta to fixed value
  } else {
    thetatest = c(0, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8)
  }

  #create table of candidate models
  model_selection=expand.grid(E=Etry, tau=tautry, theta=thetatest)
  #determine sort order
  model_selection=arrange(model_selection, tau, theta, E)

  #filter candidate models to reasonable E and tau combinations given effective ts length
  if(length(Etry)!=1 & length(tautry)!=1) {
    #if E or tau are both free
    model_selection=filter(model_selection, E^2<=serlen & E*tau/serlen<=0.2)
  } else if(length(Etry)!=1 | length(tautry)!=1) {
    #if one of E or tau are free, the other fixed
    model_selection=filter(model_selection, E^2<=serlen & E*tau/serlen<=0.2)
          #above line changed so that model sets are comparable when E is fixed to 1
  } else {
    #if both fixed, make sure E isn't too large
    while(unique(model_selection$E)^2>serlen) {
        model_selection$E=model_selection$E-1
    }
  }
  model_selection$rmse=NA
  model_selection$npred=NA
  model_selection$error1=NA

  #generate lags
  ser_lags=make_block(ser, max_lag = max(model_selection$E*model_selection$tau)+1+ifelse(!is.null(h),h,0), tau=1)[,-1]
  #ser_lags=na.omit(ser_lags) #if you wish to standardize targets across models

  #generate seasonal predictors (if relevant)
  if(seaspred) {
    ser_lags$S1=-cos(data$Month*pi/6)/0.7385489
    ser_lags$S2=-cos((data$Month-3)*pi/6)/0.7385489
  }

  #run s-map for all parameter sets
  for(i in 1:nrow(model_selection)) {
    #set h if it does not exist
    mstau=model_selection$tau[i]
    msE=model_selection$E[i]
    if(is.null(h)) {
      msh=mstau
    } else {
      msh=h
    }
    #if using first difference or growth rate as response, substitute it
    if(pgr=="fd") {
      ser_lags$col1=ser-lag(ser,msh)
    }
    if(pgr=="gr") {
      ser_lags$col1=ser_log-lag(ser_log,msh)
    }
    if(seaspred) {
      simptemp=block_lnlp(ser_lags, tp = 0, method = "s-map",
                          columns = c(paste0("col1_",mstau*(1:msE)-mstau+msh,"S1","S2")),
                          target_column = "col1", theta=model_selection$theta[i],
                          first_column_time = F, silent = T, stats_only = F)
    } else {
      simptemp=block_lnlp(ser_lags, tp = 0, method = "s-map",
                          columns = paste0("col1_",mstau*(1:msE)-mstau+msh),
                          target_column = "col1", theta=model_selection$theta[i],
                          first_column_time = F, silent = T, stats_only = F)
    }
    model_selection$rmse[i]=simptemp$rmse
    model_selection$npred[i]=simptemp$num_pred
    resultsdf=simptemp$model_output[[1]]
    resultsdf$Obs_abund=ser_or
    #get error based on abundance
    #requires backtransformation if using first difference or growth rate
    if(pgr=="none" & ylog==F) {
      resultsdf$Pred_abund=resultsdf$pred
    }
    if(pgr=="none" & ylog==T) {
      resultsdf$Pred_abund=exp(resultsdf$pred)
    }
    if(pgr=="fd") {
      resultsdf$Pred_abund=lag(ser_or,msh)+resultsdf$pred
    }
    if(pgr=="gr") {
      resultsdf$Pred_abund=lag(ser_or,msh)*exp(resultsdf$pred)
    }
    model_selection$error1[i]=-getR2(resultsdf, obse="Obs_abund", prede="Pred_abund")
    if(R2trans) { #override R2 ranking to use transformed y
      model_selection$error1[i]=-getR2(resultsdf)
    }
  }

  #calculate error accounting for df
  model_selection$error=round(model_selection$error1,2)
  bestE=model_selection$E[which.min(model_selection$error)]
  bestTau=model_selection$tau[which.min(model_selection$error)]
  bestTheta=model_selection$theta[which.min(model_selection$error)]

  #output
  if(returntable) {
    return(model_selection)
  } else {
    return(data.frame(bestE=bestE, bestTau=bestTau, bestTheta=bestTheta, Emax=max(model_selection$E), taumax=max(model_selection$tau)))
  }
}

#Get s-map model output
#This uses the above function to get the best hyperparameters.
smap_model=function(
  data=NULL, #a data frame containing the time series
  hypars=NULL, #optional output of 'besthyper' containing pre-optimized hyperparameters
  y, #either a time series vector, or the name of the time series variable in 'data'
  ylog, #whether or not to log the data (T/F)
  pgr="none", #whether response variable should be untransformed "none", first difference "fd", or growth rate "gr"
  seaspred=F, #include seasonal predictors
  Efix=NULL, #option to fix tau to a particular value
  taufix=NULL, #option to fix tau to a particular value
  thetafix=NULL, #option to fix theta to a particular value
  taumax=12, #max value of tau to try
  Emax=6, #max value of E to try
  R2trans=F, #do model selection based on transformed y rather than backtransformed y
  h=NULL #alternative steps-ahead model
  ) {

  #get hyper parameters if not already supplied
  #obviously, if you fix E and tau, this will return only a single model
  if(!is.null(hypars)) {
    hyperpars=hypars
  } else {
    hyperpars=besthyper(data=data, y=y, ylog=ylog, pgr=pgr, seaspred=seaspred, Efix=Efix, taufix=taufix, thetafix=thetafix, taumax=taumax, Emax=Emax, R2trans=R2trans, h=h)
  }

  #if providing data frame containing y, give y as character (column name)
  if(!is.null(data)) {
    ser_or=data[,y]
  } else {
    ser_or=y
  }

  #log transform
  if(ylog==T) ser=log(ser_or) else ser=ser_or
  if(pgr=="gr") ser_log=log(ser_or)

  if(is.null(h)) {
    h=hyperpars$bestTau
  }

  #store smap output
  ser_lags=make_block(ser, max_lag = hyperpars$bestE*hyperpars$bestTau+1-hyperpars$bestTau+h, tau=1)[,-1]
  if(pgr=="fd") {
    ser_lags$col1=ser-lag(ser,h)
  }
  if(pgr=="gr") {
    ser_lags$col1=ser_log-lag(ser_log,h)
  }
  if(seaspred) {
    ser_lags$S1=-cos(data$Month*pi/6)/0.7385489
    ser_lags$S2=-cos((data$Month-3)*pi/6)/0.7385489
  }
  if(seaspred) {
    smap_results=block_lnlp(ser_lags, tp = 0, method = "s-map",
                            columns = c(paste0("col1_",hyperpars$bestTau*(1:hyperpars$bestE)-hyperpars$bestTau+h),"S1","S2"),
                            target_column = "col1", theta=hyperpars$bestTheta,
                            first_column_time = F, silent = T,
                            stats_only = F, save_smap_coefficients = T)
  } else {
    smap_results=block_lnlp(ser_lags, tp = 0, method = "s-map",
                            columns = paste0("col1_",hyperpars$bestTau*(1:hyperpars$bestE)-hyperpars$bestTau+h),
                            target_column = "col1", theta=hyperpars$bestTheta,
                            first_column_time = F, silent = T,
                            stats_only = F, save_smap_coefficients = T)

  }
  resultsdf=cbind(smap_results$model_output[[1]], smap_results$smap_coefficients[[1]])

  #untransformed abundance
  resultsdf$Obs_abund=ser_or
  #get untranformed abundance predictions based on model type
  if(pgr=="none" & ylog==F) {
    resultsdf$Pred_abund=resultsdf$pred
    form="ut-ut"
  }
  if(pgr=="none" & ylog==T) {
    resultsdf$Pred_abund=exp(resultsdf$pred)
    form="log-log"
  }
  if(pgr=="fd") {
    resultsdf$Pred_abund=lag(resultsdf$Obs_abund,h)+resultsdf$pred
    form="fd-ut"
  }
  if(pgr=="gr") {
    resultsdf$Pred_abund=lag(resultsdf$Obs_abund,h)*exp(resultsdf$pred)
    #ylog is not used to get predictions, but used for jacobian construction
    if(ylog==T) form="gr-log"
    if(ylog==F) form="gr-ut"
  }

  #R2 for model
  modelR2=getR2(resultsdf)
  #R2 for untransformed abundance
  modelR2_abund=getR2(resultsdf, obse="Obs_abund", prede="Pred_abund")
  modelstats=data.frame(E=hyperpars$bestE, tau=hyperpars$bestTau, theta=hyperpars$bestTheta, h=h,
                        Emax=hyperpars$Emax, taumax=hyperpars$taumax, num_pred=smap_results$num_pred,
                        R2model=modelR2, R2abund=modelR2_abund, rho=smap_results$rho)
  return(list(modelstats=modelstats, resultsdf=resultsdf, form=form))
}

#Construct Jacobian matrices from smap coefficients
getJacobians=function(
  modelresults #output from smap_model or smap_model_options
  ){

  #jacobians are derived in terms of untransformed abudance, so will vary depending on the model and data transform.
  #the LE is (in theory) not affected by the particular form of the model or Jacobian, but this ensures that the local LEs are correct.
  form=modelresults$form

  ndim=modelresults$modelstats$E #dimension
  tau=modelresults$modelstats$tau #tau
  len=nrow(modelresults$resultsdf)
  coefs=modelresults$resultsdf[,paste0("c_",1:ndim)]

  #jacobians are stored in a 3d array
  jacobians=array(dim = c(ndim,ndim,len))

  if(form=="ut-ut") {
    if(ndim==1) {
      jacobians[,,]=coefs
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1=matrix(nrow = ndim, ncol = ndim, 0)
          J1[1,]=as.numeric(coefs[i,])
          J1[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1
        }
      }
    }
  }

  if(form=="fd-ut") {
    if(ndim==1) {
      jacobians[,,]=coefs+1
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1=matrix(nrow = ndim, ncol = ndim, 0)
          J1[1,]=as.numeric(coefs[i,])+c(1,rep(0,ndim-1))
          J1[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1
        }
      }
    }
  }

  if(form=="log-log") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs)*r_x/lag(x_obs, tau)
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-J3<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i], x_obs[i-(1:(ndim-1))*tau])
          #c(r_x[i], lag(x_obs)[i], lag(x_obs,2)[i], lag(x_obs,3)[i], lag(x_obs,4)[i])
          J2[1,]=as.numeric(coefs[i,])
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          diag(J3)<-1/x_obs[i-(1:(ndim))*tau]
          #c(1/(lag(x_obs)[i]), 1/(lag(x_obs,2)[i]), 1/(lag(x_obs,3)[i]), 1/(lag(x_obs,4)[i]), 1/(lag(x_obs,5)[i]))
          jacobians[,,i]=J1%*%J2%*%J3
        }
      }
    }
  }

  if(form=="gr-ut") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs*lag(x_obs,tau)+1)*r_x
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i], rep(1,ndim-1))
          #c(r_x[i], 1, 1, 1)
          J2[1,]=as.numeric(coefs[i,])*x_obs[i-tau]+c(1,rep(0,ndim-1))
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          jacobians[,,i]=J1%*%J2
        }
      }
    }
  }

  if(form=="gr-log") {
    x_obs=modelresults$resultsdf$Obs_abund
    r_x=exp(modelresults$resultsdf$pred)
    if(ndim==1) {
      jacobians[,,]=(coefs+1)*r_x
    } else {
      for(i in 1:len) {
        if(all(!is.na(coefs[i,]))) {
          J1<-J2<-J3<-matrix(nrow = ndim, ncol = ndim, 0)
          diag(J1)<-c(r_x[i]*x_obs[i-tau], x_obs[i-(1:(ndim-1))*tau])
          #c(r_x[i]*lag(x_obs)[i], lag(x_obs)[i], lag(x_obs,2)[i], lag(x_obs,3)[i], lag(x_obs,4)[i])
          J2[1,]=as.numeric(coefs[i,])+c(1,rep(0,ndim-1))
          J2[2:ndim,1:(ndim-1)]=diag(ndim-1)
          diag(J3)<-1/x_obs[i-(1:(ndim))*tau]
          #c(1/(lag(x_obs)[i]), 1/(lag(x_obs,2)[i]), 1/(lag(x_obs,3)[i]), 1/(lag(x_obs,4)[i]), 1/(lag(x_obs,5)[i]))
          jacobians[,,i]=J1%*%J2%*%J3
        }
      }
    }
  }

  return(jacobians)
}

#Get LE (lower confidence bound) from jacobian matrices
#This is done by averaging LEs computed over several long subsegments of the data
LEshift=function(
  modelresults, #output from smap_model or smap_model_options
  jacobians #output from getJacobians
  ) {

  #modelresults is just used to get tau
  tau=modelresults$modelstats$tau

  len=dim(jacobians)[3] #time series length
  ndim=dim(jacobians)[1] #E

  #remove leading NAs
  if(ndim==1) {
    jacobians2=jacobians[(ndim*tau+1):len]
    runlen=rle(!is.na(jacobians2))
    len2=length(jacobians2)
  } else {
    jacobians2=jacobians[,,(ndim*tau+1):len]
    runlen=rle(!is.na(jacobians2[1,1,]))
    len2=dim(jacobians2)[3]
  }

  #effective ts length (longest string of non-NA values)
  serlen=max(runlen$lengths[runlen$values==TRUE])

  Tminus=3:6 #segment lengths range in length from T-3 to T-6
  LEseg=data.frame(SegLen=(serlen-max(Tminus)):(serlen-min(Tminus)), le_mean=NA, le_sd=NA, le_ci=NA, le_n=NA) %>%
    filter(SegLen>0)

  for(i in 1:nrow(LEseg)) {
    SegLen=LEseg$SegLen[i]
    LEtemp=numeric(len2-SegLen+1)
    for(j in 1:length(LEtemp)) {
      if(ndim==1) {
        Jacs1=jacobians2[j:(j+SegLen-1)]
      } else {
        Jacs1=jacobians2[,,j:(j+SegLen-1)]
      }
      LEtemp2=numeric(tau)
      if(any(is.na(Jacs1))) {
        LEtemp2=NA
      } else {
        for(a in 1:tau) {
          indices=seq(from=a, to=ifelse(ndim==1,length(Jacs1),dim(Jacs1)[3]), by=tau)
          if(ndim==1) {
            Jacs=Jacs1[indices]
            LEtemp2[a]=mean(log(abs(Jacs)), na.rm=T)/tau
          } else {
            Jacs=Jacs1[,,indices]
            nk=dim(Jacs)[3]
            Jk=Jacs[,,1]
            QR=qr(Jk)
            R=qr.R(QR)
            Q=qr.Q(QR)
            Rcum=R
            for(k in 2:nk) {
              Jk=Jacs[,,k]
              QR=qr(Jk%*%Q)
              R=qr.R(QR)
              Q=qr.Q(QR)
              Rcum=R%*%Rcum
            }
            LEtemp2[a]=1/nk*log(max(abs(diag(Rcum))))/tau
          }
        }
      }
      LEtemp[j]=mean(LEtemp2, na.rm=F) #if tau>1, will result in more than i+1 segments being averaged unless you set na.rm=F
    }
    LEseg$le_n[i]=length(which(!is.na(LEtemp)))
    LEseg$le_mean[i]=mean(LEtemp, na.rm=T)
    LEseg$le_sd[i]=sd(LEtemp, na.rm=T)
    LEseg$le_ci[i]=LEseg$le_sd[i]/sqrt(LEseg$le_n[i])*qt(p=0.95, df=LEseg$le_n[i]-1)
  }

  minmean=min(LEseg$le_mean) #this is not necessarily the mean with lowest CI
  minci=min(LEseg$le_mean-LEseg$le_ci) #lowest of the lower confidence bounds on LE

  return(list(LEseg=LEseg, minmean=minmean, minci=minci))
}

#computes local eigenvalue and det
LElocal = function(modelresults, jacobians) {

  #modelresults is just used to get tau
  tau=modelresults$modelstats$tau

  len=dim(jacobians)[3] #time series length
  ndim=dim(jacobians)[1] #E

  lle=numeric(length = len)
  det=numeric(length = len)
  # lle3=numeric(length = len)

  if(ndim==1) {
    lle=log(abs(drop(jacobians)))/tau
    det[]=lle
  }
  else {
    for(i in 1:len) {
      Jac=jacobians[,,i]
      if(any(is.na(Jac))) {
        lle[i]=NA
        det[i]=NA
      } else {
        lle[i]=log(max(abs(eigen(Jac, only.values = T)$values)))/tau
        det[i]=log(abs(det(Jac)))/tau

        #lle[i]=log(max(svd(Jac)$d))
        #alleigen[i,]=eigen(Jac, only.values = T)$values

        # if(i>(len-2)) {
        #   lle3[i]=NA
        # } else {
        #   Jac1=jacobians[,,i+1]
        #   Jac2=jacobians[,,i+2]
        #   if(any(is.na(Jac1)) | any(is.na(Jac2))) {
        #     lle3[i]=NA
        #   } else {
        #     svd[i,]=log(svd(Jac2%*%Jac1%*%Jac)$d)/3
        #     lle3[i]=max(svd[i,])
        #   }
        # }

      }
    }
  }
  stability=data.frame(lle=lle,det=det)
  return(stability)
}

#LE by calendar year
LEannual=function(data, modelresults, jacobians) {

  #modelresults is just used to get tau
  tau=modelresults$modelstats$tau

  ndim=dim(jacobians)[1] #E

  Year=data$Year
  years=unique(Year)

  # LE=numeric(length = length(years))
  # for(i in 1:length(years)) {
  #   Jacs=jacobians[,,which(Year==years[i])]
  #   if(any(is.na(Jacs))) {
  #     LE[i]=NA
  #   } else {
  #     nk=dim(Jacs)[3]
  #     Jk=Jacs[,,1]
  #     QR=qr(Jk)
  #     R=qr.R(QR)
  #     Q=qr.Q(QR)
  #     Rcum=R
  #     for(k in 2:nk) {
  #       Jk=Jacs[,,k]
  #       QR=qr(Jk%*%Q)
  #       R=qr.R(QR)
  #       Q=qr.Q(QR)
  #       Rcum=R%*%Rcum
  #     }
  #     LE[i]=1/nk*log(max(abs(diag(Rcum))))
  #   }
  # }
  # return(data.frame(Year=years, AnnualLE=LE))


  LE=numeric(length = length(years))
  for(i in 1:length(years)) {
    if(ndim==1) {
      Jacs1=jacobians[which(Year==years[i])]
    } else {
      Jacs1=jacobians[,,which(Year==years[i])]
    }
    njacs=length(which(Year==years[i]))
    LEtemp2=numeric(tau)
    if(any(is.na(Jacs1)) | njacs<12 | tau>6) {
      LEtemp2=NA
    } else {
      for(a in 1:tau) {
        indices=seq(from=a, to=njacs, by=tau)
        if(ndim==1) {
          Jacs=Jacs1[indices]
          LEtemp2[a]=mean(log(abs(Jacs)), na.rm=T)/tau
        } else {
          Jacs=Jacs1[,,indices]
          nk=dim(Jacs)[3]
          Jk=Jacs[,,1]
          QR=qr(Jk)
          R=qr.R(QR)
          Q=qr.Q(QR)
          Rcum=R
          for(k in 2:nk) {
            Jk=Jacs[,,k]
            QR=qr(Jk%*%Q)
            R=qr.R(QR)
            Q=qr.Q(QR)
            Rcum=R%*%Rcum
          }
          LEtemp2[a]=1/nk*log(max(abs(diag(Rcum))))/tau
        }
      }
    }
    LE[i]=mean(LEtemp2, na.rm=F) #if tau>1, will result in more than i+1 segments being averaged unless you set na.rm=F
  }
  return(data.frame(Year=years, AnnualLE=LE))

}

#Power spectrum
getspectrum=function(dataframe, ts, lam=0.01, scale=T, trim=T) {

  xts=dataframe[,ts]

  #approximate spectrum calculated using ridge regression
  if(all(is.na(xts))) {
    return(data.frame(Frequency=NA, Period=NA, Power=NA))
  }

  #standardize data
  if(scale==T) {
    xtss=(xts-mean(xts, na.rm=T))/sd(xts, na.rm=T)
  } else {
    xtss=xts-mean(xts, na.rm=T)
  }
  #acf used in earlier version, no longer used
  #zooacf=acf(xtss, na.action = na.pass, lag.max = length(xtss), plot = F)$acf[,1,1] #compute acf, max lag is ts length
  #zooacf=zooacf[!is.na(zooacf)] #exclude NAs at end to get actual ts length

  #get spectrum
  if(trim) {  # ts length, exclude NAs on ends
    tsl=length(which(!is.na(xtss))[1]:which(!is.na(xtss))[length(which(!is.na(xtss)))])
  } else { # do not trim NAs on ends
    tsl=length(xtss)
  }
  tsle=floor(tsl/2)*2 #ts length rounded down to even number (if odd)
  fi=(1:(tsle/2))/tsle #frequencies (cycles per timestep)
  per=1/fi #periods (timesteps per cycle)
  wi=2*pi*fi #frequencies in radians per timestep

  times=1:length(xtss)
  cosbf=cos(outer(times,wi))
  sinbf=sin(outer(times,wi))
  allbf=cbind(cosbf,sinbf) #all basis functions

  y=xtss[complete.cases(xtss)] #remove missing timepoints
  X=allbf[complete.cases(xtss),] #remove missing timepoints
  coefsdir=solve(t(X)%*%X + lam*diag(ncol(X)))%*%t(X)%*%y
  lmr=sqrt(coefsdir[1:(length(coefsdir)/2),]^2+coefsdir[(length(coefsdir)/2+1):length(coefsdir),]^2)

  return(data.frame(Frequency=fi, Period=per, Power=lmr))
}

getspectrum12=function(dataframe, ts, scale=T, phase=F) {
  #just for period 12

  xts=dataframe[,ts]

  if(all(is.na(xts))) {
    return(data.frame(cos=NA, sin=NA, Power=NA))
  }
  #standardize data
  if(scale==T) {
    xtss=(xts-mean(xts, na.rm=T))/sd(xts, na.rm=T)
  } else {
    xtss=xts-mean(xts, na.rm=T)
  }

  fi=1/12 #frequencies (cycles per timestep)
  wi=2*pi*fi #frequencies in radians per timestep

  times=1:length(xtss)
  cosbf=cos(outer(times,wi))
  sinbf=sin(outer(times,wi))
  allbf=cbind(cosbf,sinbf)

  coefs=lm(xtss~allbf)$coefficients[2:3]
  names(coefs)=NULL
  lmr=sqrt(coefs[1]^2+coefs[2]^2)
  #ph=atan(-coefs[2]/coefs[1])

  ms=1:12
  y=coefs[1]*cos(pi*ms/6)+coefs[2]*sin(pi*ms/6)
  maxm=ms[which.max(y)]

  #return(data.frame(cos=coefs[1], sin=coefs[2], Power=lmr))
  if(phase) {
    return(maxm)
  } else {
    return(lmr)
  }
}

#best fitting Smap model
bestsmap=function(data, taumax=12, Efix=NULL, taufix=NULL, y, h=NULL) {
  modelresults3=smap_model_options(data=data, y=y, model=3, Efix=Efix, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  modelresults4=smap_model_options(data=data, y=y, model=4, Efix=Efix, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  modelresults5=smap_model_options(data=data, y=y, model=5, Efix=Efix, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  R3a=modelresults3$modelstats$R2abund
  R4a=modelresults4$modelstats$R2abund
  R5a=modelresults5$modelstats$R2abund
  bestmodel=which.max(c(R3a,R4a,R5a))
  modelresults=list(modelresults3,modelresults4,modelresults5)[[bestmodel]]
  return(modelresults)
}

#best fitting 1-d model
best1dsmap=function(data, taumax=12, taufix=NULL, y, h=NULL) {
  modelresults3=smap_model_options(data=data, y=y, model=3, Efix=1, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  modelresults4=smap_model_options(data=data, y=y, model=4, Efix=1, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  modelresults5=smap_model_options(data=data, y=y, model=5, Efix=1, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  R3a=modelresults3$modelstats$R2abund
  R4a=modelresults4$modelstats$R2abund
  R5a=modelresults5$modelstats$R2abund
  bestmodel=which.max(c(R3a,R4a,R5a))
  modelresults=list(modelresults3,modelresults4,modelresults5)[[bestmodel]]
  return(modelresults)
}

#best fitting AR model
bestAR=function(data, taumax=12, taufix=NULL, y, h=NULL) {
  modelresults3=smap_model_options(data=data, y=y, model=3, thetafix=0, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  modelresults4=smap_model_options(data=data, y=y, model=4, thetafix=0, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  modelresults5=smap_model_options(data=data, y=y, model=5, thetafix=0, taumax=taumax, taufix=taufix, seaspred=F, h=h)
  R3a=modelresults3$modelstats$R2abund
  R4a=modelresults4$modelstats$R2abund
  R5a=modelresults5$modelstats$R2abund
  bestmodel=which.max(c(R3a,R4a,R5a))
  modelresults=list(modelresults3,modelresults4,modelresults5)[[bestmodel]]
  return(modelresults)
}

#log linear
fit_loglinear=function(data, var, exrad=0) {
  x=data[,var]
  lnx=log(x); lnx_1=lag(lnx)
  df=data.frame(lnx=lnx, lnx_1=lnx_1)
  # m=lm(lnx~lnx_1, data = df)
  # df$insamp=NA
  # df$insamp[!is.na(lnx) & !is.na(lnx_1)]=predict(m)
  df$loo=NA
  for(i in 1:length(lnx)) {
    if(!is.na(lnx[i]) & !is.na(lnx_1[i])) {
      exclude=max(1,i-exrad):min(nrow(df),i+exrad)
      dfi=df[-exclude,]
      mi=lm(lnx~lnx_1, data=dfi)
      df$loo[i]=predict(mi,newdata=df[i,])
    }
  }
  out=data.frame(obs=x,pred=exp(df$loo))
  return(out)
}

#ricker
fit_ricker=function(data, var, exrad=0) {
  x=data[,var]
  x_1=lag(x); gr=log(x/lag(x))
  df=data.frame(x_1=x_1, gr=gr)
  # m=lm(gr~x_1, data = df)
  # df$insamp=NA
  # df$insamp[!is.na(gr) & !is.na(x_1)]=predict(m)
  df$loo=NA
  for(i in 1:length(gr)) {
    if(!is.na(!is.na(gr[i]) & !is.na(x_1[i]))) {
      exclude=max(1,i-exrad):min(nrow(df),i+exrad)
      dfi=df[-exclude,]
      mi=lm(gr~x_1, data=dfi)
      df$loo[i]=predict(mi,newdata=df[i,])
    }
  }
  out=data.frame(obs=x,pred=x_1*exp(df$loo))
  return(out)
}

memorystats=function(jacobians) {
  if(dim(jacobians)[1]==1) {
    return(data.frame(meanlag=1,varlag=0,entropylag=0,edimlag=1))
  }
  firstrow=t(jacobians[1,,])
  pis=apply(firstrow, 1, function(x) abs(x)/sum(abs(x)))
  meanlag=mean(apply(pis,2,function(x) sum((1:length(x))*x)),na.rm=T)
  varlag=mean(apply(pis,2,function(x) {
    Es=1:length(x)
    m=sum(Es*x)
    v=sum((Es-m)^2*x)
    return(v)
    }),na.rm=T)
  entropylag=mean(apply(pis,2,function(x) -sum(x*log(x))),na.rm=T)
  edimlag=mean(apply(pis,2,function(x) 1/sum(x^2)),na.rm=T)
  return(data.frame(meanlag=meanlag,varlag=varlag,entropylag=entropylag,edimlag=edimlag))
}
