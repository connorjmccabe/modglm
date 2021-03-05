
margplot<-function(model, vars,foc, mod, modlevels, data, modnames=NULL,hyps="means",sm=T){

  if(!all(between(modlevels,min(data[,mod]),max(data[,mod])))){warning("Note: at least one hypothetical moderator level in 'modlevels' is outside the range of your moderator variable. Please check your plotted moderator values.")}

  (int.var <- paste(vars, collapse = ":"))

  (b <- model$coef)

  if(model$call[1]=="gee()"){
    dftemp<-na.omit(data[,which(names(data) %in% names(model$coefficients))])
    X<-as.data.frame(cbind(rep(1,nrow(dftemp)),dftemp))
  }else(X<-as.data.frame(cbind(rep(1,nrow(model$model)),model$model[,-1])))

  colnames(X)[1]<-"(Intercept)"

  if (int.var %in% names(model$coefficients)){X[,int.var]<-X[,vars[1]]*X[,vars[2]]}

  X<-as.matrix(X)

  #Step 1: extract variance-covariance and point estimates that we will use later in our simulation.
  vcov<-vcov(model)
  pes<-coef(model)

  #step 2: set up counterfactual (i.e. hypothetical) scenarios.

  #2a: define a hypothetical range for your focal (x-axis) variable.
  #Below, I'm defining an object representing a sequence of numbers for the focal variable (focal.seq)
  #ranging from the minimum to maximum values of the focal predictor, in increments of .1.
  #I've also included 1 SD below the mean, the mean, and 1 SD above the mean as default values.

  focal.seq <- sort(c(-1*sd(data[,foc],na.rm=TRUE),
                      mean(data[,foc],na.rm=TRUE),
                      1*sd(data[,foc],na.rm=TRUE),
                      seq(min(data[,foc], na.rm = TRUE),
                          max(data[,foc], na.rm = TRUE),
                          .1*sd(data[,foc],na.rm=TRUE))))

  #2b: set up an empty data frame for storing hypotheticals.
  xim<-as.data.frame(matrix(NA,nrow=length(focal.seq),ncol=length(pes)))
  colnames(xim)<-colnames(X)
  for(i in 1:ncol(xim)){
    if(hyps[1]=="means")(xim[,i]<-mean(X[,i],na.rm=T))
    else(xim[,i]<-(xim[,i]<-hyps[i]))
    # (xim[,i]<-hyps[i]))
  }

  xim.modlvs<<-list()
  for(i in 1:length(modlevels)){
    xim.modlvs[[i]]<-xim
    xim.modlvs[[i]][,foc]<-focal.seq
    xim.modlvs[[i]][,mod]<-modlevels[i]
    xim.modlvs[[i]][,int.var]<-xim.modlvs[[i]][,foc]*xim.modlvs[[i]][,mod]
  }
  library(data.table)
  ximall<-rbindlist(xim.modlvs)

  #Step 4: SIMULATION.
  #define number of simulations. We will put this as 10,000 sims.
  sims <- 10000
  require(MASS) #needed for the mvrnorm() function.

  #4a: simulate parameters. The basis for why this works is that we assume asymtotic multivariate normality of
  #the parameter estimates when estimating our GLM. So, we can simulate these from the estimated
  #vcov and point estimates to generate approximate confidence regions.
  simparams<-mvrnorm(sims,pes,vcov)

  #Now, I have a set of 10000 coefficients with which to matrix multiply my hypothetical scenarios generated above.
  #Because we're doing matrix algebra, i'll need to convert these to a data matrix:

  matxi<-data.matrix(ximall)


  #IWe next multiply each ROW of my hypothetical matrix by these parameters, then transform all of these into the metric of interest and sort the results.

  #set up objects for storing results:
  lower<-pe<-upper<-rep(NA,nrow(matxi))
  resm.model<-list()

  #Loop through each ROW of our hypothetical scenarios and compute lower, median, and upper limits for our point estimates.
  for(i in 1:nrow(matxi)){
    simmu <- simparams%*%matxi[i,]

    if(model$family$link == "logit"){simy <- 1/(1 + exp(-simmu))}
    else if (model$family$link == "log"){simy<-exp(simmu)}
    else if (model$family$link == "identity"){simy <- simmu}

    simy<-sort(simy)


    length.simy <- length(simy)
    low <- up <- NULL

    low <- c( low,simy[trunc((1-.95)/2*length.simy)] )
    up  <- c( up, simy[trunc((1-(1-.95)/2)*length.simy)])

    resm.model$lower <- rbind(resm.model$low,low)
    resm.model$pe <- c(resm.model$pe,mean(simy))
    resm.model$upper <- rbind(resm.model$up,up)
  }

  #Put into a dataframe for plotting in ggplot:
  dfplot.res<-as.data.frame(cbind(matxi[,foc],matxi[,mod],resm.model$pe,resm.model$lower,resm.model$upper),row.names=FALSE)
  colnames(dfplot.res)<-c(foc,mod,"pe","lower","upper")

  if(!is.null(modnames)){for(i in 1:length(modlevels)){dfplot.res[,mod][dfplot.res[,mod]==modlevels[i]]<-modnames[i]}}

  #Now, we can plot the results:
  require(ggplot2)

  (margplot<-ggplot(data=dfplot.res,aes(x=x1,y=pe,group=as.factor(dfplot.res[,mod]),fill=as.factor(dfplot.res[,mod]))) +
      geom_line(aes(color=as.factor(dfplot.res[,mod]))) +
      geom_ribbon(aes(x=x1, ymin = lower, ymax = upper), alpha=.5) +
      scale_fill_brewer(palette = "Pastel1")+scale_color_brewer(palette = "Set1")+
      # ylim(0,1) +
      ylab("DV") +
      # facet_grid(~dfplot.res[,mod]) +
      theme_bw() +
      theme(panel.grid=element_blank(),
            text = element_text(size=10),
            plot.title = element_text(hjust=0.5)) +
      guides(fill=guide_legend(title=mod),
             color=guide_legend(title=mod)))

  if(sm){margplot<-margplot +facet_grid(~dfplot.res[,mod])}

  if(model$family$link == "logit"){margplot<-margplot + ylim(0,1) + ylab("Probability")}

  invisible(margplot)
}
