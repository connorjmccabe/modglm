#' Plotting marginal effects for interactions in nonlinear probability and count models
#'
#' @param model The estmiated model object.
#' Currently, this may take logit or Poisson model objects
#' estimated using the 'stats' package or negative binomial
#' model objects tested using the 'MASS" package.
#' @param vars The interacting variable names.
#' @param data The data frame on which the estimated model is based.
#' @param hyps User-specified levels of the predictor variables
#' for displaying marginal effects. By default,
#' this is specified at the mean values of all included covariates.
#' @param foc  Name of the variable that will be the focal (i.e., x-axis) variable
#' @param mod Name of hte variable that will be the moderator (i.e., legend) variable
#' @param modlevels A vector of hypothetical levels of the moderator
#' at which to display marginal effects.
#' @param modnames If the user wants to call the moderator levels something other
#' Than what is provided in the data itself, the user can supply a vector of names here.
#' @param sm Stands for "small multiple". If TRUE, marginal effects will be displayed
#' as small multiple (i.e., trellis) plots. if FALSE, all marginal effects will appear
#' on a single plot.
#'
#' @return
#' @export
#'
#' @examples
#' #Building from the modglm example, this outputs a ggplot object plotting the effect of x1 on
#' #the count of y across males and females:
#' set.seed(1678)
#' library(ggplot2)
#' b0 <- -3.8 ##Intercept
#' b1 <- .35 ###Sensation Seeking Effect
#' b2 <- .9 #Premeditation  Effect
#' b3 <- 1.1 #Sex covariate effect
#' b13<- .2 #product term coefficient
#' n<-1000 #Sample Size
#' mu<-rep(0,2) #Specify means
#' S<-matrix(c(1,.5,.5,1),nrow=2,ncol=2) #Specify covariance matrix
#' sigma <- 1 #Level 1 error
#'
#' rawvars<-MASS::mvrnorm(n=n, mu=mu, Sigma=S) #simulates our
#' #continuous predictors from a multivariate normal distribution
#' cat<-rbinom(n=n,1,.5)
#' id<-seq(1:n)
#' eij <- rep(rnorm(id, 0, sigma))
#' xb<-  (b0) + (b1) * (rawvars[,1]) + (b2) * (rawvars[,2]) + (b3)*cat + b13*cat*(rawvars[,1]) + eij
# (b3) * (rawvars[,1]*rawvars[,2]) +
#'
#Generate Poisson data
#' ct <- exp(xb)
#' y <- rpois(n,ct)
#'
#' df <- data.frame(y=y,senseek=rawvars[,1],premed=rawvars[,2],male=cat)
#'
#Estimate a Poisson model regressing y on sensation seeking,
#'#premeditation, sex, and the interaction between sensation seeking and sex:
#'
#' pois<-glm(y ~ senseek + premed + male + senseek:male, data=df,family="poisson")
#'
#' p<-margplot(model=pois, vars=c("senseek","male"),foc="senseek", mod="male",
#'modlevels=c(0,1),data=df, hyps="means",sm=FALSE,modnames=c("Female","Male"))
#'
#'#Because this is a ggplot object, we can modify anything further by adding elements to this plot. E.g.:
#'p + ylab("Count")
margplot<-function(model, vars, data,hyps="means", foc, mod, modlevels,modnames=NULL,sm=T){


  if(!all(data.table::between(modlevels,min(data[,mod]),max(data[,mod])))){warning("Note: at least one hypothetical moderator level in 'modlevels' is outside the range of your moderator variable. Please check your plotted moderator values.")}

  (int.varpossible <- c(paste(vars, collapse = ":"),paste(rev(vars), collapse = ":")))
  (int.var<-int.varpossible[(int.varpossible %in% names(model$coefficients))])
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
  ximall<-data.table::rbindlist(xim.modlvs)

  #Step 4: SIMULATION.
  #define number of simulations. We will put this as 10,000 sims.
  sims <- 10000

  #4a: simulate parameters. The basis for why this works is that we assume asymtotic multivariate normality of
  #the parameter estimates when estimating our GLM. So, we can simulate these from the estimated
  #vcov and point estimates to generate approximate confidence regions.
  simparams<-MASS::mvrnorm(sims,pes,vcov)

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

  if(!is.null(modnames)){
    for(i in 1:length(modlevels)){
      dfplot.res[,mod][dfplot.res[,mod]==modlevels[i]]<-modnames[i]
    }
    dfplot.res[,mod]<-factor(dfplot.res[,mod],levels=modnames)
    }

  #Now, we can plot the results:

  (margplot<-ggplot2::ggplot(data=dfplot.res,ggplot2::aes(x=dfplot.res[,foc],y=pe,group=as.factor(dfplot.res[,mod]),fill=as.factor(dfplot.res[,mod]))) +
      ggplot2::geom_line(ggplot2::aes(color=as.factor(dfplot.res[,mod]))) +
      ggplot2::geom_ribbon(ggplot2::aes(x=dfplot.res[,foc], ymin = lower, ymax = upper), alpha=.5) +
      ggplot2::scale_fill_brewer(palette = "Pastel1")+ggplot2::scale_color_brewer(palette = "Set1")+
      # ylim(0,1) +
      ggplot2::ylab("DV") +
      # facet_grid(~dfplot.res[,mod]) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid=element_blank(),
            text = element_text(size=10),
            plot.title = element_text(hjust=0.5)) +
      ggplot2::guides(fill=guide_legend(title=mod),
             color=guide_legend(title=mod)))

  if(sm){margplot<-margplot +ggplot2::facet_grid(~dfplot.res[,mod])}

  if(model$family$link == "logit"){margplot<-margplot + ggplot2::ylim(0,1) + ggplot2::ylab("Probability")}

  return(margplot)
}
