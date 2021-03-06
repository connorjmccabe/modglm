#' Computing Interactions in nonlinear probability and count models
#'
#' @param model The estmiated model object.
#' #Currently, this may ntake logit or Poisson model objects
#' #estimated using the 'stats' package or negative binomial
#' #model objects ested using the 'MASS" package.
#' @param vars The interacting variable names.
#' @param data The data frame on which the estimated model is based.
#' @param hyps User-specified levels of the predictor variables
#' #for evaluating the interaction function. By default,
#' #this is specified at the mean values of all included covariates.
#' @param plotby An option to view the interaction effect
#' #estimates plot across categories of another variable. Default is NULL.
#' @param type The type of interaction being estimated.
#' #Options are ""cpd", "fd", and "dd". "cpd" indicates
#' #cross-partial derivative, and should be used when both
#' #variables are continuous. "fd" indicates finite
#' #difference, and should be used when one variable
#' #is continuous and the other is discrete.
#' #"dd" indicates double (finite) difference, and should be
#' #used if when both variables are categorical.
#'
#' @return
#' @export
#'
#' @examples
#'
#' #Simulate a dataset
#' set.seed(1678)
#'require(ggplot2)
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
#' #Evaluated the interaction function (i.e., computes
#' #finite difference values) for sensation seeking and sex
#' pois.ints<-modglm(model=pois, vars=c("senseek","male"), data=df, type="fd", hyps="means")
#'
#' names(pois.ints)
#'
#'#`obints` provides the interaction effect conditioned on each observation in the data. E.g:
#' head(pois.ints$obints)
#'
#'#`inthyp` provides the results of the hypothetical condition
#'# specified by 'hyps' above. In this case, this is
#'#represented at the mean of all predictors:
#' pois.ints$inthyp
#'
#'#`aie` refers to the average interaction effect.
#'#This is computed as the mean of all interaction
#'#effects in the observed data:
#' pois.ints$aie
#'
#'#`desc` provides several other helpful descriptors
#'#of the interaction effect that researchers may
#'#wish to report:
#' pois.ints$desc
#'
#'#`intsplot` provides a graphical depiction of the interaction
#'#point estimates computed observation-wise, plotted against
#'#the model-predicted outcome (see also Ai & Norton, 2003)
#'
#' pois.ints$intsplot
#'
modglm<-function(model, vars, data, hyps="means", plotby=NULL,type="cpd")
{
  if(length(vars)>2){stop("You have selected more than two interaction predictors. Please select only two predictors.")}
  if(length(vars)<2){stop("You have not selected enough predictors. Please select two predictors.")}

  if(2 %in% lapply(data[,vars],function(x){length(unique(x))}) & type=="cpd"){warning("Note: one or more of your interaction variables may be dichotomous, but continuous variable interaction ('cpd') has been requested. Please check the variables and consider finite differences ('fd') or discrete finite differences ('dd') as the interaction 'type' if these variables are dichotomous.")}

  if(FALSE %in% (lapply(data[,vars],function(x){length(unique(x))})>2) & type=="dd"){warning("Note: one or more of your interaction variables may be continuous, but discrete variable interaction ('dd') has been specified. Please check the variables and consider finite differences ('fd') or cross partial derivatives ('cpd') as the interaction 'type' if one or more of these variables is continuous.")}

  if((length(unique(data[,vars[1]]))>2)!=(length(unique(data[,vars[2]]))>2) & type!="fd"){warning("Note: one of your variables appears continous and the other discrete. Finite differences ('fd') may be most appropriate for the interaction 'type' given these variables.")}
  ints<<-list()
  #This defines a string for the interaction term
  (int.varpossible <- c(paste(vars, collapse = ":"),paste(rev(vars), collapse = ":")))
  if(any(int.varpossible %in% names(model$coefficients)))(int.var<-int.varpossible[(int.varpossible %in% names(model$coefficients))])
  else(int.var<-"No product term specified")
  jacs<<-list()

  (b <- model$coef)

  if(model$call[1]=="gee()"){
    dftemp<-na.omit(data[,which(names(data) %in% names(model$coefficients))])
    X<-as.data.frame(cbind(rep(1,nrow(dftemp)),dftemp))
  }

  else(X<-as.data.frame(cbind(rep(1,nrow(model$model)),model$model[,-1])))

  if(hyps[1]=="means"){cfs <- matrix(colMeans(X),nrow=1)}
  else(cfs<-hyps)

  Xcfs<-as.data.frame(matrix(rep(cfs,each=nrow(X)),nrow=nrow(X)))
  colnames(X)[1]<-"(Intercept)"
  colnames(Xcfs)<-colnames(X)

  xmats<<-list()
  xmats$obs<-X
  xmats$means<-Xcfs

  ints<<-list()

  for(i in 1:length(xmats)){

    X<-xmats[[i]]

    if (int.var %in% names(model$coefficients)){
      X[,int.var]<-X[,vars[1]]*X[,vars[2]]

    }

    X<-as.matrix(X)

    (xb <- (X %*% b)[,,drop=F])

    if(type=="cpd"){
      if (model$family$link == "logit")
      {
        hat<-1/(1+ exp(-xb))
        #computes var(Yi)
        deriv1 <- exp(-xb)*(1+exp(-xb))^(-2)
        deriv2 <- ((exp(-xb)-1)*exp(-xb))/((exp(-xb)+1)^3)
        deriv3 <- (exp(-xb)*(exp(2*-xb)-4*exp(-xb)+1))/((exp(-xb)+1)^4)
      }

      else if (model$family$link == "log")
      {
        hat<-exp(xb)
        deriv1 <- exp(xb)
        deriv2 <- exp(xb)
        deriv3 <- exp(xb)
      }
      else if (model$family$link == "identity")
      {
        hat<-xb
        deriv1 <- 1
        deriv2 <- 0
        deriv3 <- 0

      }
    }
    else if(type=="fd"){
      if(model$call[1]=="gee()"){dum <- vars[which(sapply(apply(dftemp[, vars], 2, table), length) ==2)]}
      else(dum <- vars[which(sapply(apply(model$model[, vars], 2, table), length) ==2)])
      cont <- vars[which(vars != dum)]
      X1 <- X2 <- as.data.frame(X)
      X1[, dum] <- 1
      X2[, dum] <- 0

      if (int.var %in% names(model$coefficients)==T){

        X1[, int.var] <- X1[, cont] * X1[, dum]
        X2[, int.var] <- X2[, cont] * X2[, dum]
        bint<-b[int.var]
      }
      else{
        bint<-0
      }
      X1<-as.matrix(X1)
      X2<-as.matrix(X2)

      x1b<-X1 %*% b
      x2b<-X2 %*% b

      if(model$family$link == "logit"){
        hat <- 1/(1+ exp(-xb))

        hat1 <- 1/(1+ exp(-x1b))
        d1f1 <- exp(-x1b)*(1+exp(-x1b))^(-2)#phi1
        d2f1 <- ((exp(-x1b)-1)*exp(-x1b))/((exp(-x1b)+1)^3) #d2f1

        hat2 <- 1/(1+ exp(-x2b))
        d1f2 <- exp(-x2b)*(1+exp(-x2b))^(-2) #phi2
        d2f2 <- ((exp(-x2b)-1)*exp(-x2b))/((exp(-x2b)+1)^3) #d2f2

      }

      else if (model$family$link == "log"){
        hat<-exp(X %*% b)

        hat1 <- exp(X1 %*% b)
        d1f1 <- exp(X1 %*% b)
        d2f1 <- exp(X1 %*% b)
        hat2 <- exp(X2 %*% b)
        d1f2 <- exp(X2 %*% b)
        d2f2 <- exp(X2 %*% b)
      }
      else if (model$family$link == "identity"){
        hat<-X %*% b

        hat1 <- X1 %*% b
        d1f1 <- 1
        d2f1 <- 0
        hat2 <- X2 %*% b
        d1f2 <- 1
        d2f2 <- 0
      }

      ie1 <- (b[cont] + bint) * d1f1
      ie2 <- b[cont] * d1f2
      int.est <- ie1 - ie2

      deriv1 <- d1f1 - d1f2 + b[cont] * X[, cont] * (d2f1 - d2f2) +
        bint * X[, cont] * d2f1
      deriv2 <- (b[cont] + bint) * d2f1
      deriv3 <- d1f1 + (b[cont] + bint) * d2f1 * X[, cont]
      deriv0 <- (b[cont] + bint) * d2f1 - b[cont] * d2f2
    }

    else if(type=="dd"){

      X00 <- X01 <- X10 <- X11<-as.data.frame(X)

      X00[, vars[1]] <- 0
      X00[, vars[2]] <- 0
      X01[, vars[1]] <- 0
      X01[, vars[2]] <- 1
      X10[, vars[1]] <- 1
      X10[, vars[2]] <- 0
      X11[, vars[1]] <- 1
      X11[, vars[2]] <- 1

      if (int.var %in% names(model$coefficients)==T){

        X00[, int.var] <- X00[, vars[1]] *X00[, vars[2]]
        X01[, int.var] <- X01[, vars[1]] *X01[, vars[2]]
        X10[, int.var] <- X10[, vars[1]] *X10[, vars[2]]
        X11[, int.var] <- X11[, vars[1]] *X11[, vars[2]]

        bint<-b[int.var]
      }
      else{
        bint<-0
      }
      X00<-as.matrix(X00)
      X01<-as.matrix(X01)
      X10<-as.matrix(X10)
      X11<-as.matrix(X11)

      x00b<-X00 %*% b
      x01b<-X01 %*% b
      x10b<-X10 %*% b
      x11b<-X11 %*% b

      if(model$family$link == "logit"){
        hat <- 1/(1+ exp(-xb))

        hat00 <- 1/(1+ exp(-x00b))
        d1f00 <- exp(-x00b)*(1+exp(-x00b))^(-2)

        hat01 <- 1/(1+ exp(-x01b))
        d1f01 <- exp(-x01b)*(1+exp(-x01b))^(-2)

        hat10 <- 1/(1+ exp(-x10b))
        d1f10 <- exp(-x10b)*(1+exp(-x10b))^(-2)

        hat11 <- 1/(1+ exp(-x11b))
        d1f11 <- exp(-x11b)*(1+exp(-x11b))^(-2)

      }

      else if (model$family$link == "log"){
        hat<-exp(X %*% b)

        hat00 <- exp(x00b)
        d1f00 <- exp(x00b)

        hat01 <- exp(x01b)
        d1f01 <- exp(x01b)

        hat10 <- exp(x10b)
        d1f10 <- exp(x10b)

        hat11 <- exp(x11b)
        d1f11 <- exp(x11b)
      }
      else if (model$family$link == "identity"){
        hat<-X %*% b

        hat00 <- x00b
        d1f00 <- x00b

        hat01 <- x01b
        d1f01 <- x01b

        hat10 <- x10b
        d1f10 <- x10b

        hat11 <- x11b
        d1f11 <- x11b
      }

      int.est <- (hat11-hat10)-(hat01-hat00)

      deriv1 <- d1f11-d1f10
      deriv2 <- d1f11-d1f01
      deriv3 <- d1f11
      deriv0 <- (d1f11-d1f01)-(d1f10-d1f00)
    }

    if(type=="cpd"){


      if (int.var %in% names(model$coefficients)==F){bint<-0}
      else(bint<-b[int.var])

      int.est <- bint * deriv1 +
        (b[vars[1]] + bint * X[, vars[2]]) *
        (b[vars[2]] + bint * X[, vars[1]]) *
        deriv2

      b1b4x2 <- b[vars[1]] + bint * X[, vars[2]]
      b2b4x1 <- b[vars[2]] + bint * X[, vars[1]]

      #Taking derivative of the interaction term with respect to X1
      deriv11 <- bint * deriv2 * X[, vars[1]] +
        b2b4x1 * deriv2 +
        b1b4x2 * b2b4x1 * deriv3 * X[, vars[1]]

      #Same with respect to X2
      deriv22 <- bint * deriv2 * X[, vars[2]] +
        b1b4x2 * deriv2 +
        b1b4x2 * b2b4x1 * X[, vars[2]] * deriv3

      #With respect to X1X2
      deriv44 <- deriv1 +
        bint * deriv2 * X[, vars[1]] * X[, vars[2]] +
        X[, vars[2]] * b2b4x1 * deriv2 + X[, vars[1]] * b1b4x2 *
        deriv2 + b1b4x2 * b2b4x1 * X[, vars[1]] * X[, vars[2]] *
        deriv3

      #with respect to the intercept?
      derivcc <- bint * deriv2 + b1b4x2 * b2b4x1 * deriv3
    }

    ##NOTE: This will add covariate values to the matxi
    if (int.var %in% names(model$coefficients)==T){covars <- X[, -c(1, match(c(vars, int.var), names(b)))]}
    else {covars <- X[, -c(1, match(vars, names(b)))]}

    if (!("matrix" %in% class(covars))) {
      covars <- matrix(covars, nrow = nrow(X))
    }

    if ((int.var %in% names(model$coefficients))==T)
    {colnames(covars) <- colnames(X)[-c(1, match(c(vars, int.var),
                                                 names(b)))]}
    else{colnames(covars) <- colnames(X)[-c(1, match(vars,names(b)))]}

    if(type=="cpd"){jcovar <- apply(covars, 2, function(x) bint * deriv2 * x +b1b4x2 * b2b4x1 * x * deriv3)}
    else if(type=="fd"){jcovar <- apply(covars, 2, function(x) ((b[cont] + bint) *d2f1 - b[cont] * d2f2) * x)}
    else if(type=="dd"){jcovar <- apply(covars, 2, function(x) ((d1f11 - d1f01) - (d1f10 - d1f00)) * x)}

    jcovar <- array(jcovar, dim=dim(covars))
    dimnames(jcovar) <- dimnames(covars)

    if(type=="cpd")(jac <- cbind(deriv11, deriv22, deriv44, jcovar, derivcc)[,,drop=F])
    else if(type=="fd"| type=="dd"){jac<-cbind(deriv1, deriv2, deriv3, jcovar, deriv0)[,,drop=F]}
    colnames(jac) <- c(vars, int.var, colnames(jcovar), "(Intercept)")

    jac <- jac[, match(colnames(X), colnames(jac)), drop=F]

    if(model$call[1]=="gee()"){se <- sqrt(diag(jac %*% gee_Rap_full$robust.variance %*% t(jac)))}

    else{se <- sqrt(diag(jac %*% vcov(model) %*% t(jac)))}

    t.val <- int.est/se


    ints[[i]] <- data.frame(int.est = int.est,
                            hat = hat,
                            se.int.est = se,
                            t.val = t.val)
    jacs[[i]] <- jac
  }


  ints$aie<-data.frame(aie.est=NA,aie.se.delta=NA)
  names(ints)<-c("obints","inthyp","aie")

  ints$jac<-jacs[[1]]
  ints$vcov<-vcov(model)

  ints$aie$aie.est<-mean(ints$obints$int.est)
  if(model$call[1]=="gee()"){ints$aie$aie.se <- sqrt(as.vector(colMeans(jacs[[1]]))%*%gee_Rap_full$robust.variance %*%as.vector(t(colMeans(jacs[[1]]))))}
  else{ints$aie$aie.se.delta <- sqrt(as.vector(colMeans(jacs[[1]]))%*%vcov(model)%*%as.vector(t(colMeans(jacs[[1]]))))}

  ints$aie$aie.ll<-ints$aie$aie.est-1.96*ints$aie$aie.se.delta
  ints$aie$aie.ul<-ints$aie$aie.est+1.96*ints$aie$aie.se.delta
  ints$prop.sig<-length(which(abs(ints$obints$t.val)>1.96))/length(ints$obints$t.val)

  ints$obints$sig[abs(ints$obints$t.val)>=1.96]<-"Sig."
  ints$obints$sig[abs(ints$obints$t.val)<1.96]<-"N.S."

  ints$inthyp<-as.vector(ints$inthyp[1,])

  ints$model.summary<-summary(model)

  requireNamespace("ggplot2")

  plotdf<-ints$obints

  if(is.null(plotby)){
    ints$intsplot<-ggplot2::ggplot(data=plotdf,ggplot2::aes(x=hat,y=int.est, color=sig)) +
      ggplot2::geom_point(size=.75) +
      ggplot2::labs(x="Predicted Value",y="Interaction Effect") +
      ggplot2::theme_bw()
  }
  else{

    plotdf[,plotby]<-as.factor(data[,plotby])

    ints$intsplot<-ggplot2::ggplot(data=plotdf,ggplot2::aes(x=hat,y=betas, color=sig, fill=sig,shape=plotdf[,plotby])) +
      ggplot2::geom_point(size=1.5) +
      ggplot2::labs(x="Predicted Value",y="Interaction Effect") +
      scale_shape_manual(values=c(1,3)) +
      theme_bw()

  }

  invisible(ints)
}


