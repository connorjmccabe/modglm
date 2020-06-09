modglm<-function(model, vars, data, part=NULL, cfacs="means", plotby=NULL,dd=FALSE)
{
  ints<<-list()
#This defines a string for the interaction term
(int.var <- paste(vars, collapse = ":"))

#Define Coefficients
  if (is.list(model$coefficients)){

    # if(is.null(part)){stop("Model part must be specified for two-part models. use part= in intglm function to specify.\n")}
    # b<-model$coef$count
    if(part=="count"){b<-model$coef$count}
    if(part=="zero"){
      if(grepl("zeroinfl",rb.ssdep$call[1])==T){b<- -(model$coef$zero)}
      else{b<-model$coef$zero}
      }

    }
  else{(b <- model$coef)}
  # b <- model$coef
#Define design matrix
  X<-as.data.frame(cbind(rep(1,nrow(model$model)),model$model[,-1]))
  
  if(cfacs[1]=="means"){cfs <- matrix(colMeans(X),nrow=1)}
  else(cfs<-cfacs)
  
  Xcfs<-as.data.frame(matrix(rep(cfs,each=nrow(X)),nrow=nrow(X)))
  colnames(X)[1]<-"(Intercept)"
  colnames(Xcfs)<-colnames(X)

  xmats<<-list()
  xmats$obs<-X
  xmats$means<-Xcfs
  
  ints<<-list()
  
for(i in 1:length(xmats)){
  
  X<-xmats[[i]]
    
  if (is.list(model$coefficients)){
    if(part=="count"){
      if (int.var %in% names(model$coefficients$count)){X[,int.var]<-X[,vars[1]]*X[,vars[2]]}
    }
    if(part=="zero"){
      if (int.var %in% names(model$coefficients$zero)){X[,int.var]<-X[,vars[1]]*X[,vars[2]]}
    }
  }
  else if (int.var %in% names(model$coefficients)){
    X[,int.var]<-X[,vars[1]]*X[,vars[2]]
    
    # modvars<-colnames(vcov(model))
    
    # if((!is.na(modvars[-(which(modvars==int.var))][which(grepl(":",modvars))]))){
    #   (modvars.noint<-as.vector(na.omit(modvars[-(which(modvars==int.var))][which(grepl(":",modvars))])))
    #   
    #   otherints<-strsplit(modvars.noint,":")
    #   
    #   for(j in 1:length(otherints)){
    #     X[,modvars.noint[j]]<-X[,otherints[[j]][1]]*X[,otherints[[j]][2]]
    #   }
    # }
    
    }
  
  X<-as.matrix(X)
  
    # X<-model.matrix(model,model$model)
#matrix multiply n x p design matrix by p x 1 vector of coefficients
  (xb <- (X %*% b)[,,drop=F])
  
#Computes predicted quantities of interest
  # (phat <- plogis(q=xb))
  if (is.list(model$coefficients)){
    if (part=="zero"){
      hat<-1/(1+ exp(-xb))
      # phatmean<-1/(1+ exp(-xbmean))
      #computes var(Yi)
      deriv1 <- exp(-xb)*(1+exp(-xb))^(-2)
      deriv2 <- ((exp(-xb)-1)*exp(-xb))/((exp(-xb)+1)^3)
      deriv3 <- (exp(-xb)*(exp(2*-xb)-4*exp(-xb)+1))/((exp(-xb)+1)^4)
    }
    else if (part=="count"){
      hat<-exp(xb)
      deriv1 <- exp(xb)
      deriv2 <- exp(xb)
      deriv3 <- exp(xb)
    }
    # else{stop("Error: this two-part model is yet supported in inteff.")}
  }

  else if (model$family$link == "logit")
  {
    hat<-1/(1+ exp(-xb))
    # phatmean<-1/(1+ exp(-xbmean))
#computes var(Yi)
  deriv1 <- exp(-xb)*(1+exp(-xb))^(-2)
  deriv2 <- ((exp(-xb)-1)*exp(-xb))/((exp(-xb)+1)^3)
  deriv3 <- (exp(-xb)*(exp(2*-xb)-4*exp(-xb)+1))/((exp(-xb)+1)^4)
  }
  else if (model$family$link == "inverse")
  {
    hat<-1/(xb)
    deriv1 <- -(1/(xb^2))
    deriv2 <- 2/(xb^3)
    deriv3 <- -(6/(xb^4))
  }

  else if (model$family$link == "log")
  {
    if(dd==F){
    hat<-exp(xb)
    deriv1 <- exp(xb)
    deriv2 <- exp(xb)
    deriv3 <- exp(xb)
    }
    if(dd==T){
      dum <- vars[which(sapply(apply(model$model[, vars], 2, table), length) ==
                          2)]
      cont <- vars[which(vars != dum)]
      X1 <- X2 <- X
      X1[, dum] <- 1
      X1[, int.var] <- X1[, cont] * X1[, dum]
      phat1 <- exp(X1 %*% b)
      phi1 <- exp(X1 %*% b)
      d2f1 <- exp(X1 %*% b)
      ie1 <- (b[cont] + b[int.var]) * phi1
      X2[, dum] <- 0
      X2[, int.var] <- X2[, cont] * X2[, dum]
      phat2 <- exp(X2 %*% b)
      phi2 <- exp(X2 %*% b)
      d2f2 <- exp(X2 %*% b)
      ie2 <- b[cont] * phi2
      int.est <- ie1 - ie2
      
      deriv1 <- phi1 - phi2 + b[cont] * X[, cont] * (d2f1 - d2f2) +
        b[int.var] * X[, cont] * d2f1
      deriv2 <- (b[cont] + b[int.var]) * d2f1
      deriv3 <- phi1 + (b[cont] + b[int.var]) * d2f1 * X[, cont]
      deriv0 <- (b[cont] + b[int.var]) * d2f1 - b[cont] * d2f2
      # dum <- names(which(lapply(sapply(as.data.frame(model$model),unique),length)==2))
      # # (dum <- vars[which(sapply(apply(X[, vars], 2, table), length) ==
      # #                     2)])
      # cont <- vars[which(vars != dum)]
      # X<-as.data.frame(X)
      # X1 <- X2 <- X
      # X1[, dum] <- 1
      # # X1[, int.var] <- X1[, cont] * X1[, dum]
      # 
      # hat1 <- exp(as.matrix(X1) %*% b)
      # 
      # d1f1 <- exp(as.matrix(X1) %*% b)
      # d2f1 <- exp(as.matrix(X1) %*% b)
      # 
      # # phi1 <- phat1 * (1 - phat1)
      # # d2f1 <- phi1 * (1 - 2 * phat1)
      # ie1 <- (b[cont] + b[int.var]) * d1f1
      # X2[, dum] <- 0
      # # X2[, int.var] <- X2[, cont] * X2[, dum]
      # hat2 <- exp(as.matrix(X2) %*% b)
      # d1f2 <- exp(as.matrix(X2) %*% b)
      # d2f2 <- exp(as.matrix(X2) %*% b)
      # ie2 <- b[cont] * d2f2
      # int.est <- ie1 - ie2
      # 
      # deriv1 <- d1f1 - d1f2 + b[cont] * X[, cont] * (d2f1 - d2f2) + b[int.var] * X[, cont] * d2f1
      # deriv2 <- (b[cont] + b[int.var]) * d2f1
      # deriv3 <- d1f1 + (b[cont] + b[int.var]) * d2f1 * X[, cont]
      # deriv0 <- (b[cont] + b[int.var]) * d2f1 - b[cont] * d2f2
    }
  }

  else if (model$family$link == "identity")
  {
    hat<-xb
    deriv1 <- 1
    deriv2 <- 0
    deriv3 <- 0
    
  }
  
  bint<-b[int.var]
  
  if (is.list(model$coefficients)){
    if(part=="count"){
      if (int.var %in% names(model$coefficients$count)==F){bint<-0}
    }
    if(part=="zero"){
      if (int.var %in% names(model$coefficients$zero)==F){bint<-0}
    }
  }
  
  else if (int.var %in% names(model$coefficients)==F){bint<-0}
  
  lin <- bint * deriv1

if(dd==F){


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

  if(dd){jcovar <- apply(covars, 2, function(x) ((b[cont] + b[int.var]) *d2f1 - b[cont] * d2f2) * x)}
  else(jcovar <- apply(covars, 2, function(x) bint * deriv2 * x +b1b4x2 * b2b4x1 * x * deriv3))
  jcovar <- array(jcovar, dim=dim(covars))
  dimnames(jcovar) <- dimnames(covars)

  if(dd){jac<-cbind(deriv1, deriv2, deriv3, jcovar, deriv0)[,,drop=F]}
  else(jac <- cbind(deriv11, deriv22, deriv44, jcovar, derivcc)[,,drop=F])
  colnames(jac) <- c(vars, int.var, colnames(jcovar), "(Intercept)")

  jac <- jac[, match(colnames(X), colnames(jac)), drop=F]

  if(is.list(model$coefficients)){
    if(part=="count"){
     vcov<-vcov(model)[grepl("count_", colnames(vcov(model))),grepl("count_", colnames(vcov(model)))]}

    if(part=="zero"){
      vcov<-vcov(model)[grepl("zero_", colnames(vcov(model))),grepl("zero_", colnames(vcov(model)))]}

    se <- sqrt(diag(jac %*% vcov %*% t(jac)))
  }


  else{se <- sqrt(diag(jac %*% vcov(model) %*% t(jac)))}

  t.val <- int.est/se

  # beta <- b * sx/sy

  betas<-int.est*(sd((model$model[,vars[1]]*model$model[,vars[2]]),na.rm=T)/sd(model$model[,1],na.rm=T))

  ints[[i]] <- data.frame(int.est = int.est, 
                          betas=betas, 
                          lin = lin,
                          # hat = hat,
                    se.int.est = se, 
                    t.val = t.val)
  # if(dd){
  #   ints[[i]]$deriv1<-deriv1
  #   ints[[i]]$deriv2<-deriv2
  #   ints[[i]]$deriv3<-deriv3
  #   ints[[i]]$jcovar<-jcovar
  #   ints[[i]]$jac<-jac
  # }
}
  names(ints)<-c("obints","intcf")
  ints$prop.sig<-length(which(abs(ints$obints$t.val)>1.96))/length(ints$obints$t.val)

  ints$obints$sig[abs(ints$obints$t.val)>=1.96]<-"Sig."
  ints$obints$sig[abs(ints$obints$t.val)<1.96]<-"N.S."

  ints$intcf<-as.vector(ints$intcf[1,])

  ints$model.summary<-summary(model)

  require(ggplot2)

  plotdf<-ints$obints
  if(is.null(plotby)){
  ints$intsplot<-ggplot(data=plotdf,aes(x=hat,y=betas, color=sig)) +
    geom_point(size=.75) +
    labs(x="Predicted Value",y="Interaction Effect") +
    theme_bw()
  }
  else{

    plotdf[,plotby]<-as.factor(data[,plotby])

    ints$intsplot<-ggplot(data=plotdf,aes(x=hat,y=betas, color=sig, fill=sig,shape=plotdf[,plotby])) +
      geom_point(size=1.5) +
      labs(x="Predicted Value",y="Interaction Effect") +
      # scale_fill_manual(values=c("white","black")) +
      # scale_color_manual(values=c("red","blue")) +
      scale_shape_manual(values=c(1,3)) +
      theme_bw()
  
  }
  
invisible(ints)
}

