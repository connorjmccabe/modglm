modglm<-function(model, vars, data, part=NULL, hyps="means", plotby=NULL,type="cpd")
{
  ints<<-list()
#This defines a string for the interaction term
(int.var <- paste(vars, collapse = ":"))

#Define Coefficients
  # if (is.list(model$coefficients)){
  # 
  #   # if(is.null(part)){stop("Model part must be specified for two-part models. use part= in intglm function to specify.\n")}
  #   # b<-model$coef$count
  #   if(part=="count"){b<-model$coef$count}
  #   if(part=="zero"){
  #     if(grepl("zeroinfl",rb.ssdep$call[1])==T){b<- -(model$coef$zero)}
  #     else{b<-model$coef$zero}
  #     }
  # 
  #   }
  # else{
    (b <- model$coef)
    # }
  # b <- model$coef
#Define design matrix
  # if(model$call[1]=="gee()"){
  #   dftemp<-na.omit(data[,which(names(data) %in% names(model$coefficients))])
  #   X<-as.data.frame(cbind(rep(1,nrow(dftemp)),dftemp))
  # }
  
  # else(
    X<-as.data.frame(cbind(rep(1,nrow(model$model)),model$model[,-1]))
    # )
  
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
    
  # if (is.list(model$coefficients)){
  #   if(part=="count"){
  #     if (int.var %in% names(model$coefficients$count)){X[,int.var]<-X[,vars[1]]*X[,vars[2]]}
  #   }
  #   if(part=="zero"){
  #     if (int.var %in% names(model$coefficients$zero)){X[,int.var]<-X[,vars[1]]*X[,vars[2]]}
  #   }
  # }
  # else 
    if (int.var %in% names(model$coefficients)){
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
  
  # if (is.list(model$coefficients)){
  #   if (part=="zero"){
  #     hat<-1/(1+ exp(-xb))
  #     # phatmean<-1/(1+ exp(-xbmean))
  #     #computes var(Yi)
  #     deriv1 <- exp(-xb)*(1+exp(-xb))^(-2)
  #     deriv2 <- ((exp(-xb)-1)*exp(-xb))/((exp(-xb)+1)^3)
  #     deriv3 <- (exp(-xb)*(exp(2*-xb)-4*exp(-xb)+1))/((exp(-xb)+1)^4)
  #   }
  #   else if (part=="count"){
  #     hat<-exp(xb)
  #     deriv1 <- exp(xb)
  #     deriv2 <- exp(xb)
  #     deriv3 <- exp(xb)
  #   }
  #   # else{stop("Error: this two-part model is yet supported in inteff.")}
  # }

  # else 
if(type=="cpd"){
    if (model$family$link == "logit")
  {
    hat<-1/(1+ exp(-xb))
#computes var(Yi)
  deriv1 <- exp(-xb)*(1+exp(-xb))^(-2)
  deriv2 <- ((exp(-xb)-1)*exp(-xb))/((exp(-xb)+1)^3)
  deriv3 <- (exp(-xb)*(exp(2*-xb)-4*exp(-xb)+1))/((exp(-xb)+1)^4)
  }
  # else if (model$family$link == "inverse")
  # {
  #   hat<-1/(xb)
  #   deriv1 <- -(1/(xb^2))
  #   deriv2 <- 2/(xb^3)
  #   deriv3 <- -(6/(xb^4))
  # }

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
      # if(model$call[1]=="gee()"){dum <- vars[which(sapply(apply(dftemp[, vars], 2, table), length) ==2)]}
      # else(
        dum <- vars[which(sapply(apply(model$model[, vars], 2, table), length) ==2)]
        # )
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
    
    int.est <- (hat11-hat10)-(hat01-hat00)
    
    deriv1 <- d1f11-d1f10
    deriv2 <- d1f11-d1f01
    deriv3 <- d1f11
    deriv0 <- (d1f11-d1f01)-(d1f10-d1f00)
  }

  # if (is.list(model$coefficients)){
  #   if(part=="count"){
  #     if (int.var %in% names(model$coefficients$count)==F){bint<-0}
  #   }
  #   if(part=="zero"){
  #     if (int.var %in% names(model$coefficients$zero)==F){bint<-0}
  #   }
  # }
  
  # else

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
# 
  # if(dd)(jcovar <- apply(covars, 2, function(x) (2 + 1) *d2f1 - 2 * d2f2 * x))

  if(type=="cpd"){jcovar <- apply(covars, 2, function(x) bint * deriv2 * x +b1b4x2 * b2b4x1 * x * deriv3)}
  else if(type=="fd"){jcovar <- apply(covars, 2, function(x) ((b[cont] + bint) *d2f1 - b[cont] * d2f2) * x)}
  else if(type=="dd"){jcovar <- apply(covars, 2, function(x) ((d1f11 - d1f01) - (d1f10 - d1f00)) * x)}
  
  jcovar <- array(jcovar, dim=dim(covars))
  dimnames(jcovar) <- dimnames(covars)

  if(type=="cpd")(jac <- cbind(deriv11, deriv22, deriv44, jcovar, derivcc)[,,drop=F])
  else if(type=="fd"| type=="dd"){jac<-cbind(deriv1, deriv2, deriv3, jcovar, deriv0)[,,drop=F]}
  colnames(jac) <- c(vars, int.var, colnames(jcovar), "(Intercept)")
  
  jac <- jac[, match(colnames(X), colnames(jac)), drop=F]

  # if(is.list(model$coefficients)){
  #   if(part=="count"){
  #    vcov<-vcov(model)[grepl("count_", colnames(vcov(model))),grepl("count_", colnames(vcov(model)))]}
  # 
  #   if(part=="zero"){
  #     vcov<-vcov(model)[grepl("zero_", colnames(vcov(model))),grepl("zero_", colnames(vcov(model)))]}
  # 
  #   se <- sqrt(diag(jac %*% vcov %*% t(jac)))
  # }


  # else{
  if(model$call[1]=="gee()"){se <- sqrt(diag(jac %*% model$robust.variance %*% t(jac)))}

  else{se <- sqrt(diag(jac %*% vcov(model) %*% t(jac)))}
    # }

  t.val <- int.est/se


  ints[[i]] <- data.frame(int.est = int.est,
                          hat = hat,
                    se.int.est = se,
                    t.val = t.val)
}

  names(ints)<-c("obints","inthyp")
  ints$prop.sig<-length(which(abs(ints$obints$t.val)>1.96))/length(ints$obints$t.val)

  ints$obints$sig[abs(ints$obints$t.val)>=1.96]<-"Sig."
  ints$obints$sig[abs(ints$obints$t.val)<1.96]<-"N.S."

  ints$inthyp<-as.vector(ints$inthyp[1,])

  ints$model.summary<-summary(model)

  require(ggplot2)

  plotdf<-ints$obints
  
  if(is.null(plotby)){
  ints$intsplot<-ggplot(data=plotdf,aes(x=hat,y=int.est, color=sig)) +
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
