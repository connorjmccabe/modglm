#inteff replication

modglm<-function (model, vars, data)
{
#This defines a string for the interaction term
(int.var <- paste(vars, collapse = ":"))

#Define Coefficients
(b <- model$coef)

#Define an Xi; could be observed data OR mean-level data. lets go with observed for now
(X <- model.matrix(model, model$model))

# logit_cc <-
# function (obj = obj, int.var = int.var, vars = vars, b = b, X = X)
# {
#matrix multiply n x p matrix by p x 1 vector of coeffiicents
  (xb <- (X %*% b)[,,drop=F]) #What is drop=F?

#Computes predicted probabilities
  # (phat <- plogis(q=xb))
  if (model$family$link == "logit")
  {
    hat<-1/(1+ exp(-xb))
    # phatmean<-1/(1+ exp(-xbmean))
#computes var(Yi)
  d1f <- exp(-xb)*(1+exp(-xb))^(-2)
  d2f <- ((exp(-xb)-1)*exp(-xb))/((exp(-xb)+1)^3)
  d3f <- (exp(-xb)*(exp(2*-xb)-4*exp(-xb)+1))/((exp(-xb)+1)^4)
  }
  if (model$family$link == "log")
  {
    hat<-exp(xb)
    d1f <- exp(xb)
    d2f <- exp(xb)
    d3f <- exp(xb)
  }

  linear <- b[int.var] * d1f
#Connor's version:
  modelcc <- b[int.var] * d1f + #linear (computed above) plus...
    (b[vars[1]] + b[int.var] * X[, vars[2]]) * #coef$ss + (coef$sspre*Xpremed), times...
    (b[vars[2]] + b[int.var] * X[, vars[1]]) * #coef$pre + (coef$sspre*Xss), times...
    d2f #var(Yi) * (1-(2*predicted probabilities))
  
#UPDATE: I get this a little now. This code is designed to build the Jacobian matrix for my model. Put succinctly from https://cran.r-project.org/web/packages/modmarg/vignettes/delta-method.html#fnref6:
  
  # "So practically speaking, to get our variance, we’ll pre- and post-multiply the partial derivatives of the inverse link function by the original variance-covariance matrix from the regression"
  
  b1b4x2 <- b[vars[1]] + b[int.var] * X[, vars[2]]
  b2b4x1 <- b[vars[2]] + b[int.var] * X[, vars[1]]
  
  #Looks like we're taking the derivative of the interaction term with respect to senseek
  deriv11 <- b[int.var] * d2f * X[, vars[1]] + #int.coef. *2nd deriv * Xsenseek plus
              b2b4x1 * d2f +  #plus...
    b1b4x2 * b2b4x1 * d3f * X[, vars[1]] # #
  
  #Same with respect to premeditation
  deriv22 <- b[int.var] * d2f * X[, vars[2]] + b1b4x2 * d2f +
    b1b4x2 * b2b4x1 * X[, vars[2]] * d3f
  
  #then...with respect to the interaction?
  deriv44 <- d1f + b[int.var] * d2f * X[, vars[1]] * X[, vars[2]] +
    X[, vars[2]] * b2b4x1 * d2f + X[, vars[1]] * b1b4x2 *
    d2f + b1b4x2 * b2b4x1 * X[, vars[1]] * X[, vars[2]] *
    d3f
  
  #with respect to the intercept?
  derivcc <- b[int.var] * d2f + b1b4x2 * b2b4x1 * d3f
  
##NOTE: This will add covariate values to the matxi
  others <- X[, -c(1, match(c(vars, int.var), names(b)))]
  if (!("matrix" %in% class(others))) {
    others <- matrix(others, nrow = nrow(X))
  }
  colnames(others) <- colnames(X)[-c(1, match(c(vars, int.var),
                                              names(b)))]
#Because we take derivative with respect to some covariate, the equation simplifies to what we had for our intercept computation, but mulitplied by X.covar values.
  
  nn <- apply(others, 2, function(x) b[int.var] * d2f * x +
                b1b4x2 * b2b4x1 * x * d3f)
  nn <- array(nn, dim=dim(others))
  dimnames(nn) <- dimnames(others)
#####

  mat123 <- cbind(deriv11, deriv22, deriv44, nn, derivcc)[,,drop=F]
  colnames(mat123) <- c(vars, int.var, colnames(nn), "(Intercept)")
  
  mat123 <- mat123[, match(colnames(X), colnames(mat123)), drop=F] #reorders mat123 to be same as model
  model_se <- sqrt(diag(mat123 %*% vcov(model) %*% t(mat123))) #computes SEs (see my paper)
  model_t <- modelcc/model_se #computes t values
  out <- data.frame(int_eff = modelcc, linear = linear, hat = hat,
                    se_int_eff = model_se, zstat = model_t)
  invisible(out)
}
