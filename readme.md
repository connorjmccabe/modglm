# Modglm: Evaluating interaction effects in logit and count GLMs

The following are instructions for using the functions provided by our `modglm` code. We strongly encourage that users first read our accompanying manuscript before using this code. This manuscript is currently under review, though is available presently as a preprinted version on PsyArxiv:

(Omitted for review purposes)

This code is largely adapted from the `intEff` function `DAMisc` package for logit and probit models (Armstrong & Armstrong, 2020), which itself was an adaptation of the `inteff` command in Stata (Norton, Wang, & Ai, 2004). Citations for each of these excellent resouces are as follows:

Armstrong, D., & Armstrong, M. D. (2020). Package `DAMisc`. 
URL: ftp://cygwin.uib.no/pub/cran/web/packages/DAMisc/DAMisc.pdf 

Norton, E. C., Wang, H., & Ai, C. (2004). Computing interaction effects and standard errors in logit and probit models. The Stata Journal, 4(2), 154-167.
URL: https://journals.sagepub.com/doi/abs/10.1177/1536867X0400400206

However, we note key differences between our code and those provided by the above resources. First, `modglm` provides functionality for computing interaction effects in Poisson and negative binomial models in addition to logit models. Second, we provide additional functionality that accomodates Implication 1 of our manuscript. This includes 1.) allowing users to more flexibly estimate interaction effects conditioned on user-specified hypothetical scenarios; 2.) creating a plot that summarizes point estimates of the interaction effect in the observed data; and 3.) providing additional output that may be relevant in summarizing the results (e.g., the proportion of significant effects in the sample). Third, interaction effects can be estimated without the inclusion of a product term in these models to accomodate Implication 2 of our manuscript (i.e. that product terms are not required to estimate interaction in GLMs).
  
## Using `modglm`: Example

### Setup

The following instructions are based upon the simulated Poisson example presented in Equation 17 of our manuscript. Code for generating this data are as follows:

```
set.seed(1678)

b0 <- -3.8 ##Intercept
b1 <- .35 ###X1 Effect
b2 <- .9 #X2 Effect
b3 <- 1.1 #Sex covariate effect
b13<- .2 #product term coefficient
n<-1000 #Sample Size
mu<-rep(0,2) #Specify means
S<-matrix(c(1,.5,.5,1),nrow=2,ncol=2) #Specify covariance matrix
sigma <- 1 #Level 1 error

require(MASS)
rawvars<-mvrnorm(n=n, mu=mu, Sigma=S) #simulates our continuous predictors from a multivariate normal distribution
cat<-rbinom(n=n,1,.5)
id<-seq(1:n)
eij <- rep(rnorm(id, 0, sigma))
xb<-  (b0) + (b1) * (rawvars[,1]) + (b2) * (rawvars[,2]) + (b3)*cat + b13*cat*(rawvars[,1]) + eij
# (b3) * (rawvars[,1]*rawvars[,2]) + 

#Generate Poisson data
ct <- exp(xb)
y <- rpois(n,ct)

df <- data.frame(y=y,x1=rawvars[,1],x2=rawvars[,2],female=cat)
```

We may then use this data to estimate the model as follows:

```
pois<-glm(y ~ x1 + x2 + female + x1:female, data=df,family="poisson")
```
Finally, we may source the `modglm` using the following code:

```
require(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/Modglmtemp/Modglm/master/modglm.R", 
                          ssl.verifypeer = FALSE)))
```

### The `modglm` Function

Assume that we aim to use `modglm` to estimate the interaction between x1 and female in the above model. We estimate interaction effects using `modglm` as follows:

```
pois.ints<-modglm(model=pois, vars=c("x1","female"), data=df, type="fd", hyps="means")
```

Above, `modglm` requires at minimum four inputs, with one additional optional input.

The first is the estimated model object (e.g., `model=pois`). Currently, `model` may take logit or Poisson model objects estimated using the `stats` package and negative binomial model objects estimating using the `MASS` package. We hope to include additional GLMs in `modglm` in the near future, including zero-inflated and hurdle models and generalized estimating equation (GEE) GLMs.

The second is a 2-element vector of the variables included in the interaction term (e.g., `vars=c("x1","female")`). As noted above, `modglm` makes no requirement that a product term need be specified in the model.

The third is the data frame used in estimating the model (e.g., `data=df`).

The fourth is the type of interaction being specified (e.g., `type="fd"`). This corresponds directly with the definition of interaction being used based on the variable type of the predictors. Three options are available, which mirrior our definition of interaction in Equation 10 in our manuscript. For continuous variable interactions, specify `type="cpd"` for computing the second-order cross-partial derivative. For continuous-by-discrete variable interactions, specify `type="fd"` for computing the finite difference in the partial derivative. For discrete variable interactions, specify `type="dd"` for computing the double finite difference.

Finally, `modglm` will optionally produce an interaction point estimate at a specified hypothetical condition using the `hyp` input. By default, this is specified at the mean values of all included covariates. However, this can be modified by providing a vector of hypothetical values for the predictors involved in the model. For example, using `hyps=c(c(1,-.5,.25,0)` will provide estimates for when x1 is -0.5, x2 is 0.25, and female is 0. Note that a 1 must be provided at the beginning of this vector to carry forward the intercept value. 

### `moglm` Output

`modglm` produces a list of objects summarizing the results. Noting that we have defined our output as `pois.ints`, use (for example) `names(pois.ints)` to view the elements of this list, which provides the following:

```
> names(pois.ints)
[1] "obints"        "inthyp"        "prop.sig"      "model.summary" "intsplot" 
```

`obints` provides a data frame of values computed observation-wise in the data. In order of columns, these include the interaction point estimates in the data (`int.est`), the predicted value of the outcome (`hat`), the delta method standard error of interaction point estimate (`se.int.est`), the t-value (`t.val`), and the significance designation based on the t-value (`sig`).

`inthyp` provides the results of the hypothetical condition specified by `hyps` as described above. This has the identical format as `obints` but contains only a single row of values.

`prop.sig` provides the proportion of significant values among the observed data. This value may be helpful to users in summarizing the results in-text.

`model.summary` is strictly a reproduction of a results summary table provided by the model (i.e. `summary(pois)`).

Finally, `intsplot` provides a graphical depiction of the interaction point estimates computed observation-wise, plotted against the model-predicted outcome (see also Ai & Norton, 2003). This plot is created using `ggplot2`. This provides a snapshot summary of the interaction effects present in the data, including the significance values and the potential range in the observed interaction effects.