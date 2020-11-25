# Modglm: Evaluating interaction effects in logit and count GLMs

The following are instructions for using the functions provided by our `modglm` code. We strongly encourage that users first read our accompanying manuscript before using this code. This manuscript is currently under review, though is available presently as a preprinted version on PsyArxiv:

https://psyarxiv.com/th94c

This code is largely adapted from the `intEff` function `DAMisc` package for logit and probit models (Armstrong & Armstrong, 2020), which itself was an adaptation of the `inteff` command in Stata (Norton, Wang, & Ai, 2004). Citations for each of these excellent resouces are as follows:

Armstrong, D., & Armstrong, M. D. (2020). Package `DAMisc`. 
URL: ftp://cygwin.uib.no/pub/cran/web/packages/DAMisc/DAMisc.pdf 

Norton, E. C., Wang, H., & Ai, C. (2004). Computing interaction effects and standard errors in logit and probit models. The Stata Journal, 4(2), 154-167.
URL: https://journals.sagepub.com/doi/abs/10.1177/1536867X0400400206

However, we note key differences between our code and those provided by the above resources. First, `modglm` provides functionality for computing interaction effects in Poisson and negative binomial models in addition to logit models. Second, we provide additional functionality that accomodates Implication 1 of our manuscript. This includes 1.) allowing users to more flexibly estimate interaction effects conditioned on user-specified hypothetical scenarios; 2.) creating a plot that summarizes point estimates of the interaction effect in the observed data; 3.) providing additional output that may be relevant in summarizing the results (e.g., the proportion of significant effects in the sample); and 4.) computes the average interaction effect in the observed sample, as well as its standard error, for summary purposes.
  
## Using `modglm`: Example

### Setup

The following instructions are based upon the simulated Poisson example presented in Equation 17 of our manuscript. For illustration purposes, we will assume that we are examining the relation between sensation seeking (`senseek`; i.e. disposition toward novelty and excitement) and the count of past year alcohol use (`pyalc`) treating biological sex as a moderator (`male`, dummy coded such that 0 = female and 1 = male) in this relation. Premeditation (`premed`; i.e. thinking before acting) was included as a covariate in this example. Assume that both sensation seeking and premeditation are standardized variables. We hypothesize that the effect of higher sensation seeking on the count of past year alcohol use is stronger among males compared to females, above and beyond the influence of other dispositional risk indicators (i.e. premeditation). Code for generating this data are as follows:

```
set.seed(1678)

b0 <- -3.8 ##Intercept
b1 <- .35 ###Sensation Seeking Effect
b2 <- .9 #Premeditation  Effect
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

df <- data.frame(y=y,senseek=rawvars[,1],premed=rawvars[,2],male=cat)
```

We may then use this data to estimate the model as follows:

```
pois<-glm(y ~ senseek + premed + male + senseek:male, data=df,family="poisson")
```
We can view the model summary results as follows:

```
> summary(pois)

Call:
glm(formula = y ~ senseek + premed + male + senseek:male, family = "poisson", 
    data = df)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.5069  -0.5234  -0.3070  -0.1742   4.1944  

Coefficients:
             Estimate Std. Error z value Pr(>|z|)    
(Intercept)  -3.18884    0.20845 -15.298  < 2e-16 ***
senseek       0.40202    0.15467   2.599  0.00934 ** 
premed        1.09892    0.08650  12.704  < 2e-16 ***
male          1.08317    0.22503   4.813 1.48e-06 ***
senseek:male -0.01371    0.17503  -0.078  0.93757    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 936.43  on 999  degrees of freedom
Residual deviance: 570.70  on 995  degrees of freedom
AIC: 837.7

Number of Fisher Scoring iterations: 6
```

Finally, we may source the `modglm` using the following code:

```
require(RCurl)
eval(parse(text = getURL("https://raw.githubusercontent.com/Modglmtemp/Modglm/master/modglm.R", 
                          ssl.verifypeer = FALSE)))
```

### The `modglm` Function

Assume that we aim to use `modglm` to estimate the interaction between sensation seeking and sex in the above model. We estimate interaction effects using `modglm` as follows:

```
pois.ints<-modglm(model=pois, vars=c("senseek","male"), data=df, type="fd", hyps="means")
```

Above, `modglm` requires at minimum four inputs, with one additional optional input.

The first is the estimated model object (e.g., `model=pois`). Currently, `model` may take logit or Poisson model objects estimated using the `stats` package and negative binomial model objects estimating using the `MASS` package. We hope to include additional GLMs in `modglm` in the near future, including zero-inflated and hurdle models and generalized estimating equation (GEE) GLMs.

The second is a 2-element vector of the variables included in the interaction term (e.g., `vars=c("senseek","male")`). As noted above, `modglm` makes no requirement that a product term need be specified in the model.

The third is the data frame used in estimating the model (e.g., `data=df`).

The fourth is the type of interaction being specified (e.g., `type="fd"`). This corresponds directly with the definition of interaction being used based on the variable type of the predictors. Three options are available, which mirrior our definition of interaction in Equation 10 in our manuscript. For continuous variable interactions, specify `type="cpd"` for computing the second-order cross-partial derivative. For continuous-by-discrete variable interactions, specify `type="fd"` for computing the finite difference in the partial derivative. For discrete variable interactions, specify `type="dd"` for computing the double (finite) difference.

Finally, `modglm` will optionally produce an interaction point estimate at a specified hypothetical condition using the `hyp` input (i.e., interactions at representative values as described in Implication 1 of the manuscript). By default, this is specified at the mean values of all included covariates (i.e. interaction effect at sample means). However, this can be modified by providing a vector of hypothetical values for the predictors involved in the model. For example, using `hyps=c(c(1,-.5,.25,0)` will provide estimates for when sensation seeking is 0.5 standard deviations from its own mean, premeditation is at 0.25 standard deviation above its own mean, and sex is 0 (i.e. female). Note that a 1 must be provided at the beginning of this vector to carry forward the intercept value. 

### `modglm` Output

`modglm` produces a list of objects summarizing the results. Noting that we have defined our output as `pois.ints`, use (for example) `names(pois.ints)` to view the elements of this list, which provides the following:

```
> names(pois.ints)
[1] "obints"        "inthyp"        "aie"           "desc"          "model.summary" "intsplot"  
```

`obints` provides a data frame of values computed observation-wise in the data. In order of columns, these include the interaction point estimates (`int.est`) and the predicted value of the outcome (`hat`) in the observed data. Other values provided are the delta method standard error of each interaction point estimate (`se.int.est`), the t-value (`t.val`), and the significance designation based on the t-value (`sig`). Example ouput for the first six observations is the following:

```
> head(pois.ints$obints)
      int.est        hat  se.int.est    t.val  sig
1 0.048693740 0.19323076 0.018002356 2.704854 Sig.
2 0.012310515 0.01668499 0.005918071 2.080157 Sig.
3 0.103610260 0.14034059 0.041662718 2.486882 Sig.
4 0.025082212 0.09972100 0.010919692 2.296971 Sig.
5 0.005403181 0.02132563 0.001940828 2.783957 Sig.
6 0.216222800 0.29885428 0.119301334 1.812409 N.S.
```

Above, `int.est` refers to the interaction effect estimate for each observation, `hat` refers to the predicted count of alcohol use, and `se.int.est`, `t.val`, and `sig` refer to the delta method standard error of the estimate, the corresponding t-value, and whether or not the estimate was statistically significant at the 0.05 level, respectively.

`inthyp` provides the results of the hypothetical condition specified by `hyps` as described above. This has the identical format as `obints` but contains only a single row of values:

```
> pois.ints$inthyp
     int.est        hat se.int.est    t.val   inthyp.ll inthyp.ul
1 0.02980042 0.06981395 0.01067509 2.791584 0.008877233 0.0507236

```
Here, given the default `hyp="means"`, the `int.est` value indicates the interaction effect estimate when sensation seeking, premeditation, and biological sex are each at their respective means. We may infer from these values that the interaction effect was 0.030 (95% CI = [0.01, 0.05]) at the hypothetical mean of all predictor variables, indicating that the marginal effect of sensation seeking on the count of past year alcohol use was stronger among males than females at these levels.

`aie` refers to the average interaction effect and follows a similar formatting:
```
> pois.ints$aie
     aie.est aie.se.delta     aie.ll    aie.ul
1 0.07240055   0.03350906 0.00672279 0.1380783
```
Here, we made a similar inference that the average interaction effect was significant and positive across observations (est. = 0.072, 95% CI = [0.01, 0.14]). We note in our manuscript that this value may be a more appropriate single estimate to provide as an in-text description of the interaction effect relative to the hypothetical mean value (e.g., Hanmer & Kalkan, 2013).

`desc` provides several other helpful descriptors of the interaction effect that researchers may wish to report. For the current example:

```
> pois.ints$desc
  int.range prop.sig prop.pos prop.neg
1    0-1.86    0.858        1        0
```

`int.range` refers to the range of interaction effects observed in the data. Here, the interaction effect ranged from 0 to 1.86 across observations. The `prop.sig` value indicates the proportion of values which were significant in the sample (i.e. 85.8% of interactions were significant). `prop.pos` and `prop.neg` indicate the proportion of interactions in the sample which were of positive or negative sign, respectively. Given all lower-order terms and the product term was positive here, 100% of the interactions were positive. However, we note that when the product term is negative, and/or in the case of logistic models where the nature of the logit function may change the interaction effect, it is possible that interaction effects may be either positive or negative in sign in the same sample. Substantively, we can state here that the marginal effect of sensation seeking on the count of past year alcohol use was stronger among males than females across the entirety of the sample, and was statistically different from zero for approximately 85% of the sample. These figures speak to the robustness of the positive interaction effect we have observed on average and at the hypothetical mean of all predictors.

`model.summary` is strictly a reproduction of a results summary table provided by the model (i.e. `summary(pois)`).

Finally, `intsplot` provides a graphical depiction of the interaction point estimates computed observation-wise, plotted against the model-predicted outcome (see also Ai & Norton, 2003). This plot is created using `ggplot2`. This provides a snapshot summary of the interaction effects present in the data, including the significance values and the potential range in the observed interaction effects:

```
> pois.ints$intsplot
```
![Alt Text](https://github.com/Modglmtemp/Modglm/blob/master/pois.intsplotex.png)
