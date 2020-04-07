# Modglm: Evaluating interaction effects in logit and count GLMs

The following are instructions for using the functions provided by the modglm package. We strongly encourage that users first read our accompanying manuscript before using this code. This manuscript is currently under review, though is available presently as a preprinted version on PsyArxiv:

LINK

This code is largely adapted from the `DAMisc` package for logit and probit models (Armstrong & Armstrong, 2020), which itself was an adaptation of the inteff command in Stata (Norton, Wang, & Ai, 2004). Citations for each of these excellent resouces are as follows:

Armstrong, D., & Armstrong, M. D. (2020). Package `DAMisc`. 
URL: ftp://cygwin.uib.no/pub/cran/web/packages/DAMisc/DAMisc.pdf

Norton, E. C., Wang, H., & Ai, C. (2004). Computing interaction effects and standard errors in logit and probit models. The Stata Journal, 4(2), 154-167.
URL: https://journals.sagepub.com/doi/abs/10.1177/1536867X0400400206

However, we note key differences between our code and those provided by the above resources. First, modglm provides functionality for computing interaction effects in Poisson and negative binomial models in addition to logit models, and can accomodate models involving generalized estimating equations (GEE). Second, we provide additional functionality that accomodates Implication 1 of our manuscript that allows users to more flexibly estimate interaction effects conditioned on user-specified hypothetical scenarios; creates a plot that summarizes point estimates of the interaction effect in the observed data; and provides output that may be relevant in summarizing the results (e.g., the proportion of significant effects in the sample). Third, interaction effects can be estimated without the inclusion of a product term in these models to accomodate Implication 2 of our manuscript (i.e. that product terms are not required to estimate interaction in GLMs).
  
## Using modglm: Example

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

We may use this data to then estimate the model as follows:

```
```

### Step 1:

![alt text](https://github.com/connorjmccabe/InterActive/blob/master/image%20files/Picture1.png)

### Step 2:

The model being estimated is *Y = X + Z + XZ*.

![alt text](https://github.com/connorjmccabe/InterActive/blob/master/image%20files/Picture2.png)

Note that X and Z are centered and Y is untransformed (raw) by default, though multiple scaling options are available. No covariates or quadratic effects will be specified in the current example.

### Step 3: Click the “Run Analysis” button (depicted above) to produce a summary table of coefficients and interaction graphics.

This will generate a summary table of coefficients, conduct regions of significance analyses, and will create the interaction graphics.

### Step 4: Adjust the values of each small multiple as desired.

![alt text](https://github.com/connorjmccabe/InterActive/blob/master/image%20files/Picture3.png)

Each multiple corresponds with a particular level of the moderator. By default, these range from -2 to 2 standard deviations because this will generally be representative of the moderator range given the moderator is approximately normal. However, we encourage you to modify these values to explore this functionality. Once small multiple values are decided upon, proceed to step 5 to customize the plot.

### Step 5: Click on the “Customize Plot” tab to change axes and plot titles, select the “greyscale” option, and download the finalized plot.

![alt text](https://github.com/connorjmccabe/InterActive/blob/master/image%20files/Picture4.png)

Users may use these customization options to ready this plot for use in a manuscript. Clicking the “Download Plot” button will output a .png file.

## Additional Functions in interActive

### Examine the marginal effects plot.

Click on the “Marginal Effects Plot” tab to view a depiction of the regions of significance analyses. Note that this plot displays standardized versions of X and Z predictor variables. Customization options for this plot are provided on this tab as well, and users can download this plot using the “Download Marginal Effects Plot” button below the graphic.

### View the plot estimates and raw data.

Users may view the uploaded data by clicking on the “Raw Data” tab to help identify problems in their data (e.g., missing values, outliers, etc.). Users can also view the plot estimates directly by clicking on the “Plot Estimates” tab to view the estimates used to create the small multiples graphic. 

### Download the R objects used to generate plots and upload them directly to R/Rstudio.

Users can download an .rds file using the “Download Plot Data” button located on the upper-righthand side of the screen under the "Plot" tab. Advanced users can read this file into R or Rstudio to customize their plots or examine their models in greater detail. Example code for reading an .rds file into R or Rstudio is provided below:
  
  ```
readRDS("~/Downloads/plotdata.rds")
```

### Explore the interaction with skewed predictor variables.
