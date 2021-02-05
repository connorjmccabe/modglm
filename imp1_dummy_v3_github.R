set.seed(1678)

b0 <- -3.8 ##Intercept
b1 <- .35 ###X1 Effect
b2 <- .9 #X2 Effect
b3 <- 1.1 #Sex covariate effect
b13<- .2
n<-1000 #Sample Size
mu<-rep(0,2) #Specify means
S<-matrix(c(1,.5,.5,1),nrow=2,ncol=2) #Specify covariance matrix
sigma <- 1 #Level 1 error

require(MASS)
rawvars<-mvrnorm(n=n, mu=mu, Sigma=S) #simulates our predictors from a multivariate normal distribution
cat<-rbinom(n=n,1,.5)
id<-seq(1:n)
eij <- rep(rnorm(id, 0, sigma))
xb<-  (b0) + (b1) * (rawvars[,1]) + (b2) * (rawvars[,2]) + (b3)*cat + b13*cat*(rawvars[,1]) + eij
# (b3) * (rawvars[,1]*rawvars[,2]) + 
# linear combination with a bias

#POISSON

ct <- exp(xb)
y <- rpois(n,ct)

df <- data.frame(y=y,x1=rawvars[,1],x2=rawvars[,2],female=cat)
df.test<-rbind(c(1,mean(df$x1),quantile(df$x2,.25),1),
               c(1,mean(df$x1),quantile(df$x2,.50),1),
               c(1,mean(df$x1),quantile(df$x2,.75),1),
               c(1,mean(df$x1),quantile(df$x2,.25),0),
               c(1,mean(df$x1),quantile(df$x2,.50),0),
               c(1,mean(df$x1),quantile(df$x2,.75),0),
               df)
# logit<-glm(y ~ x1 + x2 + female + x1:female, data=df,family="binomial")
pois<-glm(y ~ x1 + x2 + female + x1:female, data=df,family="poisson")
#source("~/Documents/modglm/modglm.R")

# 0.65*exp(0.65*mean(pois$model$x1)+1.57)-0.33*exp(0.33*mean(pois$model$x1)+0.63)
x1seq<-seq(-3,3,.1)
xi.mean.cat0<-data.frame(int=rep(1,length(x1seq)),
                         x1=x1seq,
                         x2=rep(mean(df$x2),length(x1seq)),
                         # x3=rep(mean(df$x3),length(x1seq)),
                         cat=rep(0,length(x1seq)),
                         x1cat=rep(mean(df$x1*df$female),length(x1seq)))

xi.q25.cat0<-data.frame(int=rep(1,length(x1seq)),
                        x1=x1seq,
                        x2=rep(quantile(df$x2)["25%"],length(x1seq)),
                        # x3=rep(mean(df$x3),length(x1seq)),
                        cat=rep(0,length(x1seq)),
                        x1cat=rep(mean(df$x1*df$female),length(x1seq)))

xi.q75.cat0<-data.frame(int=rep(1,length(x1seq)),
                        x1=x1seq,
                        x2=rep(quantile(df$x2)["75%"],length(x1seq)),
                        cat=rep(0,length(x1seq)),
                        x1cat=rep(mean(df$x1*df$female),length(x1seq)))

xi.mean.cat1<-data.frame(int=rep(1,length(x1seq)),
                         x1=x1seq,
                         x2=rep(mean(df$x2),length(x1seq)),
                         cat=rep(1,length(x1seq)),
                         x1cat=rep(mean(df$x1*df$female),length(x1seq)))

xi.q25.cat1<-data.frame(int=rep(1,length(x1seq)),
                        x1=x1seq,
                        x2=rep(quantile(df$x2)["25%"],length(x1seq)),
                        cat=rep(1,length(x1seq)),
                        x1cat=rep(mean(df$x1*df$female),length(x1seq)))

xi.q75.cat1<-data.frame(int=rep(1,length(x1seq)),
                        x1=x1seq,
                        x2=rep(quantile(df$x2)["75%"],length(x1seq)),
                        cat=rep(1,length(x1seq)),
                        x1cat=rep(mean(df$x1*df$female),length(x1seq)))

# xi.mean.cat0$x1x2<-xi.mean.cat0$x1*xi.mean.cat0$x2
# xi.q25.cat0$x1x2<-xi.q25.cat0$x1*xi.q25.cat0$x2
# xi.q75.cat0$x1x2<-xi.q75.cat0$x1*xi.q75.cat0$x2
xi.mean.cat0<-as.matrix(xi.mean.cat0)
xi.q25.cat0<-as.matrix(xi.q25.cat0)
xi.q75.cat0<-as.matrix(xi.q75.cat0)

# xi.mean.cat1$x1x2<-xi.mean.cat1$x1*xi.mean.cat1$x2
# xi.q25.cat1$x1x2<-xi.q25.cat1$x1*xi.q25.cat1$x2
# xi.q75.cat1$x1x2<-xi.q75.cat1$x1*xi.q75.cat1$x2
xi.mean.cat1<-as.matrix(xi.mean.cat1)
xi.q25.cat1<-as.matrix(xi.q25.cat1)
xi.q75.cat1<-as.matrix(xi.q75.cat1)

predct.mean.cat0<-as.data.frame(cbind(x1seq,rep("50th Percentile",length(x1seq)),exp(xi.mean.cat0%*%pois$coefficients)))
predct.q25.cat0<-as.data.frame(cbind(x1seq,rep("25th Percentile",length(x1seq)),exp(xi.q25.cat0%*%pois$coefficients)))
predct.q75.cat0<-as.data.frame(cbind(x1seq,rep("75th Percentile",length(x1seq)),exp(xi.q75.cat0%*%pois$coefficients)))

predct.mean.cat1<-as.data.frame(cbind(x1seq,rep("50th Percentile",length(x1seq)),exp(xi.mean.cat1%*%pois$coefficients)))
predct.q25.cat1<-as.data.frame(cbind(x1seq,rep("25th Percentile",length(x1seq)),exp(xi.q25.cat1%*%pois$coefficients)))
predct.q75.cat1<-as.data.frame(cbind(x1seq,rep("75th Percentile",length(x1seq)),exp(xi.q75.cat1%*%pois$coefficients)))

dfpreds.count.cat0<-rbind(predct.q25.cat0,predct.mean.cat0,predct.q75.cat0)
dfpreds.count.cat1<-rbind(predct.q25.cat1,predct.mean.cat1,predct.q75.cat1)
colnames(dfpreds.count.cat0)<-colnames(dfpreds.count.cat1)<-c("X1","X2 Percentile","pred_count")

dfpreds.count.cat0$X1<-as.numeric(as.character(dfpreds.count.cat0$X1))
dfpreds.count.cat0$pred_count<-as.numeric(as.character(dfpreds.count.cat0$pred_count))

dfpreds.count.cat1$X1<-as.numeric(as.character(dfpreds.count.cat1$X1))
dfpreds.count.cat1$pred_count<-as.numeric(as.character(dfpreds.count.cat1$pred_count))

require(ggplot2)

dfpreds.count.cat0$Sex<-"Male"
dfpreds.count.cat1$Sex<-"Female"
dfpreds.count<-rbind(dfpreds.count.cat0,dfpreds.count.cat1)

(poisplot<-ggplot(data=dfpreds.count,aes(x=X1,y=pred_count,linetype=Sex)) + 
    geom_line() +
    # ylim(0,1) +
    # ggtitle("Male") +
    ylab(expression(paste(hat(E),"[Y|",bold(x),"]"))) +
    facet_grid(~`X2 Percentile`)+
    theme_bw() +
    xlab(expression(italic(x)[1])) +
    theme(panel.grid.major =element_blank(),
          text = element_text(size=10),
          plot.title = element_text(hjust=0.5)))


ggsave("covarex.pdf",poisplot, width = 5.8, height = 3)

#######
#The rest of this code is for making other figures in the manuscript. Kept here in case it is helpful.
#######

# source("~/Documents/modglm/modglm.R")
# # source("~/Dropbox/Postdoc/GLM Int Bias/Garbage code/inteff_cm.R")
# 
# pois.ints<-modglm(model=pois,vars=c("x1","female"),data=df,type="fd",hyps="means")
# 
# pois.ints_test<-modglm(model=pois,vars=c("x1","female"),data=df,type="fd",hyps=c(1,-.5,.25,0))
# 
# mean(pois.ints$obints$int.est)
# mean(pois.ints$obints$se.int.est)
# 
# pois.ints$obints
# pois.ints$inthyp
# pois.ints$prop.sig
# pois.ints$model.summary
# pois.ints$intsplot
# 
# pois.ints$obints$Sex<-pois$model$female
# pois.ints$obints$` `<-pois.ints$obints$sig
# pois.ints$obints$Sex[pois.ints$obints$Sex==0]<-"Male"
# pois.ints$obints$Sex[pois.ints$obints$Sex==1]<-"Female"
# 
# pois.ints$obints$x1<-pois$model$x1
# pois.ints$obints$x2<-pois$model$x2
# 
# pois.ints$obints$x2.level[pois.ints$obints$x2<median(pois.ints$obints$x2)]<-"Lower Quartile of X2"
# pois.ints$obints$x2.level[pois.ints$obints$x2>=median(pois.ints$obints$x2)]<-"Upper Quartile of X2"
# 
# colMeans(pois.ints[[1]])%*%vcov(pois)%*%t(colMeans(pois.ints[[1]]))
# sig<-pois.ints$obints[which(pois.ints$obints$sig=="Sig."),]
# 
# nrow(sig)/nrow(pois.ints$obints)
# 
# 
# ints.lower<-pois.ints$obints[which(pois.ints$obints$x2.level=="Lower Quartile of X2"),]
# ints.upper<-pois.ints$obints[which(pois.ints$obints$x2.level=="Upper Quartile of X2"),]
# sig.ints.lower<-ints.lower[which(ints.lower$sig=="Sig."),]
# sig.ints.upper<-ints.upper[which(ints.upper$sig=="Sig."),]
# nrow(sig.ints.lower)
# nrow(sig.ints.upper)
# 
# range(sig.ints.lower$int.est)
# range(sig.ints.upper$int.est)
# 
# pois.ints$intsplot
# pois.ints$prop.sig
# # pois.ints$obints$betas[which(pois.ints$obints$betas>3*sd(pois.ints$obints$betas,na.rm=T) | pois.ints$obints$betas<3*sd(pois.ints$obints$betas))]<-NA
# pois.ints$obints$x2.level<-as.factor(pois.ints$obints$x2.level)
# levels(pois.ints$obints$x2.level)<-c(expression(paste("Below the Median of ",x[2]),paste("Above the Median of ",x[2])))
# 
# (p1<-ggplot(data=pois.ints$obints,aes(x=x1,y=int.est, color=` `)) + 
#         geom_point() +
#         labs(x="x1",y="Standardized Interaction Effect") +
#         scale_color_manual(values=c("grey","black")) +
#         ylim(quantile(pois.ints$obints$betas,.025),quantile(pois.ints$obints$betas,.975)) +
#         theme_bw() +
#         facet_grid(~x2.level, labeller=label_parsed) +
#         theme(panel.grid.major=element_blank(),
#               text = element_text(size=10),
#               plot.title = element_text(hjust=0.5)))
# 
# ggsave("intsexx1.pdf",p1, width = 5.8, height = 3)
# 

