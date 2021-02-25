set.seed(62346)
b0 <- 1 ##Intercept
b1 <- .7 ###X Slope Effect
b2 <- 1.5 #Z Slope Effect
b3 <- .1 #X*Z Slope Effect
n<-10000 #Sample Size
mu<-rep(0,2) #Specify means
S<-matrix(c(1,.5,.5,1),nrow=2,ncol=2) #Specify covariance matrix
sigma <- 1 #Level 1 error

require(MASS)
rawvars<-mvrnorm(n=n, mu=mu, Sigma=S) #simulates our predictors from a multivariate normal distribution
id<-seq(1:n)
eij <- rep(rnorm(id, 0, sigma))
xb<-  (b0) + (b1) * (rawvars[,1]) + (b2) * (rawvars[,2]) + (b3) * (rawvars[,1]) * (rawvars[,2]) + eij  # linear combination with a bias

pr <- 1/(1+exp(-xb))         # pass through an inv-logit function
y <- rbinom(n,1,pr)      # bernoulli response variable

df <- data.frame(id=id,y=y,x1=rawvars[,1],x2=rawvars[,2])
#LOGIT

logit<-glm( y~x1+x2 + x1:x2,data=df,family="binomial")

confint(logit)

xi.obs<-as.matrix(data.frame(int=rep(1,nrow(df)),x1=df$x1, x2=df$x2,x1x2=df$x1*df$x2))
predpr.obs<-1/(1+exp(-(xi.obs%*%logit$coefficients)))

x1seq<-seq(-3,3,.1)
x1seq<-seq(-3,3,.1)
xi.mean<-data.frame(int=rep(1,length(x1seq)),
                    x1=x1seq,
                    x2=rep(mean(df$x2),length(x1seq)))

xi.q25<-data.frame(int=rep(1,length(x1seq)),
                   x1=x1seq,
                   x2=rep(quantile(df$x2)["25%"],length(x1seq)))

xi.q75<-data.frame(int=rep(1,length(x1seq)),
                   x1=x1seq,
                   x2=rep(quantile(df$x2)["75%"],length(x1seq)))
xi.mean$x1x2<-xi.mean$x1*xi.mean$x2
xi.q25$x1x2<-xi.q25$x1*xi.q25$x2
xi.q75$x1x2<-xi.q75$x1*xi.q75$x2

xi.mean<-as.matrix(xi.mean)
xi.q25<-as.matrix(xi.q25)
xi.q75<-as.matrix(xi.q75)

predct.mean<-as.data.frame(cbind(x1seq,rep("50%",length(x1seq)),1/(1+exp(-(xi.mean%*%logit$coefficients)))))
predct.q25<-as.data.frame(cbind(x1seq,rep("25%",length(x1seq)),1/(1+exp(-(xi.q25%*%logit$coefficients)))))
predct.q75<-as.data.frame(cbind(x1seq,rep("75%",length(x1seq)),1/(1+exp(-(xi.q75%*%logit$coefficients)))))

dfpreds.logit<-rbind(predct.q25,predct.mean,predct.q75)
colnames(dfpreds.logit)<-c("x1","x2 Percentile","pred_count")
dfpreds.logit$x1<-as.numeric(as.character(dfpreds.logit$x1))
dfpreds.logit$pred_count<-as.numeric(as.character(dfpreds.logit$pred_count))

require(ggplot2)

(logitplot<-ggplot(data=dfpreds.logit,aes(x=x1,y=pred_count,linetype=`x2 Percentile`)) + 
    geom_line() +
    ylim(0,1) +
    ylab("Pr(Y)") +
###The commented code below will add the single-unit change estimates to this plot (see Figure 5 of manuscript), but will probably not be helpful for a manuscript. Uncomment to see these
    
    # geom_segment(x=1,
    #              y=dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"],
    #              xend=2,
    #              yend=dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"], size=.25,alpha=1, linetype="solid",color="black") +
    # geom_segment(x=2,
    #              y=dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"],
    #              xend=2,
    #              yend=dfpreds.logit[which(dfpreds.logit$x1==2 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"], size=.25,alpha=1, linetype="solid", color="black") +
    # geom_segment(x=1,
    #              y=dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"],
    #              xend=2,
    #              yend=dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"], size=.25,alpha=1, linetype="solid",color="black") +
    # geom_segment(x=2,
    #              y=dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"],
    #              xend=2,
    #              yend=dfpreds.logit[which(dfpreds.logit$x1==2 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"], size=.25,alpha=1, linetype="solid", color="black") +
    # 
    # geom_segment(x=-2,
    #              y=dfpreds.logit[which(dfpreds.logit$x1==-2 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"],
    #              xend=-1,
    #              yend=dfpreds.logit[which(dfpreds.logit$x1==-2 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"], size=.25,alpha=1, linetype="solid",color="black") +
    # geom_segment(x=-1,
    #              y=dfpreds.logit[which(dfpreds.logit$x1==-2 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"],
    #              xend=-1,
    #              yend=dfpreds.logit[which(dfpreds.logit$x1==-1 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"], size=.25,alpha=1, linetype="solid", color="black") +
    # geom_segment(x=-2,
    #              y=dfpreds.logit[which(dfpreds.logit$x1==-2 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"],
    #              xend=-1,
    #              yend=dfpreds.logit[which(dfpreds.logit$x1==-2 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"], size=.25,alpha=1, linetype="solid",color="black") +
    # geom_segment(x=-1,
    #              y=dfpreds.logit[which(dfpreds.logit$x1==-2 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"],
    #              xend=-1,
    #              yend=dfpreds.logit[which(dfpreds.logit$x1==-1 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"], size=.25,alpha=1, linetype="solid", color="black") +
    # xlab(expression(italic(x)[1])) +
    # ylab(expression(paste(hat(E),"[Y|",bold(x),"]"))) +
    # scale_linetype_discrete(name=expression(paste(x[2]," Percentile")))+
  theme_bw() +
    theme(panel.grid.major=element_blank(),
          text = element_text(size=10),
          plot.title = element_text(hjust=0.5)))  # theme(legend.position="none"))

ggsave("logitex.pdf",logitplot, width = 5.8, height = 3)

#######
#The rest of this code is for making other figures in the manuscript. Kept here in case it is helpful.
#######
# dfpreds.logit[which(dfpreds.logit$x1==2 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"]- dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"]
# 
# dfpreds.logit[which(dfpreds.logit$x1==2 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"]- dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"]
# 
# dfpreds.logit[which(dfpreds.logit$x1==-1 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"]- dfpreds.logit[which(dfpreds.logit$x1==-2 & dfpreds.logit$`x2 Percentile`=="75%"),"pred_count"]
# 
# dfpreds.logit[which(dfpreds.logit$x1==-1 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"]- dfpreds.logit[which(dfpreds.logit$x1==-2 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"]
# 
# # dfpreds.logit[which(dfpreds.logit$x1==2 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"]- dfpreds.logit[which(dfpreds.logit$x1==1 & dfpreds.logit$`x2 Percentile`=="25%"),"pred_count"]
# 
# 
# -0.05231793- 1.96*0.01106469
# -0.05231793 + 1.96*0.01106469
# source("~/Documents/modglm/modglm.R")
# 
# logitglm<-modglm(model=logit,vars=c("x1","x2"),data=df)
# 
# nrow(logitglm$obints %>% filter(sig=="Sig." & int.est<0))/nrow(logitglm$obints)
# range(logitglm$obints %>% filter(sig=="Sig." & int.est<0)%>%select(int.est))
# 
# nrow(logitglm$obints %>% filter(sig=="Sig." & int.est>0))/nrow(logitglm$obints)
# range(logitglm$obints %>% filter(sig=="Sig." & int.est>0)%>%select(int.est))
# 
# logitglm$inthyp$int.est
# logitglm$inthyp$int.est-1.96*logitglm$inthyp$se.int.est
# logitglm$inthyp$int.est+1.96*logitglm$inthyp$se.int.est
# 
# x1<-logit$model[1,"x1"]
# x2<-logit$model[1,"x2"]
# 
# int<-logit$coefficients["(Intercept)"]
# b1<-logit$coefficients["x1"]
# b2<-logit$coefficients["x2"]
# eta<-int+b1*x1+b2*x2
# 
# b12<-0
#   
# #(gam12<-1.05*((exp(-(1+.7*x1+1.5*x2))-1)/(exp(2+1.4*x1+3*x2)*(1+exp(-(1+.7*x1+1.5*x2))^5))))
# d1f <- exp(-eta)*(1+exp(-eta))^(-2)
# d2f <- ((exp(-eta)-1)*exp(-eta))/((exp(-eta)+1)^3)
# linear <- b12 * d1f
# #Connor's version:
# (modelcc <- b12 * d1f + #linear (computed above) plus...
#  (b1 + b12 * x2) * #coef$ss + (coef$sspre*Xpremed), times...
#  (b2 + b12 * x1) * #coef$pre + (coef$sspre*Xss), times...
#   d2f) #var(Yi) * (1-(2*predicted probabilities))
# 
# b1*b2*d2f
# 
# tfrac<-exp(-eta)-1
# bfrac<-exp(2*eta)*(1+exp(-eta)^5)
# 
# -(int+b1*x1+b2*x2)
# 
# -int-b1*x1-b2*x2
# 
# b1*b2*(tfrac/bfrac)
# 
# b1*b2*(exp(-(eta))*((1+exp(-(eta)))^-2))*((exp(-(eta))-1)*(exp(-(eta)))/((exp(-(eta))+1)^3))
# 
# logitglm$intsplot
# logitglm$obints$sig[abs(logitglm$obints$t.val)>=1.96]<-"p<.05"
# logitglm$obints$sig[abs(logitglm$obints$t.val)<1.96]<-"p>.05"
# 
# neg<-logitglm$obints$int.est[which(logitglm$obints$int.est<0 & abs(logitglm$obints$t.val)>1.96)]
# pos<-logitglm$obints$int.est[which(logitglm$obints$int.est>0 & abs(logitglm$obints$t.val)>1.96)]
# 
# length(neg)/length(logitglm$obints$int.est)
# range(neg)
# 
# length(pos)/length(logitglm$obints$int.est)
# range(pos)
# 
# ggplot() + geom_point(aes(x=df$x2,y=logitglm$obints$betas, color=logitglm$obints$sig)) + theme_bw()
# 
# # ggplot() + geom_point(aes(x=df$x1,y=logitglm$obints$betas, color=logitglm$obints$sig)) + theme_bw()
