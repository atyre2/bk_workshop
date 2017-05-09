# Supplement  2.  Cade, B. S.  2015.  Model averaging and muddled multimodel inferences.  Ecology.
# 
# R script to generate the simulations in Appendix B for multi-part compositional predictors within a zero-truncated Poisson regression model similar to the breeding sage-grouse count model of Rice et al. (2013).

###Set up simulations for compositional data, 5 parts, one small part 
###deleted with proportions similar to those for sage-grouse breeding 
###habitat as in Rice et al. (2013), 80% sagebrush/shrub, 5% herbaceous, 
###10% forest, 4% alpine, and 1% ###urban (deleted from consideration).

set.seed(87)
pos.tot <- runif(200,min=0.8,max=1.0)
urban.tot <- pmin(runif(200,min=0.0,max=0.02),1.0 - pos.tot)

neg.tot <- (1.0 - pmin(pos.tot + urban.tot,1))


x1<- pmax(pos.tot - runif(200,min=0.05,max=0.30),0)

x3<- pmax(neg.tot - runif(200,min=0.0,max=0.10),0)
x2<- pmax(pos.tot - x1 - x3/2,0)

x4<- pmax(1 - x1 - x2 - x3 - urban.tot,0)

cor(cbind(x1,x2,x3,x4,urban.tot))
cor(cbind(x1 + x2, x3 + x4))

###Poisson model of means changing as a function of sagebrush/shrub (x1) 
###and herbaceous (x2) cover types.  Because x3 + x4 is near perfect 
###linear function of x1 + x2, means also change with exp(2.206 -9.515x3
###-5.894x4), asymptotically
 
mean.y <- exp(-5.8 + 6.3*x1 + 15.2*x2)
y <- rpois(200,mean.y)

cor(cbind(y,x1,x2,x3,x4))

###Create zero-truncated sample
y.zerotrunc <- y[y>=1]
x1.zt <- x1[y>=1]
x2.zt <- x2[y>=1]
x3.zt <- x3[y>=1]
x4.zt <- x4[y>=1]

library(VGAM)

zero.trunc.x1 <- vglm(y.zerotrunc~x1.zt,family=pospoisson())
zero.trunc.x1x2 <- vglm(y.zerotrunc~x1.zt+x2.zt,family=pospoisson())
zero.trunc.x1x3 <- vglm(y.zerotrunc~x1.zt+x3.zt,family=pospoisson())
zero.trunc.x1x4 <- vglm(y.zerotrunc~x1.zt + x4.zt,family=pospoisson())
zero.trunc.x1x3x4 <- vglm(y.zerotrunc~x1.zt + x3.zt + x4.zt,family=pospoisson())
zero.trunc.x1x2x4 <- vglm(y.zerotrunc~x1.zt + x2.zt + x4.zt,family=pospoisson())
zero.trunc.x1x2x3 <- vglm(y.zerotrunc~x1.zt + x2.zt + x3.zt,family=pospoisson())
zero.trunc.x1x2x3x4 <- vglm(y.zerotrunc~x1.zt + x2.zt + x3.zt + x4.zt,family=pospoisson())

###To examine negative relationship
zero.trunc.x3x4 <- vglm(y.zerotrunc~ x3.zt + x4.zt,family=pospoisson())


###Compute AIC and weights.

AIC.x1 <- AIC(zero.trunc.x1)
AIC.x1x2 <-AIC(zero.trunc.x1x2)
AIC.x1x3 <-AIC(zero.trunc.x1x3)
AIC.x1x4 <-AIC(zero.trunc.x1x4)
AIC.x1x3x4 <-AIC(zero.trunc.x1x3x4)
AIC.x1x2x4 <-AIC(zero.trunc.x1x2x4)
AIC.x1x2x3 <-AIC(zero.trunc.x1x2x3)
AIC.x1x2x3x4 <-AIC(zero.trunc.x1x2x3x4)

minAIC <- min(AIC.x1,AIC.x1x2,AIC.x1x3,AIC.x1x4, AIC.x1x3x4,AIC.x1x2x4,AIC.x1x2x3,AIC.x1x2x3x4)

d.AIC.x1 <- AIC.x1 - minAIC
d.AIC.x1x2 <- AIC.x1x2 - minAIC
d.AIC.x1x3 <- AIC.x1x3 - minAIC
d.AIC.x1x4 <- AIC.x1x4 - minAIC
d.AIC.x1x2x4 <- AIC.x1x2x4 - minAIC
d.AIC.x1x3x4 <- AIC.x1x3x4 - minAIC
d.AIC.x1x2x3 <- AIC.x1x2x3 - minAIC
d.AIC.x1x2x3x4 <- AIC.x1x2x3x4 - minAIC

sum.wts <- ((exp(-0.5 * d.AIC.x1)) + (exp(-0.5 * d.AIC.x1x2)) + (exp(-0.5 * d.AIC.x1x3)) + (exp(-0.5 * d.AIC.x1x4)) + (exp(-0.5 * d.AIC.x1x2x4)) + (exp(-0.5 * d.AIC.x1x3x4)) + (exp(-0.5 * d.AIC.x1x2x3)) + (exp(-0.5 * d.AIC.x1x2x3x4)))   

AIC.x1.wt <- (exp(-0.5 * d.AIC.x1))/sum.wts   
AIC.x1x2.wt <- (exp(-0.5 * d.AIC.x1x2))/sum.wts   
AIC.x1x3.wt <- (exp(-0.5 * d.AIC.x1x3))/sum.wts   
AIC.x1x4.wt <- (exp(-0.5 * d.AIC.x1x4))/sum.wts   
AIC.x1x2x3.wt <- (exp(-0.5 * d.AIC.x1x2x3))/sum.wts   
AIC.x1x3x4.wt <- (exp(-0.5 * d.AIC.x1x3x4)) /sum.wts   
AIC.x1x2x4.wt <- (exp(-0.5 * d.AIC.x1x2x4)) /sum.wts   
AIC.x1x2x3x4.wt <- (exp(-0.5 * d.AIC.x1x2x3x4)) /sum.wts 

###Estimates, se, and model-averaged estimates for B1.

comp.B1.x1 <- coef(zero.trunc.x1)[2] 
comp.B1.x1x2 <- coef(zero.trunc.x1x2)[2] 
comp.B1.x1x3 <- coef(zero.trunc.x1x3)[2] 
comp.B1.x1x4 <- coef(zero.trunc.x1x4)[2] 
comp.B1.x1x2x3 <- coef(zero.trunc.x1x2x3)[2] 
comp.B1.x1x2x4 <- coef(zero.trunc.x1x2x4)[2] 
comp.B1.x1x3x4 <- coef(zero.trunc.x1x3x4)[2] 
comp.B1.x1x2x3x4 <- coef(zero.trunc.x1x2x3x4)[2] 

comp.seB1.x1 <- coef(summary(zero.trunc.x1))[2,2] 
comp.seB1.x1x2 <- coef(summary(zero.trunc.x1x2))[2,2] 
comp.seB1.x1x3 <- coef(summary(zero.trunc.x1x3))[2,2] 
comp.seB1.x1x4 <- coef(summary(zero.trunc.x1x4))[2,2] 
comp.seB1.x1x2x3 <- coef(summary(zero.trunc.x1x2x3))[2,2] 
comp.seB1.x1x2x4 <- coef(summary(zero.trunc.x1x2x4))[2,2] 
comp.seB1.x1x3x4 <- coef(summary(zero.trunc.x1x3x4))[2,2] 
comp.seB1.x1x2x3x4 <- coef(summary(zero.trunc.x1x2x3x4))[2,2] 


comp.B1.wt <- ((comp.B1.x1 * AIC.x1.wt) + (comp.B1.x1x2 * AIC.x1x2.wt) + (comp.B1.x1x3 * AIC.x1x3.wt) +(comp.B1.x1x4 * AIC.x1x4.wt) + (comp.B1.x1x2x3 * AIC.x1x2x3.wt) + (comp.B1.x1x3x4 * AIC.x1x3x4.wt) + (comp.B1.x1x2x4 * AIC.x1x2x4.wt) + (comp.B1.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

comp.B1.wt.se <- (((((comp.B1.x1 - comp.B1.wt)^2 + comp.seB1.x1^2)^0.5) * AIC.x1.wt) + ((((comp.B1.x1x2 - comp.B1.wt)^2 + comp.seB1.x1x2^2)^0.5) * AIC.x1x2.wt) + ((((comp.B1.x1x3 - comp.B1.wt)^2 + comp.seB1.x1x3^2)^0.5) * AIC.x1x3.wt) + ((((comp.B1.x1x4 - comp.B1.wt)^2 + comp.seB1.x1x4^2)^0.5) * AIC.x1x4.wt) + ((((comp.B1.x1x2x3 - comp.B1.wt)^2 + comp.seB1.x1x2x3^2)^0.5) * AIC.x1x2x3.wt) + ((((comp.B1.x1x3x4 - comp.B1.wt)^2 + comp.seB1.x1x3x4^2)^0.5) * AIC.x1x3x4.wt) + ((((comp.B1.x1x2x4 - comp.B1.wt)^2 + comp.seB1.x1x2x4^2)^0.5) * AIC.x1x2x4.wt) + ((((comp.B1.x1x2x3x4 - comp.B1.wt)^2 + comp.seB1.x1x2x3x4^2)^0.5) * AIC.x1x2x3x4.wt)) 

###Estimate partial standard deviations of predictors.

###To get VIF shift to glm 
pois.x1 <- glm(y.zerotrunc~x1.zt,family=poisson)
pois.x1x2 <- glm(y.zerotrunc~x1.zt+x2.zt,family=poisson)
pois.x1x3 <- glm(y.zerotrunc~x1.zt+x3.zt,family=poisson)
pois.x1x4 <- glm(y.zerotrunc~x1.zt + x4.zt,family=poisson)
pois.x1x3x4 <- glm(y.zerotrunc~x1.zt + x3.zt + x4.zt,family=poisson)
pois.x1x2x4 <- glm(y.zerotrunc~x1.zt + x2.zt + x4.zt,family=poisson)
pois.x1x2x3 <- glm(y.zerotrunc~x1.zt + x2.zt + x3.zt,family=poisson)
pois.x1x2x3x4 <- glm(y.zerotrunc~x1.zt + x2.zt + x3.zt + x4.zt, family=poisson)


library(car)
###inverse of VIF is proportion of variance in x1 that is not linearly 
####related to other predictors in model


comp.ivif.x1x2 <- (vif(pois.x1x2)[1])^-1
comp.ivif.x1x3 <- (vif(pois.x1x3)[1])^-1
comp.ivif.x1x4 <- (vif(pois.x1x4)[1])^-1
comp.ivif.x1x2x3 <- (vif(pois.x1x2x3)[1])^-1
comp.ivif.x1x2x4 <- (vif(pois.x1x2x4)[1])^-1
comp.ivif.x1x3x4 <- (vif(pois.x1x3x4)[1])^-1
comp.ivif.x1x2x3x4 <- (vif(pois.x1x2x3x4)[1])^-1

comp.partSD.x1 <- sd(x1.zt)

comp.partSD.x1.x1x2 <- comp.ivif.x1x2^0.5 * ((nobs(zero.trunc.x1x2)-1)/(nobs(zero.trunc.x1x2) -(length(coef(zero.trunc.x1x2))-1)))^0.5 * sd(x1.zt) 

comp.partSD.x1.x1x3 <- comp.ivif.x1x3^0.5 * ((nobs(zero.trunc.x1x3)-1)/(nobs(zero.trunc.x1x3) -(length(coef(zero.trunc.x1x3))-1)))^0.5 * sd(x1.zt)
 
comp.partSD.x1.x1x4 <- comp.ivif.x1x4^0.5 * ((nobs(zero.trunc.x1x4)-1)/(nobs(zero.trunc.x1x4) -(length(coef(zero.trunc.x1x4))-1)))^0.5 * sd(x1.zt)
 
comp.partSD.x1.x1x2x3 <- comp.ivif.x1x2x3^0.5 * ((nobs(zero.trunc.x1x2x3)-1)/(nobs(zero.trunc.x1x2x3) -(length(coef(zero.trunc.x1x2x3))-1)))^0.5 * sd(x1.zt)
 
comp.partSD.x1.x1x2x4 <- comp.ivif.x1x2x4^0.5 * ((nobs (zero.trunc.x1x2x4)-1)/(nobs(zero.trunc.x1x2x4) -(length(coef(zero.trunc.x1x2x4))-1)))^0.5 * sd(x1.zt)
 
comp.partSD.x1.x1x3x4 <- comp.ivif.x1x3x4^0.5 * ((nobs(zero.trunc.x1x3x4)-1)/(nobs(zero.trunc.x1x3x4) -(length(coef(zero.trunc.x1x3x4))-1)))^0.5 * sd(x1.zt)
 
comp.partSD.x1.x1x2x3x4 <- comp.ivif.x1x2x3x4^0.5 * ((nobs(zero.trunc.x1x2x3x4)-1)/(nobs(zero.trunc.x1x2x3x4) -(length(coef(zero.trunc.x1x2x3x4))-1)))^0.5 * sd(x1.zt)

###inverse of VIF is proportion of variance in x2 that is not linearly 
####related to other predictors in model


comp.ivif.x1x2.x2 <- (vif(pois.x1x2)[2])^-1
comp.ivif.x1x3.x2 <- 1
comp.ivif.x1x4.x2 <- 1
comp.ivif.x1x2x3.x2 <- (vif(pois.x1x2x3)[2])^-1
comp.ivif.x1x2x4.x2 <- (vif(pois.x1x2x4)[2])^-1
comp.ivif.x1x3x4.x2 <- 1
comp.ivif.x1x2x3x4.x2 <- (vif(pois.x1x2x3x4)[2])^-1

comp.partSD.x2 <- sd(x2.zt)

comp.partSD.x2.x1x2 <- comp.ivif.x1x2.x2^0.5 * ((nobs(zero.trunc.x1x2)-1)/(nobs(zero.trunc.x1x2) -(length(coef(zero.trunc.x1x2))-1)))^0.5 * sd(x2.zt) 

comp.partSD.x2.x1x2x3 <- comp.ivif.x1x2x3.x2^0.5 * ((nobs(zero.trunc.x1x2x3)-1)/(nobs(zero.trunc.x1x2x3) -(length(coef(zero.trunc.x1x2x3))-1)))^0.5 * sd(x2.zt)
 
comp.partSD.x2.x1x2x4 <- comp.ivif.x1x2x4.x2^0.5 * ((nobs (zero.trunc.x1x2x4)-1)/(nobs(zero.trunc.x1x2x4) -(length(coef(zero.trunc.x1x2x4))-1)))^0.5 * sd(x2.zt)
 
comp.partSD.x2.x1x2x3x4 <- comp.ivif.x1x2x3x4.x2^0.5 * ((nobs(zero.trunc.x1x2x3x4)-1)/(nobs(zero.trunc.x1x2x3x4) -(length(coef(zero.trunc.x1x2x3x4))-1)))^0.5 * sd(x2.zt) 

###inverse of VIF is proportion of variance in x3 that is not linearly 
###related to other predictors in model


comp.ivif.x1x2.x3 <- 1
comp.ivif.x1x3.x3 <- (vif(pois.x1x3)[2])^-1
comp.ivif.x1x4.x3 <- 1
comp.ivif.x1x2x3.x3 <- (vif(pois.x1x2x3)[3])^-1
comp.ivif.x1x2x4.x3 <- 1
comp.ivif.x1x3x4.x3 <- (vif(pois.x1x3x4)[2])^-1
comp.ivif.x1x2x3x4.x3 <- (vif(pois.x1x2x3x4)[3])^-1

comp.partSD.x3 <- sd(x3.zt)

comp.partSD.x3.x1x3 <- comp.ivif.x1x3.x3^0.5 * ((nobs(zero.trunc.x1x3)-1)/(nobs(zero.trunc.x1x3) -(length(coef(zero.trunc.x1x3))-1)))^0.5 * sd(x3.zt) 
 
comp.partSD.x3.x1x2x3 <- comp.ivif.x1x2x3.x3^0.5 * ((nobs(zero.trunc.x1x2x3)-1)/(nobs(zero.trunc.x1x2x3) -(length(coef(zero.trunc.x1x2x3))-1)))^0.5 * sd(x3.zt)
 
comp.partSD.x3.x1x3x4 <- comp.ivif.x1x3x4.x3^0.5 * ((nobs(zero.trunc.x1x3x4)-1)/(nobs(zero.trunc.x1x3x4) -(length(coef(zero.trunc.x1x3x4))-1)))^0.5 * sd(x3.zt)
 
comp.partSD.x3.x1x2x3x4 <- comp.ivif.x1x2x3x4.x3^0.5 * ((nobs(zero.trunc.x1x2x3x4)-1)/(nobs(zero.trunc.x1x2x3x4) -(length(coef(zero.trunc.x1x2x3x4))-1)))^0.5 * sd(x3.zt) 

###inverse of VIF is proportion of variance in x4 that is not linearly 
###related to other predictors in model


comp.ivif.x1x2.x4 <- 1
comp.ivif.x1x3.x4 <- 1
comp.ivif.x1x4.x4 <- (vif(pois.x1x4)[2])^-1
comp.ivif.x1x2x3.x4 <- 1
comp.ivif.x1x2x4.x4 <- (vif(pois.x1x2x4)[3])^-1
comp.ivif.x1x3x4.x4 <- (vif(pois.x1x3x4)[3])^-1
comp.ivif.x1x2x3x4.x4 <- (vif(pois.x1x2x3x4)[4])^-1

comp.partSD.x4 <- sd(x4.zt)

comp.partSD.x4.x1x4 <- comp.ivif.x1x4.x4^0.5 * ((nobs(zero.trunc.x1x4)-1)/(nobs(zero.trunc.x1x4) -(length(coef(zero.trunc.x1x4))-1)))^0.5 * sd(x4.zt) 
 
comp.partSD.x4.x1x2x4 <- comp.ivif.x1x2x4.x4^0.5 * ((nobs(zero.trunc.x1x2x4)-1)/(nobs(zero.trunc.x1x2x4) -(length(coef(zero.trunc.x1x2x4))-1)))^0.5 * sd(x4.zt)
 
comp.partSD.x4.x1x3x4 <- comp.ivif.x1x3x4.x4^0.5 * ((nobs(zero.trunc.x1x3x4)-1)/(nobs(zero.trunc.x1x3x4) -(length(coef(zero.trunc.x1x3x4))-1)))^0.5 * sd(x4.zt)
 
comp.partSD.x4.x1x2x3x4 <- comp.ivif.x1x2x3x4.x4^0.5 * ((nobs(zero.trunc.x1x2x3x4)-1)/(nobs(zero.trunc.x1x2x3x4) -(length(coef(zero.trunc.x1x2x3x4))-1)))^0.5 * sd(x4.zt) 


###Estimate standardized estimates by transforming design matrix.
std.zero.trunc.x1 <- vglm(y.zerotrunc~ I(x1.zt/comp.partSD.x1),family=pospoisson())
std.zero.trunc.x1x2 <- vglm(y.zerotrunc~I(x1.zt/comp.partSD.x1.x1x2) + I(x2.zt/comp.partSD.x2.x1x2),family=pospoisson())
std.zero.trunc.x1x3 <- vglm(y.zerotrunc~I(x1.zt/comp.partSD.x1.x1x3)+ I(x3.zt/comp.partSD.x3.x1x3),family=pospoisson())
std.zero.trunc.x1x4 <- vglm(y.zerotrunc~I(x1.zt/comp.partSD.x1.x1x4) + I(x4.zt/comp.partSD.x4.x1x4),family=pospoisson())
std.zero.trunc.x1x3x4 <- vglm(y.zerotrunc~I(x1.zt/comp.partSD.x1.x1x3x4) + I(x3.zt/comp.partSD.x3.x1x3x4) + I(x4.zt/comp.partSD.x4.x1x3x4),family=pospoisson())
std.zero.trunc.x1x2x4 <- vglm(y.zerotrunc~I(x1.zt/comp.partSD.x1.x1x2x4) + I(x2.zt/comp.partSD.x2.x1x2x4) + I(x4.zt/comp.partSD.x2.x1x2x4),family=pospoisson())
std.zero.trunc.x1x2x3 <- vglm(y.zerotrunc~I(x1.zt/comp.partSD.x1.x1x2x3) + I(x2.zt/comp.partSD.x2.x1x2x3) + I(x3.zt/comp.partSD.x3.x1x2x3),family=pospoisson())
std.zero.trunc.x1x2x3x4 <- vglm(y.zerotrunc~I(x1.zt/comp.partSD.x1.x1x2x3x4) + I(x2.zt/comp.partSD.x2.x1x2x3x4)  + I(x3.zt/comp.partSD.x3.x1x2x3x4)  + I(x4.zt/comp.partSD.x4.x1x2x3x4) ,family=pospoisson())

###Standardized estimates,se, and model-averaged estimates for B1.

comp.stdB1.x1 <- coef(std.zero.trunc.x1)[2] 
comp.stdB1.x1x2 <- coef(std.zero.trunc.x1x2)[2] 
comp.stdB1.x1x3 <- coef(std.zero.trunc.x1x3)[2] 
comp.stdB1.x1x4 <- coef(std.zero.trunc.x1x4)[2] 
comp.stdB1.x1x2x3 <- coef(std.zero.trunc.x1x2x3)[2] 
comp.stdB1.x1x2x4 <- coef(std.zero.trunc.x1x2x4)[2] 
comp.stdB1.x1x3x4 <- coef(std.zero.trunc.x1x3x4)[2] 
comp.stdB1.x1x2x3x4 <- coef(std.zero.trunc.x1x2x3x4)[2] 

comp.sestdB1.x1 <- coef(summary(std.zero.trunc.x1))[2,2] 
comp.sestdB1.x1x2 <- coef(summary(std.zero.trunc.x1x2))[2,2] 
comp.sestdB1.x1x3 <- coef(summary(std.zero.trunc.x1x3))[2,2] 
comp.sestdB1.x1x4 <- coef(summary(std.zero.trunc.x1x4))[2,2] 
comp.sestdB1.x1x2x3 <- coef(summary(std.zero.trunc.x1x2x3))[2,2] 
comp.sestdB1.x1x2x4 <- coef(summary(std.zero.trunc.x1x2x4))[2,2] 
comp.sestdB1.x1x3x4 <- coef(summary(std.zero.trunc.x1x3x4))[2,2] 
comp.sestdB1.x1x2x3x4 <- coef(summary(std.zero.trunc.x1x2x3x4))[2,2] 


comp.stdB1.wt <- ((comp.stdB1.x1 * AIC.x1.wt) + (comp.stdB1.x1x2 * AIC.x1x2.wt) + (comp.stdB1.x1x3 * AIC.x1x3.wt) +(comp.stdB1.x1x4 * AIC.x1x4.wt) + (comp.stdB1.x1x2x3 * AIC.x1x2x3.wt) + (comp.stdB1.x1x3x4 * AIC.x1x3x4.wt) + (comp.stdB1.x1x2x4 * AIC.x1x2x4.wt) + (comp.stdB1.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

comp.stdB1.wt.se <- (((((comp.stdB1.x1 - comp.stdB1.wt)^2 + comp.sestdB1.x1^2)^0.5) * AIC.x1.wt) + ((((comp.stdB1.x1x2 - comp.stdB1.wt)^2 + comp.sestdB1.x1x2^2)^0.5) * AIC.x1x2.wt) + ((((comp.stdB1.x1x3 - comp.stdB1.wt)^2 + comp.sestdB1.x1x3^2)^0.5) * AIC.x1x3.wt) + ((((comp.stdB1.x1x4 - comp.stdB1.wt)^2 + comp.sestdB1.x1x4^2)^0.5) * AIC.x1x4.wt) + ((((comp.stdB1.x1x2x3 - comp.stdB1.wt)^2 + comp.sestdB1.x1x2x3^2)^0.5) * AIC.x1x2x3.wt) + ((((comp.stdB1.x1x3x4 - comp.stdB1.wt)^2 + comp.sestdB1.x1x3x4^2)^0.5) * AIC.x1x3x4.wt) + ((((comp.stdB1.x1x2x4 - comp.stdB1.wt)^2 + comp.sestdB1.x1x2x4^2)^0.5) * AIC.x1x2x4.wt) + ((((comp.stdB1.x1x2x3x4 - comp.stdB1.wt)^2 + comp.sestdB1.x1x2x3x4^2)^0.5) * AIC.x1x2x3x4.wt)) 


###To convert model-averaged standardized estimates back to given 
###covariance structure.

###Correlation among x1 and x2

comp.stdB1.wt*((comp.partSD.x1.x1x2)^-1)
((comp.stdB1.wt.se^2)*((comp.partSD.x1.x1x2)^-2))^0.5
 
comp.stdB2.wt*((comp.partSD.x2.x1x2)^-1)
((comp.stdB2.wt.se^2)*((comp.partSD.x2.x1x2)^-2))^0.5

###Correlation among x1,x2,x3, and x4

comp.stdB1.wt*((comp.partSD.x1.x1x2x3x4)^-1)
((comp.stdB1.wt.se^2)*((comp.partSD.x1.x1x2x3x4)^-2))^0.5

comp.stdB2.wt*((comp.partSD.x2.x1x2x3x4)^-1)
((comp.stdB2.wt.se^2)*((comp.partSD.x2.x1x2x3x4)^-2))^0.5

###Other predictors
###X2
###Estimates, se, and model-averaged estimates for B2.

comp.B2.x1 <- 0.0 
comp.B2.x1x2 <- coef(zero.trunc.x1x2)[3] 
comp.B2.x1x3 <- 0.0 
comp.B2.x1x4 <- 0.0 
comp.B2.x1x2x3 <- coef(zero.trunc.x1x2x3)[3] 
comp.B2.x1x2x4 <- coef(zero.trunc.x1x2x4)[3] 
comp.B2.x1x3x4 <- 0.0 
comp.B2.x1x2x3x4 <- coef(zero.trunc.x1x2x3x4)[3] 

comp.seB2.x1 <- 0.0 
comp.seB2.x1x2 <- coef(summary(zero.trunc.x1x2))[3,2] 
comp.seB2.x1x3 <- 0.0 
comp.seB2.x1x4 <- 0.0 
comp.seB2.x1x2x3 <- coef(summary(zero.trunc.x1x2x3))[3,2] 
comp.seB2.x1x2x4 <- coef(summary(zero.trunc.x1x2x4))[3,2] 
comp.seB2.x1x3x4 <- 0.0 
comp.seB2.x1x2x3x4 <- coef(summary(zero.trunc.x1x2x3x4))[3,2] 


comp.B2.wt <- ((comp.B2.x1 * AIC.x1.wt) + (comp.B2.x1x2 * AIC.x1x2.wt) + (comp.B2.x1x3 * AIC.x1x3.wt) +(comp.B2.x1x4 * AIC.x1x4.wt) + (comp.B2.x1x2x3 * AIC.x1x2x3.wt) + (comp.B2.x1x3x4 * AIC.x1x3x4.wt) + (comp.B2.x1x2x4 * AIC.x1x2x4.wt) + (comp.B2.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

comp.B2.wt.se <- (((((comp.B2.x1 - comp.B2.wt)^2 + comp.seB2.x1^2)^0.5) * AIC.x1.wt) + ((((comp.B2.x1x2 - comp.B2.wt)^2 + comp.seB2.x1x2^2)^0.5) * AIC.x1x2.wt) + ((((comp.B2.x1x3 - comp.B2.wt)^2 + comp.seB2.x1x3^2)^0.5) * AIC.x1x3.wt) + ((((comp.B2.x1x4 - comp.B2.wt)^2 + comp.seB2.x1x4^2)^0.5) * AIC.x1x4.wt) + ((((comp.B2.x1x2x3 - comp.B2.wt)^2 + comp.seB2.x1x2x3^2)^0.5) * AIC.x1x2x3.wt) + ((((comp.B2.x1x3x4 - comp.B2.wt)^2 + comp.seB2.x1x3x4^2)^0.5) * AIC.x1x3x4.wt) + ((((comp.B2.x1x2x4 - comp.B2.wt)^2 + comp.seB2.x1x2x4^2)^0.5) * AIC.x1x2x4.wt) + ((((comp.B2.x1x2x3x4 - comp.B2.wt)^2 + comp.seB2.x1x2x3x4^2)^0.5) * AIC.x1x2x3x4.wt)) 
  
###Standardized estimates, se, and model-averaged estimates for B2.

comp.stdB2.x1 <- 0.0 
comp.stdB2.x1x2 <- coef(std.zero.trunc.x1x2)[3] 
comp.stdB2.x1x3 <- 0.0 
comp.stdB2.x1x4 <- 0.0 
comp.stdB2.x1x2x3 <- coef(std.zero.trunc.x1x2x3)[3] 
comp.stdB2.x1x2x4 <- coef(std.zero.trunc.x1x2x4)[3] 
comp.stdB2.x1x3x4 <- 0.0
comp.stdB2.x1x2x3x4 <- coef(std.zero.trunc.x1x2x3x4)[3] 

comp.sestdB2.x1 <- 0.0
comp.sestdB2.x1x2 <- coef(summary(std.zero.trunc.x1x2))[3,2] 
comp.sestdB2.x1x3 <- 0.0 
comp.sestdB2.x1x4 <- 0.0 
comp.sestdB2.x1x2x3 <- coef(summary(std.zero.trunc.x1x2x3))[3,2] 
comp.sestdB2.x1x2x4 <- coef(summary(std.zero.trunc.x1x2x4))[3,2] 
comp.sestdB2.x1x3x4 <- 0.0 
comp.sestdB2.x1x2x3x4 <- coef(summary(std.zero.trunc.x1x2x3x4))[3,2] 


comp.stdB2.wt <- ((comp.stdB2.x1 * AIC.x1.wt) + (comp.stdB2.x1x2 * AIC.x1x2.wt) + (comp.stdB2.x1x3 * AIC.x1x3.wt) +(comp.stdB2.x1x4 * AIC.x1x4.wt) + (comp.stdB2.x1x2x3 * AIC.x1x2x3.wt) + (comp.stdB2.x1x3x4 * AIC.x1x3x4.wt) + (comp.stdB2.x1x2x4 * AIC.x1x2x4.wt) + (comp.stdB2.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

comp.stdB2.wt.se <- (((((comp.stdB2.x1 - comp.stdB2.wt)^2 + comp.sestdB2.x1^2)^0.5) * AIC.x1.wt) + ((((comp.stdB2.x1x2 - comp.stdB2.wt)^2 + comp.sestdB2.x1x2^2)^0.5) * AIC.x1x2.wt) + ((((comp.stdB2.x1x3 - comp.stdB2.wt)^2 + comp.sestdB2.x1x3^2)^0.5) * AIC.x1x3.wt) + ((((comp.stdB2.x1x4 - comp.stdB2.wt)^2 + comp.sestdB2.x1x4^2)^0.5) * AIC.x1x4.wt) + ((((comp.stdB2.x1x2x3 - comp.stdB2.wt)^2 + comp.sestdB2.x1x2x3^2)^0.5) * AIC.x1x2x3.wt) + ((((comp.stdB2.x1x3x4 - comp.stdB2.wt)^2 + comp.sestdB2.x1x3x4^2)^0.5) * AIC.x1x3x4.wt) + ((((comp.stdB2.x1x2x4 - comp.stdB2.wt)^2 + comp.sestdB2.x1x2x4^2)^0.5) * AIC.x1x2x4.wt) + ((((comp.stdB2.x1x2x3x4 - comp.stdB2.wt)^2 + comp.sestdB2.x1x2x3x4^2)^0.5) * AIC.x1x2x3x4.wt)) 


###X3
###Estimates, se, and model-averaged estimates for B3.

comp.B3.x1 <- 0.0 
comp.B3.x1x2 <- 0.0 
comp.B3.x1x3 <- coef(zero.trunc.x1x3)[3] 
comp.B3.x1x4 <- 0.0 
comp.B3.x1x2x3 <- coef(zero.trunc.x1x2x3)[4] 
comp.B3.x1x2x4 <- 0.0 
comp.B3.x1x3x4 <- coef(zero.trunc.x1x3x4)[3]  
comp.B3.x1x2x3x4 <- coef(zero.trunc.x1x2x3x4)[4] 

comp.seB3.x1 <- 0.0 
comp.seB3.x1x2 <- 0.0
comp.seB3.x1x3 <- coef(summary(zero.trunc.x1x3))[3,2]  
comp.seB3.x1x4 <- 0.0 
comp.seB3.x1x2x3 <- coef(summary(zero.trunc.x1x2x3))[4,2] 
comp.seB3.x1x2x4 <- 0.0 
comp.seB3.x1x3x4 <- coef(summary(zero.trunc.x1x3x4))[3,2]  
comp.seB3.x1x2x3x4 <- coef(summary(zero.trunc.x1x2x3x4))[4,2] 


comp.B3.wt <- ((comp.B3.x1 * AIC.x1.wt) + (comp.B3.x1x2 * AIC.x1x2.wt) + (comp.B3.x1x3 * AIC.x1x3.wt) +(comp.B3.x1x4 * AIC.x1x4.wt) + (comp.B3.x1x2x3 * AIC.x1x2x3.wt) + (comp.B3.x1x3x4 * AIC.x1x3x4.wt) + (comp.B3.x1x2x4 * AIC.x1x2x4.wt) + (comp.B3.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

comp.B3.wt.se <- (((((comp.B3.x1 - comp.B3.wt)^2 + comp.seB3.x1^2)^0.5) * AIC.x1.wt) + ((((comp.B3.x1x2 - comp.B3.wt)^2 + comp.seB3.x1x2^2)^0.5) * AIC.x1x2.wt) + ((((comp.B3.x1x3 - comp.B3.wt)^2 + comp.seB3.x1x3^2)^0.5) * AIC.x1x3.wt) + ((((comp.B3.x1x4 - comp.B3.wt)^2 + comp.seB3.x1x4^2)^0.5) * AIC.x1x4.wt) + ((((comp.B3.x1x2x3 - comp.B3.wt)^2 + comp.seB3.x1x2x3^2)^0.5) * AIC.x1x2x3.wt) + ((((comp.B3.x1x3x4 - comp.B3.wt)^2 + comp.seB3.x1x3x4^2)^0.5) * AIC.x1x3x4.wt) + ((((comp.B3.x1x2x4 - comp.B3.wt)^2 + comp.seB3.x1x2x4^2)^0.5) * AIC.x1x2x4.wt) + ((((comp.B3.x1x2x3x4 - comp.B3.wt)^2 + comp.seB3.x1x2x3x4^2)^0.5) * AIC.x1x2x3x4.wt)) 

###Standardized estimates, se, and model-averaged estimates for B3.

comp.stdB3.x1 <- 0.0 
comp.stdB3.x1x2 <- 0.0 
comp.stdB3.x1x3 <- coef(std.zero.trunc.x1x3)[3] 
comp.stdB3.x1x4 <- 0.0 
comp.stdB3.x1x2x3 <- coef(std.zero.trunc.x1x2x3)[4] 
comp.stdB3.x1x2x4 <- 0.0 
comp.stdB3.x1x3x4 <- coef(std.zero.trunc.x1x3x4)[3] 
comp.stdB3.x1x2x3x4 <- coef(std.zero.trunc.x1x2x3x4)[4] 

comp.sestdB3.x1 <- 0.0 
comp.sestdB3.x1x2 <- 0.0 
comp.sestdB3.x1x3 <- coef(summary(std.zero.trunc.x1x3))[3,2] 
comp.sestdB3.x1x4 <- 0.0 
comp.sestdB3.x1x2x3 <- coef(summary(std.zero.trunc.x1x2x3))[4,2] 
comp.sestdB3.x1x2x4 <- 0.0 
comp.sestdB3.x1x3x4 <- coef(summary(std.zero.trunc.x1x3x4))[3,2] 
comp.sestdB3.x1x2x3x4 <- coef(summary(std.zero.trunc.x1x2x3x4))[4,2] 


comp.stdB3.wt <- ((comp.stdB3.x1 * AIC.x1.wt) + (comp.stdB3.x1x2 * AIC.x1x2.wt) + (comp.stdB3.x1x3 * AIC.x1x3.wt) +(comp.stdB3.x1x4 * AIC.x1x4.wt) + (comp.stdB3.x1x2x3 * AIC.x1x2x3.wt) + (comp.stdB3.x1x3x4 * AIC.x1x3x4.wt) + (comp.stdB3.x1x2x4 * AIC.x1x2x4.wt) + (comp.stdB3.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

comp.stdB3.wt.se <- (((((comp.stdB3.x1 - comp.stdB3.wt)^2 + comp.sestdB3.x1^2)^0.5) * AIC.x1.wt) + ((((comp.stdB3.x1x2 - comp.stdB3.wt)^2 + comp.sestdB3.x1x2^2)^0.5) * AIC.x1x2.wt) + ((((comp.stdB3.x1x3 - comp.stdB3.wt)^2 + comp.sestdB3.x1x3^2)^0.5) * AIC.x1x3.wt) + ((((comp.stdB3.x1x4 - comp.stdB3.wt)^2 + comp.sestdB3.x1x4^2)^0.5) * AIC.x1x4.wt) + ((((comp.stdB3.x1x2x3 - comp.stdB3.wt)^2 + comp.sestdB3.x1x2x3^2)^0.5) * AIC.x1x2x3.wt) + ((((comp.stdB3.x1x3x4 - comp.stdB3.wt)^2 + comp.sestdB3.x1x3x4^2)^0.5) * AIC.x1x3x4.wt) + ((((comp.stdB3.x1x2x4 - comp.stdB3.wt)^2 + comp.sestdB3.x1x2x4^2)^0.5) * AIC.x1x2x4.wt) + ((((comp.stdB3.x1x2x3x4 - comp.stdB3.wt)^2 + comp.sestdB3.x1x2x3x4^2)^0.5) * AIC.x1x2x3x4.wt)) 

###X4
###Estimates, se, and model-averaged estimates for B4.


comp.B4.x1 <- 0.0 
comp.B4.x1x2 <- 0.0 
comp.B4.x1x3 <- 0.0 
comp.B4.x1x4 <- coef(zero.trunc.x1x4)[3] 
comp.B4.x1x2x3 <- 0.0 
comp.B4.x1x2x4 <- coef(zero.trunc.x1x2x4)[4]  
comp.B4.x1x3x4 <- coef(zero.trunc.x1x3x4)[4]  
comp.B4.x1x2x3x4 <- coef(zero.trunc.x1x2x3x4)[5] 

comp.seB4.x1 <- 0.0 
comp.seB4.x1x2 <- 0.0
comp.seB4.x1x3 <- 0.0  
comp.seB4.x1x4 <- coef(summary(zero.trunc.x1x4))[3,2]   
comp.seB4.x1x2x3 <- 0.0
comp.seB4.x1x2x4 <- coef(summary(zero.trunc.x1x2x4))[4,2]   
comp.seB4.x1x3x4 <- coef(summary(zero.trunc.x1x3x4))[4,2]  
comp.seB4.x1x2x3x4 <- coef(summary(zero.trunc.x1x2x3x4))[5,2] 

comp.B4.wt <- ((comp.B4.x1 * AIC.x1.wt) + (comp.B4.x1x2 * AIC.x1x2.wt) + (comp.B4.x1x3 * AIC.x1x3.wt) +(comp.B4.x1x4 * AIC.x1x4.wt) + (comp.B4.x1x2x3 * AIC.x1x2x3.wt) + (comp.B4.x1x3x4 * AIC.x1x3x4.wt) + (comp.B4.x1x2x4 * AIC.x1x2x4.wt) + (comp.B4.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

comp.B4.wt.se <- (((((comp.B4.x1 - comp.B4.wt)^2 + comp.seB4.x1^2)^0.5) * AIC.x1.wt) + ((((comp.B4.x1x2 - comp.B4.wt)^2 + comp.seB4.x1x2^2)^0.5) * AIC.x1x2.wt) + ((((comp.B4.x1x3 - comp.B4.wt)^2 + comp.seB4.x1x3^2)^0.5) * AIC.x1x3.wt) + ((((comp.B4.x1x4 - comp.B4.wt)^2 + comp.seB4.x1x4^2)^0.5) * AIC.x1x4.wt) + ((((comp.B4.x1x2x3 - comp.B4.wt)^2 + comp.seB4.x1x2x3^2)^0.5) * AIC.x1x2x3.wt) + ((((comp.B4.x1x3x4 - comp.B4.wt)^2 + comp.seB4.x1x3x4^2)^0.5) * AIC.x1x3x4.wt) + ((((comp.B4.x1x2x4 - comp.B4.wt)^2 + comp.seB4.x1x2x4^2)^0.5) * AIC.x1x2x4.wt) + ((((comp.B4.x1x2x3x4 - comp.B4.wt)^2 + comp.seB4.x1x2x3x4^2)^0.5) * AIC.x1x2x3x4.wt)) 

###Standardized estimates, se, and model-averaged estimates for B4.

comp.stdB4.x1 <- 0.0
comp.stdB4.x1x2 <- 0.0 
comp.stdB4.x1x3 <- 0.0 
comp.stdB4.x1x4 <- coef(std.zero.trunc.x1x4)[3] 
comp.stdB4.x1x2x3 <- 0.0 
comp.stdB4.x1x2x4 <- coef(std.zero.trunc.x1x2x4)[4] 
comp.stdB4.x1x3x4 <- coef(std.zero.trunc.x1x3x4)[4] 
comp.stdB4.x1x2x3x4 <- coef(std.zero.trunc.x1x2x3x4)[5] 

comp.sestdB4.x1 <- 0.0 
comp.sestdB4.x1x2 <- 0.0 
comp.sestdB4.x1x3 <- 0.0 
comp.sestdB4.x1x4 <- coef(summary(std.zero.trunc.x1x4))[3,2] 
comp.sestdB4.x1x2x3 <- 0.0 
comp.sestdB4.x1x2x4 <- coef(summary(std.zero.trunc.x1x2x4))[4,2] 
comp.sestdB4.x1x3x4 <- coef(summary(std.zero.trunc.x1x3x4))[4,2] 
comp.sestdB4.x1x2x3x4 <- coef(summary(std.zero.trunc.x1x2x3x4))[5,2] 


comp.stdB4.wt <- ((comp.stdB4.x1 * AIC.x1.wt) + (comp.stdB4.x1x2 * AIC.x1x2.wt) + (comp.stdB4.x1x3 * AIC.x1x3.wt) +(comp.stdB4.x1x4 * AIC.x1x4.wt) + (comp.stdB4.x1x2x3 * AIC.x1x2x3.wt) + (comp.stdB4.x1x3x4 * AIC.x1x3x4.wt) + (comp.stdB4.x1x2x4 * AIC.x1x2x4.wt) + (comp.stdB4.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

comp.stdB4.wt.se <- (((((comp.stdB4.x1 - comp.stdB4.wt)^2 + comp.sestdB4.x1^2)^0.5) * AIC.x1.wt) + ((((comp.stdB4.x1x2 - comp.stdB4.wt)^2 + comp.sestdB4.x1x2^2)^0.5) * AIC.x1x2.wt) + ((((comp.stdB4.x1x3 - comp.stdB4.wt)^2 + comp.sestdB4.x1x3^2)^0.5) * AIC.x1x3.wt) + ((((comp.stdB4.x1x4 - comp.stdB4.wt)^2 + comp.sestdB4.x1x4^2)^0.5) * AIC.x1x4.wt) + ((((comp.stdB4.x1x2x3 - comp.stdB4.wt)^2 + comp.sestdB4.x1x2x3^2)^0.5) * AIC.x1x2x3.wt) + ((((comp.stdB4.x1x3x4 - comp.stdB4.wt)^2 + comp.sestdB4.x1x3x4^2)^0.5) * AIC.x1x3x4.wt) + ((((comp.stdB4.x1x2x4 - comp.stdB4.wt)^2 + comp.sestdB4.x1x2x4^2)^0.5) * AIC.x1x2x4.wt) + ((((comp.stdB4.x1x2x3x4 - comp.stdB4.wt)^2 + comp.sestdB4.x1x2x3x4^2)^0.5) * AIC.x1x2x3x4.wt)) 

###Compare model-averaged predictions versus predictions from model-
###averaged estimates.

correct.modavg.pred <- (predictvglm(zero.trunc.x1,type="response"))*AIC.x1.wt + 
(predictvglm(zero.trunc.x1x2,type="response"))*AIC.x1x2.wt + 
(predictvglm(zero.trunc.x1x3,type="response"))*AIC.x1x3.wt + 
(predictvglm(zero.trunc.x1x4,type="response"))*AIC.x1x4.wt + 
(predictvglm(zero.trunc.x1x3x4,type="response"))*AIC.x1x3x4.wt + 
(predictvglm(zero.trunc.x1x2x4,type="response"))*AIC.x1x2x4.wt + 
(predictvglm(zero.trunc.x1x2x3,type="response"))*AIC.x1x2x3.wt + 
(predictvglm(zero.trunc.x1x2x3x4,type="response"))*AIC.x1x2x3x4.wt
 
###Need model-averaged estimates for B0

comp.B0.x1 <- coef(zero.trunc.x1)[1] 
comp.B0.x1x2 <- coef(zero.trunc.x1x2)[1] 
comp.B0.x1x3 <- coef(zero.trunc.x1x3)[1] 
comp.B0.x1x4 <- coef(zero.trunc.x1x4)[1] 
comp.B0.x1x2x3 <- coef(zero.trunc.x1x2x3)[1] 
comp.B0.x1x2x4 <- coef(zero.trunc.x1x2x4)[1] 
comp.B0.x1x3x4 <- coef(zero.trunc.x1x3x4)[1] 
comp.B0.x1x2x3x4 <- coef(zero.trunc.x1x2x3x4)[1] 

comp.seB0.x1 <- coef(summary(zero.trunc.x1))[1,2] 
comp.seB0.x1x2 <- coef(summary(zero.trunc.x1x2))[1,2] 
comp.seB0.x1x3 <- coef(summary(zero.trunc.x1x3))[1,2] 
comp.seB0.x1x4 <- coef(summary(zero.trunc.x1x4))[1,2] 
comp.seB0.x1x2x3 <- coef(summary(zero.trunc.x1x2x3))[1,2] 
comp.seB0.x1x2x4 <- coef(summary(zero.trunc.x1x2x4))[1,2] 
comp.seB0.x1x3x4 <- coef(summary(zero.trunc.x1x3x4))[1,2] 
comp.seB0.x1x2x3x4 <- coef(summary(zero.trunc.x1x2x3x4))[1,2] 


comp.B0.wt <- ((comp.B0.x1 * AIC.x1.wt) + (comp.B0.x1x2 * AIC.x1x2.wt) + (comp.B0.x1x3 * AIC.x1x3.wt) +(comp.B0.x1x4 * AIC.x1x4.wt) + (comp.B0.x1x2x3 * AIC.x1x2x3.wt) + (comp.B0.x1x3x4 * AIC.x1x3x4.wt) + (comp.B0.x1x2x4 * AIC.x1x2x4.wt) + (comp.B0.x1x2x3x4 * AIC.x1x2x3x4.wt)) 

logmu <- comp.B0.wt + comp.B1.wt*x1.zt + comp.B2.wt*x2.zt + comp.B3.wt*x3.zt + comp.B4.wt*x4.zt 
 
incorrect.modavg.pred <- (exp(logmu))/(1.0 - exp(- exp(logmu)))

summary(correct.modavg.pred - incorrect.modavg.pred)

###To see how changing magnitude of estimates by altering AIC weights 
###impacts differences in prediction approaches.

correct.modavg.pred.eqwt <- (predictvglm(zero.trunc.x1,type="response"))*(1/8) + 
(predictvglm(zero.trunc.x1x2,type="response"))*(1/8) + 
(predictvglm(zero.trunc.x1x3,type="response"))*(1/8) + 
(predictvglm(zero.trunc.x1x4,type="response"))*(1/8) + 
(predictvglm(zero.trunc.x1x3x4,type="response"))*(1/8) + 
(predictvglm(zero.trunc.x1x2x4,type="response"))*(1/8) + 
(predictvglm(zero.trunc.x1x2x3,type="response"))*(1/8) + 
(predictvglm(zero.trunc.x1x2x3x4,type="response"))*(1/8)
 

comp.B0.eqwt <- ((comp.B0.x1 * (1/8)) + (comp.B0.x1x2 * (1/8)) + (comp.B0.x1x3 * (1/8)) +(comp.B0.x1x4 * (1/8)) + (comp.B0.x1x2x3 * (1/8)) + (comp.B0.x1x3x4 * (1/8)) + (comp.B0.x1x2x4 * (1/8)) + (comp.B0.x1x2x3x4 * (1/8))) 

comp.B1.eqwt <- ((comp.B1.x1 * (1/8)) + (comp.B1.x1x2 * (1/8)) + (comp.B1.x1x3 * (1/8)) +(comp.B1.x1x4 * (1/8)) + (comp.B1.x1x2x3 * (1/8)) + (comp.B1.x1x3x4 * (1/8)) + (comp.B1.x1x2x4 * (1/8)) + (comp.B1.x1x2x3x4 * (1/8))) 

comp.B2.eqwt <- ((comp.B2.x1 * (1/8)) + (comp.B2.x1x2 * (1/8)) + (comp.B2.x1x3 * (1/8)) +(comp.B2.x1x4 * (1/8)) + (comp.B2.x1x2x3 * (1/8)) + (comp.B2.x1x3x4 * (1/8)) + (comp.B2.x1x2x4 * (1/8)) + (comp.B2.x1x2x3x4 * (1/8))) 

comp.B3.eqwt <- ((comp.B3.x1 * (1/8)) + (comp.B3.x1x2 * (1/8)) + (comp.B3.x1x3 * (1/8)) +(comp.B3.x1x4 * (1/8)) + (comp.B3.x1x2x3 * (1/8)) + (comp.B3.x1x3x4 * (1/8)) + (comp.B3.x1x2x4 * (1/8)) + (comp.B3.x1x2x3x4 * (1/8))) 

comp.B4.eqwt <- ((comp.B4.x1 * (1/8)) + (comp.B4.x1x2 * (1/8)) + (comp.B4.x1x3 * (1/8)) +(comp.B4.x1x4 * (1/8)) + (comp.B4.x1x2x3 * (1/8)) + (comp.B4.x1x3x4 * (1/8)) + (comp.B4.x1x2x4 * (1/8)) + (comp.B4.x1x2x3x4 * (1/8))) 


logmu.eqwt <- comp.B0.eqwt + comp.B1.eqwt*x1.zt + comp.B2.eqwt*x2.zt + comp.B3.eqwt*x3.zt + comp.B4.eqwt*x4.zt 
 
incorrect.modavg.pred.eqwt <- (exp(logmu.eqwt))/(1.0 - exp(- exp(logmu.eqwt)))

summary(correct.modavg.pred.eqwt - incorrect.modavg.pred.eqwt)

###Predictions at extremes
newData <- cbind(0.90,0.08,0.001,0.001)
colnames(newData) <- c("x1.zt","x2.zt","x3.zt","x4.zt")
newData <- as.data.frame(newData)


correct.modavg.pred.extreme <- (predictvglm(zero.trunc.x1,newdata=newData,type="response"))*AIC.x1.wt + 
(predictvglm(zero.trunc.x1x2,newdata=newData,type="response"))*AIC.x1x2.wt + 
(predictvglm(zero.trunc.x1x3,newdata=newData,type="response"))*AIC.x1x3.wt + 
(predictvglm(zero.trunc.x1x4,newdata=newData,type="response"))*AIC.x1x4.wt + 
(predictvglm(zero.trunc.x1x3x4,newdata=newData,type="response"))*AIC.x1x3x4.wt + 
(predictvglm(zero.trunc.x1x2x4,newdata=newData,type="response"))*AIC.x1x2x4.wt + 
(predictvglm(zero.trunc.x1x2x3,newdata=newData,type="response"))*AIC.x1x2x3.wt + 
(predictvglm(zero.trunc.x1x2x3x4,newdata=newData,type="response"))*AIC.x1x2x3x4.wt

logmu.extreme <- comp.B0.wt + comp.B1.wt*newData[1,1] + comp.B2.wt*newData[1,2] + comp.B3.wt*newData[1,3] + comp.B4.wt*newData[1,4] 
 
incorrect.modavg.pred.extreme <- (exp(logmu.extreme))/(1.0 - exp(- exp(logmu.extreme)))

summary(correct.modavg.pred.extreme - incorrect.modavg.pred.extreme)

###Test to compare lack of equivalence between model-averaged 
###predictions and predictions from model-averaged coefficients for 
###models that are nonlinear functions of parameters.

x1.new <- 0.80
x2.new <- 0.05
B01 <- -5
B02 <- -4
B03 <- -6

B11 <- 6.0
B12 <- 7.0
B13 <- 5.5

B21 <- 16.0
B22 <- 17.0
B23 <- 13.0

modavg.pred <- mean(c((exp(B01 + B11*x1.new + B21*x2.new)),(exp(B02 + B12*x1.new + B22*x2.new)),(exp(B03 + B13*x1.new + B23*x2.new))))

avgB0 <- mean(c(B01,B02,B03))
avgB1 <- mean(c(B11,B12,B13))
avgB2 <- mean(c(B21,B22,B23))


pred.modavg <- exp(avgB0 + avgB1*x1.new + avgB2*x2.new)

modavg.pred - pred.modavg

