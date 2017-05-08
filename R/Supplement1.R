# Supplement 1.  Cade, B.S. 2015.  Model averaging and muddled multimodel inferences.  Ecology.
# 
# 
# The following R script reads in the college GPA example data from Burnham and Anderson 2002:226, estimates the 16 linear least squares regression models, computes AICc, AICc weights, variance inflation factors, partial standard deviations of predictors, standardizes estimates by partial standard deviations, computes model-averaged standardized estimates and their standard errors, and computes the model-averaged ratio of t-statistics for unstandardized estimates (equivalent to model-averaged ratio of standardized estimates).  The code is written to be transparent with respect to the mathematical operations rather than for efficiency.
# 
# Literature Cited
# Burnham, K.  P., and D. R. Anderson.  2002.  Model selection and multimodel inference: a practical information-theoretic approach.  Springer-Verlag, New York, NY, USA.

###Create GPA data set in B&A 2002(226)


gpa.y <- c(1.97,2.74,2.19,2.60,2.98,1.65,1.89,2.38,2.66,1.96,3.14,1.96,2.20,3.90,2.02,3.61,3.07,2.63,3.11,3.20)

sat.math.x1 <- c(321,718,358,403,640,237,270,418,443,359,669,409,582,750,451,645,791,521,594,653)

sat.verb.x2 <- c(247,436,578,447,563,342,472,356,327,385,664,518,364,632,435,704,341,483,665,606)

hs.math.x3 <- c(2.30,3.80,2.98,3.58,3.38,1.48,1.67,3.73,3.09,1.54,3.21,2.77,1.47,3.14,1.54,3.50,3.20,3.59,3.42,3.69)

hs.engl.x4 <- c(2.63,3.57,2.57,2.21,3.48,2.14,2.64,2.52,3.20,3.46,3.37,2.60,2.90,3.49,3.20,3.74,2.93,3.32,2.70,3.52)

cor(cbind(gpa.y,sat.math.x1,sat.verb.x2,hs.math.x3,hs.engl.x4))


gpa.01 <- lm(gpa.y ~ sat.math.x1)
gpa.02 <- lm(gpa.y ~ sat.verb.x2)
gpa.03 <- lm(gpa.y ~ hs.math.x3)
gpa.04 <- lm(gpa.y ~ hs.engl.x4)
gpa.05 <- lm(gpa.y ~ sat.math.x1 + sat.verb.x2)
gpa.06 <- lm(gpa.y ~ sat.math.x1 + hs.math.x3)
gpa.07 <- lm(gpa.y ~ sat.verb.x2 + hs.math.x3)
gpa.08 <- lm(gpa.y ~ sat.math.x1 + hs.engl.x4)
gpa.09 <- lm(gpa.y ~ sat.verb.x2 + hs.engl.x4)
gpa.10 <- lm(gpa.y ~ hs.math.x3 + hs.engl.x4)
gpa.11 <- lm(gpa.y ~ sat.math.x1 + sat.verb.x2 + hs.math.x3)
gpa.12 <- lm(gpa.y ~ sat.math.x1 + sat.verb.x2 + hs.engl.x4)
gpa.13 <- lm(gpa.y ~ sat.math.x1 + hs.math.x3 + hs.engl.x4)
gpa.14 <- lm(gpa.y ~ sat.verb.x2 + hs.math.x3 + hs.engl.x4)
gpa.15 <- lm(gpa.y ~ sat.math.x1 + sat.verb.x2 + hs.math.x3 + hs.engl.x4)
gpa.16 <- lm(gpa.y ~ 1)


###Compute AICc for each of the 16 models

gpa.11.AICc <- AIC(gpa.11) + ((2*(length(coef(gpa.11))+1))*((length(coef(gpa.11))+1) + 1))/(nobs(gpa.11) - (length(coef(gpa.11))+1) - 1)

gpa.05.AICc <- AIC(gpa.05) + ((2*(length(coef(gpa.05))+1))*((length(coef(gpa.05))+1) + 1))/(nobs(gpa.05) - (length(coef(gpa.05))+1) - 1)

gpa.06.AICc <- AIC(gpa.06) + ((2*(length(coef(gpa.06))+1))*((length(coef(gpa.06))+1) + 1))/(nobs(gpa.06) - (length(coef(gpa.06))+1) - 1)

gpa.15.AICc <- AIC(gpa.15) + ((2*(length(coef(gpa.15))+1))*((length(coef(gpa.15))+1) + 1))/(nobs(gpa.15) - (length(coef(gpa.15))+1) - 1)

gpa.12.AICc <- AIC(gpa.12) + ((2*(length(coef(gpa.12))+1))*((length(coef(gpa.12))+1) + 1))/(nobs(gpa.12) - (length(coef(gpa.12))+1) - 1)

gpa.13.AICc <- AIC(gpa.13) + ((2*(length(coef(gpa.13))+1))*((length(coef(gpa.13))+1) + 1))/(nobs(gpa.13) - (length(coef(gpa.13))+1) - 1)

gpa.01.AICc <- AIC(gpa.01) + ((2*(length(coef(gpa.01))+1))*((length(coef(gpa.01))+1) + 1))/(nobs(gpa.01) - (length(coef(gpa.01))+1) - 1)

gpa.08.AICc <- AIC(gpa.08) + ((2*(length(coef(gpa.08))+1))*((length(coef(gpa.08))+1) + 1))/(nobs(gpa.08) - (length(coef(gpa.08))+1) - 1)

gpa.02.AICc <- AIC(gpa.02) + ((2*(length(coef(gpa.02))+1))*((length(coef(gpa.02))+1) + 1))/(nobs(gpa.02) - (length(coef(gpa.02))+1) - 1)

gpa.03.AICc <- AIC(gpa.03) + ((2*(length(coef(gpa.03))+1))*((length(coef(gpa.03))+1) + 1))/(nobs(gpa.03) - (length(coef(gpa.03))+1) - 1)

gpa.04.AICc <- AIC(gpa.04) + ((2*(length(coef(gpa.04))+1))*((length(coef(gpa.04))+1) + 1))/(nobs(gpa.04) - (length(coef(gpa.04))+1) - 1)

gpa.07.AICc <- AIC(gpa.07) + ((2*(length(coef(gpa.07))+1))*((length(coef(gpa.07))+1) + 1))/(nobs(gpa.07) - (length(coef(gpa.07))+1) - 1)

gpa.09.AICc <- AIC(gpa.09) + ((2*(length(coef(gpa.09))+1))*((length(coef(gpa.09))+1) + 1))/(nobs(gpa.09) - (length(coef(gpa.09))+1) - 1)

gpa.10.AICc <- AIC(gpa.10) + ((2*(length(coef(gpa.10))+1))*((length(coef(gpa.10))+1) + 1))/(nobs(gpa.10) - (length(coef(gpa.10))+1) - 1)

gpa.14.AICc <- AIC(gpa.14) + ((2*(length(coef(gpa.14))+1))*((length(coef(gpa.14))+1) + 1))/(nobs(gpa.14) - (length(coef(gpa.14))+1) - 1)

gpa.16.AICc <- AIC(gpa.16) + ((2*(length(coef(gpa.16))+1))*((length(coef(gpa.16))+1) + 1))/(nobs(gpa.16) - (length(coef(gpa.16))+1) - 1)


###Compute minimum AICc

minAICc <- min(c(gpa.11.AICc,gpa.05.AICc,gpa.06.AICc, gpa.15.AICc,gpa.12.AICc,gpa.13.AICc,gpa.01.AICc,gpa.08.AICc,gpa.02.AICc,gpa.03.AICc,gpa.04.AICc,gpa.07.AICc,gpa.09.AICc,gpa.10.AICc,gpa.14.AICc,gpa.16.AICc)) 

###Compute delta AICc

d.gpa.11.AICc <- gpa.11.AICc - minAICc
d.gpa.05.AICc <- gpa.05.AICc - minAICc
d.gpa.06.AICc <- gpa.06.AICc - minAICc
d.gpa.15.AICc <- gpa.15.AICc - minAICc
d.gpa.12.AICc <- gpa.12.AICc - minAICc
d.gpa.13.AICc <- gpa.13.AICc - minAICc
d.gpa.01.AICc <- gpa.01.AICc - minAICc
d.gpa.08.AICc <- gpa.08.AICc - minAICc
d.gpa.02.AICc <- gpa.02.AICc - minAICc
d.gpa.03.AICc <- gpa.03.AICc - minAICc
d.gpa.04.AICc <- gpa.04.AICc - minAICc
d.gpa.07.AICc <- gpa.07.AICc - minAICc
d.gpa.09.AICc <- gpa.09.AICc - minAICc
d.gpa.10.AICc <- gpa.10.AICc - minAICc
d.gpa.14.AICc <- gpa.14.AICc - minAICc
d.gpa.16.AICc <- gpa.16.AICc - minAICc


###Compute sum of the AICc weights

sumAICcWts <- exp(-0.5*d.gpa.11.AICc) + exp(-0.5*d.gpa.05.AICc)+ exp(-0.5*d.gpa.06.AICc) + exp(-0.5*d.gpa.15.AICc) + exp(-0.5*d.gpa.12.AICc) + exp(-0.5*d.gpa.13.AICc) + exp(-0.5*d.gpa.01.AICc) + exp(-0.5*d.gpa.08.AICc) + exp(-0.5*d.gpa.02.AICc) + exp(-0.5*d.gpa.03.AICc) + exp(-0.5*d.gpa.04.AICc) + exp(-0.5*d.gpa.07.AICc) + exp(-0.5*d.gpa.09.AICc) + exp(-0.5*d.gpa.10.AICc) + exp(-0.5*d.gpa.14.AICc) + exp(-0.5*d.gpa.16.AICc) 

###AICc weights for each model

gpa.11.wt <- exp(-0.5*d.gpa.11.AICc)/sumAICcWts
gpa.05.wt <- exp(-0.5*d.gpa.05.AICc)/sumAICcWts
gpa.06.wt <- exp(-0.5*d.gpa.06.AICc)/sumAICcWts
gpa.15.wt <- exp(-0.5*d.gpa.15.AICc)/sumAICcWts
gpa.12.wt <- exp(-0.5*d.gpa.12.AICc)/sumAICcWts
gpa.13.wt <- exp(-0.5*d.gpa.13.AICc)/sumAICcWts
gpa.01.wt <- exp(-0.5*d.gpa.01.AICc)/sumAICcWts
gpa.08.wt <- exp(-0.5*d.gpa.08.AICc)/sumAICcWts
gpa.02.wt <- exp(-0.5*d.gpa.02.AICc)/sumAICcWts
gpa.03.wt <- exp(-0.5*d.gpa.03.AICc)/sumAICcWts
gpa.04.wt <- exp(-0.5*d.gpa.04.AICc)/sumAICcWts
gpa.07.wt <- exp(-0.5*d.gpa.07.AICc)/sumAICcWts
gpa.09.wt <- exp(-0.5*d.gpa.09.AICc)/sumAICcWts
gpa.10.wt <- exp(-0.5*d.gpa.10.AICc)/sumAICcWts
gpa.14.wt <- exp(-0.5*d.gpa.14.AICc)/sumAICcWts
gpa.16.wt <- exp(-0.5*d.gpa.16.AICc)/sumAICcWts


###Variance inflation factors (vif) by model
###Inverse of VIF is proportion of variance in x_i that is not linearly ###related to other predictors in model

library(car)

###Models with sat.math.x1

ivif.x1.11 <- (vif(gpa.11)[1])^-1
ivif.x1.05 <- (vif(gpa.05)[1])^-1
ivif.x1.06 <- (vif(gpa.06)[1])^-1
ivif.x1.15 <- (vif(gpa.15)[1])^-1
ivif.x1.12 <- (vif(gpa.12)[1])^-1
ivif.x1.13 <- (vif(gpa.13)[1])^-1
ivif.x1.08 <- (vif(gpa.08)[1])^-1

###Partial standard deviations following Bring (1994).

partSD.x1.11 <- ivif.x1.11^0.5 * ((nobs(gpa.11)-1)/(nobs(gpa.11)- (length(coef(gpa.11)) -1)))^0.5 * sd(sat.math.x1)

partSD.x1.05 <- ivif.x1.05^0.5 * ((nobs(gpa.05)-1)/(nobs(gpa.05)- (length(coef(gpa.05)) -1)))^0.5 * sd(sat.math.x1)

partSD.x1.06 <- ivif.x1.06^0.5 * ((nobs(gpa.06)-1)/(nobs(gpa.06)- (length(coef(gpa.06)) -1)))^0.5 * sd(sat.math.x1)

partSD.x1.15 <- ivif.x1.15^0.5 * ((nobs(gpa.15)-1)/(nobs(gpa.15)- (length(coef(gpa.15)) -1)))^0.5 * sd(sat.math.x1)

partSD.x1.12 <- ivif.x1.12^0.5 * ((nobs(gpa.12)-1)/(nobs(gpa.12)- (length(coef(gpa.12)) -1)))^0.5 * sd(sat.math.x1)

partSD.x1.13 <- ivif.x1.13^0.5 * ((nobs(gpa.13)-1)/(nobs(gpa.13)- (length(coef(gpa.13)) -1)))^0.5 * sd(sat.math.x1)

partSD.x1.08 <- ivif.x1.08^0.5 * ((nobs(gpa.08)-1)/(nobs(gpa.08)- (length(coef(gpa.08)) -1)))^0.5 * sd(sat.math.x1)

###Estimates of standardized estimates for sat.math.x1 using partial ###SDs. Note that partial SD when using single predictor is equivalent ###to SD of that predictor.

stdB1.11 <- coef(gpa.11)[2] * partSD.x1.11
stdB1.05 <- coef(gpa.05)[2] * partSD.x1.05
stdB1.06 <- coef(gpa.06)[2] * partSD.x1.06
stdB1.15 <- coef(gpa.15)[2] * partSD.x1.15
stdB1.12 <- coef(gpa.12)[2] * partSD.x1.12
stdB1.13 <- coef(gpa.13)[2] * partSD.x1.13
stdB1.01 <- coef(gpa.01)[2] * sd(sat.math.x1)
stdB1.08 <- coef(gpa.08)[2] * partSD.x1.08

###Estimates of standard errors for standardized estimates of ###sat.math.x1.

stdB1.11.se <- ((coef(summary((gpa.11)))[2,2])^2 * partSD.x1.11^2)^0.5
stdB1.05.se <- ((coef(summary((gpa.05)))[2,2])^2 * partSD.x1.05^2)^0.5
stdB1.06.se <- ((coef(summary((gpa.06)))[2,2])^2 * partSD.x1.06^2)^0.5
stdB1.15.se <- ((coef(summary((gpa.15)))[2,2])^2 * partSD.x1.15^2)^0.5
stdB1.12.se <- ((coef(summary((gpa.12)))[2,2])^2 * partSD.x1.12^2)^0.5
stdB1.13.se <- ((coef(summary((gpa.13)))[2,2])^2 * partSD.x1.13^2)^0.5
stdB1.01.se <- ((coef(summary((gpa.01)))[2,2])^2 * (sd(sat.math.x1))^2)^0.5
stdB1.08.se <- ((coef(summary((gpa.08)))[2,2])^2 * partSD.x1.08^2)^0.5

###Renormalize weights for models with sat.math.x1

sumB1.wt<- sum(gpa.11.wt,gpa.05.wt,gpa.06.wt,gpa.15.wt,gpa.12.wt,gpa.13.wt,gpa.01.wt,gpa.08.wt)

###Model-averaged standardized estimates

stdB1.wt <- ((stdB1.11 * gpa.11.wt) + (stdB1.05 * gpa.05.wt) + (stdB1.06 * gpa.06.wt) + (stdB1.15 * gpa.15.wt) + (stdB1.12 * gpa.12.wt) + (stdB1.13 * gpa.13.wt) + (stdB1.01 * gpa.01.wt) + (stdB1.08 * gpa.08.wt))/sumB1.wt 

###Standard errors for model-averaged standardized estimates

stdB1.wt.se <- (((((stdB1.11 - stdB1.wt)^2 + stdB1.11.se^2)^0.5) * gpa.11.wt/sumB1.wt) + ((((stdB1.05 - stdB1.wt)^2 + stdB1.05.se^2)^0.5) * gpa.05.wt/sumB1.wt) + ((((stdB1.06 - stdB1.wt)^2 + stdB1.06.se^2)^0.5) * gpa.06.wt/sumB1.wt) + ((((stdB1.15 - stdB1.wt)^2 + stdB1.15.se^2)^0.5) * gpa.15.wt/sumB1.wt) + ((((stdB1.12 - stdB1.wt)^2 + stdB1.12.se^2)^0.5) * gpa.12.wt/sumB1.wt) + ((((stdB1.13 - stdB1.wt)^2 + stdB1.13.se^2)^0.5) * gpa.13.wt/sumB1.wt) + ((((stdB1.01 - stdB1.wt)^2 + stdB1.01.se^2)^0.5) * gpa.01.wt/sumB1.wt) + ((((stdB1.08 - stdB1.wt)^2 + stdB1.08.se^2)^0.5) * gpa.08.wt/sumB1.wt)) 


###Models with sat.verb.x2

###Variance inflation factors

library(car)

ivif.x2.11 <- (vif(gpa.11)[2])^-1
ivif.x2.05 <- (vif(gpa.05)[2])^-1
ivif.x2.14 <- (vif(gpa.14)[1])^-1
ivif.x2.15 <- (vif(gpa.15)[2])^-1
ivif.x2.12 <- (vif(gpa.12)[2])^-1
ivif.x2.07 <- (vif(gpa.07)[1])^-1
ivif.x2.09 <- (vif(gpa.09)[1])^-1

###Partial standard deviations following Bring (1994).

partSD.x2.11 <- ivif.x2.11^0.5 * ((nobs(gpa.11)-1)/(nobs(gpa.11)- (length(coef(gpa.11)) -1)))^0.5 * sd(sat.verb.x2)

partSD.x2.05 <- ivif.x2.05^0.5 * ((nobs(gpa.05)-1)/(nobs(gpa.05)- (length(coef(gpa.05)) -1)))^0.5 * sd(sat.verb.x2)

partSD.x2.14 <- ivif.x2.14^0.5 * ((nobs(gpa.14)-1)/(nobs(gpa.14)- (length(coef(gpa.14)) -1)))^0.5 * sd(sat.verb.x2)

partSD.x2.15 <- ivif.x2.15^0.5 * ((nobs(gpa.15)-1)/(nobs(gpa.15)- (length(coef(gpa.15)) -1)))^0.5 * sd(sat.verb.x2)

partSD.x2.12 <- ivif.x2.12^0.5 * ((nobs(gpa.12)-1)/(nobs(gpa.12)- (length(coef(gpa.12)) -1)))^0.5 * sd(sat.verb.x2)

partSD.x2.07 <- ivif.x2.07^0.5 * ((nobs(gpa.07)-1)/(nobs(gpa.07)- (length(coef(gpa.07)) -1)))^0.5 * sd(sat.verb.x2)

partSD.x2.09 <- ivif.x2.09^0.5 * ((nobs(gpa.09)-1)/(nobs(gpa.09)- (length(coef(gpa.09)) -1)))^0.5 * sd(sat.verb.x2)

###Estimates of standardized estimates for sat.verb.x2 using partial SDs

stdB2.11 <- coef(gpa.11)[3] * partSD.x2.11
stdB2.05 <- coef(gpa.05)[3] * partSD.x2.05
stdB2.14 <- coef(gpa.14)[2] * partSD.x2.14
stdB2.15 <- coef(gpa.15)[3] * partSD.x2.15
stdB2.12 <- coef(gpa.12)[3] * partSD.x2.12
stdB2.07 <- coef(gpa.07)[2] * partSD.x2.07
stdB2.09 <- coef(gpa.09)[2] * partSD.x2.09
stdB2.02 <- coef(gpa.02)[2] * sd(sat.verb.x2)


###Estimates of standard errors for standardized estimates of ###sat.verb.x2.

stdB2.11.se <- ((coef(summary((gpa.11)))[3,2])^2 * partSD.x2.11^2)^0.5
stdB2.05.se <- ((coef(summary((gpa.05)))[3,2])^2 * partSD.x2.05^2)^0.5
stdB2.14.se <- ((coef(summary((gpa.14)))[2,2])^2 * partSD.x2.14^2)^0.5
stdB2.15.se <- ((coef(summary((gpa.15)))[3,2])^2 * partSD.x2.15^2)^0.5
stdB2.12.se <- ((coef(summary((gpa.12)))[3,2])^2 * partSD.x2.12^2)^0.5
stdB2.07.se <- ((coef(summary((gpa.07)))[2,2])^2 * partSD.x2.07^2)^0.5
stdB2.09.se <- ((coef(summary((gpa.09)))[2,2])^2 * partSD.x2.09^2)^0.5
stdB2.02.se <- ((coef(summary((gpa.02)))[2,2])^2 * (sd(sat.verb.x2))^2)^0.5

###Renormalize weights for models with sat.verb.x2.
 
sumB2.wt<- sum(gpa.11.wt,gpa.05.wt,gpa.14.wt,gpa.15.wt,gpa.12.wt,gpa.07.wt,gpa.09.wt,gpa.02.wt)

###Model-averaged standardized estimates.

stdB2.wt <- ((stdB2.11 * gpa.11.wt) + (stdB2.05 * gpa.05.wt) + (stdB2.14 * gpa.14.wt) + (stdB2.15 * gpa.15.wt) + (stdB2.12 * gpa.12.wt) + (stdB2.07 * gpa.07.wt) + (stdB2.09 * gpa.09.wt) + (stdB2.02 * gpa.02.wt))/sumB2.wt 

###Standard errors for model-averaged standardized estimates.

stdB2.wt.se <- (((((stdB2.11 - stdB2.wt)^2 + stdB2.11.se^2)^0.5) * gpa.11.wt/sumB2.wt) + ((((stdB2.05 - stdB2.wt)^2 + stdB2.05.se^2)^0.5) * gpa.05.wt/sumB2.wt) + ((((stdB2.14 - stdB2.wt)^2 + stdB2.14.se^2)^0.5) * gpa.14.wt/sumB2.wt) + ((((stdB2.15 - stdB2.wt)^2 + stdB2.15.se^2)^0.5) * gpa.15.wt/sumB2.wt) + ((((stdB2.12 - stdB2.wt)^2 + stdB2.12.se^2)^0.5) * gpa.12.wt/sumB2.wt) + ((((stdB2.07 - stdB2.wt)^2 + stdB2.07.se^2)^0.5) * gpa.07.wt/sumB2.wt) + ((((stdB2.09 - stdB2.wt)^2 + stdB2.09.se^2)^0.5) * gpa.09.wt/sumB2.wt) + ((((stdB2.02 - stdB2.wt)^2 + stdB2.02.se^2)^0.5) * gpa.02.wt/sumB2.wt)) 

###Models with hs.math.x3

###Variance inflation factors by model.

library(car)

ivif.x3.11 <- (vif(gpa.11)[3])^-1
ivif.x3.07 <- (vif(gpa.07)[2])^-1
ivif.x3.06 <- (vif(gpa.06)[2])^-1
ivif.x3.15 <- (vif(gpa.15)[3])^-1
ivif.x3.13 <- (vif(gpa.13)[2])^-1
ivif.x3.14 <- (vif(gpa.14)[2])^-1
ivif.x3.10 <- (vif(gpa.10)[1])^-1

###Partial standard deviations following Bring (1994).

partSD.x3.11 <- ivif.x3.11^0.5 * ((nobs(gpa.11)-1)/(nobs(gpa.11)- (length(coef(gpa.11)) -1)))^0.5 * sd(hs.math.x3)

partSD.x3.07 <- ivif.x3.07^0.5 * ((nobs(gpa.07)-1)/(nobs(gpa.07)- (length(coef(gpa.07)) -1)))^0.5 * sd(hs.math.x3)

partSD.x3.06 <- ivif.x3.06^0.5 * ((nobs(gpa.06)-1)/(nobs(gpa.06)- (length(coef(gpa.06)) -1)))^0.5 * sd(hs.math.x3)

partSD.x3.15 <- ivif.x3.15^0.5 * ((nobs(gpa.15)-1)/(nobs(gpa.15)- (length(coef(gpa.15)) -1)))^0.5 * sd(hs.math.x3)

partSD.x3.13 <- ivif.x3.13^0.5 * ((nobs(gpa.13)-1)/(nobs(gpa.13)- (length(coef(gpa.13)) -1)))^0.5 * sd(hs.math.x3)

partSD.x3.14 <- ivif.x3.14^0.5 * ((nobs(gpa.14)-1)/(nobs(gpa.14)- (length(coef(gpa.14)) -1)))^0.5 * sd(hs.math.x3)

partSD.x3.10 <- ivif.x3.10^0.5 * ((nobs(gpa.10)-1)/(nobs(gpa.10)- (length(coef(gpa.10)) -1)))^0.5 * sd(hs.math.x3)


###Standardized estimates for hs.math.x3 using partial SDs.

stdB3.11 <- coef(gpa.11)[4] * partSD.x3.11
stdB3.07 <- coef(gpa.07)[3] * partSD.x3.07
stdB3.06 <- coef(gpa.06)[3] * partSD.x3.06
stdB3.15 <- coef(gpa.15)[4] * partSD.x3.15
stdB3.03 <- coef(gpa.03)[2] * sd(hs.math.x3)
stdB3.13 <- coef(gpa.13)[3] * partSD.x3.13
stdB3.14 <- coef(gpa.14)[3] * partSD.x3.14
stdB3.10 <- coef(gpa.10)[2] * partSD.x3.10


###Estimates of standard errors for standardized estimates of ###hs.math.x3.

stdB3.11.se <- ((coef(summary((gpa.11)))[4,2])^2 * partSD.x3.11^2)^0.5
stdB3.07.se <- ((coef(summary((gpa.07)))[3,2])^2 * partSD.x3.07^2)^0.5
stdB3.06.se <- ((coef(summary((gpa.06)))[3,2])^2 * partSD.x3.06^2)^0.5
stdB3.15.se <- ((coef(summary((gpa.15)))[4,2])^2 * partSD.x3.15^2)^0.5
stdB3.03.se <- ((coef(summary((gpa.03)))[2,2])^2 * (sd(hs.math.x3))^2)^0.5
stdB3.13.se <- ((coef(summary((gpa.13)))[3,2])^2 * partSD.x3.13^2)^0.5
stdB3.14.se <- ((coef(summary((gpa.14)))[3,2])^2 * partSD.x3.14^2)^0.5
stdB3.10.se <- ((coef(summary((gpa.10)))[2,2])^2 * partSD.x3.10^2)^0.5


###Renormalize weights for models with hs.math.x3

sumB3.wt<- sum(gpa.11.wt,gpa.07.wt,gpa.06.wt,gpa.15.wt,gpa.03.wt,gpa.13.wt,gpa.14.wt,gpa.10.wt)

###Model-averaged standardized estimates.

stdB3.wt <- ((stdB3.11 * gpa.11.wt) + (stdB3.07 * gpa.07.wt) + (stdB3.06 * gpa.06.wt) + (stdB3.15 * gpa.15.wt) + (stdB3.03 * gpa.03.wt) + (stdB3.13 * gpa.13.wt) + (stdB3.14 * gpa.14.wt) + (stdB3.10 * gpa.10.wt))/sumB3.wt 

###Standard errors for model-averaged standardized estimates.

stdB3.wt.se <- (((((stdB3.11 - stdB3.wt)^2 + stdB3.11.se^2)^0.5) * gpa.11.wt/sumB3.wt) + ((((stdB3.07 - stdB3.wt)^2 + stdB3.07.se^2)^0.5) * gpa.07.wt/sumB3.wt) + ((((stdB3.06 - stdB3.wt)^2 + stdB3.06.se^2)^0.5) * gpa.06.wt/sumB3.wt) + ((((stdB3.15 - stdB3.wt)^2 + stdB3.15.se^2)^0.5) * gpa.15.wt/sumB3.wt) + ((((stdB3.03 - stdB3.wt)^2 + stdB3.03.se^2)^0.5) * gpa.03.wt/sumB3.wt) + ((((stdB3.13 - stdB3.wt)^2 + stdB3.13.se^2)^0.5) * gpa.13.wt/sumB3.wt) + ((((stdB3.14 - stdB3.wt)^2 + stdB3.14.se^2)^0.5) * gpa.14.wt/sumB3.wt) + ((((stdB3.10 - stdB3.wt)^2 + stdB3.10.se^2)^0.5) * gpa.10.wt/sumB3.wt)) 

###Models with hs.engl.x4

###Variance inflation factors.

library(car)

ivif.x4.15 <- (vif(gpa.15)[4])^-1
ivif.x4.12 <- (vif(gpa.12)[3])^-1
ivif.x4.13 <- (vif(gpa.13)[3])^-1
ivif.x4.08 <- (vif(gpa.08)[2])^-1
ivif.x4.14 <- (vif(gpa.14)[3])^-1
ivif.x4.10 <- (vif(gpa.10)[2])^-1
ivif.x4.09 <- (vif(gpa.09)[2])^-1

###Partial standard deviations following Bring (1994).

partSD.x4.15 <- ivif.x4.15^0.5 * ((nobs(gpa.15)-1)/(nobs(gpa.15)- (length(coef(gpa.15)) -1)))^0.5 * sd(hs.engl.x4)

partSD.x4.12 <- ivif.x4.12^0.5 * ((nobs(gpa.12)-1)/(nobs(gpa.12)- (length(coef(gpa.12)) -1)))^0.5 * sd(hs.engl.x4)

partSD.x4.13 <- ivif.x4.13^0.5 * ((nobs(gpa.13)-1)/(nobs(gpa.13)- (length(coef(gpa.13)) -1)))^0.5 * sd(hs.engl.x4)

partSD.x4.08 <- ivif.x4.08^0.5 * ((nobs(gpa.08)-1)/(nobs(gpa.08)- (length(coef(gpa.08)) -1)))^0.5 * sd(hs.engl.x4)

partSD.x4.14 <- ivif.x4.14^0.5 * ((nobs(gpa.14)-1)/(nobs(gpa.14)- (length(coef(gpa.14)) -1)))^0.5 * sd(hs.engl.x4)

partSD.x4.10 <- ivif.x4.10^0.5 * ((nobs(gpa.10)-1)/(nobs(gpa.10)- (length(coef(gpa.10)) -1)))^0.5 * sd(hs.engl.x4)

partSD.x4.09 <- ivif.x4.09^0.5 * ((nobs(gpa.09)-1)/(nobs(gpa.09)- (length(coef(gpa.09)) -1)))^0.5 * sd(hs.engl.x4)

###Standardized estimates for hs.engl.x4 using partial SDs.

stdB4.15 <- coef(gpa.15)[5] * partSD.x4.15
stdB4.12 <- coef(gpa.12)[4] * partSD.x4.12
stdB4.13 <- coef(gpa.13)[4] * partSD.x4.13
stdB4.08 <- coef(gpa.08)[3] * partSD.x4.08
stdB4.14 <- coef(gpa.14)[4] * partSD.x4.14
stdB4.10 <- coef(gpa.10)[3] * partSD.x4.10
stdB4.09 <- coef(gpa.09)[3] * partSD.x4.09
stdB4.04 <- coef(gpa.04)[2] * sd(hs.engl.x4)


###Estimates of standard errors for standardized estimates of ###hs.engl.x4.

stdB4.15.se <- ((coef(summary((gpa.15)))[5,2])^2 * partSD.x4.15^2)^0.5
stdB4.12.se <- ((coef(summary((gpa.12)))[4,2])^2 * partSD.x4.12^2)^0.5
stdB4.13.se <- ((coef(summary((gpa.13)))[4,2])^2 * partSD.x4.13^2)^0.5
stdB4.08.se <- ((coef(summary((gpa.08)))[3,2])^2 * partSD.x4.08^2)^0.5
stdB4.14.se <- ((coef(summary((gpa.14)))[4,2])^2 * partSD.x4.14^2)^0.5
stdB4.10.se <- ((coef(summary((gpa.10)))[3,2])^2 * partSD.x4.10^2)^0.5
stdB4.09.se <- ((coef(summary((gpa.09)))[3,2])^2 * partSD.x4.09^2)^0.5
stdB4.04.se <- ((coef(summary((gpa.04)))[2,2])^2 * (sd(hs.engl.x4))^2)^0.5


###Renormalize weights for models with hs.engl.x4.

sumB4.wt<- sum(gpa.15.wt,gpa.12.wt,gpa.13.wt,gpa.08.wt,gpa.14.wt,gpa.10.wt,gpa.09.wt,gpa.04.wt)

###Model-averaged standardized estimate.

stdB4.wt <- ((stdB4.15 * gpa.15.wt) + (stdB4.12 * gpa.12.wt) + (stdB4.13 * gpa.13.wt) + (stdB4.08 * gpa.08.wt) + (stdB4.14 * gpa.14.wt) + (stdB4.10 * gpa.10.wt) + (stdB4.09 * gpa.09.wt) + (stdB4.04 * gpa.04.wt))/sumB4.wt 

###Standard error for model-averaged standardized estimate.

stdB4.wt.se <- (((((stdB4.15 - stdB4.wt)^2 + stdB4.15.se^2)^0.5) * gpa.15.wt/sumB4.wt)  + ((((stdB4.12 - stdB4.wt)^2 + stdB4.12.se^2)^0.5) * gpa.12.wt/sumB4.wt) + ((((stdB4.13 - stdB4.wt)^2 + stdB4.13.se^2)^0.5) * gpa.13.wt/sumB4.wt) + ((((stdB4.08 - stdB4.wt)^2 + stdB4.08.se^2)^0.5) * gpa.08.wt/sumB4.wt) + ((((stdB4.14 - stdB4.wt)^2 + stdB4.14.se^2)^0.5) * gpa.14.wt/sumB4.wt) + ((((stdB4.10 - stdB4.wt)^2 + stdB4.10.se^2)^0.5) * gpa.10.wt/sumB4.wt) + ((((stdB4.09 - stdB4.wt)^2 + stdB4.09.se^2)^0.5) * gpa.09.wt/sumB4.wt) + ((((stdB4.04 - stdB4.wt)^2 + stdB4.04.se^2)^0.5) * gpa.04.wt/sumB4.wt)) 


###To convert model-averaged standardized estimates back to given ###covariance structure. Here correlation among all 4 predictors in ###model 15. 

stdB1.wt*((partSD.x1.15)^-1)
((stdB1.wt.se^2)*((partSD.x1.15)^-2))^0.5

stdB2.wt*((partSD.x2.15)^-1)
((stdB2.wt.se^2)*((partSD.x2.15)^-2))^0.5

stdB3.wt*((partSD.x3.15)^-1)
((stdB3.wt.se^2)*((partSD.x3.15)^-2))^0.5

stdB4.wt*((partSD.x4.15)^-1)
((stdB4.wt.se^2)*((partSD.x4.15)^-2))^0.5

###Ratio of standardized model-averaged estimates

stdB2.wt/stdB1.wt
stdB3.wt/stdB1.wt
stdB4.wt/stdB1.wt

###AICc weighted average ratios of t-statistics for sat.math.x1.

RI.x1.11 <- abs(coef(summary(gpa.11))[2,3])/max(abs(coef(summary(gpa.11))[-1,3]))
RI.x1.05 <- abs(coef(summary(gpa.05))[2,3])/max(abs(coef(summary(gpa.05))[-1,3]))
RI.x1.06 <- abs(coef(summary(gpa.06))[2,3])/max(abs(coef(summary(gpa.06))[-1,3]))
RI.x1.15 <- abs(coef(summary(gpa.15))[2,3])/max(abs(coef(summary(gpa.15))[-1,3]))
RI.x1.12 <- abs(coef(summary(gpa.12))[2,3])/max(abs(coef(summary(gpa.12))[-1,3]))
RI.x1.13 <- abs(coef(summary(gpa.13))[2,3])/max(abs(coef(summary(gpa.13))[-1,3]))
RI.x1.01 <- abs(coef(summary(gpa.01))[2,3])/max(abs(coef(summary(gpa.01))[-1,3]))
RI.x1.08 <- abs(coef(summary(gpa.08))[2,3])/max(abs(coef(summary(gpa.08))[-1,3]))

avgRI.x1 <- (RI.x1.11*gpa.11.wt + RI.x1.05*gpa.05.wt + RI.x1.06*gpa.06.wt + RI.x1.15*gpa.15.wt + RI.x1.12*gpa.12.wt + RI.x1.13*gpa.13.wt + RI.x1.01*gpa.01.wt + RI.x1.08*gpa.08.wt)/sumB1.wt

###AICc weighted average ratios of t-statistics for sat.verb.x2

RI.x2.11 <- abs(coef(summary(gpa.11))[3,3])/max(abs(coef(summary(gpa.11))[-1,3]))
RI.x2.05 <- abs(coef(summary(gpa.05))[3,3])/max(abs(coef(summary(gpa.05))[-1,3]))
RI.x2.14 <- abs(coef(summary(gpa.14))[2,3])/max(abs(coef(summary(gpa.14))[-1,3]))
RI.x2.15 <- abs(coef(summary(gpa.15))[3,3])/max(abs(coef(summary(gpa.15))[-1,3]))
RI.x2.12 <- abs(coef(summary(gpa.12))[3,3])/max(abs(coef(summary(gpa.12))[-1,3]))
RI.x2.07 <- abs(coef(summary(gpa.07))[2,3])/max(abs(coef(summary(gpa.07))[-1,3]))
RI.x2.09 <- abs(coef(summary(gpa.09))[2,3])/max(abs(coef(summary(gpa.09))[-1,3]))
RI.x2.02 <- abs(coef(summary(gpa.02))[2,3])/max(abs(coef(summary(gpa.02))[-1,3]))

avgRI.x2 <- (RI.x2.11*gpa.11.wt + RI.x2.05*gpa.05.wt + RI.x2.14*gpa.14.wt + RI.x2.15*gpa.15.wt + RI.x2.12*gpa.12.wt + RI.x2.07*gpa.07.wt + RI.x2.09*gpa.09.wt + RI.x2.02*gpa.02.wt)/sumB2.wt


###AICc weighted average ratios of t-statistics for hs.math.x3

RI.x3.11 <- abs(coef(summary(gpa.11))[4,3])/max(abs(coef(summary(gpa.11))[-1,3]))
RI.x3.07 <- abs(coef(summary(gpa.07))[3,3])/max(abs(coef(summary(gpa.07))[-1,3]))
RI.x3.06 <- abs(coef(summary(gpa.06))[3,3])/max(abs(coef(summary(gpa.06))[-1,3]))
RI.x3.15 <- abs(coef(summary(gpa.15))[4,3])/max(abs(coef(summary(gpa.15))[-1,3]))
RI.x3.03 <- abs(coef(summary(gpa.03))[2,3])/max(abs(coef(summary(gpa.03))[-1,3]))
RI.x3.13 <- abs(coef(summary(gpa.13))[3,3])/max(abs(coef(summary(gpa.13))[-1,3]))
RI.x3.14 <- abs(coef(summary(gpa.14))[3,3])/max(abs(coef(summary(gpa.14))[-1,3]))
RI.x3.10 <- abs(coef(summary(gpa.10))[2,3])/max(abs(coef(summary(gpa.10))[-1,3]))

avgRI.x3 <- (RI.x3.11*gpa.11.wt + RI.x3.07*gpa.07.wt + RI.x3.06*gpa.06.wt + RI.x3.15*gpa.15.wt + RI.x3.03*gpa.03.wt + RI.x3.13*gpa.13.wt + RI.x3.14*gpa.14.wt + RI.x3.10*gpa.10.wt)/sumB3.wt

###AICc weighted average ratios of t-statistics for hs.engl.x4

RI.x4.15 <- abs(coef(summary(gpa.15))[5,3])/max(abs(coef(summary(gpa.15))[-1,3]))
RI.x4.12 <- abs(coef(summary(gpa.12))[4,3])/max(abs(coef(summary(gpa.12))[-1,3]))
RI.x4.13 <- abs(coef(summary(gpa.13))[4,3])/max(abs(coef(summary(gpa.13))[-1,3]))
RI.x4.08 <- abs(coef(summary(gpa.08))[3,3])/max(abs(coef(summary(gpa.08))[-1,3]))
RI.x4.14 <- abs(coef(summary(gpa.14))[4,3])/max(abs(coef(summary(gpa.14))[-1,3]))
RI.x4.10 <- abs(coef(summary(gpa.10))[3,3])/max(abs(coef(summary(gpa.10))[-1,3]))
RI.x4.09 <- abs(coef(summary(gpa.09))[3,3])/max(abs(coef(summary(gpa.09))[-1,3]))
RI.x4.04 <- abs(coef(summary(gpa.04))[2,3])/max(abs(coef(summary(gpa.04))[-1,3]))

avgRI.x4 <- (RI.x4.15*gpa.15.wt + RI.x4.12*gpa.12.wt + RI.x4.13*gpa.13.wt + RI.x4.08*gpa.08.wt + RI.x4.14*gpa.14.wt + RI.x4.10*gpa.10.wt + RI.x4.09*gpa.09.wt + RI.x4.04*gpa.04.wt)/sumB4.wt


###To demonstrate that units of sat.math.x1 correlated with 3 predictors ###are not the same units as sat.math.x1 correlated with 2 predictors ###for figure.

resid.x1.x2x3 <- lm(sat.math.x1~sat.verb.x2 + hs.math.x3)$resid
resid.x1.x2x3x4 <- lm(sat.math.x1~sat.verb.x2 + hs.math.x3 + hs.engl.x4)$resid

plot(resid.x1.x2x3x4,gpa.14$resid,xlim=c(-200,300),ylim=c(-0.9,0.9),pch=1,cex=1.5)
abline(a=0, b=0.002010)

plot(resid.x1.x2x3,gpa.07$resid,xlim=c(-200,300),ylim=c(-0.9,0.9),pch=1,cex=1.5)
abline(a=0,b=0.002185)





	 
