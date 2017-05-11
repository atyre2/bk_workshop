# compositional covariates
library(tidyverse)
library(GGally)
set.seed(87)
df <- data_frame(
  pos.tot = runif(200,min=0.8,max=1.0),
  urban.tot = pmin(runif(200,min=0.0,max=0.02),1.0 - pos.tot),
  neg.tot = (1.0 - pmin(pos.tot + urban.tot,1)),
  x1= pmax(pos.tot - runif(200,min=0.05,max=0.30),0),
  x3= pmax(neg.tot - runif(200,min=0.0,max=0.10),0),
  x2= pmax(pos.tot - x1 - x3/2,0),
  x4= pmax(1 - x1 - x2 - x3 - urban.tot,0))

ggpairs(df)
mean.y <- exp(-5.8 + 6.3*df$x1 + 15.2*df$x2)
df$y <- rpois(200,mean.y)

library(MuMIn)
globalMod <- glm(y ~ x1 + x2 + x3 + x4, 
                 data = df, na.action = na.fail, 
                 family=poisson)
fits <- dredge(globalMod)
model.sel(fits)
summary(model.avg(fits))

df2 <- df[, 4:7]
grouse_pca <- princomp(df2)
plot(grouse_pca)
biplot(grouse_pca)

df_pca <- as_data_frame(grouse_pca$scores)
df_pca$y <- df$y
globalMod_pca <- glm(y ~ Comp.1 + Comp.2 + Comp.3 + Comp.4, 
                     data = df_pca, na.action = na.fail, 
                     family=poisson())
fits <- dredge(globalMod_pca)
model.sel(fits)

str(grouse_pca)
grouse_pca$loadings
nd <- data_frame(x1 = seq(0.52, 0.85, 0.01),
                 x2 = 0.15,
                 pos = x1 + x2,
                 x3 = pmax(1-pos,0),
                 x4 = pmax(1-pos-x3,0))
nd_pca <- as_data_frame(predict(grouse_pca, nd[,-3]))
predmodel <- get.models(fits, subset = 2)
preds <- cbind(nd,
               meancount = predict(predmodel[[1]], 
                                   newdata=nd_pca, 
                                   type = "response"))
ggplot(preds, aes(x = x1, y=meancount)) + geom_line() + 
  geom_point(aes(x = x1, y = y, color = x2), data=df)

library(robCompositions)
library(dtplyr)

pivotCoord(as.data.frame(df2))
df2_replaced <- imputeBDLs(as.data.frame(df2), 
                           method="lm", 
                           dl = rep(0.005,4))
rowSums(df2_replaced$x)

df2_ilr <- pivotCoord(df2_replaced$x)
df2_ilr$y <- df$y
globalMod_ilr <- glm(y ~ `x1_x3-x2-x4` +  `x3_x2-x4`  + `x2_x4`, 
                     data = df2_ilr, na.action = na.fail, 
                     family=poisson())
fits <- dredge(globalMod_ilr)
model.sel(fits)
ggpairs(df2_ilr)

# using strategy from age-cancer example in Hron etal 2012

df2_ilr_x1 <- pivotCoord(df2_replaced$x, pivot = 1)
df2_ilr_x2 <- pivotCoord(df2_replaced$x, pivot = 3)
df2_ilr_x3 <- pivotCoord(df2_replaced$x, pivot = 2)
df2_ilr_x4 <- pivotCoord(df2_replaced$x, pivot = 4)

df2_ilr2 <- data.frame(x1i = df2_ilr_x1[,1],
                       x2i = df2_ilr_x2[,1],
                       x3i = df2_ilr_x3[,1],
                       x4i = df2_ilr_x4[,1])
ggpairs(df2_ilr2)

df2_ilr2$y <- df$y
globalMod_ilr2 <- glm(y ~ x1i + x2i + x3i + x4i, 
                     data = df2_ilr2, na.action = na.fail, 
                     family=poisson())
fits <- dredge(globalMod_ilr2)
model.sel(fits)
get.models(fits, subset = 4)

names(df2_ilr_x1) <- c("x1","x2","x3")
df2_ilr_x1$y <- df$y
gM1 <- glm(y ~ x1 + x2 + x3, 
          data = df2_ilr_x1, na.action = na.fail, 
          family=poisson()) 

df2_ilr_x2$y <- df$y
gM2 <- glm(y ~ `x2_x1-x3-x4` + `x1_x3-x4` + `x3_x4`, 
           data = df2_ilr_x2, na.action = na.fail, 
           family=poisson()) 

df2_ilr_x3$y <- df$y
gM3 <- glm(y ~ `x3_x1-x2-x4` + `x1_x2-x4` + `x2_x4`, 
           data = df2_ilr_x3, na.action = na.fail, 
           family=poisson()) 

df2_ilr_x4$y <- df$y
gM4 <- glm(y ~ `x4_x1-x3-x2` + `x1_x3-x2` + `x3_x2`, 
           data = df2_ilr_x4, na.action = na.fail, 
           family=poisson()) 

sgM4 <- summary(gM4)
sgM3 <- summary(gM3)
sgM2 <- summary(gM2)
sgM1 <- summary(gM1)
str(sgM1)

# table similar to Hron etal 2012
rbind(sgM1$coefficients[1:2,],
      sgM2$coefficients[2,],
      sgM3$coefficients[2,],
      sgM4$coefficients[2,])

imputeBDLs(nd[,c(1,4,2,5)],
           method="lm",
           dl = rep(0.005,4))

# doesn't work pivotCoord(nd[,-3])
nd_replaced <- nd %>%
  select(-pos) %>%
  filter(x3 != 0) %>%
  mutate(x1 = x1 - x1*0.005,
         x3 = x3 - x3*0.005,
         x2 = x2 - x2*0.005,
         x4 = 0.005) %>%
  select(1,3,2,4)

nd_ilr <- pivotCoord(as.data.frame(nd_replaced))

preds <- predict(globalMod_ilr, newdata = nd_ilr,
                 type = "response")
preds <- cbind(nd_replaced, meangrouse = preds)
ggplot(preds, aes(x = x1, y=meangrouse)) + 
  geom_line() + 
  geom_point(aes(x = x1, y=y), data=df)

preds <- predict(globalMod_ilr, type = "response")
preds <- data.frame(y = df$y,
           y_fit = preds)
ggplot(preds, aes(x = y, y=y_fit)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0, linetype = 2) + 
  geom_smooth(method = "lm")

## try to replicate their example
data(cancer)
lmCoDaX(cancer$discharges,
        as.matrix(cancer[,3:5]))
#fails

y <- as.numeric(apply(expendituresEU,1,sum))
lmCoDaX(y, expendituresEU, method="classical")
lmCoDaX(y, expendituresEU, method="robust")

## example of standardizing etc etc
ggpairs(iris)
iris.lm <- lm(Sepal.Width ~ Sepal.Length * Species, 
              data = iris)
iris.scaled <- MuMIn:::stdize(iris)
iris.scaled$Sepal.Width <- iris$Sepal.Width
iris.lm.sc <- lm(Sepal.Width ~ z.Sepal.Length * Species, 
              data = iris.scaled)
contrasts(iris.scaled$Species) <- contr.sum(3)
iris.lm.scz <- lm(Sepal.Width ~ z.Sepal.Length * Species, 
                 data = iris.scaled)
summary(iris.lm.scz)
