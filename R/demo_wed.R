set.seed(9377938)
p = 3 # this many variables
# this is the mean vector - set to 0 for now
Means <- rep(0, 3)
# this is the covariance matrix
Sigma <- matrix(c(1, 0, .8, 
                  0, 1, 0,
                  .8, 0, 1), nrow= p, ncol = p)
library(MASS)
library(tidyverse)
## Loading tidyverse: ggplot2
## Loading tidyverse: tibble
## Loading tidyverse: tidyr
## Loading tidyverse: readr
## Loading tidyverse: purrr
## Loading tidyverse: dplyr
## Conflicts with tidy packages ----------------------------------------------
## filter(): dplyr, stats
## lag():    dplyr, stats
## select(): dplyr, MASS
library(GGally)
## 
## Attaching package: 'GGally'
## The following object is masked from 'package:dplyr':
## 
##     nasa
df <- as_data_frame( mvrnorm(n = 100, mu= Means, Sigma=Sigma))
ggpairs(df)

df <- mutate(df, y = V1 + rnorm(n()))
ggplot(df, aes(x = V1, y=y)) + geom_point() + 
  geom_smooth(method="lm")
models <- list(y~1,
               y~V1,
               y~V2,
               y~V3,
               y~V1 + V2,
               y~V1 + V3,
               y~V2 + V3,
               y~V1 + V2 + V3)
fits <- map(models, lm, data=df)
library(broom)
coefs <- map_df(fits, tidy, .id = "model")
ggplot(coefs, aes(x = model, y=estimate)) + 
  geom_point() + 
  facet_wrap(~term) + 
  geom_errorbar(aes(ymin = estimate - 1.96*std.error,
                    ymax = estimate + 1.96*std.error))
library(car)
## 
## Attaching package: 'car'
## The following object is masked from 'package:dplyr':
## 
##     recode
## The following object is masked from 'package:purrr':
## 
##     some
vif(fits[[8]])
##        V1        V2        V3 
## 41.113603  1.006754 41.076494
# What do those numbers actually mean? We can calculate the ratio of the variances of β3
# in models 7 and 8.

filter(coefs, term == "V3",
       model > 6) %>%
  mutate(variance = std.error^2,
         vif_ = lead(variance)/variance) %>%
  select(vif_)

library(MuMIn)
model.sel(fits)
summary(model.avg(fits))

makeUpData <- function(r = 0, n = 100){
  p = 3 # this many variables
  # this is the mean vector - set to 0 for now
  Means <- rep(0, 3)
  # this is the covariance matrix
  Sigma <- matrix(c(1, 0, r, 
                    0, 1, 0,
                    r, 0, 1), nrow= p, ncol = p)
  df <- as_data_frame( mvrnorm(n = n, mu= Means, Sigma=Sigma))

  df <- mutate(df, y = V1 + rnorm(n()))
  models <- list(y~1,
                 y~V1,
                 y~V2,
                 y~V3,
                 y~V1 + V2,
                 y~V1 + V3,
                 y~V2 + V3,
                 y~V1 + V2 + V3)
  fits <- map(models, lm, data=df)
  modselTable <- as.data.frame(model.sel(fits))
  modavg <- model.avg(fits)
  fullavg <- coefTable(modavg, full=TRUE)
  result <- data_frame(topmodel = as.numeric(row.names(modselTable))[1],
                       estimate = fullavg[2,1],
                       se = fullavg[2, 2])
  return(result)
}

results <- data.frame(topmodel = rep(NA, 1000),
                      estimate  = rep(NA, 1000),
                      se = rep(NA, 1000))
for (i in 1:1000){
  results[i,] <- makeUpData(r = 0.6, n=50)
}
head(results)
table(results$topmodel)

ggplot(results, aes(x=estimate)) + 
  geom_histogram() + 
  facet_wrap(~topmodel)

nd <- data_frame(V1 = seq(-2, 2, 0.1),
                 V2 = 0,
                 V3 = 0)
ma_fits <- model.avg(fits)
pred <- cbind(nd, 
              y = 1,
              y_0 = predict(ma_fits, newdata = nd),
              y_1 = predict(ma_fits, newdata = nd, se.fit = TRUE))
ggplot(df, aes(x = V1, y = y)) + geom_point() + 
  geom_line(data = pred, aes(y = y_0)) + 
  geom_ribbon(data = pred, aes(
    ymin = y_0 - 1.96*y_1.se.fit,
    ymax = y_0 + 1.96*y_1.se.fit),
    alpha = 0.4)
 