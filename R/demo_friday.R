# playing with Anne's data

library(tidyverse)
library(GGally)
nuthatch <- read.csv("data/R_glm.csv")

# compositional covariates
ggpairs(nuthatch[,8:11])

glm1 <- glm(detec~agri + BL + etc + pine,
            data=nuthatch, family=binomial)
summary(glm1)

glm2 <- glm(detec~agri + etc + pine,
            data=nuthatch, family=binomial)
summary(glm2)

nh_pca <- princomp(nuthatch[,8:11])
biplot(nh_pca, choices = c(1,2))
summary(nh_pca)

# try just a 2 dimensional PCA
test <- princomp(nuthatch[,c(9,11)])
biplot(test)
nh_2 <- as.data.frame(nh_pca$scores)
nh_2$detec <- nuthatch[,"detec"]

glm3 <- glm(detec ~ Comp.1 + I(Comp.1^2),
            data=nh_2,
            family=binomial)
summary(glm3)

nd <- data_frame(agri = 0, 
                 BL = seq(0,1,0.01),
                 etc = 0, 
                 pine = 1-BL)
nd_scores <- predict(nh_pca, nd)
preds <- cbind(nd,
               predict(glm3, as.data.frame(nd_scores), se.fit = TRUE))

# save inverse link function to use in backtransforming predictions
inv.logit <- binomial()$linkinv

# summarize nuthatch data to reduce overplotting
nh_summary <- nuthatch %>%
  group_by(detec, BL) %>%
  summarize(count = n())

ggplot() + 
  geom_point(data = nh_summary,
              aes(x = BL,y = detec, size = count)) +
  geom_line(data = preds, 
            mapping = aes(x = BL, 
                          y = inv.logit(fit))) +
  geom_ribbon(data = preds, 
              mapping = aes(x = BL, 
                            ymin = inv.logit(fit - 2*se.fit),
                            ymax = inv.logit(fit + 2*se.fit)),
              alpha = 0.5)

# play data
play <- read_csv("data/adult howler play_CSV.csv")
names(play)[12:15] <- c("FE","TR","RE","PLAY")
ggpairs(play[,12:15])
play <- filter(play, PLAY > 0)

# strategy 1 - ignore the problem
library(betareg)
br1 <- betareg(PLAY~FE + TR + RE,
               data=play)
mutate(play, total = FE + TR + RE + PLAY) %>%
  select(total)
# add up to 1, so can't trust br1 numbers

# strategy 2 - LOO
br2 <- betareg(PLAY~FE + TR, data=play)
summary(br2)

# strategy 3 - PCA
play_pca <- princomp(play[,12:14])
summary(play_pca)
biplot(play_pca)
play_scores <- as.data.frame(play_pca$scores)
play_scores$PLAY <- play$PLAY
br3 <- betareg(PLAY ~ Comp.1 + Comp.2,
               data=play_scores)
summary(br3)
play_pca$loadings
plot(br3)

ggplot(play) + 
  geom_point(aes(x=TR, y=PLAY)) + 
  facet_wrap(~SITE) + 
  geom_smooth(aes(x=TR, y=PLAY), method="lm")

# spatial autocorrelation
coord <- read.table("data/coord.txt")
names(coord) <- c("x","y")
envir <- read.table("data/envir_response.txt", 
           header = TRUE)
names(envir) <- paste0("X",1:5)
spp <- read.table("data/species_response.txt")

df <- cbind(coord, envir, spp)

ggplot(df, aes(x=X1, y=V1)) + 
  geom_point() + 
  geom_smooth()

ggplot(df, aes(x = x, y=y, color=X1)) + 
  geom_point()

library(spatial)
fake.kr <- surf.ls(0, x=df$x, y=df$y, z=df$X1)
correlogram(fake.kr, nint = 10)

fake.kr <- surf.ls(2, x=df$x, y=df$y, z=df$X1)
correlogram(fake.kr, nint = 10)
fake.grid <- trmat(fake.kr, 0, 100,0,100, 100)
contour(fake.grid)

fake.kr <- surf.ls(3, x=df$x, y=df$y, z=df$X1)
correlogram(fake.kr, nint = 10)
fake.grid <- trmat(fake.kr, 0, 100,0,100, 100)
contour(fake.grid)

spp.kr <- surf.ls(0, x=df$x, y=df$y, z=df$V6)
correlogram(spp.kr, nint = 10)
spp.lm <- lm(V6~X1, data=df)
library(broom)
spp.aug <- augment(spp.lm)
spp.aug <- cbind(df[,c("x","y")],
                 spp.aug)
spp.kr <- surf.ls(0, x=spp.aug$x, y=spp.aug$y,
                  z = spp.aug$.resid)
spp.crg <- correlogram(spp.kr, 10)
spp.crg$t <- with(spp.crg, qt(0.975, df=cnt))
spp.crg$r <- with(spp.crg, t / sqrt(cnt-2+t^2))
lines(spp.crg$x, spp.crg$r, lty=2, col="blue")
lines(spp.crg$x, -spp.crg$r, lty=2, col="blue")


correlog(x = spp.aug$x,
         y = spp.aug$y,
         z = spp.aug$.resid,
         increment = 10)
xydist <- dist(as.matrix(coord))
gearymoran(xydist, spp.aug$.resid)

variogram(spp.kr, 10)





