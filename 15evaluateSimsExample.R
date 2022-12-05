# Model whether CO contains true prop depending on prop

setwd("~/git_repos/vagrantDnaSim/exampleOut/")

lines <- readLines("grepOut")
lines

dat <- data.frame(propFact = as.numeric(substr(lines, 6,6)),
           contained= ifelse(substr(lines, 13,13) == "Y", T, F),
           GS=as.numeric(substr(lines, 1, 3)),
           repl=as.numeric(substr(lines, 7,7))
)
head(dat)
summary(dat)
dat$GS[dat$GS==1] <- 1000


dat$prop <- (4e-5 - 5e-6*dat$propFact) * 15000000 * 1000 / dat$GS / 1000000
# 5 columns:
#  propFact - indicates insertion rate
#  contained - does confidence interval contain the true value?
#  GS - genome size in Mbp
#  repl - number of replicate
#  prop - simulated "average" genome proportion of vagrant DNA
head(dat)

# some runs fail to produce an outfile
# For all combinations, there were 10 replicates run.
# List these
expectedRuns <- data.frame(propFact=rep(0:6, each=10),
                      repl=0:9,
                      GS=rep(c(1000, 250, 500), each=7*10)
)
expectedRuns$prop <- (4e-5 - 5e-6*expectedRuns$propFact) * 15000000 * 1000 / expectedRuns$GS / 1000000
names(dat)
names(expectedRuns)

# merge observed and expected run data in order to count "hit and misses"
#  for a binomial GLM (no actually needed as we use "weights=10" later)
datWithMissing <- merge(dat, expectedRuns, by=c("propFact", "repl", "GS"), all=T)
head(datWithMissing)
datWithMissing$contained[is.na(datWithMissing$contained)] <- F
names(datWithMissing)[6] <- "prop"
head(datWithMissing)
datWithMissing$numtDep <- datWithMissing$prop * datWithMissing$GS *
  1000000/ 16000 * 0.1


datAgg <- aggregate(contained ~ prop+ GS + numtDep, data=datWithMissing, function(x) sum(x)/length(x))

library(ggplot2)
depPlot <- ggplot(datAgg, aes(numtDep, contained, fill=factor(GS))) +
  geom_point(position=position_jitter(w=0.1, h=0.01),
             col="white", shape=21, size=2)
depPlot
propPlot <- ggplot(datAgg, aes(prop, contained, fill=factor(GS))) +
  geom_point(position=position_jitter(w=0.00003, h=0.01),
             col="white", shape=21, size=2)
propPlot
#modelling

par(mfrow=c(2,2))

# predict accurate fit with genome size and numt depth
glm03 <- glm(contained ~ GS +  numtDep, weights=rep(10, nrow(datAgg)),
             data=datAgg,
             family=binomial)
glm03l <- glm(contained ~ log(GS) +  numtDep, weights=rep(10, nrow(datAgg)),
              data=datAgg,
              family=binomial)
glm03ll <- glm(contained ~ log(GS) +  log(numtDep), weights=rep(10, nrow(datAgg)),
               data=datAgg,
               family=binomial)


# fits look alright
par(mfrow=c(2,2))
plot(glm03)
plot(glm03l)
plot(glm03ll)
AIC(glm03, glm03l, glm03ll)
summary(glm03)
summary(glm03l)
# both predictors log transformed
summary(glm03ll)

# make dataframe for prediction from model fit
predDatDep <- data.frame(GS=rep(c(1000,500,250), each=length(seq(1, 4,by=0.05))),
                      numtDep=c(seq(1, 4,by=0.05))
)
p3 <- plogis(predict(glm03, predDatDep))
p3l <- plogis(predict(glm03l, predDatDep))
p3ll <- plogis(predict(glm03ll, predDatDep))

predDatDep$p3 <- p3
predDatDep$p3l <- p3l
predDatDep$p3ll <- p3ll

# best model
#pdf("../Accuracy.pdf", width = 7,height=4)
#png("../Accuracy.png", res = 300, width = 7,height=4, units="in")
ggplot(datAgg, aes(numtDep, contained, fill=factor(GS))) +
  geom_line(data=predDatDep, aes(x=numtDep, y=p3ll, col=factor(GS))) +
  geom_point(position=position_jitter(w=0.05, h=0.01),
             col="white", shape=21, size=2) +
  labs(title="Accuracy depends on mapping depth",
          subtitle = "Based on simulations",
       x="Expected mapping depth of nuclear reads to the extranuclear reference",
       y="Accuracy") +
  scale_fill_manual(name="Genome size (Mbp)",aesthetics=c("colour", "fill"),labels=c(250, 500, 1000), values=c(2,3,4))
#dev.off()




coef(glm03ll)
# acc = 3.6 - 1.5*log(gsMbp) + 6.8*log(dep)

#0.95 = c0 + c1*log(gsMbp) + c2*log(dep)

# -c1*log(gsMbp) = c0-0.95 + c2*log(dep)
#
# -c1*log(gsMbp) = c0-0.95 + c2*log(gsMbp * prop *ndep / 0.016)
# -c1*log(gsMbp) = c0-0.95 + c2*(log(gsMbp)+log(ndep)+log(prop)-log(0.016))
# -c1*log(gsMbp) = c0-0.95 + c2*log(gsMbp)+c2*log(prop)+c2*log(ndep)-c2*log(0.016)
# -c1*log(gsMbp)- c2*log(gsMbp) = c0-0.95 +c2*log(prop)+c2*log(ndep)-c2*log(0.016)
# -(c1+c2) * log(gsMbp) =  c0-0.95 +c2*log(prop)+c2*log(ndep)-c2*log(0.016)
# log(gsMbp) =  (c0-0.95 +c2*log(prop)+c2*log(ndep)-c2*log(0.016))/-(c1+c2)

# dep = gs * prop * ndep /extrGS

gs <- function(prop, c0, c1, c2, ndep){
  exp((c0-0.95 +c2*log(prop)+c2*log(ndep)-c2*log(0.016))/-(c1+c2)) * 1000000
}
curve(gs(x, coef(glm03ll)[1], coef(glm03ll)[2], coef(glm03ll)[3], ndep=0.1), 0, 0.01)
par(mfrow=c(1,1))

#pdf("../DepReq.pdf", height=5.5, width=7)
#png("../DepReq.png", height=5.5, width=7, units = "in", res=150)
curve(gs(x, coef(glm03ll)[1], coef(glm03ll)[2], coef(glm03ll)[3], ndep=0.1),
      1e-5,
      0.001,
      log="xy",
      xlab="Nuclear proportion of vagrant DNA",
      ylab="Genome Size",
      main="Sequencing depth required",
      lwd=2,
      ylim=c(100000000, 20000000000),
      xaxt = 'n',
      yaxt = 'n')
xat <- c(0.00001, 0.00002, 0.00005, 0.0001,0.0002, 0.0005, 0.001)
axis(1, at = xat,
     labels = expression("10"^-5, "2×10"^-5, "5×10"^-5, "10"^-4, "2×10"^-4, "5×10"^-4, "10"^-3)
)
yat <- c(1e8, 2e8, 5e8, 1e9, 2e9, 5e9, 1e10, 2e10)
axis(2,
     at = yat,
     labels = c("100 M", "200 M", "500 M", "1 G", "2 G", "5 G", "10 G", "20 G"),
     las=2
     )
curve(gs(x, coef(glm03ll)[1], coef(glm03ll)[2], coef(glm03ll)[3], ndep=0.2),
      1e-5,
      0.001,
      lty=2,
      lwd=2,
      add=T)

curve(gs(x, coef(glm03ll)[1], coef(glm03ll)[2], coef(glm03ll)[3], ndep=0.5),
      1e-5,
      0.01,
      lty=3,
      lwd=2,
      add=T)

curve(gs(x, coef(glm03ll)[1], coef(glm03ll)[2], coef(glm03ll)[3], ndep=1),
      1e-5,
      0.01,
      lty=4,
      lwd=2,
      add=T)

abline(v=xat,
       col="grey",
       lty=3)
abline(h=yat,
       col="grey",
       lty=3)

legend("bottomleft",
       lty=1:4,
       lwd=2,
       legend=c("0.1x", "0.2x", "0.5x", "1.0x"),
       title = "Sequencing depth"
       )
points(y=c(1, 3.5, 17)*1e9,
       x=c(0.000081, 0.00014, 0.00056),
       pch=c("P", "H", "G"))
#dev.off()
