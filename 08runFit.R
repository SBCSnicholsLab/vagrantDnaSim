library(vagrantDNA)

tDir <- commandArgs(trailingOnly = T)[1]
outPrf <- commandArgs(trailingOnly = T)[2]

trueVal <- as.numeric(readLines(paste0(tDir, outPrf, ".log")))
dat <- read.table(paste0(tDir, "transformedData.csv"), header = T)

#head(dat)
a <- rainbowPlot(dat)


cat(
  c(a$intercepts, a$depth.est),
  file=paste0(tDir, outPrf, ".log"),
  append=T
)

cat(
  ifelse(trueVal < a$intercepts[3] & trueVal > a$intercepts[2],
         "\nYES, contained\n",
         "\nNOT contained\n"),
  file=paste0(tDir, outPrf, ".log"),
  append=T
)

