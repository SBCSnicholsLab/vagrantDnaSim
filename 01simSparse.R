library(Matrix)
# Generate data
nInd <- as.numeric(commandArgs(trailingOnly = T)[1])
tDir <- commandArgs(trailingOnly = T)[2]
#nsize <- 500000000 # 500MBp
nsize <- as.numeric(commandArgs(trailingOnly = T)[3])
insRate <- as.numeric(commandArgs(trailingOnly = T)[4])
outPrf <- commandArgs(trailingOnly = T)[5]
meanSeqDep <- as.numeric(commandArgs(trailingOnly = T)[6])
cat(paste0("R: temp dir is ", tDir, "\n"))
# print(paste0("Working in directory ", tDir))

#set.seed(12344567)



#nInd <- 20
#insRate <- 0.00004
#tDir <- "~"
cat(paste0("SIM.R###################################\n"))
cat(paste0("Generating ref genome and data for ", nInd, " Samples...\n"))
cat(paste0("Insertion rate is ", insRate, ".\n"))

# A T G C - 1 2 3 4

# nuclear genome size

extraSize<- 16000
cat("Generating insert sequences...\n")
# original extranuclear reference
sampleBP <- function(l=extraSize, gc=0.3){
  sample(x=c("A", "T", "G", "C"),
         size=l,
         replace = T,
         prob=c((1-gc)/2, (1-gc)/2, gc/2, gc/2)
         )
}

sampleBPnum <- function(l=extraSize, gc=0.3){
  sample(x=1:4,
         size=l,
         replace = T,
         prob=c((1-gc)/2, (1-gc)/2, gc/2, gc/2)
  )
}

extraNuc <- sampleBPnum() # integer vector, not nucleotides!
write(paste0(">refOrig\n", paste(extraNuc, collapse = "")),
      paste0(tDir,"extraNucOrig.fa")
      )

# a function to generate a numt sequence from start, length, and mito sequence
getInsertSeq <- function(pos, len, ref=extraNuc){
  if((pos+len) > length(ref)){
    lref <- rep(ref, ((pos+len) %% length(ref)) + 1)
    #print(length(lref))
    return(lref[pos:(pos+len-1)])
  }
  ref[pos:(pos+len-1)]
}
#getInsertSeq(15990, 11)
times <- numeric()
pos <- numeric()
seqs <- list()
starts <- numeric()
# insertion loop
k <- 0 # counter of insertions
for(i in 1:15000000){ # 15M years

  if(runif(1) < extraSize*0.0115*4/3/1000000){ # test if extranuclear DNA has a substitution this year
    #print("Subst!") # for debug
    p <- sample.int(extraSize, 1)
    extraNuc[p] <- sample(4, 1)
  }
  # test if insertion happens this year
  if(runif(1) < insRate){
    k <- k + 1
    times[k] <- i
    l <- rnbinom(1, size=2, mu=1000) # insert length
    st <- sample.int(extraSize, 1) # insert start site
    starts[k] <- st
    seqs[[k]] <- getInsertSeq(st, l, ref = extraNuc)
    pos[k] <- sample.int(nsize, 1)
  }
}

# plot depth
dep=rep(0,extraSize)

cat("Insert proportion is",
       sum(sapply(seqs, length, USE.NAMES = F))/nsize * 100,
       "%.\n")
#plot(function(x) dnbinom(x, size=2, mu=1000), 0, 20000)
write(sum(sapply(seqs, length, USE.NAMES = F))/nsize,
      paste0(tDir, outPrf, ".log")
      )

cat("Generating genome sequence...\n")
genome <- sparseVector(0, 1, nsize)
#genome <- sampleBP(nsize, gc=0.4)
# write(paste0(">nucOrig\n", paste(genome, collapse = "")),
#       paste0(tDir, "genomeOrig.fa")
# )

cat("Adding inserts...\n")
#print(length(times))
for(i in 1:length(times)){
  
  e <- pos[i] + length(seqs[[i]])
  
  
  dep[c(starts[i]:(starts[i] + length(seqs[[i]])))%%extraSize] <- dep[c(starts[i]:(starts[i] + length(seqs[[i]])))%%extraSize] +1
  
  if(e > nsize){
    #print("DEBUG: Insert too long!")
    seqs[i] <- seqs[i][1:(length(seqs[i]) - e + nsize - 1)]
  }
  
  genome[pos[i]:(pos[i]+length(seqs[[i]])-1)] <- seqs[[i]]
  

}
# Plot deps along the extranuclear genome
pdf(paste0(tDir,"Deps.pdf"))
plot(dep)
grid()
dev.off()
write.table(data.frame(pos=1:extraSize,
                       dep=dep),
            paste0(tDir, "Deps.csv"),
            col.names = T,
            row.names = F,
            quote = F,
            sep=",")


cat("Writing extranuclear sequence...\n")
# which has been mutated
extraNuc[extraNuc==1] <- "A"
extraNuc[extraNuc=="2"] <- "T"
extraNuc[extraNuc=="3"] <- "G"
extraNuc[extraNuc=="4"] <- "C"
write(paste0(">ref\n", paste(extraNuc, collapse = "")),
      paste0(tDir, "extraNuc.fa")
)
cat("Writing nuclear sequence...\n")
heapSize <- 1000000
nHeaps <- ceiling(nsize / heapSize)
write(">nuc",
      paste0(tDir, "genome.fa")
      )
for(i in 1:nHeaps){
  heapEnd <- ifelse(i*heapSize <= nsize, i*heapSize, nsize)
  curHeap <- as.character(as.vector(genome[(1+((i-1)*heapSize)):heapEnd]))
  curHeap[curHeap=="1"] <- "A"
  curHeap[curHeap=="2"] <- "T"
  curHeap[curHeap=="3"] <- "G"
  curHeap[curHeap=="4"] <- "C"
  nToFill <- sum(curHeap=="0")
  curHeap[curHeap=="0"] <- sampleBP(nToFill, gc=0.4)
  write(paste(curHeap, collapse = ""),
        paste0(tDir, "genome.fa"),
        append = T
        )
}

cat("Generating random sequencing depths and numbers of extranuclear DNA copies...\n")
# random sequencing depths around 5% # or whateve is supplied
deps <- rnorm(nInd, 0.1, meanSeqDep)
deps[deps < 0]  <- -deps[deps < 0]
write.table(as.integer(round(deps * nsize/152)),
            paste0(tDir, "pairCounts"),
            col.names = F,
            row.names = F,
            quote = F)

#how much extraNuc DNA? around 3%
#enDeps <- rnorm(20, 0.03, 0.01)

#enDeps <- runif(nInd, 0, 0.02)
enDeps <- rexp(nInd, 1/0.01)
#enDeps <- 2^(log2(0.02) - (0:(nInd-1))/4)
enDeps[enDeps < 0]  <- -enDeps[enDeps < 0]
#hist(enDeps)
write.table(round(enDeps /(extraSize/nsize)),
            paste0(tDir, "enCopies"),
            col.names = F,
            row.names = F,
            quote = F)

cat("Done. sim.R\n")