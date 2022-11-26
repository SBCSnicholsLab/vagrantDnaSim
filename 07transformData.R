######################################
# Prepare data for general method ####
######################################

# set constants
nSamp = commandArgs(trailingOnly = T)[1]
tDir = commandArgs(trailingOnly = T)[2]
extraNucLen <- 16000


# a utility function to rename the mapping depth vals
getNums <- function(x){
  a <- x[,3]
  names(a) <- sapply(x[,1], function(y){
    substr(strsplit(y, "/")[[1]][[4]],2,4)
  })
  return(a)
}
nMapped <- read.table(paste0(tDir, "bpMapped"),
                      sep = "\t",
                      stringsAsFactors = F)
#head(nMapped)

nTot <- read.table(paste0(tDir, "bpTotal"),
                   sep = "\t",
                   stringsAsFactors = F)
#head(nTot)

nMapped <- getNums(nMapped)
nTot <- getNums(nTot)

#mappingProp <- nMapped / nTot
##hist(mappingProp, breaks=20)


#################
# Genotype calls

gts <- read.table(paste0(tDir, "genotypes.csv"),
                  header = T,
                  check.names = F)
#head(gts)
# counts of each allele, some are missing (-1)
#table(as.vector(unlist(gts[,1:nSamp])))

# remove lowest coverage samples
gtsHC <- gts # none removed

# visualise missingness
#image(as.matrix(gtsHC[,1:nSamp]) == -1)

# remove individuals with missing data
#indMissing <- apply(gtsHC, 2, function(x) -1 %in% x)
#gtsHC <- gtsHC[,!indMissing]
#head(gtsHC)

# # remove any site with missing data
# posMissing <- gtsHC$POS[apply(gtsHC, 1, function(x) -1 %in% x)]
# gtsHC <- gtsHC[,!gtsHC$POS %in% posMissing]

nindFilt <- ncol(gtsHC)-1
###############
# Allele counts

allCounts <- read.table(paste0(tDir, "alleleCounts.csv"))
#head(allCounts)
allCountsHC <- allCounts # opportunity to remove individuals if desired

# remove sites with missing alleles identified above
#allCountsHC <- allCountsHC[,!c(indMissing, F)]
#head(allCountsHC)
# split allele counts into separate DFs, one per allele
gtList <- lapply(1:4, function(x) {
  allCountsHC[allCountsHC$allele==x,]
})

#lapply(gtList, dim)
# 12442*4 = 49768

# flatten into four vectors
#rm(c) # In case somebody named a variable c, which would shadow the c function.
gtVectors <- lapply(gtList, function(x){
  do.call(c,x[,1:nindFilt])
})
#index of samples
sampleIndex <- rep(colnames(gtsHC)[1:nindFilt], each=nrow(gtList[[1]]))
#str(gtVectors)

# long(er) data frame for allele counts with sample index
#allDF <- data.frame(sample=sampleIndex, pos = paste0("S",rep(sprintf("%05d",gtsHC$POS), nindFilt)),
#                    ref=gtVectors[[1]], alt=gtVectors[[2]] + gtVectors[[3]] + gtVectors[[4]])
allDF <- data.frame(sample=sampleIndex, pos = rep(gtsHC$POS, nindFilt),
                    ref=gtVectors[[1]], alt=gtVectors[[2]] + gtVectors[[3]] + gtVectors[[4]])
#head(allDF)
#str(allDF)

# Depth of the alt alleles (summed)
siteDeps <- rowSums(allDF[,3:4])
altRatio <- allDF$alt / siteDeps



# DF that with all relevant data for the analysis:
mainDF <- data.frame(Sample = as.character(allDF$sample),
                     Position=allDF$pos,
                     AltProp=altRatio,
                     DP=allDF$ref+allDF$alt,
                     stringsAsFactors = F, row.names = NULL)
#head(mainDF)
# Add mapping rate, make a DF so we can use "merge"
#mappingProp

propDF <- data.frame(Sample = names(nMapped),
                     #mappingrate=mappingProp,
                     nMapped,
                     nTot,
                     stringsAsFactors = F)
#head(propDF)
mainDF <- merge(mainDF, propDF, by="Sample")

#head(mainDF)
mainDF$ylog <- log(mainDF$AltProp)
mainDF$ylog[is.infinite(mainDF$ylog)] <- NA # replace -Inf by NA, makes plotting easier
#mainDF$xnqlogis <- -qlogis(mainDF$mappingrate)
#mainDF$xlog <- log(mainDF$mappingrate)

#head(mainDF)

# remove NA lines
mainDF <- mainDF[!is.na(mainDF$AltProp),]
mainDF <- mainDF[!is.na(mainDF$ylog),]
#mainDF <- mainDF[!is.na(mainDF$xnqlogis),]

# Write out data
write.table(mainDF,
            paste0(tDir, "transformedData.csv")
            )


