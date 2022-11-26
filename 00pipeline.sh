# USAGE: 00pipeline.sh outPref GS insRate
set -e
eval "$(conda shell.bash hook)"
wDir="$(mktemp -d)/" # add slash (which does not need to be escaped it seems)
outPrf=$1
GS=$2
IR=$3
SD=0.05

#GS=500000000


echo "Temp dir is $wDir"
echo "Ins rate is $IR"

conda activate vagrantSim
nInd=20
Rscript 01simSparse.R $nInd $wDir $GS $IR $outPrf $SD && bash 02simReads.sh $wDir && bash 03map.sh $nInd $wDir && bash 04sortIndexCall.sh $wDir && bash 06mappingNumbers.sh $wDir
#Rscript 01simSparseBeta.R $nInd $wDir $GS $IR $outPrf $SD && bash 02simReads.sh $wDir && bash 03map.sh $nInd $wDir && bash 04sortIndexCall.sh $wDir && bash 06mappingNumbers.sh $wDir


echo "Extracting variants and allele counts from VCF..."
python 05processVCF.py $wDir

echo "Generating CSV for rainbow plot..."

Rscript 07transformData.R $nInd $wDir

conda deactivate

echo "Getting output data and saving to "$outPrf"transfomredData.csv..."
cp $wDir/transformedData.csv "$outPrf"transformedData.csv


echo "Running fit..."
Rscript 08runFit.R $wDir $outPrf


cp $wDir$outPrf.log .

#echo "Removing $wDir..."
#rm -rf $wDir
echo "ALL DONE."
