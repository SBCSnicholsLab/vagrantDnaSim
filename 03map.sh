
echo "MAPPING#################################################"
# $1 comes in from 00pipiline script. It is the number of samples.
# $2 is the temp dir
bwa-mem2 index "$2"extraNuc.fa
for (( j=0; j<=$1; j++ ))
do
  i=$(printf '%03d' $j)
  echo "Mapping s$i..."
  bwa-mem2 mem -t 2 -R '@RG\tID:'$i'\tSM:'$i "$2"extraNuc.fa "$2"s$i.fw.fq "$2"s$i.rw.fq | samtools view -b -o "$2"s"$i"r.bam
done

echo "Done. map"
echo ""
