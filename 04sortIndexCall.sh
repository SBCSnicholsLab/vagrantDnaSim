
echo "SORTINDEXCALL###################################"
for i in "$1"*r.bam
do
  echo "Sorting $i..."
  samtools sort -o ${i%.*}s.bam $i
done

echo "Indexing BAM files..."
for i in "$1"*rs.bam; do samtools index $i; done

echo "Calling variants..."
freebayes -f "$1"extraNuc.fa -iXu --haplotype-length -1 --min-alternate-fraction 0.001 --min-alternate-count 1 --pooled-continuous -p 1  "$1"*s.bam> "$1"sim.vcf

for i in "$1"*s.bam; do echo "Processing $i..."; samtools stats $i > $i.stats; done

# echo "Removing FA and FQ files..."
# rm *fa *fq
# echo "Removing FQ files..."
# rm *fq
echo "Done. sortIndexCall"
echo ""
