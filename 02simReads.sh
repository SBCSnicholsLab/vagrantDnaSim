
#rl=150 # read length
rl=76
i=0
echo "SIMREADS##################################"
echo "Generating genome FASTA files (containig extranuclear DNA and inserts)..."
cat $1enCopies | while read line
do
  i=$((i+1))
  if test -f $(printf "$1"ref%03d.fa $i); then
   rm  $(printf "$1"ref%03d.fa $i)
  fi

  cat $1genome.fa > $(printf "$1"ref%03d.fa $i)
  #for c in {1..$line} #does not work
  for (( c=1; c<=$line; c++ ))
  do
    cat $1extraNuc.fa | tail -n 1 >> $(printf "$1"ref%03d.fa $i)
  done
done

i=0
cat "$1"pairCounts | while read line
do
  i=$((i+1))
  echo "Number of reads is $line"
  ~/git_repos/wgsim/wgsim -1 $rl -2 $rl -r 0 -R 0 -e 0.001 -N $line $(printf "$1"ref%03d.fa $i) $(printf "$1"s%03d.fw.fq $i) $(printf "$1"s%03d.rw.fq $i)
done
echo "Done. simReads"
echo ""
