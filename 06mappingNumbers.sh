echo "Getting mapping numbers..."
grep "total length:" $1*.stats > "$1"bpTotal
grep "bases mapped (cigar):" $1*.stats > "$1"bpMapped
echo "Done."
echo ""

