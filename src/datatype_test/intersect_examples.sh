echo A.bed
cat A.bed

echo B.bed
cat B.bed
echo "intersect"
bedtools intersect -a A.bed -b B.bed

echo B2.bed
cat B2.bed
echo "intersect"
bedtools intersect -a A.bed -b B2.bed

echo B3.bed
cat B3.bed
echo intersect
bedtools intersect -a A.bed -b B3.bed

echo B4.bed
cat B4.bed
echo intersect
bedtools intersect -a A.bed -b B4.bed

echo B5.bed
cat B5.bed
echo intersect
bedtools intersect -a A.bed -b B5.bed
