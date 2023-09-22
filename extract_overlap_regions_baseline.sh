Input_dir=$1
GRNdir=$2
genome=$3
cd $Input_dir
cat ATAC.txt|cut -f 1 |sed '1d' |sed 's/:/\t/g'| sed 's/-/\t/g' > Region.bed
ml bedtools2/2.26.0-gcc/9.5.0
for i in $(seq 1 22); do
bedtools intersect -a "$GRNdir/${genome}_Peaks_chr$i.bed" -b Region.bed -wa -wb > "Region_overlap_chr$i.bed"
done
i=X
bedtools intersect -a "$GRNdir/${genome}_Peaks_chr$i.bed" -b Region.bed -wa -wb > "Region_overlap_chr$i.bed"

