Input_dir=$1
GRNdir=$2
genome=$3
cd $Input_dir
ml bedtools2/2.26.0-gcc/9.5.0
for i in $(seq 1 22); do
bedtools intersect -a "$GRNdir/${genome}_Peaks_chr$i.bed" -b Region.bed -wa -wb > "Region_overlap_chr$i.bed"
done
i=X
bedtools intersect -a "$GRNdir/${genome}_Peaks_chr$i.bed" -b Region.bed -wa -wb > "Region_overlap_chr$i.bed"

