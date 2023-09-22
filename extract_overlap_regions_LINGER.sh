Input_dir=$1
GRNdir=$2
genome=$3
cd $Input_dir
ml bedtools2/2.26.0-gcc/9.5.0
cat ATAC.txt|cut -f 1 |sed '1d' |sed 's/:/\t/g'| sed 's/-/\t/g' > Region.bed
if [ $genome = "hg38"  ]; 
then 
bedtools intersect -a Region.bed -b  $GRNdir/hg38_hg19_pair.bed -wa -wb |awk '{print $7"\t"$8"\t"$9"\t"$1"\t"$2"\t"$3}'> match_hg19_peak.bed # the first 3 are the hg19 region origianl ; the last 3 is  the query data region
fi
if [ $genome = "hg19"  ]
then
bedtools intersect -a Region.bed -b  $GRNdir/hg19_hg38_pair.bed -wa -wb |awk '{print $7"\t"$8"\t"$9"\t"$1"\t"$2"\t"$3}'> match_hg19_peak.bed # the first 3 are the hg19 region origianl ; the last 3 is  the query data region
fi

bedtools intersect -a match_hg19_peak.bed -b $GRNdir/RE_gene_corr_hg19.bed -wa -wb |awk '$2 == $8 && $3 == $9'|awk '{print $1":"$2"-"$3"\t"$4":"$5"-"$6"\t"$10}'|sort|uniq > hg19_Peak_hg19_gene_u.txt # first is hg19  original and second is the overlaped  qurey rengion ; 3 is the target gene 
# the first column is hg the second col
## motif enrhchment mapping
ml bedtools2/2.26.0-gcc/9.5.0
for i in $(seq 1 22)
do
bedtools intersect -a  match_hg19_peak.bed -b $GRNdir/MotifTarget_matrix_chr$i.bed -wa -wb  |awk '$1==$7 && $2==$8 && $3==$9 { print }' |awk '{print $4":"$5"-"$6"\t"$7":"$8"-"$9}' > MotifTarget_hg19_hg38_chr$i.txt
done
i=X
bedtools intersect -a  match_hg19_peak.bed -b $GRNdir/MotifTarget_matrix_chr$i.bed -wa -wb  |awk '$1==$7 && $2==$8 && $3==$9 { print }' |awk '{print $4":"$5"-"$6"\t"$7":"$8"-"$9}' > MotifTarget_hg19_hg38_chr$i.txt
for i in $(seq 1 22); do
bedtools intersect -a "$GRNdir/${genome}_Peaks_chr$i.bed" -b Region.bed -wa -wb > "Region_overlap_chr$i.bed"
done
i=X
bedtools intersect -a "$GRNdir/${genome}_Peaks_chr$i.bed" -b Region.bed -wa -wb > "Region_overlap_chr$i.bed"