#!/bin/sh
#PBS -l nodes=1:ppn=8,mem=50gb,vmem=52gb
#PBS -l walltime=60:00:00
#PBS -N genome
#PBS -m a
#PBS -j oe
#PBS -M X@Xmail.com

### By MinTu, 2020
### Waksman Genomics Core
### Rutgers University

### Function---1) Alternative 3'UTR length; 2) Filters; 3) novel 3'UTR fasta 4)3'UTR RNA motif 

# Set your workzone
PRE='/ingens/home/mintu/'

cd $PRE
mkdir genome results

cd $PRE/genome

# transform genome.gff3 to bed file. 
$PRE/softwares/ToGenePred/gff3ToGenePred $PRE/inputs/genome.gff3 ./genome.bed

# locates genome 3'UTR 
awk '{print $1,$12,$3,$4,$5,$6,$7,$9,$10}' genome.bed | awk '{split($8,a,","); split($9,b,","); size=""; for(i=1;i<length(a)-1;i++) {size=size""(a[i+1]-b[i])","} print$1,$2,$3,$4,$5,$6,$7,$8,$9,size }' | awk '{OFS="\t";if($10==""){$10=0} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > genome.intron.txt
cat genome.intron.txt | awk '{OFS="\t"; if($3=="+") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > genome+.intron.txt
cat genome.intron.txt | awk '{OFS="\t"; if($3=="-") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > genome-.intron.txt
cat genome+.intron.txt | awk '{split($8,a,",");split($9,b,",");A=0;B=0; for (i=1;i<length(a);i++) if($7>=a[i]&&$7<=b[i]){A=i;B=length(a)-1-i} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,A,B}'|awk '{OFS="\t"; split($10,a,","); A=$5-$7; for(i=$11;i<length(a);i++) {A-=a[i]} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,A}'> genome+.3UTR.txt
cat genome-.intron.txt |awk '{split($8,a,","); split($9,b,",");A=0;B=0; for (i=1;i<length(a);i++) if($6>=a[i]&&$6<=b[i]){A=i;B=i-1} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,A,B}'| awk '{OFS="\t"; split($10,a,",");A=$6-$4; for(i=$12;i>=0;i--){A-=a[i]} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,A}' > genome-.3UTR.txt
cat genome+.3UTR.txt genome-.3UTR.txt > genome.3UTR.txt

# genome's intronless.txt has 0/1 intron in non-3UTR region; annotation_3UTR.txt has 2 or more introns in non-3UTR region.
cat genome.3UTR.txt | awk '{OFS="\t"; split($10,a,","); if(length(a)-1-$12>=2) print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > annotation_3UTR.txt
cat genome.3UTR.txt | awk '{OFS="\t"; split($10,a,","); if(length(a)-1-$12<2) print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > intronless.txt

#  locates merged-3utr 
cat annotation_3UTR.txt | awk '{OFS="\t"; print$1,$2,$4,$5}' | sort > annotation_3UTR.sort.txt
sort genome.3UTR.txt > genome.3UTR.sort.txt
sort annotation_3UTR.txt > annotation_3UTR.sort.txt
cat genome.bed| awk '{OFS="\t"; print $1, $6, $7}' | sort > cds.sort


