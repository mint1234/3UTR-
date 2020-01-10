#!/bin/sh
#PBS -l nodes=1:ppn=8,mem=50gb,vmem=52gb
#PBS -l walltime=60:00:00
#PBS -N 3r-${SN} 
#PBS -m a
#PBS -j oe
#PBS -M X@Xmail.com

### By MinTu, 2020
### Waksman Genomics Core 
### Rutgers University

# Function:(1)find alternative 3'UTR and caculate the length; (2) Filters for false positive; (3) output novel 3'UTR fasta 4)predict secondary structures of noval 3'UTR 

# Before running, please set a "inputs" folder with the following files: 1) genome.gff3; 2)genome.fa; 3)bam files; 4)isoforms_fpkm.tracking from cuffdiff results; 4)merged.gtf files from cuffmerge results.

qsub -I -v interactive

# default settings:  working zone and parameters

PRE='/ingens/home/mintu/'

A=1  # replecates RPKM value
B=30 # 3'UTR reads value

#(1) find alternative 3'UTR and caculate the length

cd $PRE
mkdir ${SN}-outputs 2nd-${SN}  

cd $PRE/${SN}-outputs

# locates each samples' 3'UTR based on SampleName-merged.gtf
$PRE/ToGenePred/gtfToGenePred $PRE/inputs/${SN}-merged.gtf ./${SN}-merged.bed

# Initial filter: fliter out overlapping transcripts and one contig correspond to several transcripts

# Choose classcode "=,j" and fpkm from each sample's isoforms.fpkm_tracking
cat $PRE/inputs/${SN}-isoforms.fpkm_tracking | awk '{OFS="\t"; if($2=="=" || $2=="j") print$1,$2,$3,$5,$7,$10,$13,$14,$17,$18,$21}' > ./${SN}-classcode-fpkm.txt
sort ${SN}-merged.bed > ./${SN}-merged.sorted.txt
sort ${SN}-classcode-fpkm.txt > ./${SN}-classcode-fpkm.sorted.txt
join -o 2.1 2.3 2.4 2.2 1.3 2.5 2.6 2.7 2.8 2.9 1.4 1.5 1.9 1.10 2.10 2.11 ${SN}-merged.sorted.txt ${SN}-classcode-fpkm.sorted.txt > ./${SN}-merged-classcode.txt
cat ${SN}-merged-classcode.txt | awk '{OFS="\t"; if($4=="=") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > ./${SN}-merged=.txt
cat ${SN}-merged-classcode.txt | awk '{OFS="\t"; if($4=="j") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' > ./${SN}-mergedj.txt
cat ${SN}-merged=.txt | awk '{OFS="\t"; print $2,$3,$1,$4,$11,$12}' | sort  > ./${SN}-merged=.sort.txt
join -o 1.3 1.1 1.2 1.4 1.5 1.6 2.3 2.4 ${SN}-merged=.sort.txt ../genome/annotation_3UTR.sort.txt > ./${SN}-merged=-annotation.txt

# If classcode= but trans-start and trans-end do not match, join in mergej.txt
cat ${SN}-merged=-annotation.txt | awk '{OFS="\t"; if($5!=$7 || $6!=$8) print$1,$2,$3,$4,$5,$6,$7,$8}' > ./${SN}-merged=-annotation-selected.txt
sort ${SN}-merged=-annotation-selected.txt > ./${SN}-merged=-annotation-selected.sort.txt
sort ${SN}-merged=.txt > ./${SN}-merged=.sort.txt
join -o 1.1 1.2 1.3 1.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12 2.13 2.14 2.15 2.16 ${SN}-merged=-annotation-selected.sort.txt ${SN}-merged=.sort.txt > ./${SN}-merge=.selected.txt
awk '{OFS="\t"; print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' ${SN}-merge=.selected.txt > ./${SN}-merge=.selected-OFS.txt
cat ${SN}-mergedj.txt ${SN}-merge=.selected-OFS.txt > ./${SN}-merged-selected.txt

# Add cds from genome file
awk '{OFS="\t"; print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' ${SN}-merged-selected.txt| sort > ./${SN}-merged-selected.sort.txt
join -o 1.2 1.1 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 2.6 2.7 1.13 1.14 1.15 1.16 ${SN}-merged-selected.sort.txt ../genome/genome.3UTR.sort.txt > ./${SN}-selected-cds.txt
 
# 3UTR length
cat ${SN}-selected-cds.txt | awk '{OFS="\t"; split($15,a,","); split($16,b,","); size=""; for(i=1;i<length(a)-1;i++) {size=size""(a[i+1]-b[i])","} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,size,$17,$18}' > ./${SN}-selected-cds.introns.txt
 
# Delete no intron
cat ${SN}-selected-cds.introns.txt | awk '{OFS="\t"; if($17!="") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' > ./${SN}-cds-introns-selected.txt
cat ${SN}-cds-introns-selected.txt | awk '{OFS="\t"; if($5=="+") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' > ./${SN}-cds-introns-selected+.txt
cat ${SN}-cds-introns-selected.txt | awk '{OFS="\t"; if($5=="-") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' > ./${SN}-cds-introns-selected-.txt

# Delete those genes not in annotation.3UTR.txt
cat ${SN}-cds-introns-selected.txt | awk '{OFS="\t"; print$2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}'|sort > ./${SN}-cds-introns-selected.sort.txt
join -o 1.2 2.1 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 1.17 2.3 2.10 2.11 2.12 2.13 1.18 1.19 ${SN}-cds-introns-selected.sort.txt ../genome/annotation_3UTR.sort.txt > ./${SN}-annotation_3UTR.txt
 
# Delete 3UTR intron in annotation 
cat ${SN}-annotation_3UTR.txt | awk '{OFS="\t";if($18=="+") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}' > ./${SN}-annotation+_3UTR.txt
cat ${SN}-annotation_3UTR.txt | awk '{OFS="\t";if($18=="-") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}' > ./${SN}-annotation-_3UTR.txt
cat ${SN}-annotation+_3UTR.txt | awk '{OFS="\t"; split($19,a,","); size=""; for(i=1;i<$20;i++) {size=size""(a[i])","} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,size,$23,$24 }' | awk '{OFS="\t"; if($23=="," || $23=="0,") {$23=0} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' | sed 's/,,/,/g'  > ${SN}-annotation+_3UTR-newintron.txt
cat ${SN}-annotation-_3UTR.txt | awk '{OFS="\t"; split($19,a,","); size=""; for(i=$21;i<length(a);i++) {size=size""(a[i+1])","} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,size,$23,$24 }' | awk '{OFS="\t"; if($23=="," || $23=="0,") {$23=0} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' | sed 's/,,/,/g'  >  ${SN}-annotation-_3UTR-newintron.txt
cat ${SN}-annotation+_3UTR-newintron.txt ${SN}-annotation-_3UTR-newintron.txt > ${SN}-annotation_3UTR-newintrons.txt

# Remove merged introns length = 0
cat ${SN}-annotation_3UTR-newintrons.txt | awk '{if ($17!=0) print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' > ${SN}-annotation_3UTR-newintrons-selected.txt

# merged 3UTR length
cat ${SN}-annotation_3UTR-newintrons.txt | awk '{OFS="\t"; if($5=="+") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' > ${SN}-annotation_3UTR-newintrons+.txt
cat ${SN}-annotation_3UTR-newintrons.txt | awk '{OFS="\t"; if($5=="-") print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25}' > ${SN}-annotation_3UTR-newintrons-.txt
cat ${SN}-annotation_3UTR-newintrons+.txt | awk '{split($15,a,",");split($16,b,",");A=0;B=0; for (i=1;i<length(a);i++) if($14>=a[i]&&$14<=b[i]){A=i;B=length(a)-1-i}print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,A,B,$18,$19,$20,$21,$22,$23,$24,$25}'| awk '{OFS="\t"; split($17,a,","); A=$12-$14; for(i=$18;i<length(a);i++) {A-=a[i]} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,A,$20,$21,$22,$23,$24,$25,$26,$27}' > ${SN}+-3UTR-annotation_3UTR-newintrons.txt
cat ${SN}-annotation_3UTR-newintrons-.txt | awk '{split($15,a,",");split($16,b,",");A=0;B=0; for (i=1;i<length(a);i++) if($13>=a[i]&&$13<=b[i]){A=i;B=i-1} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,A,B,$18,$19,$20,$21,$22,$23,$24,$25}'| awk '{OFS="\t"; split($17,a,",");A=$13-$11; for(i=$19;i>=0;i--){A-=a[i]} print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,A,$20,$21,$22,$23,$24,$25,$26,$27}' >${SN}--3UTR-annotation_3UTR-newintrons.txt

# Delete merged 3UTR <0
cat ${SN}+-3UTR-annotation_3UTR-newintrons.txt | awk '{OFS="\t"; if($20>=0) print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}' > ${SN}+-3UTR-annotation_3UTR-newintrons.revised.txt
cat ${SN}--3UTR-annotation_3UTR-newintrons.txt | awk '{OFS="\t"; if($20>=0) print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}' > ${SN}--3UTR-annotation_3UTR-newintrons.revised.txt

# Compare intron chain
awk '{OFS="\t"; split($17,a,","); split($26,b,",");j=length(b)-1; for (i=length(a)-1;i>=1;i--) if(a[i]==b[j] && a[i-1]==b[j-1] && a[i-2]==b[j-2]) print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,i,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}' ${SN}+-3UTR-annotation_3UTR-newintrons.revised.txt> ${SN}+.result.txt
awk '{OFS="\t"; split($17,a,","); split($26,b,",");j=1; for (i=1;i<length(a);i++) if(a[i]==b[j] && a[i+1]==b[j+1] && a[i+2]==b[j+2])print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,i,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}' ${SN}--3UTR-annotation_3UTR-newintrons.revised.txt > ${SN}-.result.txt

# Selected result without comparing 3UTR and fpkm status
cat ${SN}+.result.txt ${SN}-.result.txt > ${SN}.result.txt

# Extract 3UTR length from result
cat ${SN}.result.txt | awk '{OFS="\t"; print$1,$2,$3,$4,$7,$8,$9,$10,$21,$26,$28,$29}' > ${SN}.result.fpkm-3UTR.txt

# Selected by fpkm status and value,delete both=0 
sed '/NOTEST/d' ${SN}.result.fpkm-3UTR.txt |awk '{OFS="\t"; if ($5!=0 && $7!=0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > ${SN}.3UTR.result.txt

# Selected by fpkm mean, stand deviation and coefficient of variation
cat ${SN}.3UTR.result.txt | awk '{OFS="\t"; mean=($5+$7+$11)/3; sd=sqrt((($5-mean)**2+($7-mean)**2+($11-mean)**2)/3); cv=sd/mean ; print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,mean,sd,cv,$11,$12}' > ${SN}.3utr.cv.txt
 
# Expression filter, each replicate RPKM >=1, A=1
cat ${SN}.3utr.cv.txt | awk '{OFS="\t"; if ($5>=A && $7>=A && $14>=A) print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' > ${SN}.result-fpkm-1.txt


# Selected only 1 tcons corresponding transcripts' ID

cat ${SN}.result-fpkm-1.txt | awk '{print $2}' > ${SN}.transcripts.txt
cat ${SN}.transcripts.txt | sort | uniq -c > ${SN}.tem.txt
cat ${SN}.tem.txt | awk '{OFS="\t"; if ($1==1) print $2}' > ${SN}.ID.txt
cat ${SN}.tem.txt | awk '{OFS="\t"; if ($1!=1) print $2}' > ${SN}.ID.delete.txt
cat ${SN}.result-fpkm-1.txt | awk '{OFS="\t"; print$2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' | sort > ${SN}.sort.txt
sort ${SN}.ID.txt > ${SN}.ID.sort.txt
join -o 2.2 1.1 2.3 2.4 2.5 2.6 2.7 2.8 2.9 2.10 2.11 2.12 2.13 2.14 2.15 ${SN}.ID.sort.txt ${SN}.sort.txt > ${SN}.tcons1-1transcripts.txt
awk '{OFS="\t"; print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' ${SN}.tcons1-1transcripts.txt > ${SN}.tcons1-1Transcripts.txt

# ORF prediction filter : filter out the dismatches with coding sequence

# generate alignment gff3 formatted Output
$PRE/softwares/TransDecoder-2.0.1/util/cufflinks_gtf_to_alignment_gff3.pl $PRE/inputs/${SN}-merged.gtf > ${SN}.gff3

# generate transcripts fasta file
$PRE/softwares/TransDecoder-2.0.1/util/cufflinks_gtf_genome_to_cdna_fasta.pl $PRE/inputs/${SN}-merged.gtf $PRE/inputs/genome.fa > ${SN}.fasta

# Extract the long ORFs
$PRE/softwares/TransDecoder-2.0.1/TransDecoder.LongOrfs -t ${SN}.fasta

# TransDecoder 3'UTR length
sed '/three_prime_UTR/!d' ${SN}.fasta.transdecoder_dir/longest_orfs.gff3 > ${SN}-3utr.txt
cat ${SN}-3utr.txt | awk '{OFS="\t"; print $1,$2,$3,$4,$5,$7,$5-$4+1}' > ${SN}-3utrlength.txt

# compare novel 3'UTR length
sort ${SN}.tcons1-1Transcripts.txt > ${SN}.tcons1-1Transcripts.sort.txt
sort ${SN}-3utrlength.txt > ${SN}-3utrlength.sort.txt
join -o 1.1 1.2 1.4 1.9 1.10 2.2 2.3 2.4 2.5 2.6 2.7 ${SN}.tcons1-1Transcripts.sort.txt ${SN}-3utrlength.sort.txt > ${SN}-result.txt
awk '{OFS="\t"; print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' ${SN}-result.txt > ${SN}-results.txt
cat ${SN}-results.txt | awk '{OFS="\t"; if($4==$11) print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' > ${SN}-compare.txt
cat ${SN}-compare.txt| awk '{OFS="\t"; print$1,$2,$3,$4,$5,$6,$7,$10,$11}'|sort | uniq -c > ORF${SN}.txt
awk '{OFS="\t"; print $2,$3,$4,$5,$6,$7,$9,$10}' ORF${SN}.txt | sort > ORF${SN}.sort.txt
join -o 1.1 1.2 2.3 1.3 1.4 1.5 1.6 1.7 1.8 2.5 2.7 2.14 2.13 ORF${SN}.sort.txt ${SN}.tcons1-1Transcripts.sort.txt | awk '{OFS="\t"; print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > ${SN}-ORFselected.txt
cat ${SN}-ORFselected.txt| awk '{OFS="\t"; if(($5-$6)!=0)print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > ${SN}.change.txt

cp ${SN}.change.txt $PRE/results
# line annotation: 1)TCONS 2)transcripts 3) gene ID 4) class code 5) 3UTR length 6)annotated 3UTR length  7) transdecoder 8)+/- transdecoder 9)3utr length transdecoder 10)FPKM-1 11)FPKM-2 12)FPKM-3 

# 3' UTR reads filter >=B,B=30 
$PRE/softwares/samtools sort $PRE/inputs/${SN}_R1-unique.bam -o ${SN}_R1-unique.sorted.bam
$PRE/softwares/samtools index ${SN}_R1-unique.sorted.bam
$PRE/softwares/samtools sort $PRE/inputs/${SN}_R2-unique.bam -o ${SN}_R2-unique.sorted.bam
$PRE/softwares/samtools index ${SN}_R2-unique.sorted.bam
$PRE/softwares/samtools sort $PRE/inputs/${SN}_R3-unique.bam -o ${SN}_R3-unique.sorted.bam
$PRE/softwares/samtools index ${SN}_R3-unique.sorted.bam

sort ${SN}.change.txt > ${SN}.change.sort
sort ${SN}-selected-cds.txt > ${SN}-selected-cds.sort
join -o 2.6 1.1 1.2 1.3 2.5 2.11 2.12 2.13 2.14 ${SN}.change.sort ${SN}-selected-cds.sort > ${SN}.bed
cat ${SN}.bed |awk '{OFS="\t"; if($5=="+") print$1,$9,$7,$3,0,$5;else if($5=="-") print$1,$6,$8,$3,0,$5 }'> ${SN}-3utr.bed

# TCONS id matches TRANSCRIPT id
awk '{OFS="\t"; print $4,$7}' ${SN}-3utr.bed > ${SN}-tconsID-transID.txt
sed 's/:/\t/g' ${SN}-3utr.bed > ${SN}-3UTR.bed
cat ${SN}-3UTR.bed | awk '{OFS="\t"; print $1,$3,$4,$6,$7,$8}' > ${SN}-3utr-revise.bed
$PRE/softwares/bedtools multicov -bams ${SN}_R1-unique.sorted.bam -bed ${SN}-3utr-revise.bed > ${SN}_R1-3utr-counts.txt
$PRE/softwares/bedtools multicov -bams ${SN}_R2-unique.sorted.bam -bed ${SN}-3utr-revise.bed > ${SN}_R2-3utr-counts.txt
$PRE/softwares/bedtools multicov -bams ${SN}_R3-unique.sorted.bam -bed ${SN}-3utr-revise.bed > ${SN}_R3-3utr-counts.txt
awk '{OFS="\t"; if($7>=B) print $1,$2,$3,$4,$6,$7}' ${SN}_R1-3utr-counts.txt > ${SN}_R1-3utr-counts-30.txt
awk '{OFS="\t"; if($7>=B) print $1,$2,$3,$4,$6,$7}' ${SN}_R2-3utr-counts.txt > ${SN}_R2-3utr-counts-30.txt
awk '{OFS="\t"; if($7>=B) print $1,$2,$3,$4,$6,$7}' ${SN}_R3-3utr-counts.txt > ${SN}_R3-3utr-counts-30.txt

cp ${SN}_R*-3utr-counts-30.txt $PRE/results
# line annotation: 1-chr 2-3utr_start 3-3utr_end 4-transcript 5-+/-  6-3'UTR reads

# Output alternative 3' UTR fasta 
cat $PRE/inputs/${SN}-isoforms.fpkm_tracking | awk '{OFS="\t"; if($2=="=" || $2=="j") print$1,$2,$3,$5,$7,$10,$13,$14,$17}' > ${SN}-ID.txt
cat ${SN}-ID.txt | awk '{OFS="\t"; print $1,$3}' |  sort > ${SN}-ID.sort
join -o 2.2 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 ${SN}-merged.sorted.txt  ${SN}-ID.sort > ${SN}-ID.bed
cat $PRE/genome/genome.bed| awk '{OFS="\t"; print $1, $6, $7}' | sort > cds.sort
sort ${SN}-ID.bed > ${SN}-ID.sort
join -o 1.1 1.2 1.3 1.4 1.5 1.6 2.2 2.3 1.9 1.10 1.11 ${SN}-ID.sort cds.sort > ${SN}-ID-cds.bed
cat ${SN}-ID-cds.bed | awk '{OFS="\t"; split($10,a,",");split($11,b,",");size="";for(i=1;i<=length(a)-1;i++) {size=size""(b[i]-a[i])","} print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,size}' > ${SN}-exon-chain.bed
cat ${SN}-exon-chain.bed | awk '{OFS="\t"; split($10,a,",");split($11,b,","); for (i=1;i<=length(a);i++) if ($4=="+" && $8>=a[i] && $8<=b[i])print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,i; else if ($4=="-" && $7>=a[i] && $7<=b[i]) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,i}' >${SN}-whichexon.bed
cat ${SN}-whichexon.bed | awk '{OFS="\t"; split($12,a,",");A=0;for(i=1;i<=length(a)-1;i++){A+=a[i]} print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,A}'> ${SN}-totalexon.bed
sed '/+/!d' ${SN}-totalexon.bed > ${SN}-totalexon+.bed
sed '/+/d' ${SN}-totalexon.bed > ${SN}-totalexon-.bed
cat ${SN}-totalexon+.bed | awk '{OFS="\t";split($10,a,",");split($12,c,","); U=0;A=0; for (i=1;i<$13;i++) if($4=="+") {A+=c[i];U=A+$8-a[$13]} print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,U}'> ${SN}-3utr+.bed
cat ${SN}-totalexon-.bed | awk '{OFS="\t";split($11,b,",");split($12,c,","); U=0;A=0; for (i=length(b);i>$13;i--) {A+=c[i];U=$14-A-(b[$13]-$7)} print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,U}'> ${SN}-3utr-.bed
cat ${SN}-3utr+.bed ${SN}-3utr-.bed | awk '{OFS="\t";print $2,$1}'|sort > ${SN}-3utr-ID.bed
cat ${SN}-3utr+.bed ${SN}-3utr-.bed | awk '{OFS="\t";print $2,$4,$15,$14}'|sort > ${SN}-3utr-1.bed
cat  ${SN}-3utr-1.bed | awk '{OFS="\t"; if($2=="+")print $1,$3,$4; else if($2=="-") print$1,$4-$3,$4,$5}'>${SN}-3utr.bed

# catch 3utr fasta
$PRE/softwares/gffread -w ${SN}.fa -g $PRE/inputs/genome.fa $PRE/inputs/${SN}-merged.gtf
$PRE/softwares/bedtools getfasta -fi ${SN}.fa -bed ${SN}-3utr.bed -fo ${SN}-3utr.fasta

# Output novel 3'UTR fasta 
awk '{print $4}' ${SN}_R1-3utr-counts-30.txt  |sort > ${SN}_R1-reads-30-ID.sort  
awk '{print $4}' ${SN}_R2-3utr-counts-30.txt  |sort > ${SN}_R2-reads-30-ID.sort
awk '{print $4}' ${SN}_R3-3utr-counts-30.txt  |sort > ${SN}_R3-reads-30-ID.sort
comm -1 -2 ${SN}_R1-reads-30-ID.sort ${SN}_R2-reads-30-ID.sort > ${SN}-reads-30-ID-R1R2.txt
comm -1 -2 ${SN}-reads-30-ID-R1R2.txt ${SN}_R3-reads-30-ID.sort > ${SN}-reads-30-ID.txt
sort ${SN}-reads-30-ID.txt > ${SN}-reads-30-ID.sort
awk '{OFS="\t";print $2, $1}' ${SN}.change.txt > ${SN}.change-id.txt  
sed 's/transcript://g' ${SN}.change-id.txt | sort > ${SN}.change-ID.sort 
join -o 1.2 ${SN}.change-ID.sort ${SN}-reads-30-ID.sort > ${SN}.list 
sed 's/:/\t/g' ${SN}-3utr.fasta | awk '{print $1}' > ${SN}.fasta
perl $PRE/softwares/extractFromFasta.pl ${SN}.fasta list ${SN}.list > ${SN}-alter-3UTR.fasta
mv ${SN}-alter-3UTR.fasta $PRE/results

# predict secondary structures of noval 3'UTR with minimum free enery
cd $PRE/2nd-${SN}
sed 's/TCONS_//g' $PRE/results/${SN}-alter-3UTR.fasta > ${SN}-alter-3UTR-noTCONS.fasta
$PRE/softwares/RNAfold -p < ${SN}-alter-3UTR-noTCONS.fasta > ${SN}-alter





                                                                                                                                                                                                   





