# OR_probe_design_workflow
## This repo includes scripts and details for each step. (script is crude though :( )
### 47 cattle assemblies are selected across Bos species and their OR regions are predicted by G2OR. Coordinates of gene boundary are used as OR regions in this design.
### Step 1. Build OR senquence pool for each assembly
```
for subject in `cat genome_list`  ; do
mkdir ${subject} && cd ${subject}
ln -s /home/wsj/workplace/OR/g2OR/part1/results/Mammalia/${subject}/ORannotation_itera2_final_func_dna_ORs.fasta ./
ln -s /home/wsj/workplace/OR/g2OR/part1/results/Mammalia/${subject}/ORannotation_itera2_final_low-quality.fasta ./
ln -s /home/wsj/workplace/OR/g2OR/part1/results/Mammalia/${subject}/ORannotation_itera2_final_pseu_ORs.fasta ./
ln -s /home/wsj/workplace/OR/g2OR/part1/ref_unzip/${subject}_genomic.fna ./
samtools faidx ORannotation_itera2_final_func_dna_ORs.fasta
samtools faidx ORannotation_itera2_final_low-quality.fasta 
samtools faidx ORannotation_itera2_final_pseu_ORs.fasta 
grep -v hit ORannotation_itera2_final_pseu_ORs.fasta.fai |awk -F'_' 'BEGIN{OFS="\t"} {print $1,$2,$3,"psedo"NR}'  |sort -k1,1 -k2,2n > ORannotation_itera2_final_pseu_ORs.bed1
grep hit ORannotation_itera2_final_pseu_ORs.fasta.fai |awk -F'_' 'BEGIN{OFS="\t"} {print $2,$3,$4,"psedo"NR}' |sort -k1,1 -k2,2n > ORannotation_itera2_final_pseu_ORs.bed2
awk -F'_' 'BEGIN{OFS="\t"} {print $1,$2,$3,"low"NR}' ORannotation_itera2_final_low-quality.fasta.fai |sort -k1,1 -k2,2n > ORannotation_itera2_final_low_quality.bed
awk -F'_' 'BEGIN{OFS="\t"} {print $2,$3,$4,$1}' ORannotation_itera2_final_func_dna_ORs.fasta.fai |sort -k1,1 -k2,2n > ORannotation_itera2_final_func_dna_ORs.bed
cat ORannotation_itera2_final_pseu_ORs.bed1 ORannotation_itera2_final_pseu_ORs.bed2 ORannotation_itera2_final_low_quality.bed ORannotation_itera2_final_func_dna_ORs.bed > gene_body.pool.bed
sh ../modify_bed_coor.sh gene_body.pool.bed gene_body.pool.reorder.bed
sort -k1,1 -k2,2n gene_body.pool.reorder.bed > gene_body.pool.reorder.sort.bed 
bedtools merge -i gene_body.pool.reorder.sort.bed  -c 4 -o collapse -delim "-" > gene_body.pool.merge.reorder.sort.bed
awk 'BEGIN{OFS="\t"} {print $1,$2-2000,$2-1,$4"-up"}' gene_body.pool.merge.reorder.sort.bed |awk 'BEGIN{OFS="\t"} $2<0 {$2 = 0} {print}' > gene_body.pool.merge.reorder.sort.up.bed
awk 'BEGIN{OFS="\t"} {print $1,$3+1,$3+2000,$4"-down"}' gene_body.pool.merge.reorder.sort.bed > gene_body.pool.merge.reorder.sort.down.bed
cat gene_body.pool.merge.reorder.sort.up.bed gene_body.pool.merge.reorder.sort.down.bed |sort -k1,1 -k2,2n > up_down.pool.bed
bedtools getfasta -fi ${subject}_genomic.fna -bed up_down.pool.bed -name > up_down.pool.fasta
bedtools getfasta -fi ${subject}_genomic.fna -bed gene_body.pool.merge.reorder.sort.bed -name > gene_body.pool.fasta
cd ../
done
```

