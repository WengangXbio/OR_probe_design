# OR_probe_design_workflow
## This repo includes scripts and details for each step. (script is crude though :( )
### 47 cattle assemblies are selected across Bos species and their OR regions are predicted by G2OR. Coordinates of gene boundary are used as OR regions in this design.
### Step 1. Build OR senquence pool for each assembly
```
cd /home/wgzhang/probe_design/initial
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
#### merge sequence pool together
```
cd /home/wgzhang/probe_design/initial
for subject in `cat genome_list`  ; do
cd ${subject}
cp gene_body.pool.fasta ../fasta_collection/${subject}.gene_body.pool.fasta
cp up_down.pool.fasta ../fasta_collection/${subject}.up_down.pool.fasta
cd ../
done 
```

### Step 2. Rename fasta files and build chain
```
cd /home/wgzhang/probe_design/initial/fasta_collection
for subject in `cat genome_list`  ; do
ass=$(grep ${subject} genome.chain|cut -f2)
samtools faidx ${subject}.gene_body.pool.fasta 
awk -v a=${ass} '/^>/{print ">" a"_r" sprintf("%04d", NR); next}{print}'  ${subject}.gene_body.pool.fasta > ${subject}.gene_body.pool.rename.fasta 
samtools faidx ${subject}.gene_body.pool.rename.fasta 
paste -d ' ' <(awk '{print $1}' ${subject}.gene_body.pool.fasta.fai) <(awk '{print $1}' ${subject}.gene_body.pool.rename.fasta.fai) > ${subject}.gene_body.read.chain
samtools faidx ${subject}.up_down.pool.fasta 
awk -v a=${ass} '/^>/{print ">" a"_r" sprintf("%04d", NR); next}{print}'  ${subject}.up_down.pool.fasta > ${subject}.up_down.pool.rename.fasta 
samtools faidx ${subject}.up_down.pool.rename.fasta 
paste -d ' ' <(awk '{print $1}' ${subject}.up_down.pool.fasta.fai) <(awk '{print $1}' ${subject}.up_down.pool.rename.fasta.fai) > ${subject}.up_down.read.chain
done
cat fasta_collection/*.gene_body.pool.rename.fasta > all.gene_body.pool.rename.fasta
cat fasta_collection/*.up_down.pool.rename.fasta > all.up_down.pool.rename.fasta
```
### Step 3. Remove redunant sequence in the total pools
```
cd /home/wgzhang/probe_design/initial
cd-hit -i all.gene_body.pool.rename.fasta -o all.gene_body.pool.rename.fasta.cdhit95.fa -c 0.95 -n 5 -M 16000 -T 8
cd-hit -i all.up_down.pool.rename.fasta -o all.up_down.pool.rename.fasta.cdhit95.fa -c 0.95 -n 5 -M 16000 -T 8
```
### Step 4. Build reads.cluster.chain
```
cd /home/wgzhang/probe_design/initial
temp_file="all.reads.gene_body.cluster"
echo > "$temp_file"
while IFS= read -r line; do
    if [[ $line == ">Cluster"* ]]; then
printf "%s\n" >> "$temp_file"
    else
        extracted=$(echo "$line" | awk -F'[>.]' '{printf "%s ", $2}')
        printf "%s" "$extracted" >> "$temp_file"
    fi
done < all.gene_body.pool.rename.fasta.cdhit95.fa.clstr

temp_file="all.reads.up_down.cluster"
echo > "$temp_file"
while IFS= read -r line; do
    if [[ $line == ">Cluster"* ]]; then
printf "%s\n" >> "$temp_file"
    else
        extracted=$(echo "$line" | awk -F'[>.]' '{printf "%s ", $2}')
        printf "%s" "$extracted" >> "$temp_file"
    fi
done < all.up_down.pool.rename.fasta.cdhit95.fa.clstr
```
### Step 5. Run blockParse to output all possible probes
```
cd /home/wgzhang/probe_design/initial
echo > all.gene_body_design.fastq
for read in `cut -f1 all.gene_body.pool.rename.fasta.cdhit95.fa.fai` ; do
samtools faidx all.gene_body.pool.rename.fasta.cdhit95.fa ${read} > test.fa
python /home/wgzhang/tools/OligoMiner-master/blockParse.py -l 120 -L 120 -t 40 -T 200 -O \
-f test.fa -o test
cat test.fastq >> all.gene_body_design.fastq
echo "" >> all.gene_body_design.fastq
done

echo > all.up_down_design.fastq
for read in `cut -f1 all.up_down.pool.rename.fasta.cdhit95.fa.fai` ; do
samtools faidx all.up_down.pool.rename.fasta.cdhit95.fa ${read} > test1.fa
python /home/wgzhang/tools/OligoMiner-master/blockParse.py -l 120 -L 120 -t 40 -T 200 -O \
-f test1.fa -o test1
cat test1.fastq >> all.up_down_design.fastq
echo "" >> all.up_down_design.fastq
done
```
### Step 6. Repeatmask detection and mask repeat probes
```
cd /home/wgzhang/probe_design/initial
module load seqkit/0.11.0
seqtk seq -a all.gene_body_design.fastq >all.gene_body_design.fasta
seqtk seq -a all.up_down_design.fastq >all.up_down_design.fasta
seqkit rmdup -n < all.gene_body_design.fasta > all.gene_body_design.rmdup.fasta
seqkit rmdup -n < all.up_down_design.fasta > all.up_down_design.rmdup.fasta
#########################################################
cd /home/wgzhang/probe_design/repeatmask/gene_body
module load RepeatMasker/4.1.0
module load seqkit/0.11.0
fasta_file="all.gene_body_design.rmdup.fasta"
cut -f1 all.gene_body_design.rmdup.fasta.fai > total_reads
batch_size=10000
total_batches=$(($(wc -l < total_reads) / batch_size))
for ((i = 0; i <= total_batches; i++)); do
    start=$((i * batch_size + 1))
    end=$((start + batch_size - 1))
    sed -n "${start},${end}p" total_reads > temp_reads_${i}
seqkit grep -n -f temp_reads_${i} all.gene_body_design.rmdup.fasta > temp_reads_${i}.fa
bsub -J blast -n 4 -R "span[hosts=1] rusage[mem=4GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal \
'
RepeatMasker temp_reads_'"${i}"'.fa -species "bos taurus" -pa 4 -dir gene_body_'"${i}"'
'
done

cd /home/wgzhang/probe_design/repeatmask/up_down
fasta_file="all.up_down_design.rmdup.fasta"
cut -f1 all.up_down_design.rmdup.fasta.fai > total_reads
batch_size=10000
total_batches=$(($(wc -l < total_reads) / batch_size))
for ((i = 0; i <= total_batches; i++)); do
    start=$((i * batch_size + 1))
    end=$((start + batch_size - 1))
    sed -n "${start},${end}p" total_reads > temp_reads_${i}
seqkit grep -n -f temp_reads_${i} all.up_down_design.rmdup.fasta > temp_reads_${i}.fa
bsub -J blast -n 4 -R "span[hosts=1] rusage[mem=4GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal \
'
RepeatMasker temp_reads_'"${i}"'.fa -species "bos taurus" -pa 4 -dir up_down_'"${i}"'
'
done
```
#### Extract non-repeatmask probes
```
cd /home/wgzhang/probe_design/repeatmask/gene_body
for ((i = 0; i <= 217; i++)); do
awk '{print $0}' ./up_down_${i}/temp_reads_${i}.fa.out | grep -v Low_complexity |grep -v Simple_repeat |grep -v position |grep -v begin |grep -v There
done |awk '{print $5}' |grep "a" > repeatmask_find.list
seqkit grep -v -n -f repeatmask_find.list all.gene_body_design.rmdup.fasta > all.gene_body_design.rmdup.repeatmask.fasta

cd /home/wgzhang/probe_design/repeatmask/up_down
for ((i = 0; i <= 1464; i++)); do
awk '{print $0}' ./up_down_${i}/temp_reads_${i}.fa.out | grep -v Low_complexity |grep -v Simple_repeat |grep -v position |grep -v begin |grep -v There
done |awk '{print $5}' |grep "a" > repeatmask_find.list
seqkit grep -v -n -f repeatmask_find.list all.up_down_design.rmdup.fasta > all.up_down_design.rmdup.repeatmask.fasta
```
### Step 7. Blast potential probes in batch
```
cd /home/wgzhang/probe_design/round1
module load BLAST+/2.9.0
for subject in `cut -f1 genome.chain`  ; do
bsub -J blast -n 4 -R "span[hosts=1] rusage[mem=6GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal \
'
blastn -db blast_index/'"${subject}"'_genomic -query all.gene_body_design.rmdup.repeatmask.fasta -out blast_result/'"${subject}"'.gene_body.blast -evalue 1e-30 -outfmt 6 -num_threads 10
'
done

cd /home/wgzhang/probe_design/round1
module load BLAST+/2.9.0
for subject in `cut -f1 genome.chain`  ; do
bsub -J blast -n 4 -R "span[hosts=1] rusage[mem=6GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal \
'
blastn -db blast_index/'"${subject}"'_genomic -query all.up_down_design.rmdup.repeatmask.fasta -out blast_result/'"${subject}"'.up_down.blast -evalue 1e-30 -outfmt 6 -num_threads 10
'
done
```
### Step 8. deciphering blast results
```
cd /home/wgzhang/probe_design/round1
module load BEDTools/2.27
for subject in `cut -f1 genome.chain`  ; do
awk '$4 > 90 && $3 > 90' blast_result/${subject}.gene_body.blast |awk 'BEGIN{OFS="\t"}{print $2,$9,$10,$1}' > blast_off_target/${subject}.gene_body.bed
sh modify_bed_file.sh blast_off_target/${subject}.gene_body.bed blast_off_target/${subject}.gene_body.modify.bed
bedtools sort -i blast_off_target/${subject}.gene_body.modify.bed > blast_off_target/${subject}.gene_body.modify.sort.bed
bedtools intersect -v -a blast_off_target/${subject}.gene_body.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_off_target/${subject}.gene_body.modify.sort.intersect.bed
cut -f4 blast_off_target/${subject}.gene_body.modify.sort.intersect.bed |sort |uniq |awk -F':' '{print $1,$0}' > blast_off_target/${subject}.gene_body.offtarget.probe
cut -f1 -d " " blast_off_target/${subject}.gene_body.offtarget.probe |uniq > blast_off_target/${subject}.gene_body.offtarget.probe.host
for seqhost in `cat blast_off_target/${subject}.gene_body.offtarget.probe.host`  ; do
host=$(grep ${subject} genome.chain|cut -f2)
line=$(grep ${seqhost} all.reads.gene_body.cluster |grep ${host} |wc -l)
if [ "$line" -eq 0 ]; then
    :
else
    echo "$seqhost"
fi
done > blast_off_target/${subject}.gene_body.offtarget.probe.seq 
awk -f vlookup.awk blast_off_target/${subject}.gene_body.offtarget.probe.seq blast_off_target/${subject}.gene_body.offtarget.probe |awk '$3!= "NA" {print $2}' > blast_off_target/${subject}.gene_body.offtarget_real.probe
done

for subject in `cut -f1 genome.chain`  ; do
awk '$4 > 90 && $3 > 90' blast_result/${subject}.up_down.blast |awk 'BEGIN{OFS="\t"}{print $2,$9,$10,$1}' > blast_off_target/${subject}.up_down.bed
sh modify_bed_file.sh blast_off_target/${subject}.up_down.bed blast_off_target/${subject}.up_down.modify.bed
bedtools sort -i blast_off_target/${subject}.up_down.modify.bed > blast_off_target/${subject}.up_down.modify.sort.bed
bedtools intersect -v -a blast_off_target/${subject}.up_down.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_off_target/${subject}.up_down.modify.sort.intersect.bed
cut -f4 blast_off_target/${subject}.up_down.modify.sort.intersect.bed |sort |uniq |awk -F':' '{print $1,$0}' > blast_off_target/${subject}.up_down.offtarget.probe
cut -f1 -d " " blast_off_target/${subject}.up_down.offtarget.probe |uniq > blast_off_target/${subject}.up_down.offtarget.probe.host
for seqhost in `cat blast_off_target/${subject}.up_down.offtarget.probe.host`  ; do
host=$(grep ${subject} genome.chain|cut -f2)
line=$(grep ${seqhost} all.reads.up_down.cluster |grep ${host} |wc -l)
if [ "$line" -eq 0 ]; then
    :
else
    echo "$seqhost"
fi
done > blast_off_target/${subject}.up_down.offtarget.probe.seq 
awk -f vlookup.awk blast_off_target/${subject}.up_down.offtarget.probe.seq blast_off_target/${subject}.up_down.offtarget.probe |awk '$3!= "NA"' > blast_off_target/${subject}.up_down.offtarget_real.probe
done
```
### Step 9. deciphering blast results
```



