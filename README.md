# OR_probe_design_workflow
## This repo includes scripts and details for each step. (scripts are very crude though :( )
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
### Step 9. Design round1
Probe design priority: GC(40-60) > GC(30-70) > GC(20-80) > off-target (1 genome) > off-target (2 genomes) > off-target (3 genomes)

```
cd /home/wgzhang/probe_design/round1/design
awk '{print $2}' gene_body.offtarget_real.probe  | seqkit grep -v -n -f - all.gene_body_design.rmdup.repeatmask.fasta > gene_body.ontarget_probe.fasta
awk '{print $2}' up_down.offtarget_real.probe  | seqkit grep -v -n -f - all.up_down_design.rmdup.repeatmask.fasta > up_down.ontarget_probe.fasta
sed 's/:/_/g' up_down.ontarget_probe.fasta > up_down.ontarget_probe.rename.fasta
sed 's/:/_/g' gene_body.ontarget_probe.fasta > gene_body.ontarget_probe.rename.fasta
cat up_down.ontarget_probe.rename.fasta | infoseq -auto -only -name -length -pgc stdin > up_down.gc
cat gene_body.ontarget_probe.rename.fasta | infoseq -auto -only -name -length -pgc stdin > gene_body.gc
#gene body40-60 gc probes 
awk  '$3>=40 && $3<=60' gene_body.gc > gene_body.4060.gc
awk -F'_' '{print $1"_"$2}' gene_body.4060.gc |uniq > gene_body.4060.select.seqname
for seqname in `cat gene_body.4060.select.seqname`  ; do
grep $seqname gene_body.4060.gc |shuf |head -n1 |cut -d' ' -f1
done > gene_body.4060.select
seqkit grep -n -f gene_body.4060.select gene_body.ontarget_probe.rename.fasta > gene_body.4060.select.fasta 
#gene body30-70 gc probes 
grep -v -f gene_body.4060.select.seqname gene_body.gc  |awk  '$3>=30 && $3<=70' |awk -F'_' '{print $1"_"$2}'  |uniq > gene_body.3070.select.seqname
grep -v -f gene_body.4060.select.seqname gene_body.gc  |awk  '$3>=30 && $3<=70' > gene_body.3070.gc
for seqname in `cat gene_body.3070.select.seqname`  ; do
grep $seqname gene_body.3070.gc |shuf |head -n1 |cut -d' ' -f1
done > gene_body.3070.select
seqkit grep -n -f gene_body.3070.select gene_body.ontarget_probe.rename.fasta > gene_body.3070.select.fasta 
#gene body20-80 gc probes 
cat gene_body.3070.select.seqname gene_body.4060.select.seqname |grep -v -f - gene_body.gc  |awk  '$3>=20 && $3<=80' |awk -F'_' '{print $1"_"$2}'  |uniq > gene_body.2080.select.seqname
cat gene_body.3070.select.seqname gene_body.4060.select.seqname |grep -v -f - gene_body.gc  |awk  '$3>=20 && $3<=80' > gene_body.2080.gc
for seqname in `cat gene_body.2080.select.seqname`  ; do
grep $seqname gene_body.2080.gc |shuf |head -n1 |cut -d' ' -f1
done > gene_body.2080.select
seqkit grep -n -f gene_body.2080.select gene_body.ontarget_probe.rename.fasta > gene_body.2080.select.fasta 
#up down 40-60 gc probes 
awk  '$3>=40 && $3<=60' up_down.gc > up_down.4060.gc
awk -F'_' '{print $1"_"$2}' up_down.4060.gc |uniq > up_down.4060.select.seqname
for seqname in `cat up_down.4060.select.seqname`  ; do
grep $seqname up_down.4060.gc |shuf |head -n1 |cut -d' ' -f1
done > up_down.4060.select
seqkit grep -n -f up_down.4060.select up_down.ontarget_probe.rename.fasta > up_down.4060.select.fasta 
#up down 30-70 gc probes 
grep -v -f up_down.4060.select.seqname up_down.gc  |awk  '$3>=30 && $3<=70' |awk -F'_' '{print $1"_"$2}'  |uniq > up_down.3070.select.seqname
grep -v -f up_down.4060.select.seqname up_down.gc  |awk  '$3>=30 && $3<=70' > up_down.3070.gc
for seqname in `cat up_down.3070.select.seqname`  ; do
grep $seqname up_down.3070.gc |shuf |head -n1 |cut -d' ' -f1
done > up_down.3070.select
seqkit grep -n -f up_down.3070.select up_down.ontarget_probe.rename.fasta > up_down.3070.select.fasta 
#up down 20-80 gc probes 
cat up_down.3070.select.seqname up_down.4060.select.seqname |grep -v -f - up_down.gc  |awk  '$3>=20 && $3<=80' |awk -F'_' '{print $1"_"$2}'  |uniq > up_down.2080.select.seqname
cat up_down.3070.select.seqname up_down.4060.select.seqname |grep -v -f - up_down.gc  |awk  '$3>=20 && $3<=80' > up_down.2080.gc
for seqname in `cat up_down.2080.select.seqname`  ; do
grep $seqname up_down.2080.gc |shuf |head -n1 |cut -d' ' -f1
done > up_down.2080.select
seqkit grep -n -f up_down.2080.select up_down.ontarget_probe.rename.fasta > up_down.2080.select.fasta 
#allow only one genome off target 
cat gene_body.4060.select.seqname gene_body.3070.select.seqname gene_body.2080.select.seqname |grep -v -f - all.gene_body.pool.rename.fasta.cdhit95.fa.fai |cut -f1 > gene_body.off-target.seqname
for seqname in `cat gene_body.off-target.seqname`  ; do
grep $seqname gene_body.offtarget_real.probe.uniq | awk '$1==1' |shuf |head -n1 |awk '{print $2}'
done > gene_body.off-target.g1.select
awk -F':' '{print $1}' gene_body.off-target.g1.select > gene_body.off-target.g1.select.seqname
#allow only two genome off target 
grep -v -f  gene_body.off-target.g1.select.seqname gene_body.off-target.seqname  > gene_body.off-target.seqname2
for seqname in `cat gene_body.off-target.seqname2`  ; do
grep $seqname gene_body.offtarget_real.probe.uniq | awk '$1==2' |shuf |head -n1 |awk '{print $2}'
done > gene_body.off-target.g2.select
awk -F':' '{print $1}' gene_body.off-target.g2.select > gene_body.off-target.g2.select.seqname
#allow three genome off target 
grep -v -f  gene_body.off-target.g2.select.seqname gene_body.off-target.seqname2  > gene_body.off-target.seqname3
for seqname in `cat gene_body.off-target.seqname3`  ; do
grep $seqname gene_body.offtarget_real.probe.uniq | awk '$1==3' |shuf |head -n1 |awk '{print $2}'
done > gene_body.off-target.g3.select
awk -F':' '{print $1}' gene_body.off-target.g3.select > gene_body.off-target.g3.select.seqname
cat gene_body.off-target.g1.select gene_body.off-target.g2.select gene_body.off-target.g3.select |seqkit grep -n -f - all.gene_body_design.rmdup.repeatmask.fasta >gene_body.off-target.select.fasta
#allow only one genome off target 
cat up_down.4060.select.seqname up_down.3070.select.seqname up_down.2080.select.seqname |grep -v -f - all.up_down.pool.rename.fasta.cdhit95.fa.fai |cut -f1 > up_down.off-target.seqname
for seqname in `cat up_down.off-target.seqname`  ; do
grep $seqname up_down.offtarget_real.probe.uniq | awk '$1==1' |shuf |head -n1 |awk '{print $2}'
done > up_down.off-target.g1.select
awk -F':' '{print $1}' up_down.off-target.g1.select > up_down.off-target.g1.select.seqname
#allow only two genome off target 
grep -v -f  up_down.off-target.g1.select.seqname up_down.off-target.seqname  > up_down.off-target.seqname2
for seqname in `cat up_down.off-target.seqname2`  ; do
grep $seqname up_down.offtarget_real.probe.uniq | awk '$1==2' |shuf |head -n1 |awk '{print $2}'
done > up_down.off-target.g2.select
awk -F':' '{print $1}' up_down.off-target.g2.select > up_down.off-target.g2.select.seqname
#allow three genome off target 
grep -v -f  up_down.off-target.g2.select.seqname up_down.off-target.seqname2  > up_down.off-target.seqname3
for seqname in `cat up_down.off-target.seqname3`  ; do
grep $seqname up_down.offtarget_real.probe.uniq | awk '$1==3' |shuf |head -n1 |awk '{print $2}'
done > up_down.off-target.g3.select
awk -F':' '{print $1}' up_down.off-target.g3.select > up_down.off-target.g3.select.seqname
cat up_down.off-target.g1.select up_down.off-target.g2.select up_down.off-target.g3.select |seqkit grep -n -f - all.up_down_design.rmdup.repeatmask.fasta > up_down.off-target.select.fasta
```
#### Merge all designed probes
```
cat gene_body.2080.select.fasta gene_body.4060.select.fasta gene_body.3070.select.fasta gene_body.off-target.select.fasta > ../first_round_design/gene_body.design.round1.fasta
cat up_down.2080.select.fasta up_down.4060.select.fasta up_down.3070.select.fasta up_down.off-target.select.fasta > ../first_round_design/up_down.design.round1.fasta
```

### Step 10. Blast first round probes
```
cd /home/wgzhang/probe_design/round1/
module load BLAST+/2.9.0
for subject in `cut -f1 genome.chain`  ; do
bsub -J blast -n 4 -R "span[hosts=1] rusage[mem=6GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal \
'
blastn -db blast_index/'"${subject}"'_genomic -query design/gene_body.design.round1.fasta -out design_validate/'"${subject}"'.gene_body.blast -evalue 1e-30 -outfmt 6 -num_threads 10
'
done

module load BLAST+/2.9.0
for subject in `cut -f1 genome.chain`  ; do
bsub -J blast -n 4 -R "span[hosts=1] rusage[mem=6GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal \
'
blastn -db blast_index/'"${subject}"'_genomic -query design/up_down.design.round1.fasta -out design_validate/'"${subject}"'.up_down.blast -evalue 1e-30 -outfmt 6 -num_threads 10
'
done
```
### Step 11. Summary first round probes
```
cd /home/wgzhang/probe_design/round1
for subject in `cut -f1 genome.chain`  ; do
awk '$4 > 90 && $3 > 90' design_validate/${subject}.gene_body.blast |awk 'BEGIN{OFS="\t"}{print $2,$9,$10,$1}' > blast_result_round1_validate/${subject}.gene_body.bed
sh modify_bed_file.sh blast_result_round1_validate/${subject}.gene_body.bed blast_result_round1_validate/${subject}.gene_body.modify.bed
bedtools sort -i blast_result_round1_validate/${subject}.gene_body.modify.bed > blast_result_round1_validate/${subject}.gene_body.modify.sort.bed
bedtools intersect -wa -wb -a blast_result_round1_validate/${subject}.gene_body.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_result_round1_validate/${subject}.gene_body.on-target.bed
awk -F '\t' '{print $5"_"$6"_"$7}' blast_result_round1_validate/${subject}.gene_body.on-target.bed |sort |uniq -c > blast_result_round1_validate/${subject}.gene_body.on-target.bed.stat
bedtools intersect -v -a blast_result_round1_validate/${subject}.gene_body.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_result_round1_validate/${subject}.gene_body.off-target.bed
cut -f4 blast_result_round1_validate/${subject}.gene_body.off-target.bed |awk -F '[_:]' '{print $1"_"$2}' > blast_result_round1_validate/${subject}.gene_body.off-target.seqnames
host=$(grep ${subject} genome.chain|cut -f2)
probe_hits=$(cat blast_result_round1_validate/${subject}.gene_body.bed |wc -l)
offtargetnarrow=$(grep -f blast_result_round1_validate/${subject}.gene_body.off-target.seqnames all.reads.gene_body.cluster |grep ${host} |wc -l)
offtargetgeneral=$(cat blast_result_round1_validate/${subject}.gene_body.off-target.bed |wc -l)
total_OR=$(cat target_pool/${subject}.target_regions.modify.bed |wc -l)
target1_OR=$(awk '$1 == 1' blast_result_round1_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target2_OR=$(awk '$1 == 2' blast_result_round1_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target3_OR=$(awk '$1 == 3' blast_result_round1_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target4_OR=$(awk '$1 == 4' blast_result_round1_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target5_OR=$(awk '$1 >= 5' blast_result_round1_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target_OR=$(expr $target1_OR + $target2_OR + $target3_OR + $target4_OR + $target5_OR)
ungrep=$(expr $total_OR - $target_OR)
echo -e "$subject\t$probe_hits\t$total_OR\t$target1_OR\t$target2_OR\t$target3_OR\t$target4_OR\t$target5_OR\t$ungrep\t$target_OR\t$offtargetgeneral\t$offtargetnarrow"
done > gene_body.round1.design.summary

for subject in `cut -f1 genome.chain`  ; do
awk '$4 > 90 && $3 > 90' design_validate/${subject}.up_down.blast |awk 'BEGIN{OFS="\t"}{print $2,$9,$10,$1}' > blast_result_round1_validate/${subject}.up_down.bed
sh modify_bed_file.sh blast_result_round1_validate/${subject}.up_down.bed blast_result_round1_validate/${subject}.up_down.modify.bed
bedtools sort -i blast_result_round1_validate/${subject}.up_down.modify.bed > blast_result_round1_validate/${subject}.up_down.modify.sort.bed
bedtools intersect -wa -wb -a blast_result_round1_validate/${subject}.up_down.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_result_round1_validate/${subject}.up_down.on-target.bed
awk -F '\t' '{print $5"_"$6"_"$7}' blast_result_round1_validate/${subject}.up_down.on-target.bed |sort |uniq -c > blast_result_round1_validate/${subject}.up_down.on-target.bed.stat
bedtools intersect -v -a blast_result_round1_validate/${subject}.up_down.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_result_round1_validate/${subject}.up_down.off-target.bed
cut -f4 blast_result_round1_validate/${subject}.up_down.off-target.bed |awk -F '[_:]' '{print $1"_"$2}' > blast_result_round1_validate/${subject}.up_down.off-target.seqnames
host=$(grep ${subject} genome.chain|cut -f2)
probe_hits=$(cat blast_result_round1_validate/${subject}.up_down.bed |wc -l)
offtargetnarrow=$(grep -f blast_result_round1_validate/${subject}.up_down.off-target.seqnames all.reads.up_down.cluster |grep ${host} |wc -l)
offtargetgeneral=$(cat blast_result_round1_validate/${subject}.up_down.off-target.bed |wc -l)
total_OR=$(cat target_pool/${subject}.target_regions.modify.bed |wc -l)
target1_OR=$(awk '$1 == 1' blast_result_round1_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target2_OR=$(awk '$1 == 2' blast_result_round1_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target3_OR=$(awk '$1 == 3' blast_result_round1_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target4_OR=$(awk '$1 == 4' blast_result_round1_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target5_OR=$(awk '$1 >= 5' blast_result_round1_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target_OR=$(expr $target1_OR + $target2_OR + $target3_OR + $target4_OR + $target5_OR)
ungrep=$(expr $total_OR - $target_OR)
echo -e "$subject\t$probe_hits\t$total_OR\t$target1_OR\t$target2_OR\t$target3_OR\t$target4_OR\t$target5_OR\t$ungrep\t$target_OR\t$offtargetgeneral\t$offtargetnarrow"
done > up_down.round1.design.summary
```
### Step 12. Round2 probe design 
#### Calulate round1 probe enrichment
```
cd /home/wgzhang/probe_design/round2/
for subject in `cut -f1 genome.chain`  ; do
cat on-target/${subject}.* |cut -f4 |sort |uniq -c |awk '{print $2,1/$1}' > summary_tmp/${subject}.probe.unique
cat on-target/${subject}.* |awk '{print $1,$2,$3,$4,$5"-"$6"-"$7}' > summary_tmp/${subject}.on-target.tab
file1="${subject}.probe.unique"
file2="${subject}.on-target.tab"
while IFS= read -r line; do
    query=$(echo "$line" | awk '{print $4}')
    match=$(grep -F "$query" "$file1" | awk '{print $2}')
    if [ -n "$match" ]; then
        echo "$line $match"
    else
        echo "$line"
    fi
done < "$file2" > summary_tmp/${subject}.on-target.tab.counts
for OR in `cut -f5 -d ' ' summary_tmp/${subject}.on-target.tab.counts |sort |uniq`  ; do
uniq=$(grep ${OR} summary_tmp/${subject}.on-target.tab.counts |awk -F' ' '$6==1' |wc -l)
score=$(grep ${OR} summary_tmp/${subject}.on-target.tab.counts |awk '{sum += $6} END {print a,sum}')
up=$(grep ${OR} summary_tmp/${subject}.on-target.tab.counts |awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $1,$2,$4+2500,$5-2500}' |awk '$1<$3'|wc -l)
down=$(grep ${OR} summary_tmp/${subject}.on-target.tab.counts |awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $1,$2,$4+2500,$5-2500}' |awk '$4<$2'|wc -l)
echo -e "$OR\t$score\t$uniq\t$up\t$down"
done > summary_tmp/${subject}.on-target.tab.counts.enrichment
for OR in `cut -f5 -d ' ' summary_tmp/${subject}.on-target.tab.counts |sort |uniq`  ; do
score=$(grep ${OR} summary_tmp/${subject}.on-target.tab.counts |awk '{sum += $6} END {print a,sum}')
grep ${OR} summary_tmp/${subject}.on-target.tab.counts > ana.hits
while IFS= read -r line; do
deph1=$(echo "$line" | awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $1}')
deph2=$(echo "$line" | awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $2}')
deph3=$(echo "$line" | awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $4+2500}')
deph4=$(echo "$line" | awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $5-2500}')
if [ "$deph3" -gt "$deph1" ]; then
echo -e "$line\tup"
elif [ "$deph2" -gt "$deph4" ]; then
echo -e "$line\tdown"
else
echo -e "$line\tcoding"
fi
done < ana.hits |awk  -F' ' -v a="${score}" -v b="${OR}" '$6==1 {print $4,b,a,$7}' 
done > summary_tmp/${subject}.on-target.tab.counts.probe_unique
done
```

#### Remove unique probe corresponding to high score OR
```
library(UpSetR)
library(reshape2)
library(ggplot2)
working_directory <- getwd()
files <- list.files(working_directory)
target_files <- files[grep(".on-target.tab.counts.probe_unique$", files)]
combined_data=c()
for(i in target_files){
df=read.table(i,head=F)
df$name=substr(i,1,15)
combined_data=rbind(combined_data,df)
}
combined_data=combined_data[which(combined_data$V3>=8),]
result <- dcast(combined_data,  V1 ~ name, fun.aggregate = length)
rownames(result)=result[,1]
result=result[,-1]
jpeg(file="gene_body.upset.jpeg",width=2000,height=2000)
upset(result,sets=colnames(result),order.by="freq")
dev.off()
rowsums =rowSums(result)
remove=names(rowsums)[which(rowsums>42)]
write.table(remove, file = "removelist_score8_common42.txt", row.names = FALSE, col.names = FALSE)
```
```
sed -i "s/\"//g" removelist_score8_common42.txt
```
```
cd /home/wgzhang/probe_design/round2
for OR in `cut -d' ' -f2 summary_tmp/GCA_905123885.1_ROSLIN_BTI_ANK1.on-target.tab.counts.probe_unique |sort |uniq`; do
grep ${OR} summary_tmp/GCA_905123885.1_ROSLIN_BTI_ANK1.on-target.tab.counts.probe_unique |awk -F ' ' '$4=="up"' |shuf |sed '1d' > summary_tmp/potential_up
grep ${OR} summary_tmp/GCA_905123885.1_ROSLIN_BTI_ANK1.on-target.tab.counts.probe_unique |awk -F ' ' '$4=="down"' |shuf |sed '1d' > summary_tmp/potential_down
cat  summary_tmp/potential_up summary_tmp/potential_down |awk -F' ' '{print $1}' | grep -f removelist_score8_common42.txt -
done > removelist_score8_common42.filter.txt
cp gene_body.design.round1.fasta gene_body.design.round2.fasta
seqkit grep -v -n -f removelist_score8_common42.filter.txt up_down.design.round1.fasta > up_down.design.round2.fasta
```
### Step 13. Blast second round probes
```
cd /home/wgzhang/probe_design/round2
module load BLAST+/2.9.0
for subject in `cut -f1 genome.chain`  ; do
bsub -J blast -n 4 -R "span[hosts=1] rusage[mem=6GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal \
'
blastn -db blast_index/'"${subject}"'_genomic -query gene_body.design.round2.fasta -out blast_result_round2/'"${subject}"'.gene_body.blast -evalue 1e-30 -outfmt 6 -num_threads 10
'
bsub -J blast -n 4 -R "span[hosts=1] rusage[mem=6GB] select[maxmem>224800]" -o %J.out -e %J.err -q normal \
'
blastn -db blast_index/'"${subject}"'_genomic -query up_down.design.round2.fasta -out blast_result_round2/'"${subject}"'.up_down.blast -evalue 1e-30 -outfmt 6 -num_threads 10
'
done
```
### Step 14. Summary second round probes
```
cd /home/wgzhang/probe_design/round2
for subject in `cut -f1 genome.chain`  ; do
awk '$4 > 90 && $3 > 90' blast_result/${subject}.gene_body.blast |awk 'BEGIN{OFS="\t"}{print $2,$9,$10,$1}' > blast_result_round2_validate/${subject}.gene_body.bed
sh modify_bed_file.sh blast_result_round2_validate/${subject}.gene_body.bed blast_result_round2_validate/${subject}.gene_body.modify.bed
bedtools sort -i blast_result_round2_validate/${subject}.gene_body.modify.bed > blast_result_round2_validate/${subject}.gene_body.modify.sort.bed
bedtools intersect -wa -wb -a blast_result_round2_validate/${subject}.gene_body.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_result_round2_validate/${subject}.gene_body.on-target.bed
awk -F '\t' '{print $5"_"$6"_"$7}' blast_result_round2_validate/${subject}.gene_body.on-target.bed |sort |uniq -c > blast_result_round2_validate/${subject}.gene_body.on-target.bed.stat
bedtools intersect -v -a blast_result_round2_validate/${subject}.gene_body.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_result_round2_validate/${subject}.gene_body.off-target.bed
cut -f4 blast_result_round2_validate/${subject}.gene_body.off-target.bed |awk -F '[_:]' '{print $1"_"$2}' > blast_result_round2_validate/${subject}.gene_body.off-target.seqnames
host=$(grep ${subject} genome.chain|cut -f2)
probe_hits=$(cat blast_result_round2_validate/${subject}.gene_body.bed |wc -l)
offtargetnarrow=$(grep -f blast_result_round2_validate/${subject}.gene_body.off-target.seqnames all.reads.gene_body.cluster |grep ${host} |wc -l)
offtargetgeneral=$(cat blast_result_round2_validate/${subject}.gene_body.off-target.bed |wc -l)
total_OR=$(cat target_pool/${subject}.target_regions.modify.bed |wc -l)
target1_OR=$(awk '$1 == 1' blast_result_round2_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target2_OR=$(awk '$1 == 2' blast_result_round2_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target3_OR=$(awk '$1 == 3' blast_result_round2_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target4_OR=$(awk '$1 == 4' blast_result_round2_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target5_OR=$(awk '$1 >= 5' blast_result_round2_validate/${subject}.gene_body.on-target.bed.stat |wc -l)
target_OR=$(expr $target1_OR + $target2_OR + $target3_OR + $target4_OR + $target5_OR)
ungrep=$(expr $total_OR - $target_OR)
echo -e "$subject\t$probe_hits\t$total_OR\t$target1_OR\t$target2_OR\t$target3_OR\t$target4_OR\t$target5_OR\t$ungrep\t$target_OR\t$offtargetgeneral\t$offtargetnarrow"
done > gene_body.round2.design.summary

for subject in `cut -f1 genome.chain`  ; do
awk '$4 > 90 && $3 > 90' blast_result/${subject}.up_down.blast |awk 'BEGIN{OFS="\t"}{print $2,$9,$10,$1}' > blast_result_round2_validate/${subject}.up_down.bed
sh modify_bed_file.sh blast_result_round2_validate/${subject}.up_down.bed blast_result_round2_validate/${subject}.up_down.modify.bed
bedtools sort -i blast_result_round2_validate/${subject}.up_down.modify.bed > blast_result_round2_validate/${subject}.up_down.modify.sort.bed
bedtools intersect -wa -wb -a blast_result_round2_validate/${subject}.up_down.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_result_round2_validate/${subject}.up_down.on-target.bed
awk -F '\t' '{print $5"_"$6"_"$7}' blast_result_round2_validate/${subject}.up_down.on-target.bed |sort |uniq -c > blast_result_round2_validate/${subject}.up_down.on-target.bed.stat
bedtools intersect -v -a blast_result_round2_validate/${subject}.up_down.modify.sort.bed -b target_pool/${subject}.target_regions.modify.bed > blast_result_round2_validate/${subject}.up_down.off-target.bed
cut -f4 blast_result_round2_validate/${subject}.up_down.off-target.bed |awk -F '[_:]' '{print $1"_"$2}' > blast_result_round2_validate/${subject}.up_down.off-target.seqnames
host=$(grep ${subject} genome.chain|cut -f2)
probe_hits=$(cat blast_result_round2_validate/${subject}.up_down.bed |wc -l)
offtargetnarrow=$(grep -f blast_result_round2_validate/${subject}.up_down.off-target.seqnames all.reads.up_down.cluster |grep ${host} |wc -l)
offtargetgeneral=$(cat blast_result_round2_validate/${subject}.up_down.off-target.bed |wc -l)
total_OR=$(cat target_pool/${subject}.target_regions.modify.bed |wc -l)
target1_OR=$(awk '$1 == 1' blast_result_round2_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target2_OR=$(awk '$1 == 2' blast_result_round2_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target3_OR=$(awk '$1 == 3' blast_result_round2_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target4_OR=$(awk '$1 == 4' blast_result_round2_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target5_OR=$(awk '$1 >= 5' blast_result_round2_validate/${subject}.up_down.on-target.bed.stat |wc -l)
target_OR=$(expr $target1_OR + $target2_OR + $target3_OR + $target4_OR + $target5_OR)
ungrep=$(expr $total_OR - $target_OR)
echo -e "$subject\t$probe_hits\t$total_OR\t$target1_OR\t$target2_OR\t$target3_OR\t$target4_OR\t$target5_OR\t$ungrep\t$target_OR\t$offtargetgeneral\t$offtargetnarrow"
done > up_down.round2.design.summary
```

### Step 15. calculate round2 probe enrichment
```
for subject in `cut -f1 genome.chain `  ; do
cat on-target-round2/${subject}.* |cut -f4 |sort |uniq -c |awk '{print $2,1/$1}' > summary_round2_tmp/${subject}.probe.unique
cat on-target-round2/${subject}.* |awk '{print $1,$2,$3,$4,$5"-"$6"-"$7}' > summary_round2_tmp/${subject}.on-target.tab
file1="summary_round2_tmp/${subject}.probe.unique"
file2="summary_round2_tmp/${subject}.on-target.tab"
while IFS= read -r line; do
    query=$(echo "$line" | awk '{print $4}')
    match=$(grep -F "$query" "$file1" | awk '{print $2}')
    if [ -n "$match" ]; then
        echo "$line $match"
    else
        echo "$line"
    fi
done < "$file2" > summary_round2_tmp/${subject}.on-target.tab.counts
for OR in `cut -f5 -d ' ' summary_round2_tmp/${subject}.on-target.tab.counts |sort |uniq`  ; do
uniq=$(grep ${OR} summary_round2_tmp/${subject}.on-target.tab.counts |awk -F' ' '$6==1' |wc -l)
score=$(grep ${OR} summary_round2_tmp/${subject}.on-target.tab.counts |awk '{sum += $6} END {print a,sum}')
up=$(grep ${OR} summary_round2_tmp/${subject}.on-target.tab.counts |awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $1,$2,$4+2500,$5-2500}' |awk '$1<$3'|wc -l)
down=$(grep ${OR} summary_round2_tmp/${subject}.on-target.tab.counts |awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $1,$2,$4+2500,$5-2500}' |awk '$4<$2'|wc -l)
coding=$(($(grep "${OR}" "summary_round2_tmp/${subject}.on-target.tab.counts" | wc -l) - $up - $down))
echo -e "$OR\t$score\t$uniq\t$up\t$down\t$coding"
done > summary_round2_tmp/${subject}.on-target.tab.counts.enrichment
for OR in `cut -f5 -d ' ' summary_round2_tmp/${subject}.on-target.tab.counts |sort |uniq`  ; do
score=$(grep ${OR} summary_round2_tmp/${subject}.on-target.tab.counts |awk '{sum += $6} END {print a,sum}')
grep ${OR} summary_round2_tmp/${subject}.on-target.tab.counts > summary_round2_tmp/ana.hits
while IFS= read -r line; do
deph1=$(echo "$line" | awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $1}')
deph2=$(echo "$line" | awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $2}')
deph3=$(echo "$line" | awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $4+2500}')
deph4=$(echo "$line" | awk '{print $2,$3,$5}' |awk -F'[- ]' '{print $5-2500}')
if [ "$deph3" -gt "$deph1" ]; then
echo -e "$line\tup"
elif [ "$deph2" -gt "$deph4" ]; then
echo -e "$line\tdown"
else
echo -e "$line\tcoding"
fi
done < summary_round2_tmp/ana.hits |awk  -F' ' -v a="${score}" -v b="${OR}" '$6==1 {print $4,b,a,$7}' 
done > summary_round2_tmp/${subject}.on-target.tab.counts.probe_unique
done

for subject in `cut -f1 summary_round2_tmp/genome.chain`  ; do
target=$(cat summary_round2_tmp/${subject}.on-target.tab.counts.enrichment|wc -l)
total=$(cat target_regions/${subject}.target_regions.modify.bed|wc -l)
up=$(awk '$4==0' summary_round2_tmp/${subject}.on-target.tab.counts.enrichment |wc -l)
down=$(awk '$5==0' summary_round2_tmp/${subject}.on-target.tab.counts.enrichment |wc -l)
coding=$(awk '$6==0' summary_round2_tmp/${subject}.on-target.tab.counts.enrichment |wc -l)
echo -e "$subject\t$total\t$target\t$up\t$down\t$coding"
done 
```

