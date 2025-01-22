### PIPELINE FOR PROCESSING RAW DATA of AZORES samples -Tagsteady protocol- ###
### COI, TELEO, and 18S-ceph primers ###

## STEP 1: DEMULTIPLEX MIX's
mix=(`cut -f 1 /share/projects/OCEAN_eDNA/AZORES/GIT/Mix_path_to_R1_R2_files_run2.txt`)
R1=(`cut -f 2 /share/projects/OCEAN_eDNA/AZORES/GIT/Mix_path_to_R1_R2_files_run2.txt`)
R2=(`cut -f 3 /share/projects/OCEAN_eDNA/AZORES/GIT/Mix_path_to_R1_R2_files_run2.txt`)
x=(`wc -l /share/projects/OCEAN_eDNA/AZORES/GIT/Mix_path_to_R1_R2_files_run2.txt`)
y=$((x-1))

echo "Mix,MixRaw,Sample,RawReads" > reads_stats_samples.txt

for i in `seq 0 $y`; do
  gzip -dc < ${R1[i]} > ${mix[i]}_R1.fastq
  gzip -dc < ${R2[i]} > ${mix[i]}_R2.fastq
  r="$(grep -c '^@A0' ${mix[i]}_R1.fastq)"
	echo "" > ${mix[i]}.log

  sample=(`cut -f 1 /share/projects/OCEAN_eDNA/AZORES/GIT/${mix[i]}_barcodes_run2.txt`)
  Tags=(`cut -f 2 /share/projects/OCEAN_eDNA/AZORES/GIT/${mix[i]}_barcodes_run2.txt`)
  a=(`wc -l /share/projects/OCEAN_eDNA/AZORES/GIT/${mix[i]}_barcodes_run2.txt`)
  b=$((a-1))

  for j in `seq 0 $b`; do
    cutadapt -g ${Tags[j]} -G ${Tags[j]} --untrimmed-output ${sample[j]}_${mix[i]}_R1_untagged.fastq --untrimmed-paired-output ${sample[j]}_${mix[i]}_R2_untagged.fastq -e 0.1 -O 6 -o ${sample[j]}_${mix[i]}_R1_temp.fastq -p ${sample[j]}_${mix[i]}_R2_temp.fastq ${mix[i]}_R1.fastq ${mix[i]}_R2.fastq >> ${mix[i]}.log
    echo -n ${mix[i]} >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    echo -n $r >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    echo -n ${sample[j]} >> reads_stats_samples.txt
    echo -n "," >> reads_stats_samples.txt
    s="$(grep -c '^@A0' ${sample[j]}_${mix[i]}_R1_temp.fastq)"
    echo -n $s >> reads_stats_samples.txt
    echo -e "" >> reads_stats_samples.txt
  done
done
