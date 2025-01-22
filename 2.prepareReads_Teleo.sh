################################################### Prepare reads for downstream analysis #######################################################

#it needs a three column file: sample name - path to R1 -  path to R2
#raw read files are in .gz
names=(`cut -f 1 /share/projects/OCEAN_eDNA/AZORES/GIT/Teleo_samples_20230703.txt`)
R1=(`cut -f 2 /share/projects/OCEAN_eDNA/AZORES/GIT/Teleo_samples_20230703.txt`)
R2=(`cut -f 3 /share/projects/OCEAN_eDNA/AZORES/GIT/Teleo_samples_20230703.txt`)
x=(`wc -l /share/projects/OCEAN_eDNA/AZORES/GIT/Teleo_samples_20230703.txt`)
y=$((x-1))
z=$((x-2))

echo -n "" > fasta_file.txt
echo -n "" > group.txt
echo -n "" > Teleo_samples_20230703.fasta
echo "Sample,raw,assembled,retained,cut_fwd,cut_rev,cut" > reads_stats.txt
echo -n "" > Samples_min_1read_fwd.txt
echo -n "" > Samples_min_1read_rev.txt

for i in `seq 0 $y`; do

	echo -n ${names[i]} >> reads_stats.txt
	echo -n ","   >> reads_stats.txt
	cp ${R1[i]} .
	cp ${R2[i]} .
	gzip -dc < ${R1[i]} > ${names[i]}_R1.fastq
	gzip -dc < ${R2[i]} > ${names[i]}_R2.fastq
	r="$(grep -c '^@A0' ${names[i]}_R1_temp.fastq)"
	echo -n $r >> reads_stats.txt
	echo "" > ${names[i]}.log

#merge pairs

	pear -f ${names[i]}_R1_temp.fastq -r ${names[i]}_R2_temp.fastq -v 30 -m 130 -n 80 -o ${names[i]} >> ${names[i]}.log
	echo -n "," >> reads_stats.txt
	r="$(grep -c '^@A0' ${names[i]}.assembled.fastq)"
	echo -n $r >> reads_stats.txt

	#remove low quality reads

	trimmomatic SE -phred33  ${names[i]}.assembled.fastq ${names[i]}.assembled_retained.fastq AVGQUAL:25 >> ${names[i]}.log
	echo -n ","   >> reads_stats.txt
	r="$(grep -c '^@A0' ${names[i]}.assembled_retained.fastq)"
	echo -n $r >> reads_stats.txt

#Remove forward and reverse primers from assembled sequences (forward and reverse sequences separately)

	cutadapt -g ACACCGCCCGTCACTCT...CATGGTAAGTGTACCGGAAG --discard-untrimmed --minimum-length 50 -e 0.2 -O 15 -o ${names[i]}_cut_fwd.fastq ${names[i]}.assembled_retained.fastq >> ${names[i]}.log
	cutadapt -g CTTCCGGTACACTTACCATG...AGAGTGACGGGCGGTGT --discard-untrimmed --minimum-length 50 -e 0.2 -O 15 -o ${names[i]}_cut_rev.fastq ${names[i]}.assembled_retained.fastq >> ${names[i]}.log
	echo -n ","   >> reads_stats.txt
	t="$(grep -c '^@A0' ${names[i]}_cut_fwd.fastq)"
	echo -n $t >> reads_stats.txt
	echo -n ","   >> reads_stats.txt
	q="$(grep -c '^@A0' ${names[i]}_cut_rev.fastq)"
  echo -n $q >> reads_stats.txt
	echo -n "," >> reads_stats.txt
  w=$(($t+$q))
  echo $w >> reads_stats.txt

	g="$(grep -c '^@A0' ${names[i]}_cut_fwd.fastq)"
	  if (($g>0)); then
	    echo -n ${names[i]} >> Samples_min_1read_fwd.txt
	  	echo -e "" >> Samples_min_1read_fwd.txt
	  fi
	 g="$(grep -c '^@A0' ${names[i]}_cut_rev.fastq)"
	  if (($g>0)); then
	    echo -n ${names[i]} >> Samples_min_1read_rev.txt
	    echo -e "" >> Samples_min_1read_rev.txt
	  fi

done

#Prepare .fasta and .groups files for fwd and rev separately + concatenate fwd and rev files
	echo -n "" > fasta_file_fwd.txt
	echo -n "" > group_fwd.txt
	echo -n "" > Samples_fwd.fasta
	echo -n "" > fasta_file_rev.txt
	echo -n "" > group_rev.txt
	echo -n "" > Samples_rev.fasta

# 3a) Prepare .fasta and .groups file for FORWARD
names2=(`cut -f 1 Samples_min_1read_fwd.txt`)
j=(`wc -l Samples_min_1read_fwd.txt`)
k=$((j-1))
m=$((j-2))

	for i in `seq 0 $k`; do
	  #convert into fasta
		mothur "#fastq.info(fastq=${names2[i]}_cut_fwd.fastq)"
		if (($i > $m)) ; then
			echo ${names2[i]}_cut_fwd.fasta >> fasta_file_fwd.txt
			echo ${names2[i]} >> group_fwd.txt
		else
			echo -n ${names2[i]}_cut_fwd.fasta"-" >> fasta_file_fwd.txt
			echo -n ${names2[i]}"-" >> group_fwd.txt
		fi
		cat ${names2[i]}_cut_fwd.fasta >> Samples_fwd.fasta
	done

file="fasta_file_fwd.txt"
fastas_fwd=$(cat "$file")
file="group_fwd.txt"
groups_fwd=$(cat "$file")
mothur "#make.group(fasta=$fastas_fwd, groups=$groups_fwd, output=Samples_fwd.groups)"
rm -f mothur*logfile

# 3b) Prepare .fasta and .groups file for REVERSE
names2=(`cut -f 1 Samples_min_1read_rev.txt`)
j=(`wc -l Samples_min_1read_rev.txt`)
k=$((j-1))
m=$((j-2))

for i in `seq 0 $k`; do
  #convert into fasta
  mothur "#fastq.info(fastq=${names2[i]}_cut_rev.fastq)"
  if (($i > $m)) ; then
    echo ${names2[i]}_cut_rev.fasta >> fasta_file_rev.txt
    echo ${names2[i]} >> group_rev.txt
  else
    echo -n ${names2[i]}_cut_rev.fasta"-" >> fasta_file_rev.txt
    echo -n ${names2[i]}"-" >> group_rev.txt
  fi
  cat ${names2[i]}_cut_rev.fasta >> Samples_rev.fasta
done

file="fasta_file_rev.txt"
fastas_rev=$(cat "$file")
file="group_rev.txt"
groups_rev=$(cat "$file")
mothur "#make.group(fasta=$fastas_rev, groups=$groups_rev, output=Samples_rev.groups)"
# Reverse-complement of _rev sequences in the *_rev.fasta and rename the resulting file
mothur "#reverse.seqs(fasta=Samples_rev.fasta)"
rm -f mothur*logfile

# 3c) Concatenate fasta and groups files
cat Samples_fwd.fasta Samples_rev.rc.fasta > Samples_joined_temp.fasta
cat Samples_fwd.groups Samples_rev.groups > Samples_joined_temp.groups

# Remove duplicates
awk '/^>/ { f = !a[$0]++ } f' Samples_joined_temp.fasta > Samples_joined.fasta
awk '!seen[$0]++' Samples_joined_temp.groups > Samples_joined.groups
