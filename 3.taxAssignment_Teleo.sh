
mothur

#create logfile with all manip that are going to be done
set.logfile(name=Teleo.log)

# Analyse reads (length, ambiguities) and count number of reads per group
summary.seqs(fasta=Samples_joined.fasta, processors=12)
count.groups(group=Samples_joined.groups)

# Dereplicate
unique.seqs(fasta=Samples_joined.fasta)
summary.seqs(fasta=Samples_joined.unique.fasta, name=Samples_joined.names, processors=12)

# Align sequences against the 12S rRNA database
alignmentlocation=/share/projects/OCEAN_eDNA/DATABASES/12S_LOCAL
alignmentprefix=teleo
align.seqs(fasta=Samples_joined.unique.fasta, reference=$alignmentlocation/$alignmentprefix.align, processors=12, flip=T)
summary.seqs(fasta=Samples_joined.unique.align, name=Samples_joined.names, processors=12)

#remove sequences not covering the "teleo" region of 12S and with ambiguous bases
screen.seqs(fasta=Samples_joined.unique.align, name=Samples_joined.names, group=Samples_joined.groups, minlength=40, maxlength=80, maxambig=0, processors=12) It took 1 secs to screen 173428 sequences, removed 81034
summary.seqs(fasta=Samples_joined.unique.good.align, name=Samples_joined.good.names, processors=12)
count.groups(group=Samples_joined.good.groups)

# remove columns that contain gap characters
filter.seqs(fasta=Samples_joined.unique.good.align, vertical=T, processors=12)
unique.seqs(fasta=Samples_joined.unique.good.filter.fasta, name=Samples_joined.good.names)
summary.seqs(fasta=Samples_joined.unique.good.filter.unique.fasta,name=Samples_joined.unique.good.filter.names,processors=12)

# Remove chimeras
chimera.uchime(fasta=Samples_joined.unique.good.filter.unique.fasta, name=Samples_joined.unique.good.filter.names, group=Samples_joined.good.groups, processors=12)
remove.seqs(accnos=Samples_joined.unique.good.filter.unique.denovo.uchime.accnos, fasta=Samples_joined.unique.good.filter.unique.fasta, name=Samples_joined.unique.good.filter.names, group=Samples_joined.good.groups, dups=T)
summary.seqs(fasta=Samples_joined.unique.good.filter.unique.pick.fasta, name=Samples_joined.unique.good.filter.pick.names, processors=12)
count.groups(group=Samples_joined.good.pick.groups)

# clarity, rename files
system(cp Samples_joined.unique.good.filter.unique.pick.fasta Samples_joined_all.fasta)
system(cp Samples_joined.unique.good.filter.pick.names Samples_joined_all.names)
system(cp Samples_joined.good.pick.groups Samples_joined_all.groups)

#create count_table with the new files
make.table(name=Samples_joined_all.names, group=Samples_joined_all.groups, compress=F)

quit

#Taxonomic assignment (local DB)
localDBlocation=/share/projects/OCEAN_eDNA/AZORES/databases
localDBprefix=MZGdb_Azores_FISH_12S_V2
#assign their sequences to the taxonomy
mothur "#classify.seqs(fasta=Samples_joined_all.fasta, template=$localDBlocation/$localDBprefix.fasta, taxonomy=$localDBlocation/$localDBprefix.tax, name=Samples_joined_all.names, group=Samples_joined_all.groups, method=wang, cutoff=70, processors=12)"
mothur "#phylotype(taxonomy=Samples_joined_all.$localDBprefix.wang.taxonomy)"
mothur "#make.shared(list=Samples_joined_all.$localDBprefix.wang.tx.list, count=Samples_joined_all.count_table, label=1)"
mothur "#classify.otu(list=Samples_joined_all.$localDBprefix.wang.tx.list, count=Samples_joined_all.count_table, taxonomy=Samples_joined_all.$localDBprefix.wang.taxonomy, label=1)"
sed -E s/";[a-zA-Z]*_unclassified"/";unclassified"/g Samples_joined_all.$localDBprefix.wang.tx.1.cons.taxonomy > Samples_joined_all.$localDBprefix.wang.tx.1.cons_corr.taxonomy
#merge information and create output files
#remove label and phylo_count columns
awk '{$1=$3=""; print $0}' Samples_joined_all.MZGdb_Azores_FISH_12S.wang.tx.shared > file4
#transpose file
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < file4 > file5
paste file5 Samples_joined_all.MZGdb_Azores_FISH_12S.wang.tx.1.cons_corr.taxonomy> file6
awk '{$4=""; print $0}' file6 > output
