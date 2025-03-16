# Copy input files from previous step
mothur

# Create logfile with all manip that are going to be done
set.logfile(name=Ceph18S.log)

# Analyse reads (length, ambiguities) and count number of reads per group
summary.seqs(fasta=Samples_joined.fasta, processors=12)
count.groups(group=Samples_joined.groups)

# Dereplicate
unique.seqs(fasta=Samples_joined.fasta)
summary.seqs(fasta=Samples_joined.unique.fasta, name=Samples_joined.names, processors=12)

# Remove chimeras
chimera.uchime(fasta=Samples_joined.unique.fasta, name=Samples_joined.names, group=Samples_joined.groups, processors=12)
remove.seqs(accnos=Samples_joined.unique.denovo.uchime.accnos, fasta=Samples_joined.unique.fasta, name=Samples_joined.names, group=Samples_joined.groups, dups=T)
summary.seqs(fasta=Samples_joined.unique.pick.fasta, name=Samples_joined.pick.names, processors=12)
count.groups(group=Samples_joined.pick.groups)

# Remove sequences with N
screen.seqs(fasta=Samples_joined.unique.pick.fasta, maxambig=0, minlength=120, maxlength=200, name=Samples_joined.pick.names, group=Samples_joined.pick.groups, processors=12)
summary.seqs(fasta=Samples_joined.unique.pick.good.fasta, name=Samples_joined.pick.good.names, processors=12)
count.groups(group=Samples_joined.pick.good.groups)

# Rename files for clarity
system(cp Samples_joined.unique.pick.good.fasta Samples_joined_all.fasta)
system(cp Samples_joined.pick.good.names Samples_joined_all.names)
system(cp Samples_joined.pick.good.groups Samples_joined_all.groups)

# Create count_table with the new files
make.table(name=Samples_joined_all.names, group=Samples_joined_all.groups, compress=f)

quit()

#Taxonomic assignment (local DB)
localDBlocation=/share/projects/OCEAN_eDNA/AZORES/databases
localDBprefix=MZGdb_Azores_CEPH_18S_V2
#assign their sequences to the taxonomy
mothur "#classify.seqs(fasta=Samples_joined_all.fasta, template=$localDBlocation/$localDBprefix.fasta, taxonomy=$localDBlocation/$localDBprefix.tax, name=Samples_joined_all.names, group=Samples_joined_all.groups, method=wang, cutoff=70, processors=12)"
mothur "#phylotype(taxonomy=Samples_joined_all.$localDBprefix.wang.taxonomy)"
mothur "#make.shared(list=Samples_joined_all.$localDBprefix.wang.tx.list, count=Samples_joined_all.count_table, label=1)"
mothur "#classify.otu(list=Samples_joined_all.$localDBprefix.wang.tx.list, count=Samples_joined_all.count_table, taxonomy=Samples_joined_all.$localDBprefix.wang.taxonomy, label=1)"

awk '{$1=$3=""; print $0}' Samples_joined_all.MZGdb_Azores_CEPH_18S.wang.tx.shared > file1
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < file1 > file2
paste file2 Samples_joined_all.MZGdb_Azores_CEPH_18S.wang.tx.1.cons.taxonomy > file3
#check merge is correct
awk '{$4=""; print $0}' file3 > OTU_table_local_AZORES
rm file*
