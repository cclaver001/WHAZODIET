# Whale poo and water eDNA
This is the repository where the scripts and input files associated to the paper "Mesopelagic fish and cephalopods in the diet of rorquals (Balaenoptera spp) and sperm whales (Physeter macrocephalus) around the Azores using faecal and environmental DNA" by Cristina Claver, Leire G. de Amézaga, Iñaki Mendibil, Oriol Canals, Rui Prieto, Irma Cascão, Cláudia Oliveira, Mónica A. Silva, Naiara Rodríguez-Ezpeleta are locaed. 
*Link to the preprint:  https://doi.org/10.1101/2025.03.14.643344*

## Scripts and files
Scripts used for raw data (available at SRA under Bioproject PRJNA1244105) preparation and analysis.

*Demultiplexing step:*
- 1.script_Demultiplex_run1.sh (for cephalopod data) and  1.script_Demultiplex_run2.sh (for fish data)
- Input files: M*_barcodes_run1.sh (for cephalopod data) and M*_barcodes_run2.sh (for fish data)
- Input files: Mix_path_to_R1_R2_files_run1.txt (for cephalopod data) and  Mix_path_to_R1_R2_files_run2.txt (for fish data)

*Cleaning step:*
- 2.prepareReads_Ceph18S.sh (for cephalopod data) and 2.prepareReads_Teleo.sh (for fish data)
- Input files: Ceph18S_samples_20230113.txt (for cephalopod data) and Teleo_samples_20230703 (for fish data)

*Tax assignment step:*
- 3.taxAssignment_Ceph18S.sh (for cephalopod data) and  3.taxAssignment_Teleo.sh (for fish data)
- Input files (for cephalopod data): MZGdb_Azores_CEPH_18S_V2.tax and  MZGdb_Azores_CEPH_18S_V2.fasta
- Input files (for fish data): MZGdb_Azores_FISH_12S_V2.tax and MZGdb_Azores_FISH_12S_V2.fasta

## Reference databases
Reference databases used for taxonomic assignment are porvided. Each reference databse has a .tax file and a .fasta file. 
- Cephalopod database -> MZGdb_Azores_CEPH_18S_V2
- Fish database -> MZGdb_Azores_FISH_12S_V2

## Contact
Cristina Claver (cclaver@azti.es)
