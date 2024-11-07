
# This file contains commands to create scripts for meta-analysing the overlapping metabolites in METAL

while IFS=, read -r file_name_1 file_name_2; do echo -e "SCHEME STDERR\nCUSTOMVARIABLE N\nTRACKPOSITIONS ON\nMARKER SNP\nChromosome CHR\nPosition BP\nALLELE A2 A1\nEFFECT BETA\nSTDERR SE\nPVALUE P\nLABEL N as N\nPROCESS /scratch/prj/proitsi/sumstats/Hysi_metabolites/munged_for_MR/$file_name_1\nMARKER SNP\nChromosome CHR\nPosition BP\nALLELE A2 A1\nEFFECT BETA\nSTDERR SE\nPVALUE P\nLABEL N as N\nPROCESS /scratch/prj/proitsi/sumstats/Chen_metabolites/Chen_munged/$file_name_2\nOUTFILE /scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/${file_name_1%.*}_${file_name_2%.*}_meta .txt\nANALYZE\nQUIT" > "${file_name_1%.*}_${file_name_2%.*}_script.txt"; done < /scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/match_file/chen_hysi_match_file.txt

# bash /scratch/prj/proitsi/lachlan/PhD_project_2/scripts/submit_METAL_scripts.sh
