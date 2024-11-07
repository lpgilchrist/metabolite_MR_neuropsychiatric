#!/bin/bash
#SBATCH --job-name=chen_hysi_meta
#SBATCH --output=/scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/match_file/chen_hysi_meta.log
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition cpu
#SBATCH --time=0-48:00

# Directory where your script files are located
script_dir="/scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/match_file/chen_hysi_scripts"

# Loop through the script files in the directory
for script_file in "$script_dir"/*_script.txt; do
  # Extract the base name of the script file
  base_name=$(basename "$script_file" _munged.txt_script.txt)

  # Remove "_munged.txt" from the base name
  clean_name=$(echo "$base_name" | sed 's/_munged.txt//')


  
  # Check if the output file already exists
  if [ -f "/scratch/prj/proitsi/sumstats/metabolites_meta/chen_hysi_meta/${clean_name}_meta1.txt" ]; then
    echo "Output file ${clean_name}_meta1.txt already exists. Skipping..."
  else
    echo "Output file ${clean_name}_meta1.txt does not exist. Running analysis!"
    # If output file doesn't exist, run the script
    /scratch/prj/ukbiobank/usr/Lachlan/polygenic_paper/METAL/METAL/build/bin/metal "$script_file"
  fi
done

# sbatch -p cpu /scratch/prj/proitsi/lachlan/PhD_project_2/scripts/submit_METAL_scripts.sh
