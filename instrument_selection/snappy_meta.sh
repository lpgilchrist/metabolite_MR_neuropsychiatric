#!/bin/bash
#SBATCH --job-name=snappy_matched_meta_Depression
#SBATCH --output=/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/snappy_out/meta_proxies/snappy_matched_meta_all.log
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --partition cpu
#SBATCH --time=0-48:00

module load plink/1.9-beta6.10-gcc-10.3.0

module load perl/5.34.1-gcc-10.3.0

SNAPPY_SCRIPT="/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/snappy"
INSTRUMENT_FILE_DIR="/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/meta_instruments"
FILE_DIR="/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/snappy_snps/snappy_snps_meta"
REF_PANEL_DIR="/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/snappy_reference_panel"
REF_PANEL_KEEP="/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/EUR.fam"
PLINK_MIN_R2=0.8
OUTPUT_DIRECTORY="/scratch/prj/proitsi/lachlan/PhD_project_2/snappy/snappy_out/meta_proxies"

# Get the list of instrument files
instrument_files=("$INSTRUMENT_FILE_DIR"/*_instruments.txt)
# List of target values
targets=("Depression" "AD" "PD" "ALS" "BP" "SCZ" "ANX" "MS")

# Loop over the instrument files
for instrument_file in "${instrument_files[@]}"; do
  echo "Processing init file: $instrument_file"

  # Extract the filename without the extension
  instrument_filename=$(basename "$instrument_file" _instruments.txt)

  # Loop over the list of target outcomes
  for target in "${targets[@]}"; do
    # Generate the SNP filename based on the current target
    snp_file="$FILE_DIR/meta_${target}_snps.txt"
    snp_filename=$(basename "$snp_file" | sed "s/^meta_\(.*\)_snps.txt$/\1/")

    # Check if the corresponding SNP file exists
    if [[ -f "$snp_file" ]]; then
      echo "Corresponding SNP file found: $snp_file"

      # Construct the output file name
      output_file="${instrument_filename}_${snp_filename}_matched_proxies.txt"

      # Check if the output file already exists
      if [[ -f "$OUTPUT_DIRECTORY/$output_file" ]]; then
        printf "The match/proxy file for ${output_file} exists, skipping...\n"
      else
        printf "\n${output_file} does not exist, starting matching and proxy IV identification...\n"
        "$SNAPPY_SCRIPT" init "$instrument_file" "$snp_file" | "$SNAPPY_SCRIPT" match | \
          "$SNAPPY_SCRIPT" proxy "$REF_PANEL_DIR" "$REF_PANEL_KEEP" "$PLINK_MIN_R2" > "$OUTPUT_DIRECTORY/$output_file"
        echo "Finished processing SNP file: $snp_file"
      fi
    else
      echo "Corresponding SNP file not found for: $instrument_file"
    fi
  done

  echo "Finished processing init file: $instrument_file"
done




# sbatch -p cpu /scratch/prj/proitsi/lachlan/PhD_project_2/scripts/snappy_meta_all.sh

