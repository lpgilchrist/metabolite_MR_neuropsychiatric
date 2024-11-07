#!/bin/bash

#SBATCH --job-name=mtCOJO_all
#SBATCH --output=/scratch/prj/proitsi/lachlan/PhD_project_2/mtCOJO/output/mtCOJO_all.log
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --partition cpu
#SBATCH --time=0-48:00

# Define paths
BASE_DIR="/scratch/prj/proitsi/lachlan/PhD_project_2/mtCOJO"
BFILE="/scratch/prj/proitsi/lachlan/PhD_project_2/reference_panel/EUR"
LD="${BASE_DIR}/ld_ref/eur_w_ld_chr/"
OUTPUT="${BASE_DIR}/output"
mtCOJO_path_files="$BASE_DIR/mtCOJO_path_files"

mtcojo_sum_files=("$mtCOJO_path_files"/*_mtcojo_input.txt)

# Loop over the files
for file in "${mtcojo_sum_files[@]}"; do
    filename=$(basename "$file" _mtcojo_input.txt)
    OUTPUT_PATH="${OUTPUT}/${filename}"

    # Check if output file already exists
    if [ -f "${OUTPUT_PATH}.mtcojo.cma" ]; then
        echo "Output file ${OUTPUT_PATH}.mtcojo.cma already exists. Skipping..."
    else
        # Run mtCOJO
        /scratch/prj/proitsi/lachlan/PhD_project_2/mtCOJO/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
        --bfile ${BFILE} \
        --mtcojo-file ${file} \
        --ref-ld-chr ${LD} \
        --w-ld-chr ${LD} \
        --threads 10 \
        --out "${OUTPUT_PATH}"
    fi
done


# submit: sbatch -p cpu /scratch/prj/proitsi/lachlan/PhD_project_2/scripts/mtCOJO_all.sh











