#!/bin/bash
SPACE_BIN='/groups/umcg-franke-scrna/tmp01/software/spaceranger-1.3.1/spaceranger'
RAW_DIR='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/spatial_transcriptonics/raw/'
INPUT_DIR=${RAW_DIR}'/sequence_data/'
JPEG_DIR_APPEND='/images/'
FASTQ_DIR_APPEND='raw_data/'
ALIGN_DIR_APPEND='/alignment/'
REF_DIR='/groups/umcg-franke-scrna/tmp01/external_datasets/refdata-gex-GRCh38-2020-A/'
OUT_DIR='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/spatial_transcriptonics/processed/alignment/spaceranger_out'
SCR_DIR='releases/blokland-2020/v1/spatial_transcriptonics/processed/alignment/spaceranger_out'
JOB_DIR='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/spatial_transcriptonics/processed/jobs/'

# here are the slides to look at
#SLIDES=('V10A20-013')
SLIDES=('V10A20-012')
#AREAS=('A1')
#AREAS=('C1')
AREAS=('B1', 'C1')

# these are the run names
RUN_NAMES=('211130_spikeinrun1_trimmed_filtered_synched_sorted_cr100')

# check each run
for run in ${RUN_NAMES[*]}
  do
    # create the job directory
    mkdir -p ${JOB_DIR}'/'${run}'/'
    # check each slide
    for slide in ${SLIDES[*]}
      do
        # create the
        for area in ${AREAS[*]}
          do
            # /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/spatial_transcriptonics/raw/sequence_data/211201_spikeinrun2/fastq/
            fastq=${INPUT_DIR}${run}'/'${FASTQ_DIR_APPEND}
            # V10A20-014_A1
            sample_id=${slide}'_'${area}
            # V10A20-014_A1
            sample=${slide}'_'${area}
            # /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/spatial_transcriptonics/raw/sequence_data/211201_spikeinrun2/images/V10A20-014_A1.jpg
            jpeg=${INPUT_DIR}'/'${run}'/'${JPEG_DIR_APPEND}'/'${sample}'.jpg'
            # V10A20-014
            slide=${slide}
            # A1
            area=${area}
            space_out_full=${OUT_DIR}${run}'/'${sample_id}
            mkdir -p ${space_out_full}
            # /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/spatial_transcriptonics/raw/sequence_data/211201_spikeinrun2/alignment/V10A20-014_A1.json
            loupe_dir=${INPUT_DIR}'/'${run}
            loupe_full=${loupe_dir}'/'${ALIGN_DIR_APPEND}'/'${sample}'.json'
            # where the directories and err/out files should end up
            out_sbatch=${JOB_DIR}'/'${run}'/run_'${sample}'_SBATCH.sh'
            out=${JOB_DIR}'/'${run}'/run_'${sample}'.out'
            err=${JOB_DIR}'/'${run}'/run_'${sample}'.err'
            job_name='map_'${sample}
            # paste it all together
            echo '#!/usr/bin/env bash' > ${out_sbatch}
            echo '#SBATCH --job-name=map_'${job_name} >> ${out_sbatch}
            echo '#SBATCH --output='${out} >> ${out_sbatch}
            echo '#SBATCH --error='${err} >> ${out_sbatch}
            echo '#SBATCH --time=23:59:59'  >> ${out_sbatch}
            echo '#SBATCH --cpus-per-task=12' >> ${out_sbatch}
            echo '#SBATCH --mem=128gb' >> ${out_sbatch}
            echo '#SBATCH --nodes=1'  >> ${out_sbatch}
            echo '#SBATCH --open-mode=append' >> ${out_sbatch}
            echo '#SBATCH --export=NONE'  >> ${out_sbatch}
            echo '#SBATCH --get-user-env=L' >> ${out_sbatch}
            echo '#SBATCH --tmp=4096mb' >> ${out_sbatch}
            # create the scratch folder
            echo 'mkdir -p ${TMPDIR}/'${SCR_DIR} >> ${out_sbatch}
            # go to the scratch folder
            echo 'cd ${TMPDIR}/'${SCR_DIR} >> ${out_sbatch}
            echo ${SPACE_BIN}' count --id='${sample_id} '\' >> ${out_sbatch}
            echo '--fastqs='${fastq}' \'  >> ${out_sbatch}
            echo '--transcriptome='${REF_DIR}' \' >> ${out_sbatch}
            echo '--sample='${sample}' \' >> ${out_sbatch}
            echo '--image='${jpeg}' \' >> ${out_sbatch}
            echo '--slide='${slide}' \' >> ${out_sbatch}
            echo '--area='${area}' \' >> ${out_sbatch}
            echo '--loupe-alignment='${loupe_full} >> ${out_sbatch}
            # copy the result to tmp
            echo 'cp -r '${sample_id}' '${OUT_DIR}'/'${run}'/' >> ${out_sbatch}
            # tell people what we did
            echo 'wrote job script to '${out_sbatch}
          done
      done
  done
