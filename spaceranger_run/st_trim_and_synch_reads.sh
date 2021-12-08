PAIR_PREPEND='V10A20-012_B1_S21_L001_R'
PAIR_PREPEND='V10A20-012_B1_S21_L002_R'
PAIR_PREPEND='V10A20-012_C1_S22_L001_R'
PAIR_PREPEND='V10A20-012_C1_S22_L002_R'
PAIR_PREPEND='V10A20-013_A1_S23_L001_R'
PAIR_PREPEND='V10A20-013_A1_S23_L002_R'
PAIR_PREPEND='V10A20-013_C1_S24_L001_R'
PAIR_PREPEND='V10A20-013_C1_S24_L002_R'
PAIR_APPEND='_001'
BASE_DIR='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/spatial_transcriptonics/raw/sequence_data/211130_spikeinrun1'
ORIG_DIR=${BASE_DIR}'/'
ORIG_RAWS=${ORIG_DIR}'/raw_data/'
TRIMMED_DIR=${BASE_DIR}'_trimmed/'
TRIMMED_RAWS=${TRIMMED_DIR}'/raw_data/'
FILTERED_DIR=${BASE_DIR}'_trimmed_filtered/'
FILTERED_RAWS=${FILTERED_DIR}'/raw_data/'
SYNCHED_DIR=${BASE_DIR}'_trimmed_filtered_synched/'
SYNCHED_RAWS=${SYNCHED_DIR}'/fastq/'
SORTED_DIR=${BASE_DIR}'_trimmed_filtered_synched_sorted/'
SORTED_RAWS=${SORTED_DIR}'/raw_data/'

READS_R2_LOC=${FILTERED_DIR}${PAIR_PREPEND}'2'${PAIR_APPEND}'.reads'
IDS_FROM_R2=${FILTERED_DIR}${PAIR_PREPEND}'2'${PAIR_APPEND}'.ids'
IDS_BOTH=${FILTERED_DIR}${PAIR_PREPEND}'X'${PAIR_APPEND}'.ids'

CONDA_TRIM_ENV='visiumtrim_env'
TRIMMER_LOC='/groups/umcg-franke-scrna/tmp01/software/VisiumTrim/TSO_polyA_trimming.sh'
CONDA_FILTER_ENV='seqtk'
SEQKIT_LOC='/groups/umcg-franke-scrna/tmp01/software/seqkit_2.1.0/seqkit'

# load trimming environment
conda activate ${CONDA_TRIM_ENV}

# grab read two and trim it
READ2_ORIG=${ORIG_RAWS}${PAIR_PREPEND}'2'${PAIR_APPEND}'.fastq.gz'
READ2_TRIMMED=${TRIMMED_RAWS}${PAIR_PREPEND}'2'${PAIR_APPEND}'.fastq.gz'
READ2_FILTERED=${FILTERED_RAWS}${PAIR_PREPEND}'2'${PAIR_APPEND}'.fastq.gz'

${TRIMMER_LOC} \
${READ2_ORIG} \
-o ${READ2_TRIMMED}

# deactivate this environment
conda deactivate


# load filtering environment
conda activate ${CONDA_FILTER_ENV}

# filter the R2 reads on length
for ii in $( ls ${TRIMMED_RAWS}${PAIR_PREPEND}2* -1 ); do zcat ${ii} | ${SEQKIT_LOC} seq ${ii} -m 85 > ${FILTERED_RAWS}/$( basename ${ii} | sed "s/\.gz//g" ); done


# now grab the IDs that are present in R2
zcat ${READ2_FILTERED} | zgrep @  > ${READS_R2_LOC}
# remove the @
sed 's/@//g' ${READS_R2_LOC} > ${IDS_FROM_R2}
# get only the first column
cut -d$' ' -f1 ${IDS_FROM_R2} > ${IDS_BOTH}

# I'm not sure if seqkit likes compressed files, so let's decompress first
READ1_ORIG=${ORIG_RAWS}${PAIR_PREPEND}'1'${PAIR_APPEND}'.fastq.gz'
READ1_UNCOMPRESSED=${SYNCHED_RAWS}${PAIR_PREPEND}'1'${PAIR_APPEND}'.fastq.original'
gunzip -c ${READ1_ORIG} > ${READ1_UNCOMPRESSED}

# this will be the filtered R1
READ1_SYNCHED=${SYNCHED_RAWS}${PAIR_PREPEND}'1'${PAIR_APPEND}'.fastq'
# then do the filtering
seqtk subseq \
${READ1_UNCOMPRESSED} \
${IDS_BOTH} > \
${READ1_SYNCHED}

READ1_SYNCHED_COMPRESSED=${READ1_SYNCHED}'.gz'

# comparessed is better
gzip -c ${READ1_SYNCHED} > ${READ1_SYNCHED_COMPRESSED}
# deactivate this environment
conda deactivate

# create a symlink in the synch direcgtory, that points to the R2 sequences with the same filtering
ln -s ${READ2_TRIMMED} ${SYNCHED_RAWS}'/'
# clean up
rm ${READ1_UNCOMPRESSED}

# last thing to do is sorting
READ2_SORTED=${SORTED_RAWS}${PAIR_PREPEND}'2'${PAIR_APPEND}'.fastq.gz'
READ1_SORTED=${SORTED_RAWS}${PAIR_PREPEND}'1'${PAIR_APPEND}'.fastq.gz'

# output to the new folder
zcat ${READ1_SYNCHED_COMPRESSED} \
 | paste - - - - \
  | sort -k1,1 -S 3G \
  | tr '\t' '\n' \
   | gzip > ${READ1_SORTED}
zcat ${READ2_FILTERED} \
 | paste - - - - \
  | sort -k1,1 -S 3G \
   | tr '\t' '\n' \
    | gzip > ${READ2_SORTED}
