BEDDIR=2019_11_07_FallopianTube_WGBS_Kit_Comparison/qc_assets/seascape
ISLD=${BEDDIR}/cpg_islands.bed
SEAS=${BEDDIR}/cpg_open_seas.bed
SHLV=${BEDDIR}/cpg_shelves.bed
SHOR=${BEDDIR}/cpg_shores.bed

DIRLOC=2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/align

mkdir -p seashore_bed_files

for FILE in `ls ${DIRLOC}/*.cg.sorted.mergecg.bed.gz`; do
    base="$(basename -- $FILE)"
    samp=${base%.cg.sorted.mergecg.bed.gz}

    bedtools intersect -a ${FILE} -b ${ISLD} -wa | gzip > ${samp}.island.bed.gz
    bedtools intersect -a ${FILE} -b ${SEAS} -wa | gzip > ${samp}.open_seas.bed.gz
    bedtools intersect -a ${FILE} -b ${SHLV} -wa | gzip > ${samp}.shelves.bed.gz
    bedtools intersect -a ${FILE} -b ${SHOR} -wa | gzip > ${samp}.shores.bed.gz
done
