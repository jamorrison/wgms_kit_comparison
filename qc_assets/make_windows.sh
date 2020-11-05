REFLOC="references/biscuit_gencode_rel32_p13_all"

bedtools makewindows -w 100000 \
    -g ${REFLOC}/GRCh38.p13.genome.fa.fai | \
sort -k1,1 -k2,2n > 100kb_hg38_windows.bed
