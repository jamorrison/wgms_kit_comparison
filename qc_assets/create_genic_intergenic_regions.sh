# Reference used
REF=references/biscuit_gencode_rel32_p13_all/GRCh38.p13.genome.fa.fai

# Create merged BED file of genic regions from genes.bed.gz
zcat hg38/genes.bed.gz | \
awk '{ if ($4 == "gene") { print } }' | \
bedtools merge -i - | \
gzip -c > hg38/genic_regions.bed.gz

# Need genome file sorted in the same way as BED file
sort -k1,1 -k2,2n ${REF} | \
bedtools complement -L -i hg38/genic_regions.bed.gz -g - | \
gzip -c > hg38/intergenic_regions.bed.gz
