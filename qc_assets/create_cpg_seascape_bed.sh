# Reference location
REFLOC=references/biscuit_gencode_rel32_p13_all

# CpG Islands from UCSC on only the canonical chromosomes
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz | \
gunzip -c | \
awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4; }' | \
awk '{ if ($1 ~ /^chr[1234567890XYM]{1,2}$/) { print } }' | \
bedtools sort -i - \
> cpg_islands.bed

# Create middles for finding locations
bedtools slop \
    -i cpg_islands.bed \
    -g ${REFLOC}/GRCh38.p13.genome.fa.fai \
    -b 2000 | \
bedtools merge -i - \
> shores.tmp.bed

bedtools slop \
    -i cpg_islands.bed \
    -g ${REFLOC}/GRCh38.p13.genome.fa.fai \
    -b 4000 | \
bedtools merge -i - \
> shelves.tmp.bed

# CpG Open Seas (intervening locations)
sort -k1,1 -k2,2n ${REFLOC}/GRCh38.p13.genome.fa.fai | \
bedtools complement -L -i shelves.tmp.bed -g - \
> cpg_open_seas.bed

# CpG Shelves (shores +/- 2kb)
bedtools subtract \
    -a shelves.tmp.bed \
    -b shores.tmp.bed \
> cpg_shelves.bed

# CpG Shores (island +/- 2kb)
bedtools subtract \
    -a shores.tmp.bed \
    -b cpg_islands.bed \
> cpg_shores.bed

# Clean up temporary files
rm -f shelves.tmp.bed shores.tmp.bed
