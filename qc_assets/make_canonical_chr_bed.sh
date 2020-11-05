# Create a BED file of the locations for the canonical chromosomes in the 
# gencode release 32p13 reference

REF=references/biscuit_gencode_rel32_p13_all
awk '/^chr[a-zA-Z0-9_.-]*\t/ { printf("%s\t0\t%s\n", $1, $2); }' \
    ${REF}/GRCh38.p13.genome.fa.fai > canonical_chr.bed
