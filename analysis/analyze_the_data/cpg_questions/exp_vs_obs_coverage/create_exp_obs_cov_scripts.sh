# Create pbs directory if it doesn't exist
mkdir -p pbs_mappability

# Clean out directory if files already exist
rm -f pbs_mappability/*.pbs
rm -f pbs_mappability/submit_pbs_scripts.sh

# Create a bash script to loop through pbs scripts and submit them to the queue
cat > pbs_mappability/submit_pbs_scripts.sh <<EOF
count=0
for FILE in \`ls *.pbs\`; do
    qsub \${FILE}
    sleep 1.0s
done
EOF

DIRLOC=2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/analyze_the_data/cpg_questions/exp_vs_obs_coverage
BAMLOC=2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/align

for BAM in `ls ${BAMLOC}/*.sorted.markdup.bam`; do
    base="$(basename -- $BAM)"
    samp=${base%.sorted.markdup.bam}

cat > pbs_mappability/${samp}_expobs_cov.pbs <<EOF
#!/bin/bash
#PBS -N obsexp
#PBS -e ${samp}_expobs_cov.stderr
#PBS -o ${samp}_expobs_cov.stdout
#PBS -l nodes=1:ppn=1,mem=100gb,walltime=24:00:00

cd ${DIRLOC}

# BED file with CpG locations, only includes chromosomes 1-22/X/Y/M
CPGSBED=2019_11_07_FallopianTube_WGBS_Kit_Comparison/qc_assets/bismap/cpg_bismap.bed.gz
CGISBED=2019_11_07_FallopianTube_WGBS_Kit_Comparison/qc_assets/bismap/cgi_bismap.bed.gz
RMSKBED=2019_11_07_FallopianTube_WGBS_Kit_Comparison/qc_assets/bismap/rmsk_bismap.bed.gz
EXONBED=2019_11_07_FallopianTube_WGBS_Kit_Comparison/qc_assets/bismap/exon_bismap.bed.gz
GENEBED=2019_11_07_FallopianTube_WGBS_Kit_Comparison/qc_assets/bismap/genic_regions_bismap.bed.gz
INTRBED=2019_11_07_FallopianTube_WGBS_Kit_Comparison/qc_assets/bismap/intergenic_regions_bismap.bed.gz
BISMBED=2019_11_07_FallopianTube_WGBS_Kit_Comparison/qc_assets/bismap/k100.bismap.bedgraph.gz

# (# bases in feature) [i.e. sum of weights times element length, for specific regions]
# column 2: start of element; column 3: end of element; column 4: element mappability score
CPGSLEN=\$(zcat \${CPGSBED} | awk '{ sum += \$4*(\$3-\$2) } END { print sum }')
CGISLEN=\$(zcat \${CGISBED} | awk '{ sum += \$4*(\$3-\$2) } END { print sum }')
RMSKLEN=\$(zcat \${RMSKBED} | awk '{ sum += \$4*(\$3-\$2) } END { print sum }')
EXONLEN=\$(zcat \${EXONBED} | awk '{ sum += \$4*(\$3-\$2) } END { print sum }')
GENELEN=\$(zcat \${GENEBED} | awk '{ sum += \$4*(\$3-\$2) } END { print sum }')
INTRLEN=\$(zcat \${INTRBED} | awk '{ sum += \$4*(\$3-\$2) } END { print sum }')

# (# bases in genome) [i.e. sum of weights times element length, for all available elements]
# column 2: start of element; column 3: end of element; column 4: element mappability score
REFLEN=\$(zcat \${BISMBED} | awk '{ sum += \$4*(\$3-\$2) } END { print sum }')

# Find genomic coverage
samtools view -hb -F 0x4 -q 40 ${BAM} | \\
bedtools genomecov -bg -ibam stdin | \\
gzip -c > ${samp}.genomecov.tmp.bed.gz

# (# mapped bases) [i.e. sum of weights times coverage times overlap, for the whole genome]
# column 4: coverage, column 9: # number of overlapped bases in region
MAPBASE=\$(bedtools intersect -a ${samp}.genomecov.tmp.bed.gz -b \${BISMBED} -wo | \\
    awk '{ sum += \$9*\$4 } END { print sum }')

# (# bases mapped to feature) [i.e. sum of weights times coverage times overlap, for specific regions]
# column 4: coverage, column 9: # number of overlapped bases in region
CPGSBASE=\$(bedtools intersect -a ${samp}.genomecov.tmp.bed.gz -b \${CPGSBED} -wo | \\
    awk '{ sum += \$9*\$4 } END { print sum }')
CGISBASE=\$(bedtools intersect -a ${samp}.genomecov.tmp.bed.gz -b \${CGISBED} -wo | \\
    awk '{ sum += \$9*\$4 } END { print sum }')
RMSKBASE=\$(bedtools intersect -a ${samp}.genomecov.tmp.bed.gz -b \${RMSKBED} -wo | \\
    awk '{ sum += \$9*\$4 } END { print sum }')
EXONBASE=\$(bedtools intersect -a ${samp}.genomecov.tmp.bed.gz -b \${EXONBED} -wo | \\
    awk '{ sum += \$9*\$4 } END { print sum }')
GENEBASE=\$(bedtools intersect -a ${samp}.genomecov.tmp.bed.gz -b \${GENEBED} -wo | \\
    awk '{ sum += \$9*\$4 } END { print sum }')
INTRBASE=\$(bedtools intersect -a ${samp}.genomecov.tmp.bed.gz -b \${INTRBED} -wo | \\
    awk '{ sum += \$9*\$4 } END { print sum }')

# Observed / Expected ratio
# [(# bases mapped to feature) / (# mapped bases)] / [(# bases in feature) / (# bases in genome)]

echo -e "${samp}\t\${REFLEN}\t\${CPGSLEN}\t\${CGISLEN}\t\${RMSKLEN}\t\${EXONLEN}\t\${GENELEN}\t\${INTRLEN}\t\${MAPBASE}\t\${CPGSBASE}\t\${CGISBASE}\t\${RMSKBASE}\t\${EXONBASE}\t\${GENEBASE}\t\${INTRBASE}"
EOF
done
