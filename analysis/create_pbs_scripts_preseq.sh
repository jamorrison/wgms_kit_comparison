# Script to create PBS scripts and a submission script

# Create pbs directory if it doesn't exist
mkdir -p pbs_preseq
mkdir -p preseq

# Clean out directory if files already exist
rm -f pbs_preseq/*.pbs
rm -f pbs_preseq/submit_pbs_scripts.sh

# Create a bash script to loop through pbs scripts and submit them to the queue
cat > pbs_preseq/submit_pbs_scripts.sh <<EOF
count=0
for FILE in \`ls *.pbs\`; do
    qsub \${FILE}
    sleep 1.0s
done
EOF

ALIGNDIR=2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/align

# Use pbs template file and file names from aligned BAMs to create pbs scripts
for FILE in `ls ${ALIGNDIR}/*.sorted.markdup.bam`; do
    base="$(basename -- $FILE)"
    samp=${base%.sorted.markdup.bam}

    sed "s/QQQ/${samp}/g" pbs_template_preseq > pbs_preseq/${samp}.pbs
done
