# Script to create PBS scripts and a submission script

# Create pbs directory if it doesn't exist
mkdir -p pbs_align
mkdir -p align

# Clean out directory if files already exist
rm -f pbs_align/*.pbs
rm -f pbs_align/submit_pbs_scripts.sh

# Create a bash script to loop through pbs scripts and submit them to the queue
cat > pbs_align/submit_pbs_scripts.sh <<EOF
count=0
for FILE in \`ls *.pbs\`; do
    qsub \${FILE}
    sleep 1.0s
done
EOF

TRIMDIR=2019_11_07_FallopianTube_WGBS_Kit_Comparison/trimmed_fastq

# Use pbs template file and file names from trimmed_fastqs to create pbs scripts
for FILE in `ls ${TRIMDIR}/*_L000_R1_001_val_1.fq.gz`; do
    base="$(basename -- $FILE)"
    samp=${base%_L000_R1_001_val_1.fq.gz}

    sed "s/QQQ/${samp}/g" pbs_template_align > pbs_align/${samp}.pbs
done
