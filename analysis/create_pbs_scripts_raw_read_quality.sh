# Script to create PBS scripts for analyzing raw FASTQ files in WGMS kit
# comparison

# Create pbs directory if it doesn't exist
mkdir -p pbs_raw_read_quality
mkdir -p raw_read_quality

# Clean out directory if files already exist
rm -f pbs_raw_read_quality/*.pbs
rm -f pbs_raw_read_quality/submit_pbs_scripts.sh

# Create a bash script to loop through pbs scripts and submit them to the queue
cat > pbs_raw_read_quality/submit_pbs_scripts.sh <<EOF
count=0
for FILE in \`ls *.pbs\`; do
    qsub \${FILE}
    sleep 1.0s
done
EOF

DIRLOC=2019_11_07_FallopianTube_WGBS_Kit_Comparison

for file in `ls ${DIRLOC}/raw_data/*_L000_R1_001.fastq.gz`; do
    base=$(basename -- ${file})
    samp=${base/_L000_R1_001.fastq.gz}
    pair=${file/_L000_R1_001.fastq.gz/_L000_R2_001.fastq.gz}

cat > pbs_raw_read_quality/${samp}.pbs <<EOF
#!/bin/bash
#PBS -N ${samp}
#PBS -j oe
#PBS -o ${DIRLOC}/analysis/raw_read_quality/pbs_raw_read_quality/${samp}.log
#PBS -l nodes=1:ppn=1,mem=75gb,walltime=24:00:00

# Change to working directory
cd ${DIRLOC}/analysis/raw_read_quality

python base_quality.py ${file} ${pair}
EOF
done


