# Script to create PBS scripts for trimming FASTQ files for WGMS kit comparison

DIRLOC=2019_11_07_FallopianTube_WGBS_Kit_Comparison

for file in `ls ${DIRLOC}/raw_data/*_L000_R1_001.fastq.gz`; do
    base=$(basename -- ${file})
    samp=${base/_L000_R1_001.fastq.gz}
    pair=${file/_L000_R1_001.fastq.gz/_L000_R2_001.fastq.gz}

cat > ${samp}.pbs <<EOF
#!/bin/bash
#PBS -N ${samp}
#PBS -j oe
#PBS -o ${DIRLOC}/trimmed_fastq/${samp}.log
#PBS -l nodes=1:ppn=16,mem=75gb,walltime=24:00:00

# Change to working directory
cd ${DIRLOC}/trimmed_fastq

trim_galore \\
    --illumina \\
    --trim-n \\
    --paired \\
    --fastqc \\
    --fastqc_args "--noextract" \\
    --cores 4 \\
EOF
    
    if [[ ${samp} == "FtubeAswift10ngRep2" ]] || \
       [[ ${samp} == "FtubeAswift10ng" ]] || \
       [[ ${samp} == "FtubeAswiftRep2" ]] || \
       [[ ${samp} == "FtubeAswift" ]] || \
       [[ ${samp} == "FtubeBswift10ngRep2" ]] || \
       [[ ${samp} == "FtubeBswift10ng" ]] || \
       [[ ${samp} == "FtubeBswiftRep2" ]] || \
       [[ ${samp} == "FtubeBswift" ]]; then
        echo "    --clip_R2 14 \\" >> ${samp}.pbs
    fi

    echo "    ${file} \\" >> ${samp}.pbs
    echo "    ${pair}" >> ${samp}.pbs

done
