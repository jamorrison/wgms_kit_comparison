# Find # CpGs, Coverage distribution, and # covered CpGs from Total CpG files
function do_awk_total {
    outdirc=${1}
    disfile=${2}
    inpfile=${3}
    outfile=${4}
    readtyp=${5}
    awk \
        -v out1="${disfile}" \
        -v rtyp=${readtyp} \
    '{
        n_rows += 1;
        covg[$4] += 1;
        if ($4>0) { n_cpgs += 1 }
     } END {
        for (cov in covg) { print cov"\t"covg[cov] };
        print "TotalCpGs\t"rtyp"\t"n_rows"\t"n_cpgs >> out1
    }' ${inpfile} | \
    sort -k1,1n -T ${outdirc} >> ${outfile}
}

# Find # CpGs, Coverage distribution, and # covered CpGs from CpG files that
# need to have bedtools intersect run (exon, genic, repeat masked, cpg island)
function do_awk_intersect {
    outdirc=${1}
    disfile=${2}
    inpfile=${3}
    iscfile=${4}
    outfile=${5}
    cpgtype=${6}
    readtyp=${7}
    bedtools intersect -sorted \
        -a ${inpfile} \
        -b <(bedtools merge -i ${iscfile}) | \
    awk \
        -v out1="${disfile}" \
        -v cpgt="${cpgtype}" \
        -v rtyp=${readtyp} \
    '{
        n_rows += 1;
        covg[$4] += 1;
        if ($4>0) { n_cpgs += 1 }
     } END {
        for (cov in covg) { print cov"\t"covg[cov] };
        print cpgt"\t"rtyp"\t"n_rows"\t"n_cpgs >> out1
    }' - | \
    sort -k1,1n -T ${outdirc} >> ${outfile}
}

# Export functions for GNU parallel usage
export -f do_awk_total
export -f do_awk_intersect

# Find various CpG counts and coverage
function count_cpgs {
    echo "Started at `date`"

    mkdir -p ${outdir}

    # Find genome coverage BED file for all reads (all) and MAPQ >= 40 reads (q40)
    # Pick off CpG locations from genomecov BED file
    samtools view -hb ${in_bam} | parallel -j4 -k --tmpdir ${outdir} --pipe --tee {} ::: \
        "bedtools genomecov -bga -split -ibam stdin | grep chr | LC_ALL=C sort -k1,1 -k2,2n -T ${outdir} | bedtools intersect -wo -sorted -a ${BISCUIT_CPGS} -b stdin | bedtools groupby -g 1-3 -c 7 -o min > ${outdir}/${sample}_cpg_all.bed" \
        "samtools view -q 40 -hb | bedtools genomecov -bga -split -ibam stdin | grep chr | LC_ALL=C sort -k1,1 -k2,2n -T ${outdir} | bedtools intersect -wo -sorted -a ${BISCUIT_CPGS} -b stdin | bedtools groupby -g 1-3 -c 7 -o min > ${outdir}/${sample}_cpg_q40.bed"

    # Print header lines to output files
    echo -e "Region\tReadType\tnTotalCpGs\tnCoveredCpGs" > ${outdir}/${sample}_cpg_dist_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_total_all_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_total_q40_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_exon_all_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_exon_q40_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_rmsk_all_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_rmsk_q40_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_gene_all_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_gene_q40_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_cgis_all_table.txt
    echo -e "Coverage\tnCpGs" > ${outdir}/${sample}_cpgs_cgis_q40_table.txt

    # Create CpG count and coverage files
    parallel -j10 -k --tmpdir ${outdir} --tee {} ::: \
        "do_awk_total ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_all.bed ${outdir}/${sample}_cpgs_total_all_table.txt All" \
        "do_awk_total ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_q40.bed ${outdir}/${sample}_cpgs_total_q40_table.txt Q40" \
        "do_awk_intersect ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_all.bed ${BISCUIT_EXON} ${outdir}/${sample}_cpgs_exon_all_table.txt ExonicCpGs All" \
        "do_awk_intersect ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_q40.bed ${BISCUIT_EXON} ${outdir}/${sample}_cpgs_exon_q40_table.txt ExonicCpGs Q40" \
        "do_awk_intersect ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_all.bed ${BISCUIT_RMSK} ${outdir}/${sample}_cpgs_rmsk_all_table.txt RepeatCpGs All" \
        "do_awk_intersect ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_q40.bed ${BISCUIT_RMSK} ${outdir}/${sample}_cpgs_rmsk_q40_table.txt RepeatCpGs Q40" \
        "do_awk_intersect ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_all.bed ${BISCUIT_GENE} ${outdir}/${sample}_cpgs_gene_all_table.txt GenicCpGs All" \
        "do_awk_intersect ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_q40.bed ${BISCUIT_GENE} ${outdir}/${sample}_cpgs_gene_q40_table.txt GenicCpGs Q40" \
        "do_awk_intersect ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_all.bed ${BISCUIT_CGIS} ${outdir}/${sample}_cpgs_cgis_all_table.txt CGICpGs All" \
        "do_awk_intersect ${outdir} ${outdir}/${sample}_cpg_dist_table.txt ${outdir}/${sample}_cpg_q40.bed ${BISCUIT_CGIS} ${outdir}/${sample}_cpgs_cgis_q40_table.txt CGICpGs Q40"

    if [ `wc -l ${outdir}/${sample}_cpg_dist_table.txt | awk '{print $1}'` -eq "11" ]; then
        rm -f ${outdir}/${sample}_cpg_*.bed
    fi

    echo "Ended at `date`"
}

function usage {
    >&2 echo -e "\nUsage: cpg_counter.sh [-h,--help] [-o,--outdir] assets_directory sample_name in_bam\n"
    >&2 echo -e "Required inputs:"
    >&2 echo -e "\tassets_directory    : Path to assets directory"
    >&2 echo -e "\tsample_name         : Prefix of output files"
    >&2 echo -e "\tinput_bam           : Aligned BAM from BISCUIT\n"
    >&2 echo -e "Optional inputs:"
    >&2 echo -e "\t-h,--help           : Print help message and exit"
    >&2 echo -e "\t-o,--outdir         : Output directory [DEFAULT: cpg_counts]"
}

# Initialize default values for optional inputs
outdir="cpg_counts"

# Process command line arguments
OPTS=$(getopt \
    --options ho: \
    --long help,outdir: \
    --name "$(basename "$0")" \
    -- "$@"
)
eval set -- ${OPTS}

while true; do
    case "$1" in
        -h|--help )
            usage
            exit 0
            ;;
        -o|--outdir )
            outdir="$2"
            shift 2
            ;;
        -- )
            shift
            break
            ;;
        * )
            >&2 echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Make sure there are the correct number of inputs
if [[ $# -ne 3 ]]; then
    >&2 echo "$0: Missing inputs"
    usage
    exit 1
fi

# Fill required positional arguments
assets=$1
sample=$2
in_bam=$3

# Do some checks on the given inputs
if [[ ! -d "$assets" ]]; then
    >&2 echo "Assets directory missing: $assets"
    exit 1
fi

if [[ ! -f "${in_bam}" ]]; then
    >&2 echo "Cannot locate aligned BAM: ${in_bam}"
    >&2 echo "Please provide an existing aligned BAM."
    exit 1
fi

# Set variables for supplementary BED files
BISCUIT_CPGS="${assets}/cpg.bed.gz"
BISCUIT_CGIS="${assets}/cgi.bed.gz"
BISCUIT_RMSK="${assets}/rmsk.bed.gz"
BISCUIT_EXON="${assets}/exon.bed.gz"
BISCUIT_GENE="${assets}/genes.bed.gz"

>&2 echo "## Running CpG counter script with following configuration ##"
>&2 echo "=============="
>&2 echo "Sample Name        : ${sample}"
>&2 echo "Input BAM          : ${in_bam}"
>&2 echo "Output Directory   : ${outdir}"
>&2 echo "Assets Directory   : ${assets}"
>&2 echo "CPGS               : ${BISCUIT_CPGS}"
>&2 echo "CGIS               : ${BISCUIT_CGIS}"
>&2 echo "RMSK               : ${BISCUIT_RMSK}"
>&2 echo "EXON               : ${BISCUIT_EXON}"
>&2 echo "GENE               : ${BISCUIT_GENE}"
>&2 echo "=============="

count_cpgs
