"""Extract info from samtools stats output."""
import re
import os

def process_file(fname):
    """Parse samtools stats output.

    Inputs -- fname - filename of log file

    Returns -- tuple of sample, dictionary of samtools stats data
    """
    patterns = {
        'avg_insert': r'SN\s+insert size average\:\s+(\d*[.,]?\d*)'
    }

    data = {}
    with open(fname, 'r') as f:
        lines = f.read()

        for key, pat in patterns.items():
            m = re.search(pat, lines, re.MULTILINE)
            if m is not None:
                data[key] = float(m.group(1))

    sample = os.path.basename(fname).replace('.sorted.markdup.bam.stat', '')
    sample = sample.replace('.subsampled', '')

    return sample, data

if __name__ == '__main__':
    process_file('../align/FtubeAkapaBCrep2.sorted.markdup.bam.stat')
