"""Read data in trinucleotide context files."""
import os

def process_file(fname, ext):
    """Process trinucleotide context files.

    Inputs -- fname - filename of trinucleotide context file
              ext   - tag to remove from end of file

    Returns -- tuple (sample, dictionary with trinucleotide methylation values)
    """
    sample = os.path.basename(fname).replace(ext, '')

    data = {
        'cah_methylation_percent': -1,
        'cag_methylation_percent': -1,
        'cth_methylation_percent': -1,
        'ctg_methylation_percent': -1
    }
    with open(fname) as f:
        count = len(f.readlines())
        if count != 1:
            print(f'[trinuc_meth::process_file] ERROR: Incorrect number of lines in {fname}')

        f.seek(0,0)
        d = f.readline().strip().split('\t')
        if len(d) != 4:
            print(f'[trinuc_meth::process_file] ERROR: Incorrect number of methylation values in {fname}')

        data['cah_methylation_percent'] = float(d[0])
        data['cag_methylation_percent'] = float(d[1])
        data['cth_methylation_percent'] = float(d[2])
        data['ctg_methylation_percent'] = float(d[3])

    return sample, data

if __name__ == '__main__':
    samp, data = process_file('../cpg_questions/trinuc_methylation/FtubeAkapaBC.tsv')
    print(samp, data)

