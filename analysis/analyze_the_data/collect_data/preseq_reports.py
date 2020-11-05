"""Read data in Preseq complexity curve files."""
import os

def process_file(fname):
    """Process complexity curve file.

    Inputs -- fname - filename of complexity curve file

    Returns -- dictionary with complexity curve points
    """
    sample = os.path.basename(fname).replace('.complex.ccurve.txt','')

    data = {'complexity_curve': {}}
    with open(fname, 'r') as f:
        lines = f.read().splitlines()[1:]

        for l in lines:
            fields = l.strip().split('\t')

            data['complexity_curve'][int(fields[0])] = float(fields[1])

    return sample, data

if __name__ == '__main__':
    dirpath = '2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/preseq'
    process_file(dirpath + '/FtubeAkapaBC.complex.ccurve.txt')
