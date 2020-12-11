"""Parse various BISCUITqc files."""
import glob
import sys
import os
import re

def parse_logs_align_mapq(fname):
    """Parse _mapq_table.txt

    Inputs --  fname - filename of table

    Returns -- dictionary of aligned mapq data
    """
    output = {'opt_align': 0, 'sub_align': 0, 'not_align': 0, 'mapq_percent': {}}
    with open(fname, 'r') as f:
        file_data = f.read().splitlines()[2:]

        data = {}
        for l in file_data:
            s = l.split()
            data[s[0]] = s[1] # data[MAPQ] = number of reads

        total = sum([int(cnt) for _, cnt in data.items() if _ != 'unmapped'])
        counts = []
        for mapq in range(61):
            if str(mapq) in data.keys():
                counts.append(100.0 * float(data[str(mapq)]) / total)
            else:
                counts.append(0)
        output['mapq_percent'] = dict(zip(range(61), counts))

        for mapq, cnt in data.items():
            if mapq == 'unmapped':
                output['not_align'] += int(cnt)
            elif int(mapq) >= 40:
                output['opt_align'] += int(cnt)
            else:
                output['sub_align'] += int(cnt)

    return output

def parse_logs_align_isize(fname):
    """Parse _isize_table.txt

    Inputs: fname - filename of table

    Returns: data - dictionary of insert size data
    """
    data = {'percent': {}, 'readcnt': {}}

    with open(fname, 'r') as f:
        file_data = f.read().splitlines()[2:]
        for l in file_data:
            fields = l.split('\t')
            data['percent'][int(fields[0])] = 100.0 * float(fields[1])
            data['readcnt'][int(fields[0])] = float(fields[2])

    return data

def parse_logs_dup_report(fname):
    """Parses _dup_report.txt

    Inputs -- fname - filename of table

    Returns -- dictionary of duplicate fractions
    """
    patterns = [
        (
            r'Number of duplicate reads:\s+(\d+)',
            r'Number of reads:\s+(\d+)',
            'dup_all'
        ),
        (
            r'Number of duplicate q40-reads:\s+(\d+)',
            r'Number of q40-reads:\s+(\d+)',
            'dup_q40'
        )
    ]

    output = {}
    with open(fname, 'r') as f:
        file_data = f.read()
        for pat_dup, pat_tot, key in patterns:
            m1 = re.search(pat_dup, file_data, re.MULTILINE)
            m2 = re.search(pat_tot, file_data, re.MULTILINE)
            if m1 is not None and m2 is not None:
                output[key] = (100.0 * float(m1.group(1)) / float(m2.group(1)))
            else:
                print('dup data is missing ... :(')
                print(fname)
                sys.exit(1)

    return output

def parse_logs_qc_cv(fname):
    """Parses _cv_table.txt

    Inputs -- fname - filename of table

    Returns -- dictionary of depth uniformity measures
    """
    output = {}
    targets = ['all_base', 'all_cpg', 'q40_base', 'q40_cpg']
    with open(fname, 'r') as f:
        file_data = f.read()
        for t in targets:
            m = re.search(
                '{}\t([\d\.]+)\t([\d\.]+)\t([\d\.]+)'.format(t),
                file_data, re.MULTILINE
            )
            if m is not None:
                output[t+'_mu'] = float(m.group(1))
                output[t+'_cv'] = float(m.group(3))
            else:
                print('covg data may be missing ... :(')
                print('More likely it is a number in scientific notation')
                print(fname)
                output[t+'_mu'] = 'NA'
                output[t+'_cv'] = 'NA'

    return output

def parse_logs_covdist_all_base(fname):
    """Parses _covdist_all_base_table.txt
              _covdist_all_cpg_table.txt
              _covdist_q40_base_table.txt
              _covdist_q40_cpg_table.txt

    Inputs: fname - filename of table

    Returns: data - dictionary of coverage distributions up to 30X data
    """
    dd = {}

    with open(fname, 'r') as f:
        file_data = f.read().splitlines()[2:]
        for l in file_data:
            fields = l.split()
            dd[int(float(fields[0]))] = int(float(fields[1]))

    covs = sorted([k for k in dd])[:51]
    _ccov_cnt = sum(dd.values())

    ccov_cnts = []
    for cov in covs:
        ccov_cnts.append(_ccov_cnt/1000000.0)
        _ccov_cnt -= dd[cov]

    return dict(zip(covs, ccov_cnts))

def parse_logs_read_avg_retention_rate(fname):
    """ Parses _totalReadConversionRate.txt

    Inputs: fname - filename of table

    Returns: output - dictionary of read averaged fraction of retainied cytosines by context
    """
    output = {}
    with open(fname, 'r') as f:
        try:
            m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.read().splitlines()[2])
        except IndexError:
            print('read avg data is missing :(')
            print(fname)
            sys.exit(1)
        else:
            if m is not None:
                output['rca'] = 100.0 * float(m.group(1))
                output['rcc'] = 100.0 * float(m.group(2))
                output['rcg'] = 100.0 * float(m.group(3))
                output['rct'] = 100.0 * float(m.group(4))

    return output

def parse_logs_base_avg_retention_rate(fname):
    """Parses _totalBaseConversionRate.txt

    Inputs -- fname - filename of table

    Returns -- dictionary of base averaged fraction of retainied cytosines by context
    """
    output = {}
    with open(fname, 'r') as f:
        try:
            m = re.match(r'([\d.]+)\t([\d.]+)\t([\d.]+)\t([\d.]+)', f.read().splitlines()[2])
        except IndexError:
            print('base avg data is missing :(')
            print(fname)
            sys.exit(1)
        else:
            if m is not None:
                output['bca'] = 100.0 * float(m.group(1))
                output['bcc'] = 100.0 * float(m.group(2))
                output['bcg'] = 100.0 * float(m.group(3))
                output['bct'] = 100.0 * float(m.group(4))

    return output

def parse_logs_cpg_retention_readpos(fname):
    """ Parses _CpGRetentionByReadPos.txt _CpHRetentionByReadPos.txt

    Inputs -- fname - filename of table

    Returns -- dictionary of fraction of retained cytosines for reads 1 and 2
        in either a CpH or CpG context
    """
    r1 = {'C': {}, 'R': {}}
    r2 = {'C': {}, 'R': {}}
    with open(fname, 'r') as f:
        file_data = f.read().splitlines()[2:]
        for l in file_data:
            fields = l.strip().split('\t')

            if fields[0] not in ['1', '2'] or fields[2] not in ['C', 'R']:
                return {}
            if fields[0] == '1':
                r1[fields[2]][int(fields[1])] = int(fields[3])
            elif fields[0] == '2':
                r2[fields[2]][int(fields[1])] = int(fields[3])

    r1rate = {}
    for k in sorted(r1['C'].keys()):
        if k in r1['R']:
            r1rate[k] = 100.0 * float(r1['R'][k]) / (r1['R'][k] + r1['C'][k])

    r2rate = {}
    for k in sorted(r2['C'].keys()):
        if k in r2['R']:
            r2rate[k] = 100.0 * float(r2['R'][k]) / (r2['R'][k] + r2['C'][k])

    return {'1': r1rate, '2': r2rate}

def extract_file_type(fname, dirpath, sample):
    """Extract what BISCUITqc file is being checked

    Inputs -- fname - filename to extract type from
              dirpath - Path to BISCUITqc directory
              sample - Sample name to prepend to BISCUITqc files

    Returns -- integer of file type
    """
    used_files = {
        '_cv_table.txt': 1,
        '_dup_report.txt': 2,
        '_totalBaseConversionRate.txt': 3,
        '_totalReadConversionRate.txt': 4,
        '_mapq_table.txt': 5,
        '_CpGRetentionByReadPos.txt': 6,
        '_CpHRetentionByReadPos.txt': 7,
        '_covdist_all_base_table.txt': 8,
        '_covdist_all_cpg_table.txt': 9,
        '_covdist_q40_base_table.txt': 10,
        '_covdist_q40_cpg_table.txt': 11,
        '_isize_table.txt': 12
    }

    file_type = fname.replace(dirpath, '').replace('/' + sample, '')
    if file_type in used_files.keys():
        return used_files[file_type]
    else:
        return -1


def process_dir(dirpath, sample):
    """Process BISUITqc for sample living in dirpath.

    Inputs -- dirpath - Path to BISCUITqc directory
              sample - Sample name to prepend to BISCUITqc files

    Returns -- dictionary of QC data
    """
    files = glob.glob(dirpath + '/' + sample + '_*.txt')

    data = {}
    for f in files:
        file_type = extract_file_type(f, dirpath, sample)
        if file_type == 1:
            data['uniformity'] = parse_logs_qc_cv(f)
        elif file_type == 2:
            data['dup_report'] = parse_logs_dup_report(f)
        elif file_type == 3:
            data['base_rtn'] = parse_logs_base_avg_retention_rate(f)
        elif file_type == 4:
            data['read_rtn'] = parse_logs_read_avg_retention_rate(f)
        elif file_type == 5:
            data['aligned_reads'] = parse_logs_align_mapq(f)
        elif file_type == 6:
            data['cpg_rtn_readpos'] = parse_logs_cpg_retention_readpos(f)
        elif file_type == 7:
            data['cph_rtn_readpos'] = parse_logs_cpg_retention_readpos(f)
        elif file_type == 8:
            data['covdist_all_base'] = parse_logs_covdist_all_base(f)
        elif file_type == 9:
            data['covdist_all_cpg'] = parse_logs_covdist_all_base(f)
        elif file_type == 10:
            data['covdist_q40_base'] = parse_logs_covdist_all_base(f)
        elif file_type == 11:
            data['covdist_q40_cpg'] = parse_logs_covdist_all_base(f)
        elif file_type == 12:
            data['isize_data'] = parse_logs_align_isize(f)

    return data

if __name__ == '__main__':
    dirpath = '2019_11_07_FallopianTube_WGBS_Kit_Comparison/analysis/align/FtubeAkapaBC_QC'
    sample = 'FtubeAkapaBC'
    process_dir(dirpath, sample)
