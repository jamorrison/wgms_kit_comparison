"""Various utility functions for making statistical metric plots."""
import constants

def retrieve_single_element_data(data, key, rev):
    """Retrieve elements from JSON data structure that only have one item,
       i.e. json_dict['sample']['key'] = int (or float or string)

       Will sort data based on A_ORDER and B_ORDER dictionaries

       Inputs: data - JSON dictionary
               key  - key to search within each sample
               rev  - True/False for whether to reverse sort dictionaries

       Returns: tuple (A, B) - sorted list of tuples for each biological
                               replicate [tuple is (sample, key value)]
    """
    A = []
    B = []

    for samp, dic in data.items():
        if samp.startswith('FtubeA'):
            A.append( (samp, dic[key]) )
        else:
            B.append( (samp, dic[key]) )

    A = sorted(A, key=lambda d: constants.A_ORDER[d[0]], reverse=rev)
    B = sorted(B, key=lambda d: constants.B_ORDER[d[0]], reverse=rev)

    return A, B

def retrieve_single_element_data_from_dict(data, key1, key2, rev):
    """Retrieve elements from JSON data structure that exist in dictionary,
       i.e. json_dict['sample']['key'] = dict

       Will sort data based on A_ORDER and B_ORDER dictionaries

       Inputs: data - JSON dictionary
               key1 - key to search within each sample
               key2 - key to search within each key1 dictionary
               rev  - True/False for whether to reverse sort dictionaries

       Returns: tuple (A, B) - sorted list of tuples for each biological
                               replicate [tuple is (sample, key value)]
    """
    A = []
    B = []

    for samp, dic in data.items():
        if samp.startswith('FtubeA'):
            A.append( (samp, dic[key1][key2]) )
        else:
            B.append( (samp, dic[key1][key2]) )

    A = sorted(A, key=lambda d: constants.A_ORDER[d[0]], reverse=rev)
    B = sorted(B, key=lambda d: constants.B_ORDER[d[0]], reverse=rev)

    return A, B

def retrieve_single_element_data_one_output(data, key, rev):
    """Retrieve elements from JSON data structure that only have one item,
       i.e. json_dict['sample']['key'] = int (or float or string)

       Will sort data based on REPLICATE_ORDER dictionary

       Inputs: data - JSON dictionary
               key  - key to search within each sample
               rev  - True/False for whether to reverse sort dictionaries

       Returns: sorted list of tuples for each biological
                replicate [tuple is (sample, key value)]
    """
    C = []

    for samp, dic in data.items():
        C.append( (samp, dic[key]) )

    C = sorted(C, key=lambda d: constants.REPLICATE_ORDER[d[0]], reverse=rev)

    return C

def retrieve_single_element_data_from_dict_one_output(data, key1, key2, rev):
    """Retrieve elements from JSON data structure that exist in dictionary,
       i.e. json_dict['sample']['key'] = dict

       Will sort data based on REPLICATE_ORDER dictionary

       Inputs: data - JSON dictionary
               key1 - key to search within each sample
               key2 - key to search within each key1 dictionary
               rev  - True/False for whether to reverse sort dictionaries

       Returns: sorted list of tuples for each biological
                replicate [tuple is (sample, key value)]
    """
    C = []

    for samp, dic in data.items():
        C.append( (samp, dic[key1][key2]) )

    C = sorted(C, key=lambda d: constants.REPLICATE_ORDER[d[0]], reverse=rev)

    return C

def retrieve_aligned_reads_total_only(data, rev):
    """Retrieve total number of aligned reads.

    Will sort data based on REPLICATE_ORDER dictionary

    Inputs: data - JSON dictionary
            rev  - True/False for whether to reverse sort dictionaries

    Returns: sorted list of tuples for each biological replicate
             [tuple is (sample, total number of aligned reads)]
    """
    C = []

    for samp in data.keys():
        dic = data[samp]['aligned_reads']
        tot = sum([val for _, val in dic.items() if _ != 'mapq_percent'])
        C.append( (samp, tot) )

    C = sorted(C, key=lambda d: constants.REPLICATE_ORDER[d[0]], reverse=rev)

    return C

def retrieve_aligned_reads(data, rev):
    """Retrieve aligned_reads entries and calculate values used in plots.

    Will sort data based on A_ORDER and B_ORDER dictionaries

    Inputs: data - JSON dictionary
            rev  - True/False for whether to reverse sort dictionaries

    Returns: tuple (A, B) - sorted list of tuples for each biological
                            replicate
                            [tuple is (sample, [opt %, sub %, not %], total # reads)]
    """
    A = []
    B = []

    for samp in data.keys():
        if samp.startswith('FtubeA'):
            dic = data[samp]['aligned_reads']
            tot = sum([val for _, val in dic.items() if _ != 'mapq_percent'])
            dat = [
                100.0 * dic['opt_align'] / float(tot),
                100.0 * dic['sub_align'] / float(tot),
                100.0 * dic['not_align'] / float(tot)
            ]
            A.append( (samp, dat, tot) )
        else:
            dic = data[samp]['aligned_reads']
            tot = sum([val for _, val in dic.items() if _ != 'mapq_percent'])
            dat = [
                100.0 * dic['opt_align'] / float(tot),
                100.0 * dic['sub_align'] / float(tot),
                100.0 * dic['not_align'] / float(tot)
            ]
            B.append( (samp, dat, tot) )

    A = sorted(A, key=lambda d: constants.A_ORDER[d[0]], reverse=rev)
    B = sorted(B, key=lambda d: constants.B_ORDER[d[0]], reverse=rev)

    return A, B

def retrieve_data_points_from_dict(data, key1, rev):
    """Retrieve elements from JSON data structure that are a dictionary of data
       points, i.e. json_dict['sample']['key'] = {'pt1': 0, 'pt2': 1, ...}

       Will sort data based on A_ORDER and B_ORDER dictionaries

       Inputs: data - JSON dictionary
               key1 - key to search within each sample
               rev  - True/False for whether to reverse sort dictionaries

       Returns: tuple (A, B) - sorted list of ntuples for each biological
                               replicate
                               [tuple is (sample, xvals, yvals)]
    """
    A = []
    B = []

    for samp in data.keys():
        if samp.startswith('FtubeA'):
            dic = data[samp][key1]
            val = [(int(k), round(float(v), 2)) for k, v in dic.items()]
            x_vals = [i[0] for i in val]
            y_vals = [i[1] for i in val]
            A.append( (samp, x_vals, y_vals) )
        else:
            dic = data[samp][key1]
            val = [(int(k), round(float(v), 2)) for k, v in dic.items()]
            x_vals = [i[0] for i in val]
            y_vals = [i[1] for i in val]
            B.append( (samp, x_vals, y_vals) )

    A = sorted(A, key=lambda d: constants.A_ORDER[d[0]], reverse=rev)
    B = sorted(B, key=lambda d: constants.B_ORDER[d[0]], reverse=rev)

    return A, B

def retrieve_data_points_from_dict_in_dict(data, key1, key2, rev):
    """Retrieve elements from JSON data structure that are a dictionary of data
       points, i.e.
       json_dict['sample']['key'] = {'1': {'pt1': 0, 'pt2': 1, ...}, '2': {...}}

       Will sort data based on A_ORDER and B_ORDER dictionaries

       Inputs: data - JSON dictionary
               key1 - key to search within each sample
               key2 - key to search within each dictionary
               rev  - True/False for whether to reverse sort dictionaries

       Returns: tuple (A, B) - sorted list of ntuples for each biological
                               replicate
                               [tuple is (sample, xvals, yvals)]
    """
    A = []
    B = []

    for samp in data.keys():
        if samp.startswith('FtubeA'):
            dic = data[samp][key1][key2]
            val = [(int(k), round(float(v), 2)) for k, v in dic.items()]
            x_vals = [i[0] for i in val]
            y_vals = [i[1] for i in val]
            A.append( (samp, x_vals, y_vals) )
        else:
            dic = data[samp][key1][key2]
            val = [(int(k), round(float(v), 2)) for k, v in dic.items()]
            x_vals = [i[0] for i in val]
            y_vals = [i[1] for i in val]
            B.append( (samp, x_vals, y_vals) )

    A = sorted(A, key=lambda d: constants.A_ORDER[d[0]], reverse=rev)
    B = sorted(B, key=lambda d: constants.B_ORDER[d[0]], reverse=rev)

    return A, B

