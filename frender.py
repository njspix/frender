'''
Fastq REad Name DemultiplexER

TODO: ADD ONE SENTENCE EXPLANATION

Creator:
    Nathan Spix

Date Created:
    September 2020

Requires:
    pandas 1.2.2
    regex 2020.11.13

Inputs:
    TODO: ADD DESCRIPTIONS OF INPUTS

OUTPUT/RETURNS:
    TODO: ADD DESCRIPTION OF OUTPUT FILES
'''
import argparse
from pathlib import Path
import pandas as pd
import gzip
from itertools import zip_longest
import sys

try:
    import regex
except ModuleNotFoundError:
    print('ERROR: Install regex module')
    print('\tpip install regex\tOR\tconda install -c conda-forge regex')
    sys.exit(0)

def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks.
       Ex. grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"

    Inputs - 
        iterable  - string to group in blocks of length n
        n         - block size
        fillvalue - if block sizes of uneven length, fill missing values with
                    fillvalue
    Returns -
        iter object with fixed-length blocks
    """
    args = [iter(iterable)] * n

    return zip_longest(*args, fillvalue=fillvalue)

def fuz_match_list(pattern, set_of_strings):
    """Given a query string and a list of strings, return a list of indices
       where a match (+/- 1 substitution) is found

    Inputs -
        pattern        - pattern to search for
        set_of_strings - strings to search for pattern in
    Returns -
        list of where matches (with +/- 1 substitution) occur in set_of_strings
    """
    pattern = regex.compile("(?:"+pattern+"){s<=1}", regex.IGNORECASE)
    matches = [bool(regex.match(pattern, string)) for string in set_of_strings]
    
    return [i for i, val in enumerate(matches) if val]

def frender(barcode, fastq_1, fastq_2,
            out_dir='.', preefix='', ihopped='',
            cnflict='', undeter=''):
    """Main module function.

    Inputs -
        barcode - CSV file containing barcodes
        fastq_1 - Read 1 FASTQ file
        fastq_2 - Read 2 FASTQ file
        out_dir - Output directory name [default: '.']
        preefix - Prefix to add to output files [default: '']
        ihopped - Name of FASTQ file to save index hopped reads to
                  [default: '']
        cnflict - Name of FASTQ file to save conflicting reads to
                  [default: '']
        undeter - Name of FASTQ file to save undetermined reads to
                  [default: '']
    Returns -
        Creates demultiplexed FASTQ files from fastq_1 and fastq_2 and a summary
        of the index hopped reads. Will also create index hopped FASTQ file,
        conflicting barcode FASTQ file, and undetermined FASTQ file if supplied
        during function call.
    """
    # Not all data includes a header in barcodeAssociationTable.txt
    # Check for how to correctly load data
    cols=['AccessionID', 'ClientAccessionID', 'Index1', 'Index2']
    test = pd.read_csv(barcode, nrows=1, header=None)

    # As of 15 Sep 2020, early scWGBS data has four columns and newer data
    # has seven columns. There is one dataset with three columns, and this
    # is related to a collaborative project with human fetal DNA.
    if (len(test.columns) == 3):
        print('WARNING: This barcode file only has 3 columns.', end=' ')
        print('It likely is a single-indexed library')
        sys.exit(0)

    #TODO: Handle single index
    #TODO: Check barcodes are all the same length, match [ATCGatcg]
    if (test[0][0].startswith('Acc') or test[0][0].startswith('Client')):
        indexes = pd.read_csv(barcode, usecols=cols, index_col='AccessionID')
    else:
        indexes = pd.read_csv(barcode, header=None, names=cols,
                              index_col='AccessionID')

    all_idx1 = indexes['Index1'].tolist()
    all_idx2 = indexes['Index2'].tolist()
    ids = list(indexes.index)

    # Create output directory and prep prefix
    out_dir = out_dir.rstrip('/') + '/' # Add trailing slash
    if out_dir != './':
        Path(out_dir).mkdir(parents=True, exist_ok=True)

    if (preefix != '' and not preefix.endswith('_')):
        preefix = preefix + '_'

    # Open output files
    r1_files = []
    r2_files = []

    for id in ids:
        r1_files.append(open(f'{out_dir}{preefix}{id}_R1.fastq', 'w'))
        r2_files.append(open(f'{out_dir}{preefix}{id}_R2.fastq', 'w'))
        #r1_files.append(f'{out_dir}{preefix}{id}_R1.fastq')
        #r2_files.append(f'{out_dir}{preefix}{id}_R2.fastq')

    if (ihopped != ''):
        r1_hop = open(f'{out_dir}{preefix}{ihopped}_R1.fastq', 'w')
        r2_hop = open(f'{out_dir}{preefix}{ihopped}_R2.fastq', 'w')
    if (cnflict != ''):
        r1_con = open(f'{out_dir}{preefix}{cnflict}_R1.fastq', 'w')
        r2_con = open(f'{out_dir}{preefix}{cnflict}_R2.fastq', 'w')
    if (undeter != ''):
        r1_und = open(f'{out_dir}{preefix}{undeter}_R1.fastq', 'w')
        r2_und = open(f'{out_dir}{preefix}{undeter}_R2.fastq', 'w')

    # DataFrame to store index hopping information
    hops = pd.DataFrame()

    # Start processing data
    with gzip.open(fastq_1, 'rt') as read1, gzip.open(fastq_2, 'rt') as read2:
        reads_1 = grouper(read1, 4, '')
        reads_2 = grouper(read2, 4, '')

        barcode_dict = {}

        record_count = 0
        for record_1, record_2 in zip(reads_1, reads_2):
            assert (len(record_1) == 4) and (len(record_2) == 4)
            if (record_count%10000 == 1):
                print(f"Processed {record_count} reads")
            record_count += 1

            # Parse record header line to compare name and extract barcode
            r1_head = record_1[0].rstrip('\n').split(' ')
            r2_head = record_2[0].rstrip('\n').split(' ')

            if (r1_head[0] != r2_head[0]):
                print('WARNING: Skipping mismatched reads', end=' ')
                print(f'{r1_head[0]} != {r2_head[0]} at line', end=' ')
                print(f'{4*record_count} of {fastq_1} and {fastq_2}')
                continue
            
            # Barcode indices
            code = r1_head[1].split(':')[-1] # full index name
            idx1 = code.split('+')[0]        # index 1 in full name
            idx2 = code.split('+')[1]        # index 2 in full name

            if code in barcode_dict: # seen this combination before

                # could be an index hop:
                if (barcode_dict[code] == 'hop'):

                    # write to files if requested
                    if (ihopped != ''):
                        for line in record_1:
                            r1_hop.write(str(line))
                        for line in record_2:
                            r2_hop.write(str(line))

                    # update hop count table
                    idx1_matches = fuz_match_list(idx1, all_idx1)
                    idx2_matches = fuz_match_list(idx2, all_idx2)
                    i_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                    i_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                    hops.loc[i_idx1, i_idx2] += 1 

                # could be a conflict:
                elif (barcode_dict[code] == 'conflict') and (cnflict != ''):
                    for line in record_1:
                        r1_con.write(str(line))
                    for line in record_2:
                        r2_con.write(str(line))

                # could be unkonwn:
                elif (barcode_dict[code] == 'undetermined') and (undeter != ''):
                    for line in record_1:
                        r1_und.write(str(line))
                    for line in record_2:
                        r2_und.write(str(line))
                
                # or it could be a valid demux:
                elif barcode_dict[code] not in ['hop', 'conflict', 'undetermined']:
                    for line in record_1:
                        r1_files[barcode_dict[code]].write(str(line))
                    for line in record_2:
                        r2_files[barcode_dict[code]].write(str(line))

            else: # need to process this read
                idx1_matches = fuz_match_list(idx1, all_idx1)
                idx2_matches = fuz_match_list(idx2, all_idx2)

                if (bool(idx1_matches) and bool(idx2_matches)):
                    # Can find at least one barcode match for both indices
                    match_isec = set(idx1_matches).intersection(idx2_matches)

                    if (len(match_isec) == 0):
                        # No matches, so index hop
                        barcode_dict[code] = 'hop' # add to known barcodes

                        # This is currently imprecise. A 'hop' could match more than one barcode in either index, 
                        # idx1_matches = [4] and idx2_matches = [2, 3]. In this case, the matched indices could 
                        # be identical, or they could both be 1 substitution away from the target (bad index design, but that's not my problem here :-)
                        # This doesn't account for that second case.
                        i_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                        i_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                        try:
                            hops.loc[i_idx1, i_idx2] += 1 
                        except KeyError:
                            # Combination of indices hasn't been initialized
                            hops.loc[i_idx1, i_idx2] = 1  

                        # Write to output file (if applicable)
                        if (ihopped != ''):
                            for line in record_1:
                                r1_hop.write(str(line))
                            for line in record_2:
                                r2_hop.write(str(line))
                    elif (len(match_isec) == 1):
                        # good read; idx1 and idx2 line up in exactly one spot
                        demux_id = indexes.index[match_isec.pop()]

                        barcode_dict[code] = indexes.index.get_loc(demux_id)
                        for line in record_1:
                            r1_files[barcode_dict[code]].write(str(line))
                        for line in record_2:
                            r2_files[barcode_dict[code]].write(str(line))
                    else:
                        # Read matches to more than one possible output file
                        barcode_dict[code] = 'conflict'
                        if (cnflict != ''):
                            for line in record_1:
                                r1_con.write(str(line))
                            for line in record_2:
                                r2_con.write(str(line))
                else:
                    # Can't find a match for at least one barcode
                    barcode_dict[code] = 'undetermined'
                    if (undeter != ''):
                        for line in record_1:
                            r1_und.write(str(line))
                        for line in record_2:
                            r2_und.write(str(line))

    # Close output files
    for i in range(len(r1_files)):
        r1_files[i].close()
        r2_files[i].close()

    if (ihopped != ''):
        r1_hop.close()
        r2_hop.close()
    if (cnflict != ''):
        r1_con.close()
        r2_con.close()
    if (undeter != ''):
        r1_und.close()
        r2_und.close()

    # Write index hopping report
    hops = hops.reset_index().melt(id_vars='index',
                                   var_name='idx2',
                                   value_name= 'num_hops_observed')
    hops = hops[hops['num_hops_observed'] > 0]
    hops.columns = ['idx1', 'idx2', 'num_hops_observed']
    hops.to_csv(f'{out_dir}{preefix}barcode_hops.csv', index=False)

def frender_se(barcode, fastq,
            out_dir='.', preefix='', ihopped='',
            cnflict='', undeter=''):
    """Demultiplexes a single-end read FASTQ file based on barcodes in read name.

    Inputs -
        barcode - CSV file containing barcodes
        fastq - FASTQ file
        out_dir - Output directory name [default: '.']
        preefix - Prefix to add to output files [default: '']
        ihopped - Name of FASTQ file to save index hopped reads to
                  [default: '']
        cnflict - Name of FASTQ file to save conflicting reads to
                  [default: '']
        undeter - Name of FASTQ file to save undetermined reads to
                  [default: '']
    Returns -
        Creates demultiplexed FASTQ files and a summary
        of the index hopped reads. Will also create index hopped FASTQ file,
        conflicting barcode FASTQ file, and undetermined FASTQ file if supplied
        during function call.
    """
    # Not all data includes a header in barcodeAssociationTable.txt
    # Check for how to correctly load data
    cols=['AccessionID', 'ClientAccessionID', 'Index1', 'Index2']
    test = pd.read_csv(barcode, nrows=1, header=None)

    # As of 15 Sep 2020, early scWGBS data has four columns and newer data
    # has seven columns. There is one dataset with three columns, and this
    # is related to a collaborative project with human fetal DNA.
    if (len(test.columns) == 3):
        print('WARNING: This barcode file only has 3 columns.', end=' ')
        print('It likely is a single-indexed library')
        sys.exit(0)

    #TODO: Handle single index
    #TODO: Check barcodes are all the same length, match [ATCGatcg]
    if (test[0][0].startswith('Acc') or test[0][0].startswith('Client')):
        indexes = pd.read_csv(barcode, usecols=cols, index_col='AccessionID')
    else:
        indexes = pd.read_csv(barcode, header=None, names=cols,
                              index_col='AccessionID')

    all_idx1 = indexes['Index1'].tolist()
    all_idx2 = indexes['Index2'].tolist()
    ids = list(indexes.index)

    # Create output directory and prep prefix
    out_dir = out_dir.rstrip('/') + '/' # Add trailing slash
    if out_dir != './':
        Path(out_dir).mkdir(parents=True, exist_ok=True)

    if (preefix != '' and not preefix.endswith('_')):
        preefix = preefix + '_'

    # Open output file
    r1_files = []

    for id in ids:
        r1_files.append(open(f'{out_dir}{preefix}{id}_R1.fastq', 'w'))
        #r1_files.append(f'{out_dir}{preefix}{id}_R1.fastq')

    if (ihopped != ''):
        r1_hop = open(f'{out_dir}{preefix}{ihopped}_R1.fastq', 'w')
    if (cnflict != ''):
        r1_con = open(f'{out_dir}{preefix}{cnflict}_R1.fastq', 'w')
    if (undeter != ''):
        r1_und = open(f'{out_dir}{preefix}{undeter}_R1.fastq', 'w')

    # DataFrame to store index hopping information
    hops = pd.DataFrame()

    # Start processing data
    with gzip.open(fastq, 'rt') as read1:
        reads = grouper(read1, 4, '')

        barcode_dict = {}

        record_count = 0
        for record in reads:
            assert (len(record) == 4)
            if (record_count%10000 == 1):
                print(f"Processed {record_count} reads")
            record_count += 1

            # Parse record header line to compare name and extract barcode
            r1_head = record[0].rstrip('\n').split(' ')
            
            # Barcode indices
            code = r1_head[1].split(':')[-1] # full index name
            idx1 = code.split('+')[0]        # index 1 in full name
            idx2 = code.split('+')[1]        # index 2 in full name

            if code in barcode_dict: # seen this combination before

                # could be an index hop:
                if (barcode_dict[code] == 'hop'):

                    # write to files if requested
                    if (ihopped != ''):
                        for line in record:
                            r1_hop.write(str(line))

                    # update hop count table
                    idx1_matches = fuz_match_list(idx1, all_idx1)
                    idx2_matches = fuz_match_list(idx2, all_idx2)
                    i_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                    i_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                    hops.loc[i_idx1, i_idx2] += 1 

                # could be a conflict:
                elif (barcode_dict[code] == 'conflict') and (cnflict != ''):
                    for line in record:
                        r1_con.write(str(line))

                # could be unkonwn:
                elif (barcode_dict[code] == 'undetermined') and (undeter != ''):
                    for line in record:
                        r1_und.write(str(line))
                
                # or it could be a valid demux:
                elif barcode_dict[code] not in ['hop', 'conflict', 'undetermined']:
                    for line in record:
                        r1_files[barcode_dict[code]].write(str(line))

            else: # need to process this read
                idx1_matches = fuz_match_list(idx1, all_idx1)
                idx2_matches = fuz_match_list(idx2, all_idx2)

                if (bool(idx1_matches) and bool(idx2_matches)):
                    # Can find at least one barcode match for both indices
                    match_isec = set(idx1_matches).intersection(idx2_matches)

                    if (len(match_isec) == 0):
                        # No matches, so index hop
                        barcode_dict[code] = 'hop' # add to known barcodes

                        # This is currently imprecise. A 'hop' could match more than one barcode in either index, 
                        # idx1_matches = [4] and idx2_matches = [2, 3]. In this case, the matched indices could 
                        # be identical, or they could both be 1 substitution away from the target (bad index design, but that's not my problem here :-)
                        # This doesn't account for that second case.
                        i_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                        i_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                        try:
                            hops.loc[i_idx1, i_idx2] += 1 
                        except KeyError:
                            # Combination of indices hasn't been initialized
                            hops.loc[i_idx1, i_idx2] = 1  

                        # Write to output file (if applicable)
                        if (ihopped != ''):
                            for line in record:
                                r1_hop.write(str(line))
                    elif (len(match_isec) == 1):
                        # good read; idx1 and idx2 line up in exactly one spot
                        demux_id = indexes.index[match_isec.pop()]

                        barcode_dict[code] = indexes.index.get_loc(demux_id)
                        for line in record:
                            r1_files[barcode_dict[code]].write(str(line))
                    else:
                        # Read matches to more than one possible output file
                        barcode_dict[code] = 'conflict'
                        if (cnflict != ''):
                            for line in record:
                                r1_con.write(str(line))
                else:
                    # Can't find a match for at least one barcode
                    barcode_dict[code] = 'undetermined'
                    if (undeter != ''):
                        for line in record:
                            r1_und.write(str(line))

    # Close output files
    for i in range(len(r1_files)):
        r1_files[i].close()

    if (ihopped != ''):
        r1_hop.close()
    if (cnflict != ''):
        r1_con.close()
    if (undeter != ''):
        r1_und.close()

    # Write index hopping report
    hops = hops.reset_index().melt(id_vars='index',
                                   var_name='idx2',
                                   value_name= 'num_hops_observed')
    hops = hops[hops['num_hops_observed'] > 0]
    hops.columns = ['idx1', 'idx2', 'num_hops_observed']
    hops.to_csv(f'{out_dir}{preefix}barcode_hops.csv', index=False)

def frender_scan(barcode, fastq_1,
            out_dir='.', preefix=''):
    """Scan a single fastq file, counting exact and inexact barcode matches, conflicting barcodes, index hops, and undetermined reads. No demultiplexing is performed.

    Inputs -
        barcode - CSV file containing barcodes
        fastq_1 - Read 1 FASTQ file
        out_dir - Output directory name [default: '.']
        preefix - Prefix to add to output files [default: '']

    Returns -
        TODO: doc...
    """
    # Not all data includes a header in barcodeAssociationTable.txt
    # Check for how to correctly load data
    cols=['AccessionID', 'ClientAccessionID', 'Index1', 'Index2']
    test = pd.read_csv(barcode, nrows=1, header=None)

    # As of 15 Sep 2020, early scWGBS data has four columns and newer data
    # has seven columns. There is one dataset with three columns, and this
    # is related to a collaborative project with human fetal DNA.
    if (len(test.columns) == 3):
        print('WARNING: This barcode file only has 3 columns.', end=' ', file=sys.stderr)
        print('It likely is a single-indexed library', file=sys.stderr)
        sys.exit(0)

    #TODO: Handle single index
    #TODO: Check barcodes are all the same length, match [ATCGatcg]
    if (test[0][0].startswith('Acc') or test[0][0].startswith('Client')):
        indexes = pd.read_csv(barcode, usecols=cols, index_col='AccessionID')
    else:
        indexes = pd.read_csv(barcode, header=None, names=cols,
                              index_col='AccessionID')

    all_idx1 = indexes['Index1'].tolist()
    all_idx2 = indexes['Index2'].tolist()
    ids = list(indexes.index)

    # Create output directory and prep prefix
    out_dir = out_dir.rstrip('/') + '/' # Add trailing slash
    if out_dir != './':
        Path(out_dir).mkdir(parents=True, exist_ok=True)

    if (preefix != '' and not preefix.endswith('_')):
        preefix = preefix + '_'

    test = pd.DataFrame()
    record_count = 0

    with gzip.open('test.fastq.gz', 'rt') as read_file:
        for read_head in islice(read_file, 0, None, 4):
            if record_count%10000==1:
                print(f"processed {record_count}")
            
            code = read_head.rstrip('\n').split(' ')[1].split(':')[-1]
            idx1 = code.split('+')[0]
            idx2 = code.split('+')[1]

            try:
                test.loc[idx1, idx2] += 1
            except KeyError:
                test.loc[idx1, idx2] = 1
            record_count += 1

    test2 = test.stack().to_frame("total_reads").reindex(columns = ['total_reads', 'matched_idx1', 'matched_idx2','read_type', 'sample_name'])
    test2.index = test2.index.rename(['idx1','idx2'])
    array = test2.index.values

    for i in array:
        idx1 = i[0]
        idx2 = i[1]

        idx1_matches = fuz_match_list(idx1, all_idx1)
        idx2_matches = fuz_match_list(idx2, all_idx2)

        if (bool(idx1_matches) and bool(idx2_matches)):
            # Can find at least one barcode match for both indices
            match_isec = set(idx1_matches).intersection(idx2_matches)
            
            if (len(match_isec) == 0):
                #this is an index hop
                matched_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                matched_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                
                test2.loc[(idx1, idx2), 'matched_idx1'] = matched_idx1
                test2.loc[(idx1, idx2), 'matched_idx2'] = matched_idx2
                test2.loc[(idx1, idx2), 'read_type'] = 'index_hop'
                test2.loc[(idx1, idx2), 'sample_name'] = ''

            elif (len(match_isec) == 1):
                # this is a good read
                matched_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                matched_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                sample_name = indexes.index[match_isec.pop()]
                
                test2.loc[(idx1, idx2), 'matched_idx1'] = matched_idx1
                test2.loc[(idx1, idx2), 'matched_idx2'] = matched_idx2
                test2.loc[(idx1, idx2), 'read_type'] = 'demuxable'
                test2.loc[(idx1, idx2), 'sample_name'] = sample_name
            
            else: 
                # this is an ambiguous read
                matched_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                matched_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                
                test2.loc[(idx1, idx2), 'matched_idx1'] = matched_idx1
                test2.loc[(idx1, idx2), 'matched_idx2'] = matched_idx2
                test2.loc[(idx1, idx2), 'read_type'] = 'ambiguous'
                test2.loc[(idx1, idx2), 'sample_name'] = ''
        else:
            # this is an undetermined read
            test2.loc[(idx1, idx2), 'matched_idx1'] = ''
            test2.loc[(idx1, idx2), 'matched_idx2'] = ''
            test2.loc[(idx1, idx2), 'read_type'] = 'undetermined'
            test2.loc[(idx1, idx2), 'sample_name'] = ''

    # Write report
    print(test2.tocsv())

    #TODO: write out all information in csv format, e.g.:
    # idx_1,idx_2,class,sample_name,total_reads
    # ACTGA,CGATG,{index_hop,ambiguous,undetermined,demuxable},{sample_name,NA},131245125

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Demultiplex Undetermined FastQ files with given barcodes.'
    )
    parser.add_argument(
        '-o', metavar='output_dir',
        default='.',
        help="output files directory [default: '.']"
    )
    parser.add_argument(
        '-p', metavar='file_prefix',
        default='' ,
        help='prefix added to output files [default: '']')
    parser.add_argument(
        '-i', metavar='index_hopped.fq',
        default='',
        help='output fastq for index hopped reads [default: '']'
    )
    parser.add_argument(
        '-c', metavar='conflicting.fq',
        default='',
        help='output fastq for reads with conflicting barcodes [default: '']'
    )
    parser.add_argument(
        '-u', metavar='undetermined.fq',
        default='',
        help='output fastq for remaining undetermined reads [default: '']'
    )
    parser.add_argument(
        '-b', metavar='barcode.txt',
        help='Barcode association table, csv format',
        required=True
    )
    parser.add_argument(
        '-s', '--scan', 
        action = 'store_true',
        help='Scan a single fastq file, counting exact and inexact barcode matches, conflicting barcodes, index hops, and undetermined reads. No demultiplexing is performed.',
    )

    parser.add_argument('fastqs', nargs=argparse.REMAINDER) 

    args = parser.parse_args()
    
    if scan:
        print(f"Scanning {fastqs[0]}...", file=sys.stderr)
        frender_scan(
                args.b,
                args.fastqs[0],
                args.o,
                args.p
            )

    else:
        #single end case
        if ( len(args.fastqs) == 1 ):
            print("Only one input fastq detected. Running in single-end mode...")
            frender_se(
                args.b,
                args.fastqs[0],
                args.o,
                args.p,
                args.i,
                args.c,
                args.u
            )

        #paired-end case
        elif ( len(args.fastqs) == 2 ):
            print("Two input fastqs detected. Running in paired-end mode...")
            frender(
                args.b,
                args.fastqs[0],
                args.fastqs[1],
                args.o,
                args.p,
                args.i,
                args.c,
                args.u
            )
        else: raise TypeError("Wrong number of fastq files provided. Please specify only one or two fastq files to use as input!")

#TODO: add 'scan' mode. Rip through read1 Undetermined_* file (doesn't matter if single or paired-end...)
# Classify and count all the barcodes present, then dump into csv or some other report. 
# Have an option to proceed with standard frender() if any un-demultiplexed reads are found.