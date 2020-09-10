import argparse
import os
import pandas as pd
import gzip
from itertools import zip_longest
import regex

parser = argparse.ArgumentParser()
parser.add_argument("-c", help = "Barcode association table, csv format", required = True)
parser.add_argument("-r1", help = "Undetermined read 1", required = True)
#todo: make this optional
parser.add_argument("-r2", help = "Undetermined read 2", required = True)
parser.add_argument("-o", help = "output files directory", default = './split_files/')

args = parser.parse_args()

# create output dir
try:
    os.makedirs(args.o)
except OSError as e:
    if e.errno != 17: #17 = file exists
        raise

#read barcode association file:
indexes = pd.read_csv(args.c, header=0, names=['id', 'idx1', 'idx2']).set_index('id')
all_idx1 = indexes['idx1'].tolist() 
all_idx2 = indexes['idx2'].tolist()
#todo: handle single index
#todo: check that barcodes are all same length, match [ACTGactg]

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def fuz_match_list(pattern, set_of_strings):
    "Given a query string and a list of strings, return a list of indices where a match (+/- 1 substitution) is found"
    pattern = regex.compile("(?:"+pattern+"){s<=1}", regex.IGNORECASE)
    matches = [bool(regex.match(pattern, set_of_strings[i])) for i in range (len(set_of_strings))]
    matches = [i for i, val in enumerate(matches) if val] 
    return matches

hops = pd.DataFrame()

with gzip.open(args.r1, 'rt') as read1:
    for record in grouper(read1, 4, ''):
        assert len(record) == 4
        idx1 = record[0].split(":")[-1].split("+")[0].rstrip('\n')
        idx2 = record[0].split(":")[-1].split("+")[1].rstrip('\n')
        idx1_matches = fuz_match_list(idx1, all_idx1)
        idx2_matches = fuz_match_list(idx2, all_idx2) 
        if (bool(idx1_matches) & bool(idx2_matches)):
            #print(idx1, idx2)
            if len(set(idx1_matches).intersection(idx2_matches)) == 0: # index hop
                try:
                    hops.loc[set([all_idx1[i] for i in idx1_matches]).pop(),set([all_idx2[i] for i in idx2_matches]).pop()] += 1
                except KeyError:
                    hops.loc[set([all_idx1[i] for i in idx1_matches]).pop(),set([all_idx2[i] for i in idx2_matches]).pop()] = 1
            elif len(set(idx1_matches).intersection(idx2_matches)) == 1: # good read; idx1 and idx2 line up in one spot
                # print(indexes.iloc[set(idx1_matches).intersection(idx2_matches).pop()]) #this logic might be faulty
                pass
            else:
                raise NameError('more than one unique match for barcodes!')

print(hops)


#todo: logic to prevent duplication of effort: store 'known' barcodes in a dict as you read the file.
#then, consult the list before processing the record.


