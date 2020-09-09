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
#todo: handle single index
#todo: check that barcodes are all same length, match [ACTGactg]

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def fuz_match_set(pattern, set_of_strings):
    "Given a query string and a list/set of strings, return a list denoting whether the query string +/- 1 substitution is found in each item"
    pattern = regex.compile("(?:"+pattern+"){s<=1}", regex.IGNORECASE)
    matches = [bool(regex.match(pattern, set_of_strings[i])) for i in range (len(set_of_strings))]
    return matches

with gzip.open(args.r1, 'rt') as read1:
    for record in grouper(read1, 4, ''):
        assert len(record) == 4
        idx1 = record[0].split(":")[-1].split("+")[0].rstrip('\n')
        idx2 = record[0].split(":")[-1].split("+")[1].rstrip('\n')
        if sum(fuz_match_set(idx1, list(set(indexes['idx1'].tolist()))))==1:
            print(idx1, idx2)


