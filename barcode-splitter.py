import argparse
import os
import pandas as pd
import gzip
from itertools import zip_longest

parser = argparse.ArgumentParser()
parser.add_argument("-c", help="Barcode association table, csv format", required = True)
parser.add_argument("-r1", help="Undetermined read 1", required = True)
#todo: make this optional
parser.add_argument("-r2", help="Undetermined read 2", required = True)
parser.add_argument("-o", help = "output files directory", default = './split_files/')

args = parser.parse_args()

# create output dir
try:
    os.makedirs(args.o)
except OSError as e:
    if e.errno != 17: #17 = file exists
        raise

#read barcode association file:
df = pd.read_csv(args.c, header=0, names=['id', 'idx1', 'idx2']).set_index('id')
#todo: handle single index
#todo: check that barcodes are all same length, match [ACTGactg]
indexes = df.to_dict('index')

#like
# {'11-A12_reseq1': {'idx1': 'CTTGTACT', 'idx2': 'ATCACGAT'}, '9-H9_reseq1': {'idx1': 'GATCAGCG', 'idx2': 'ACTTGAAT'}}

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

with gzip.open(args.r1, 'rt') as read1:
    for record in grouper(read1, 4, ''):
        assert len(record) == 4
        print(record[0].split(':')[-1].split("+"))

