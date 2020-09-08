import argparse
import os
import pandas as pd
import gzip
from itertools import islice

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

with gzip.open(args.r1, 'rt') as read1:
    id,seq,plus,qual = islice(read1, 4)
    bc1, bc2 = id.split(":")[-1].split("+")
    if bc1 in 