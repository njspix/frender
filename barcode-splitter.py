import argparse
import os
import pandas as pd
import gzip
from itertools import zip_longest
import regex

parser = argparse.ArgumentParser()
parser.add_argument("-c", help = "Barcode association table, csv format", required = True)
parser.add_argument("-f", help = "Undetermined fastq", required = True)
parser.add_argument("-o", help = "output files directory", default = './split_files')

args = parser.parse_args()

#read barcode association file:
indexes = pd.read_csv(args.c, header=0, names=['id', 'idx1', 'idx2']).set_index('id')
all_idx1 = indexes['idx1'].tolist() 
all_idx2 = indexes['idx2'].tolist()
ids = list(indexes.index)
#todo: handle single index
#todo: check that barcodes are all same length, match [ACTGactg]

# create output dir
args.o = args.o.rstrip('/')+'/' #add trailing slash just to make sure
#todo: what if this doesn't start with './'?
try:
    os.makedirs(args.o)
except OSError as e:
    if e.errno != 17: #17 = file exists
        raise

#open output files:
for i in range(len(ids)):
    ids[i] = open(f'{args.o}{ids[i]}.fastq', 'w')

undeter = open(f'{args.o}undetermined.fastq', 'w')
hop = open(f'{args.o}hop.fastq', 'w')

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

with gzip.open(args.f, 'rt') as read1:
    for record in grouper(read1, 4, ''):
        assert len(record) == 4
        idx1 = record[0].split(":")[-1].split("+")[0].rstrip('\n')
        idx2 = record[0].split(":")[-1].split("+")[1].rstrip('\n')
        idx1_matches = fuz_match_list(idx1, all_idx1)
        idx2_matches = fuz_match_list(idx2, all_idx2) 
        if (bool(idx1_matches) & bool(idx2_matches)):
            #print(idx1, idx2)
            if len(set(idx1_matches).intersection(idx2_matches)) == 0: # index hop

                for line in record:
                    hop.write(str(line))

                try:
                    hops.loc[set([all_idx1[i] for i in idx1_matches]).pop(),set([all_idx2[i] for i in idx2_matches]).pop()] += 1 
                except KeyError: #this combination of indices hasn't been initialized
                    hops.loc[set([all_idx1[i] for i in idx1_matches]).pop(),set([all_idx2[i] for i in idx2_matches]).pop()] = 1
            elif len(set(idx1_matches).intersection(idx2_matches)) == 1: # good read; idx1 and idx2 line up in one spot
                demux_id = indexes.index[set(idx1_matches).intersection(idx2_matches).pop()]
                for line in record:
                    ids[indexes.index.get_loc(demux_id)].write(str(line))
            else:
                raise NameError('more than one unique match for barcodes!')
        else:
            for line in record:
                undeter.write(str(line))

#close output files:
for i in range(len(ids)):
    ids[i].close()
hop.close()
undeter.close()

hops = hops.reset_index().melt(id_vars = 'index', var_name='idx2', value_name= 'num_hops_observed')
hops = hops[hops['num_hops_observed'] > 0]
hops.columns = ['idx1', 'idx2', 'num_hops_observed']
hops.to_csv(f'{args.o}barcode_hops.csv', index = False)

#todo: logic to prevent duplication of effort: store 'known' barcodes in a dict as you read the file.
#then, consult the list before processing the record.


