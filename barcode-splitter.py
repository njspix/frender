import argparse
import os
import pandas as pd
import gzip
from itertools import zip_longest
import regex

parser = argparse.ArgumentParser()
parser.add_argument("-c", help = "Barcode association table, csv format", required = True)
parser.add_argument("-o", help = "output files directory", default = ".")
parser.add_argument("-p", "--prefix", help = "prefix to be added to demuxed output files and index hop summary")
parser.add_argument("-h", "--hopfile", help = "(optional) output fastq for index hopped reads", default = "")
parser.add_argument("-u", "--undetermined", help = "(optional) output fastq for remaining undetermined reads", default = "")
parser.add_argument("-c", "--conflict", help = "(optional) output fastq for reads with conflicting barcodes", default = "")
parser.add_argument("fastqs", nargs= 2, required = True)

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

if args.p != "":
    args.p = args.p+"_"

r1_files = []
r2_files = []

#open output files:
for i in range(len(ids)):
    r1_files.append(open(f'{args.o}{args.p}{ids[i]}_R1.fastq', 'w'))
    r2_files.append(open(f'{args.o}{args.p}{ids[i]}_R2.fastq', 'w'))

if args.u != "" :
    r1_undeter = open(f'{args.o}{args.u}', 'w')
    r2_undeter = open(f'{args.o}{args.u}', 'w')

if args.h != "" :
    r1_hop = open(f'{args.o}{args.h}', 'w')
    r2_hop = open(f'{args.o}{args.h}', 'w')

if args.c != "":
    r1_conflict = open(f'{args.o}{args.c}', 'w')
    r2_conflict = open(f'{args.o}{args.c}', 'w')


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

with gzip.open(args.fastqs[0], 'rt') as read1:
    with gzip.open(args.fastqs[1], 'rt') as read2:
        reads_1 = grouper(read1, 4, '')
        reads_2 = grouper(read2, 4, '')

        for record_1,record_2 in zip(reads_1,reads_2):
            assert len(record_1) == 4
            idx1 = record_1[0].split(":")[-1].split("+")[0].rstrip('\n')
            idx2 = record_1[0].split(":")[-1].split("+")[1].rstrip('\n')
            idx1_matches = fuz_match_list(idx1, all_idx1)
            idx2_matches = fuz_match_list(idx2, all_idx2) 
            if (bool(idx1_matches) & bool(idx2_matches)):
                #print(idx1, idx2)
                if len(set(idx1_matches).intersection(idx2_matches)) == 0: # index hop

                    for line in record_1:
                        r1_hop.write(str(line))
                    for line in record_2:
                        r2_hop.write(str(line))

                    try:
                        hops.loc[set([all_idx1[i] for i in idx1_matches]).pop(),set([all_idx2[i] for i in idx2_matches]).pop()] += 1 
                    except KeyError: #this combination of indices hasn't been initialized
                        hops.loc[set([all_idx1[i] for i in idx1_matches]).pop(),set([all_idx2[i] for i in idx2_matches]).pop()] = 1
                elif len(set(idx1_matches).intersection(idx2_matches)) == 1: # good read; idx1 and idx2 line up in one spot
                    demux_id = indexes.index[set(idx1_matches).intersection(idx2_matches).pop()]
                    for line in record_1:
                        r1_files[indexes.index.get_loc(demux_id)].write(str(line))
                    for line in record_2:
                        r2_files[indexes.index.get_loc(demux_id)].write(str(line))
                else:
                    raise NameError('more than one unique match for barcodes!')
            else:
                for line in record_1:
                    r1_undeter.write(str(line))                
                for line in record_2:
                    r2_undeter.write(str(line))

#close output files:
for i in range(len(r1_files)):
    r1_files[i].close()
for i in range(len(r2_files)):
    r2_files[i].close()

if args.u != "" :
    r1_undeter.close()
    r2_undeter.close()

if args.h != "" :
    r1_hop.close()
    r2_hop.close()

if args.c != "" :
    r1_conflict.close()
    r2_conflict.close()

hops = hops.reset_index().melt(id_vars = 'index', var_name='idx2', value_name= 'num_hops_observed')
hops = hops[hops['num_hops_observed'] > 0]
hops.columns = ['idx1', 'idx2', 'num_hops_observed']
hops.to_csv(f'{args.o}{args.p}barcode_hops.csv', index = False)

# todo: logic to prevent duplication of effort: store 'known' barcodes in a dict as you read the file.
# then, consult the list before processing the record.


