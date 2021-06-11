# FRENDER

###### Fastq REad Name DemultiplexER

---
## Overview
Paired-end FASTQ formatted files produced by Illumina's ```bcl2fastq``` tool contain index read sequences in the read header:

```@E00503:288:HG7YWCCX2:3:2104:8775:63806 1:N:0:ATCACGTT+ATCACGAT```
([Wikipedia](https://en.wikipedia.org/wiki/FASTQ_format) has a further explanation of the this format if you're interested) 

Here ```ATCACGTT``` is the read 1 index and ```ATCACGAT``` is the read 2 index. These index sequences may be of interest to you for a couple of reasons:
* You didn't supply a sample sheet and thus all your reads are in an ```Undetermined``` file, and you want to demultiplex them yourself
* You have read about [index hopping](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/index-hopping.html#:~:text=What%20is%20Index%20Hopping%3F%20Index%20hopping%20or%20index,to%20a%20different%20index%20%28in%20the%20multiplexed%20pool%29.) and are curious whether this affects your data (if you use combinatorial dual indexes, it does :wink:), and to what extent
* You use a sequencing service or vendor, and want to double-check whether your files are correctly demultiplexed, or (if not), re-demultiplex them yourself

Frender is designed to address these needs. Given a ```.fastq.gz``` file and a ```.csv``` file associating barcodes with samples, Frender will do the following:

1. Read through the ```.fastq.gz``` file  and tally all unique index combinations, noting the total number of reads that possess each combination of indexes.
2. Attempt to match each unique index combination with a sample in the provided ```.csv``` file.
3. Output a ```.csv``` file listing the assigned classification for each index combination, along with the number of reads with that index combination and other helpful information.

## Features

* Multicore support for fast processing (~2.5 minutes for a 10 GB compressed fastq file, using 8 cores; most of this time is spent reading the gzip file.)
* No dependencies (besides the Python3 standard libraries) allows for easy installation
* Fuzzy matching: allow a specified number of mismatches when matching indexes found in the ```.fastq.gz``` file with those specified in the ```.csv``` file
* Reverse complement mode: some Illumina machines (e.g. HiSeq 4000) read the reverse complement of Index 2 rather than the forward sequence (due to their chemistry). This can lead to confusion when demultiplexing. To address this, the ```-rc``` flag instructs Frender to try both the forward and reverse-complement sequence of each index 2.

---


## Usage

#### Sample command:
```python3 ./frender_scan.py -c 8 -n 1 -rc -b barcode_table.csv -f input.fastq.gz -o test.csv```

Option|Explanation
--|--
-c 8 | process with 8 cpus
-n 1 | allow 1 mismatch when trying to match indexes (applied separately to index 1 and index 2)
-rc | consider the reverse complement as well as the forward sequence of index 2
-b | barcode association table, ```.csv``` format
-f | ```.fastq.gz``` input file
-o | output file name

#### More details:

##### ```-c``` (cores):
* ```1```: default, use one core
* ```0```: autodetect number of available cpus
* ```0.01 - 0.99```: use a fraction of available cpus, e.g. enter ```0.5``` to use half the available cpus

##### ```-b``` (barcode association table):

This file must be in comma-separated value format (although the ```.csv``` extension is unnecessary) and must contain a header row with column names similar to **'index1'**, **'index2'**, and **'id'** (case-insensitive regexes ```".*index.*1.*"```, ```".*index.*2.*"```, and ```".*id.*"``` are used to identify these columns). All barcodes must be the same length, which needs to match the length of the barcodes in the ```.fastq.gz``` file. 

##### ```-o``` (output file):

Frender produces a CSV file with the following columns:
* idx1 - Index 1 in barcode
* idx2 - Index 2 in barcode
* reads - total number of reads associated with this pair of barcodes
* matched_idx1 - if index 1 matches an index 1 in the supplied table, the first match is printed here
* matched_idx2 - if index 2 matches an index 2 in the supplied table, the first match is printed here
* read_type - one of:
    - 'undetermined'
    - 'index_hop' (exactly one match found for each index 1 and index 2, but the matched indexes are associated with different samples)
    - 'ambiguous' (barcodes match more than one sample)
    - 'demuxable' (barcodes match to one sample)
* sample_name - if 'read_type' is 'demuxable', the sample name associated with this pair of indexes

If rc_mode is True, the following columns are also included:
* matched_rc_idx2 - if the reverse complment of index 2 matches an index 2 in the supplied table, the first match is printed here (in original, not reverse complemented, format)
* rc_read_type - the read type found using the *reverse complement* of the supplied index 2
* rc_sample_name - if 'rc_read_type' is 'demuxable', the sample name associated with index 1 and the *reverse complement* of the supplied index 2.

Specifying an output file name is optional. If none is specified, the output will be named ```"frender-scan-results_{fastq_name}.csv"```

---
## Development

I wrote this script partially to address our lab's specific needs and partially to develop my skillset. I do hope to add updates in the future (demultiplexing support, single index support, more options) but no guarantee!

Thanks to @jamorrison for help with refactoring and polishing.

This software is released under GPL v3 or later. 
