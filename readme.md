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


## Usage

```./frender_scan.py --help``` provides a listing of the various options:


