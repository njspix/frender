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
