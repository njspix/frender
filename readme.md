# frender

Fastq REad Name DemultiplexER

---

## Overview

Next-generation sequencing projects often pool multiple samples in the same sequencing lane to increase throughput and data-to-cost ratio. Reads from such 'multiplexed' sequencing runs must then be deconvolved in a process known as demultiplexing. Often, this will be done by Illumina's `bcl2fastq` tool using a 'sample sheet', which provides associations between samples and their index sequences. This tool produces a `fastq` file (or files, if paired-end sequencing was performed) for each sample, as well as an `Undetermined` file(s) containing reads that could not be unambiguously assigned to a sample.

Generally, these files are the starting point for subsequent analyses. However, you might desire to inspect these files for various reasons:

* You didn't supply a sample sheet and thus all your reads are in an `Undetermined` file, and you want to demultiplex them yourself
* You have read about [index hopping](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/index-hopping.html#:~:text=What%20is%20Index%20Hopping%3F%20Index%20hopping%20or%20index,to%20a%20different%20index%20%28in%20the%20multiplexed%20pool%29.) and are curious whether this affects your data (if you use combinatorial dual indexes, it does :wink: ), and to what extent
* You use a sequencing service or vendor, and want to double-check whether your files are correctly demultiplexed, or (if not), re-demultiplex them yourself
* Some of your samples are mysteriously missing or have nonsensical data, while other samples appear to be just fine
* You are dissatisfied with the proportion of reads in the `Undetermined` file and want to investigate if you could recover more data (for example, by allowing 1 or 2 mismatches when identifying indexes)

`frender` is designed to address these needs, and more!

## Features

* **Reverse complement mode:** some Illumina machines (e.g. HiSeq 4000) read the reverse complement of Index 2 rather than the forward sequence (due to their chemistry). This can be confusing, but some sequencing providers  automatically 'correct' if you submit index 2 sequences in the wrong orientation. This is wonderful when it works, but when it doesn't (e.g. low read number per sample), you may get some very confusing results! To address this, `frender` has a  `-rc` flag, which tries both the forward and reverse-complement sequence of each index 2, and chooses the orientation that makes the most sense in the context of all the sequence data provided.
* **Allowing mismatches:** `frender` allows you to specify how many mismatches to tolerate when comparing sequencing data to your sample sheet. This allows you to:

  * Check if your sequencing provider demultiplexed your data according to your specification
  * See if you can recover more data by specifying a larger number of mismatches (since most indexes are 3-4 edits away from each other, this could be reasonable in some circumstances)
  * Filter your data more stringently by specifying a smaller number of mismatches

* **Full data reporting:** `frender` provides a complete report containing all unique index combinations and their frequency in a `.csv` format. This can easily be picked up (e.g. by `R` scripts) for further analysis of index hopping, sample swaps, contamination, etc.
* No dependencies (besides the Python 3 standard libraries) to reduce installation headaches
* Multicore support for fast processing

---

## Usage

### Quick Start

1. Scan a file, list of files, or directory:

    * `python3 ./frender.py scan -n 1 -rc -b barcode_table.csv input.fastq.gz [input2.fastq.gz] [input3.fastq.gz]`
    * `python3 ./frender.py scan -n 1 -rc /path/to/sequencing/directory`

    Option|Explanation
    --|--
    -n 1 | allow 1 mismatch when trying to match indexes
    -rc | consider the reverse complement as well as the forward sequence of index 2
    -b | barcode association table, `.csv` format

2. Demultiplex based on scan results:

    * `python3 ./frender.py demux -r frender-scan-results_1-mismatches_sequencing-directory.csv input.fastq.gz [input2.fastq.gz] [input3.fastq.gz]`
    * `python3 ./frender.py demux -r frender-scan-results_1-mismatches_sequencing-directory.csv /path/to/sequencing/directory`

    Option|Explanation
    --|--
    -r | `frender` scan result file (contains information necessary for demultiplexing)

### More details

#### `scan` subcommand

##### Inputs

* `-b` (barcode association table or sample sheet)
  * Must be in csv format or contain csv formatted data after the line `[Data]` (Illumina sample sheets use this format)
  * This file must contain a header row with column names similar to **'index'**, **'index2'**, and **'id'** or **'name'**. All indexes in a column must be the same length (although index 1 and index 2 can be different lengths); this should correspond to the length of the indexes in the `fastq` files to be scanned.
  * You must specify a barcode association table unless you supply `frender` with a single directory path that already contains such a file.
    * If this is the case, the directory will be recursively searched for `.csv` or `.txt` files that match `barcode.*association` or `sample.*sheet` (case insensitive). If multiple files are found, the one with the shortest path will be selected.
    * If you specify a barcode association table with the `-b` option and also supply a directory containing such a file, the explicitly specified file takes precedence

* `input.fastq.gz` or `/path/to/dir` (input file(s)/directory)
  * If a directory is specified, it will be search recursively for `*.fastq.gz` and `*.fq.gz` files. For the purposes of the `scan` function, only files matching `_R1_` (read 1 files) are used, as the barcodes are identical in both read 1 and read 2.

##### Options

###### `-c` (cores)

(default: `1`)
* `0`: use all available cores
* `0 < c < 1`: use a fraction of available cores (e.g. `-c 0.5` would use 4 cores on an 8-core machine)
* `c ≥ 1`: use the specified number of cores 

###### `-n` (number of mismatches)

* Each barcode is allowed to have this number of mismatches. At this time you cannot specify a different number of mismatches for index 1 and index 2.
* No default; this parameter must be specified

###### `-rc` (reverse complement mode)

* If this flag is set, for each sample, `frender` will count the number of reads that could be assigned to that sample (from the provided fastq files) using the forward (supplied) index 2 sequence, and the number of reads that could be assigned with the reverse complement of that sequence. If the reverse complement would yield more demuxable reads than the forward sequence, that sample is flagged as 'reverse complement'. A report is printed at the command line (and saved to a `frender-index-2-calls.csv` file):

    Sample Name|Supplied Index 2|Reads supporting (forward)|Reverse complement Index 2|Reads supporting (rev comp)|Final call
    --|--|--|--|--|--
    FT-SA49175|ACTTGAAT|3374|ATTCAAGT|0|forward
    FT-SA49189|ATGATCTG|1|CAGATCAT|1021281|reverse complement

###### `-s` (sample)

* Optional; if specified, only the first `n` reads in the file (s) will be examined.

###### `-o` (output infix)

* Optional; if specified, this string will be added to the output file name

###### `-p` (prefix to remove)

* Optional; if specified, this string will be removed from the sample id's provided in the barcode association table (for the purpose of calling properly/improperly demuxed)

##### Outputs

###### `frender-scan-results_1-mismatches_testdir.csv` (barcode tallies)

A data file containing the following metrics for each unique combination of barcodes discovered in the scanned `fastq` files:

```txt
idx1,idx2,reads,matched_idx1,matched_idx2,read_type,sample_name,demux_ok
GATCAGCG,CTTGTAAT,1000000,GATCAGCG,CTTGTAAT,demuxable,FT-SA49186,False
NTCACGTT,NTCACGAT,59,ATCACGTT,ATCACGAT,index_hop,True
NTCCCGTT,NAGATCAT,59,,,undetermined,True
```

`reads` refers to the total number of reads that possess the specified `idx1` and `idx2`, summed over all the provided input files.

Possible `read_type`s are:

* `demuxable`: this barcode combination can be uniquely assigned to one sample file
* `index_hop`: this barcode combination matches at least one of the supplied index 1s and at least one of the supplied index 2s, but none of the matched indexes are associated with the same sample.
* `ambiguous`: this barcode combination matches the indexes assigned to more than one sample
* `undetermined`: none of the above conditions are met

If `idx1` and `idx2` can be unambiguously assigned to a sample in the provided barcode association table, `sample_name` is populated accordingly. This field is blank for all other reads.

After analyzing all supplied files and the barcodes found in them, `frender` will determine if the barcodes are correctly distributed among the provided files according to the following rules:

* `Index-hop` and `ambiguous` reads are found _only_ in files whose names match `index-hop` and `ambiguous` (case-insensitive), respectively, or in a file whose name matches `undetermined` (case-insensitive).
* `Undetermined` reads are found _only_ in a file whose name matches `undetermined` (case-insensitive).
* `Demuxable` reads are found _only_ in a file whose name matches the sample id in the barcode association file.

If `frender` detects any incorrectly demuxed reads, a warning will be issued and the `demux_ok` flag will be set to 'False' for that barcode. A list of incorrectly demuxed files will also be printed to the terminal.

###### `frender-index-2-calls_1-mismatches_testdir.csv` (Index 2 calls)

Produced only if `-rc` is specified. A data file containing the folllowing metrics for each sample in the barcode association file:

```txt
sample_name,supplied_index_2,reads_supplied_index_2,rc_index_2,reads_rc_index_2,use_rc
FT-SA49175,ACTTGAAT,3374,ATTCAAGT,0,FALSE
FT-SA49189,ATGATCTG,1,CAGATCAT,1021281,TRUE
```

`reads_supplied_index_2` and `reads_rc_index_2` refer to the number of reads that could be demuxed if the forward or reverse complement of the supplied index 2 is used, respectively. `use_rc` is `FALSE` if `frender` used the supplied index 2 sequence for that sample when writing the `scan-results` file, and `TRUE` if the reverse complement was used. The format of the `scan-results` file is unchanged whether or not the `-rc` flag was set.

#### `demux` subcommand

##### Inputs

* `input.fastq.gz` or `/path/to/dir` (input file(s)/directory)

  * If a directory is specified, it will be search recursively for `*.fastq.gz` and `*.fq.gz` files. To allow `frender` to pair files properly, Read 1 and read 2 files must have *identical* paths and names except for the single number (`1` or `2`) denoting the read.
  * If individual files are specified, you need to include both read 1 and read 2 files. `frender` will pair them up for you automatically based on their names, as above.
* `-r` (frender result file)
  * A `frender-scan-result` file is required. This file associates index combinations with read types and sample names. The content of this file is affected by the parameters of the `frender scan` command as well as the input `fastq` files provided to that command.
  * As a result, this file **must** be generated using the exact same files that are fed to `frender demux`.

##### Options

By default, `frender` produces a pair of output files (read 1 and read 2) for each sample in the sample sheet, as well as a pair of `Undetermined` files, a pair of `Ambiguous` files, and a pair of `Index-hop` files (corresponding to the read types above). You can alter this behavior with the following switches:

Flag|Action
--|--
`-a`|Don't produce separate files for ambiguous reads (these reads will be sent to the `Undetermined` file unless `-u` is used; the name of the `Undetermined` file will be altered to reflect this)
`-i`|Don't produce separate files for index hopped reads (these reads will be sent to the `Undetermined` file unless `-u` is used; the name of the `Undetermined` file will be altered to reflect this)
`-s`|Don't produce individual sample files
`-u`|Don't produce `Undetermined` files

`-o` (output infix)

* Optional; if specified, this string will be added to the output file names

`-d` (output directory)

* All output files will be collected in this directory. Defaults to `./frender-demux-output_{current date and time}`.

##### Outputs

Gzipped fastq files as noted above.

---

## Development

I wrote this script partially to address our lab's specific needs and partially to develop my skillset. I  hope to add updates in the future (speedup/optimization, single index support, more options) but no guarantees are made. I hope this may be helpful to others in similar situations.

Thanks to @jamorrison for help with refactoring and polishing.

This software is released under GPL v3 or later.
