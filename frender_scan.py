"""
Fastq REad Name DemultiplexER

TODO: ADD ONE SENTENCE EXPLANATION

Creator:
    Nathan Spix

Date Created:
    May 2021

Requires:
    package version
    package version

Inputs:
    TODO: ADD DESCRIPTIONS OF INPUTS

OUTPUT/RETURNS:
    TODO: ADD DESCRIPTION OF OUTPUT FILES
"""

import argparse
from pathlib import Path
import csv
import pandas as pd
import gzip
from itertools import zip_longest, islice, repeat
import sys
from time import perf_counter
import numpy as np
from multiprocessing import Pool

try:
    import regex
except ModuleNotFoundError:
    print("ERROR: Install regex module", file=sys.stderr)
    print(
        "\tpip install regex\tOR\tconda install -c conda-forge regex", file=sys.stderr
    )
    sys.exit(0)


def fuz_match_list(pattern, set_of_strings, num_subs):
    """Given a query string and a list of strings, return a list of indices
       where a match (+/- allowed substitutions) is found

    Inputs -
        pattern        - pattern to search for
        set_of_strings - strings to search for pattern in
        num_subs       - number of substitutions to allow
    Returns -
        list of where matches (with +/- allowed substitutions) occur in set_of_strings
    """
    if set_of_strings == []:
        return []

    else:
        pattern = regex.compile(
            "(?:" + pattern + "){s<=" + str(num_subs) + "}", regex.IGNORECASE
        )
        matches = [bool(regex.match(pattern, string)) for string in set_of_strings]

        return [i for i, val in enumerate(matches) if val]


def reverse_complement(string):
    return string.translate(str.maketrans("ATGCN", "TACGN"))[::-1]


def analyze_barcode(idx1, idx2, all_indexes, num_subs, rc_flag=False):
    all_idx1 = all_indexes["Index1"].tolist()
    all_idx2 = (
        [reverse_complement(i) for i in all_indexes["Index2"].tolist()]
        if rc_flag
        else all_indexes["Index2"].tolist()
    )

    idx1_matches = fuz_match_list(idx1, all_idx1, num_subs)
    idx2_matches = fuz_match_list(idx2, all_idx2, num_subs)

    if bool(idx1_matches) and bool(idx2_matches):
        # Can find at least one barcode match for both indices
        matched_idx1 = all_idx1[idx1_matches[0]]
        matched_idx2 = all_idx2[idx2_matches[0]]

        match_isec = set(idx1_matches).intersection(idx2_matches)

        if len(match_isec) == 0:
            read_type = "index_hop"
            sample_name = ""

        elif len(match_isec) == 1:
            # this is a good read
            sample_name = all_indexes.index[match_isec.pop()]
            read_type = "demuxable"

        else:
            # this is an ambiguous read
            read_type = "ambiguous"
            sample_name = ""

    else:
        matched_idx1 = ""
        matched_idx2 = ""
        read_type = "undetermined"
        sample_name = ""

    return {
        "matched_idx1": matched_idx1,
        "matched_idx2": matched_idx2,
        "read_type": read_type,
        "sample_name": sample_name,
        "rc_flag": rc_flag,
    }


def analyze_barcode_wrapper(barcode, num_reads, indexes):
    idx1, idx2 = barcode.split("+")

    # analyze barcode using supplied idx2
    temp = analyze_barcode(idx1, idx2, indexes, 1)

    # then if matched_idx1 is not na, update only 'matched_rc_idx2', 'rc_read_type', 'rc_sample_name' for rc idx2
    if not temp["matched_idx1"] == "":
        test2 = analyze_barcode(idx1, idx2, indexes, 1, True)
        final = {
            "idx1": idx1,
            "idx2": idx2,
            "reads": num_reads,
            "matched_idx1": temp["matched_idx1"],
            "matched_idx2": temp["matched_idx2"],
            "read_type": temp["read_type"],
            "sample_name": temp["sample_name"],
            "matched_rc_idx2": test2["matched_idx2"],
            "rc_read_type": test2["read_type"],
            "rc_sample_name": test2["sample_name"],
        }
        return final

    # but otherwise, update 'matched_idx1', 'matched_rc_idx2', 'rc_read_type', 'rc_sample_name' for rc idx2
    else:
        test2 = analyze_barcode(idx1, idx2, indexes, 1, True)
        final = {
            "idx1": idx1,
            "idx2": idx2,
            "reads": num_reads,
            "matched_idx1": test2["matched_idx1"],
            "matched_idx2": temp["matched_idx2"],
            "read_type": temp["read_type"],
            "sample_name": temp["sample_name"],
            "matched_rc_idx2": test2["matched_idx2"],
            "rc_read_type": test2["read_type"],
            "rc_sample_name": test2["sample_name"],
        }
        return final


def frender_scan(
    barcode,
    fastq,
    cores,
    num_subs,
    rc_mode,
    out_csv_name,
):
    """Scan a single fastq file, counting exact and inexact barcode matches, conflicting barcodes, index hops, and undetermined reads. No demultiplexing is performed.

    Inputs -
        barcode - CSV file containing barcodes
        fastq_1 - Read 1 FASTQ file
        out_dir - Output directory name [default: '.']
        preefix - Prefix to add to output files [default: '']
        num_subs - number of substitutions allowed

    Returns -
        TODO: doc...
    """
    # Not all data includes a header in barcodeAssociationTable.txt
    # Check for how to correctly load data
    cols = ["AccessionID", "ClientAccessionID", "Index1", "Index2"]
    test = pd.read_csv(barcode, nrows=1, header=None)

    # As of 15 Sep 2020, early scWGBS data has four columns and newer data
    # has seven columns. There is one dataset with three columns, and this
    # is related to a collaborative project with human fetal DNA.
    if len(test.columns) == 3:
        print(
            "WARNING: This barcode file only has 3 columns.", end=" ", file=sys.stderr
        )
        print("It likely is a single-indexed library", file=sys.stderr)
        sys.exit(0)

    # TODO: Handle single index
    # TODO: Check barcodes are all the same length, match [ATCGatcg]
    if test[0][0].startswith("Acc") or test[0][0].startswith("Client"):
        indexes = pd.read_csv(barcode, usecols=cols, index_col="AccessionID")
    else:
        indexes = pd.read_csv(barcode, header=None, names=cols, index_col="AccessionID")

    barcode_counter = {}
    record_count, new_count = 0, 0
    start = perf_counter()

    with gzip.open(fastq_1, "rt") as read_file:
        for read_head in islice(read_file, 0, None, 4):
            if record_count % counter_interval == 1:
                end = perf_counter()
                print(
                    f"Processed {record_count} reads in {round(end-start, 2)} sec, {round(counter_interval/(end-start),0)} reads/sec, {new_count} new barcodes found",
                    file=sys.stderr,
                    end="\r",
                )
                start, new_count = perf_counter(), 0

            code = read_head.rstrip("\n").split(" ")[1].split(":")[-1]

            try:
                barcode_counter[code] += 1
            except KeyError:
                barcode_counter[code] = 1
                new_count += 1
            record_count += 1

    print(
        f"\nScanning complete! Analyzing {len(barcode_counter)} barcodes...",
        file=sys.stderr,
    )
    # Turn the barcode_counter dict into a pd dataframe

    if cores > 1:
        with Pool(processes=cores) as pool:
            print(f"multiprocessing with {cores} cores")
            results = pool.starmap(
                analyze_barcode_wrapper,
                zip(barcode_counter, barcode_counter.values(), repeat(indexes)),
            )

    else:
        results = list(
            map(
                analyze_barcode_wrapper,
                barcode_counter,
                barcode_counter.values(),
                repeat(indexes),
            )
        )

    keys = results[0].keys()
    with open(out_csv_name, "w", newline="") as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Demultiplex Undetermined FastQ files with given barcodes."
    )
    parser.add_argument(
        ["-o", "--outfile"],
        metavar="output_csv",
        default="frender_scan_output.csv",
        help="output csv file name",
    )
    parser.add_argument(
        ["-b", "--barcode_table"],
        metavar="barcodeAssociationTable.csv",
        help="Barcode association table, csv format",
        required=True,
    )
    parser.add_argument(
        ["-n", "--num_subs"],
        metavar="num_subs",
        help="Number of substitutions allowed in barcode when matching (default = 1)",
        default=1,
        type=int,
        required=True,
    )
    parser.add_argument(
        ["-c", "--cores"],
        metavar="cores",
        help="Number of cores to use for analysis, default = 1",
        default=1,
        type=int,
    )
    parser.add_argument(
        ["-rc", "--reverse_complement"],
        action="store_true",
        metavar="reverse_complement",
        help="Also scan for reverse complement of index 2 (to check for mistakes with e.g. HiSeq 4000 and other systems)",
    )
    parser.add_argument(
        "fastq",
        nargs=1,
        help="Gzipped fastq file to be scanned. If analyzing paired-end data,only read 1 need be analyzed. ",
    )

    args = parser.parse_args()

    if not args.outfile.endswith(".csv"):
        args.outfile = args.outfile + ".csv"

    rc_mode_text = (
        "both supplied and reverse-complement index 2 sequences"
        if args.rc
        else "index 2 sequences as supplied"
    )
    print(f"Scanning {args.fastq} using {rc_mode_text}...", file=sys.stderr)
    frender_scan(
        args.barcode_table,
        args.fastq,
        args.cores,
        args.num_subs,
        args.reverse_complement,
        args.outfile,
    )
