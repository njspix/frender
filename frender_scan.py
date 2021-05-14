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
import csv
import gzip
from itertools import islice, repeat
import sys
from time import perf_counter
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


def analyze_barcode_wrapper(barcode, num_reads, indexes, num_subs, rc_mode):
    idx1, idx2 = barcode.split("+")

    # analyze barcode using supplied idx2
    temp = analyze_barcode(idx1, idx2, indexes, num_subs)

    if rc_mode:
        # then if matched_idx1 is not na, update only 'matched_rc_idx2', 'rc_read_type', 'rc_sample_name' for rc idx2
        if not temp["matched_idx1"] == "":
            test2 = analyze_barcode(idx1, idx2, indexes, num_subs, True)
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
            test2 = analyze_barcode(idx1, idx2, indexes, num_subs, True)
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
    else:
        return temp


def get_col(pattern, cols):
    a = [bool(regex.match(pattern, string, flags=regex.IGNORECASE)) for string in cols]
    return [i for i, x in enumerate(a) if x][0]


def get_indexes(barcode_assoc_table):
    with open(barcode_assoc_table, newline="") as f:
        header = next(csv.reader(f))
        id_col = get_col(".*id.*", header)
        idx1_col = get_col(".*index.*1.*", header)
        idx2_col = get_col(".*index.*2.*", header)
        interesting_cols = [header[id_col], header[idx1_col], header[idx2_col]]

        all_indexes = {}
        with open(barcode_assoc_table, newline="") as f:
            for row in csv.DictReader(f):
                for column, value in row.items():
                    if column in interesting_cols:
                        all_indexes.setdefault(column, []).append(value)
    return all_indexes


def frender_scan(
    barcode_assoc_table,
    fastq,
    cores,
    num_subs,
    rc_mode,
    out_csv_name,
):
    """Scan a single fastq file, counting exact and inexact barcode matches, conflicting barcodes, index hops, and undetermined reads. No demultiplexing is performed.

    Inputs -
        barcode_assoc_table - CSV file containing barcodes
        fastq - Read 1 FASTQ file
        cores - number of cores to use when processing barcodes
        num_subs - number of substitutions allowed
        rc_mode - consider reverse complement of index 2 in analysis?
        out_csv_name - Output csv name

    Returns -
        CSV file with the following columns:
            idx1 - Index 1 in barcode
            idx2 - Index 2 in barcode
            reads - total number of reads associated with this pair of barcodes
            matched_idx1 - if index 1 matches an index 1 in the supplied table, the first match is printed here
            matched_idx2 - if index 2 matches an index 2 in the supplied table, the first match is printed here
            read_type - one of:
                - 'undetermined'
                - 'index_hop' (exactly one match found for each index 1 and index 2, but the matched indexes are associated with different samples)
                - 'ambiguous' (barcodes match more than one sample)
                - 'demuxable' (barcodes match to one sample)
            sample_name - if 'read_type' is 'demuxable', the sample name associated with this pair of indexes
        If rc_mode is TRUE, the following columns are also included:
            matched_rc_idx2 - if the reverse complment of index 2 matches an index 2 in the supplied table, the first match is printed here (in original, not reverse complemented, format)
            rc_read_type - the read type found using the *reverse complement* of the supplied index 2
            rc_sample_name - if 'rc_read_type' is 'demuxable', the sample name associated with index 1 and the *reverse complement* of the supplied index 2
    """

    indexes = get_indexes(barcode_assoc_table)
    barcode_counter = {}

    with gzip.open(fastq, "rt") as read_file:
        for read_head in islice(read_file, 0, None, 4):
            code = read_head.rstrip("\n").split(" ")[1].split(":")[-1]

            try:
                barcode_counter[code] += 1
            except KeyError:
                barcode_counter[code] = 1

    print(
        f"\nScanning complete! Analyzing {len(barcode_counter)} barcodes...",
        file=sys.stderr,
    )

    if cores > 1:
        with Pool(processes=cores) as pool:
            print(f"multiprocessing with {cores} cores")
            results = pool.starmap(
                analyze_barcode_wrapper,
                zip(
                    barcode_counter,
                    barcode_counter.values(),
                    repeat(indexes),
                    repeat(num_subs),
                    repeat(rc_mode),
                ),
            )

    else:
        results = list(
            map(
                analyze_barcode_wrapper,
                barcode_counter,
                barcode_counter.values(),
                repeat(indexes),
                repeat(num_subs),
                repeat(rc_mode),
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
        "-i",
        nargs=1,
        help="Gzipped fastq file to be scanned. If analyzing paired-end data,only read 1 need be analyzed. ",
        required=True,
        metavar="input.fastq.gz",
    )
    parser.add_argument(
        "-b",
        help="Barcode association table, csv format",
        required=True,
        metavar="barcodeAssociationTable.csv",
    )
    parser.add_argument(
        "-n",
        help="Number of substitutions allowed in barcode when matching (default = 1)",
        default=1,
        type=int,
        required=True,
        metavar="num_subs",
    )
    parser.add_argument(
        "-rc",
        "--reverse_complement",
        action="store_true",
        help="Also scan for reverse complement of index 2 (to check for mistakes with e.g. HiSeq 4000 and other systems)",
    )
    parser.add_argument(
        "-c",
        "--cores",
        help="Number of cores to use for analysis, default = 1",
        default=1,
        type=int,
        metavar="cores",
    )
    parser.add_argument(
        "-o",
        default="frender_scan_output.csv",
        help="output csv file name",
        required=True,
        metavar="output.csv",
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
