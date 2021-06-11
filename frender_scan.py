"""
Fastq REad Name DemultiplexER

TODO: ADD ONE SENTENCE EXPLANATION

Creator:
    Nathan Spix

Date Created:
    May 2021

Inputs:
    TODO: ADD DESCRIPTIONS OF INPUTS

OUTPUT/RETURNS:
    TODO: ADD DESCRIPTION OF OUTPUT FILES
"""

import argparse
import csv
import gzip
from itertools import islice, repeat
from multiprocessing import Pool
import re
import os
from math import floor


def get_indexes_of_approx_matches(query, list_of_strings, hamming_dist):
    """Returns a list containing *indexes* of matches to query in list_of_strings within hamming_dist.
    Since all strings must be the same length, hamming_dist is equivalent to the number of substitutions/differences between strings.
    Case insensitive.
    """
    if list_of_strings == []:
        return []

    else:
        result = []
        for i in range(len(list_of_strings)):
            str1, str2 = query.lower(), list_of_strings[i].lower()
            assert len(str1) == len(
                str2
            ), f"Barcode {str1} doesn't match length of supplied barcode {str2}"
            if len([0 for a, b in zip(str1, str2) if a != b]) <= hamming_dist:
                result += [i]
            else:
                pass
        return result


def reverse_complement(string):
    return string.translate(str.maketrans("ATGCNatgcn", "TACGNtacgn"))[::-1]


def analyze_barcode(idx1, idx2, all_idx1, all_idx2, all_ids, num_subs, rc_flag=False):
    """Determine which sample a combination of indexes belongs to.
    Inputs:
        - idx1: the first barcode (read 1 index/barcode)
        - idx2: the second barcode (read 2 index/barcode)
        - all_indexes: dict of lists "idx1", "idx2", and "id", which contain, in the same order,
            all index1s, index2s, and ids from the input csv.
        - num_subs: number of substitutions to allow when matching fastq barcodes to barcodes in the input csv
        rc_flag: if True, try to match using the reverse complement of index 2 as well as the sequence in the input csv
    Returns:
        A dict containing these entries:
            - "matched_idx1": the first index 1 in the input csv that matched
            - "matched_idx2": the first index 2 in the input csv that matched
            - "read_type": one of:
                - "index_hop": matches a supplied index1 and index2, but these matches are not associated with the same sample (e.g. maches index 1 for sample 'sample_x' and index 2 for sample 'sample_y')
                - "demuxable": maches one supplied index1 and index2 - an uniquely assignable read
                - "undetermined": no matches found for index1, index2, or both
                - "ambiguous": could be assigned to more than one sample, e.g. index 1 and index 2 both match 'sample x' and 'sample y'
            - "sample_name": if read_type is "demuxable", the associated sample id in the innput csv
            - "rc_flag": was the computation done with the rc_flag True?
    """
    all_idx2 = [reverse_complement(i) for i in all_idx2] if rc_flag else all_idx2

    idx1_matches = get_indexes_of_approx_matches(idx1, all_idx1, num_subs)
    idx2_matches = get_indexes_of_approx_matches(idx2, all_idx2, num_subs)

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
            sample_name = all_ids[match_isec.pop()]
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


def analyze_barcode_wrapper(
    barcode, num_reads, all_idx1, all_idx2, all_ids, num_subs, rc_mode
):
    """Convenience function to call and format the output of analyze_barcode.
    Inputs:
        - barcode: extracted barcode from fastq header line (format [ACTG]+\+[ACTG]+)
        - num_reads: number of reads with this barcode in fastq file
        - all_idx1, all_idx2, all_ids: lists of index1/index2/sample ids from input csv file. Each list must be in the same order
        - num_subs: number of substitutions allowed when matching barcode to supplied indexes
        - rc_mode: if True, try to match using the reverse complement of index 2 as well as the sequence in the input csv

    """
    idx1, idx2 = barcode.split("+")[0:2]

    # analyze barcode using supplied idx2
    temp = analyze_barcode(
        idx1, idx2, all_idx1, all_idx2, all_ids, num_subs, rc_flag=False
    )

    result = {
        "idx1": idx1,
        "idx2": idx2,
        "reads": num_reads,
        "matched_idx1": temp["matched_idx1"],
        "matched_idx2": temp["matched_idx2"],
        "read_type": temp["read_type"],
        "sample_name": temp["sample_name"],
    }

    if rc_mode:
        rc_temp = analyze_barcode(
            idx1, idx2, all_idx1, all_idx2, all_ids, num_subs, True
        )

        # if we already have a match for idx1, don't update it
        idx1_match = (
            rc_temp["matched_idx1"]
            if temp["matched_idx1"] == ""
            else temp["matched_idx1"]
        )

        result.update(
            {
                "matched_idx1": idx1_match,
                "matched_rc_idx2": rc_temp["matched_idx2"],
                "rc_read_type": rc_temp["read_type"],
                "rc_sample_name": rc_temp["sample_name"],
            }
        )

    return result


def get_col(pattern, cols):
    """Return the index of the first entry in cols matching pattern, case insensitive"""
    a = [bool(re.match(pattern, string, flags=re.IGNORECASE)) for string in cols]
    return [i for i, x in enumerate(a) if x][0]


def get_indexes(barcode_assoc_table):
    """Get indexes and ids from a supplied csv file
    Inputs: csv format file with columns matching ".*id.*", ".*index.*1.*", and ".*index.*2.*" (case insensitive)
    Returns: dict of 3 lists (id, idx1, idx2) containing the values from their respective columns in the csv. Entries are in the same order.
    """
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
    all_indexes["id"] = all_indexes.pop(header[id_col])
    all_indexes["idx1"] = all_indexes.pop(header[idx1_col])
    all_indexes["idx2"] = all_indexes.pop(header[idx2_col])

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
        fastq - Read 1 FASTQ file (read 2, if present, should be identical, and is not used)
        cores - number of cores to use when processing barcodes:
            0: use all available cores (autodetect)
            0 < cores < 1: use fraction of available cores, e.g. cores = 0.5 will use 4 cores on an 8-core system (autodetect)
            cores >= 1: use this number of cores
        num_subs - number of substitutions allowed when matching discovered barcodes against supplied indexes. Applied independently to each index
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
        If rc_mode is True, the following columns are also included:
            matched_rc_idx2 - if the reverse complment of index 2 matches an index 2 in the supplied table, the first match is printed here (in original, not reverse complemented, format)
            rc_read_type - the read type found using the *reverse complement* of the supplied index 2
            rc_sample_name - if 'rc_read_type' is 'demuxable', the sample name associated with index 1 and the *reverse complement* of the supplied index 2
    """

    # Calculate number of cores to use
    assert cores >= 0, "Number of cores is negative... what does that mean?"
    try:
        avail_cores = len(os.sched_getaffinity(0))
    except AttributeError:
        avail_cores = os.cpu_count()

    if cores == 0:
        temp = avail_cores
    elif cores > 0 & cores < 1:
        temp = max(floor(cores * avail_cores), 1)
    else:  # cores >= 1
        temp = cores
    cores = temp

    # Count all barcodes found in fastq file
    barcode_counter = {}
    with gzip.open(fastq, "rt") as read_file:
        for read_head in islice(read_file, 0, None, 4):
            code = (
                read_head.rstrip("\n").split(" ")[1].split(":")[-1]
            )  # works for header format @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:AAAAAAAA+GGGGGGGG
            try:
                barcode_counter[code] += 1
            except KeyError:
                barcode_counter[code] = 1

    print(
        f"Scanning complete! Analyzing {len(barcode_counter)} barcodes...",
    )

    indexes = get_indexes(barcode_assoc_table)
    all_idx1 = indexes["idx1"]
    all_idx2 = indexes["idx2"]
    all_ids = indexes["id"]

    if cores > 1:
        with Pool(processes=cores) as pool:
            print(f"multiprocessing with {cores} cores")
            results = pool.starmap(
                analyze_barcode_wrapper,
                zip(
                    barcode_counter,
                    barcode_counter.values(),
                    repeat(all_idx1),
                    repeat(all_idx2),
                    repeat(all_ids),
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
                repeat(all_idx1),
                repeat(all_idx2),
                repeat(all_ids),
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

    if not args.o.endswith(".csv"):
        args.o = args.o + ".csv"

    rc_mode_text = (
        "both supplied and reverse-complement index 2 sequences"
        if args.reverse_complement
        else "index 2 sequences as supplied"
    )
    print(f"Scanning {args.i[0]} using {rc_mode_text}...")
    frender_scan(
        args.b,
        args.i[0],
        args.cores,
        args.n,
        args.reverse_complement,
        args.o,
    )
