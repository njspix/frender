"""
Fastq REad Name DemultiplexER

Scans a fastq file and parses barcode names in the read headers, given known barcodes and sample ids.

Creator:
    Nathan Spix

Date Created:
    May 2021

Inputs:
    - input.fastq.gz: Gzipped fastq file to be scanned. If analyzing paired-end data,only read 1 need be analyzed.
    - barcodeAssociationTable.csv: Barcode association table, csv format. Must include header with column names similar to 'index1', 'index2', and 'id'
    - num_subs: number of substitutions to allow when trying to match barcodes
    - rc: if set, also scan for reverse complement of index 2 (to check for mistakes with e.g. HiSeq 4000 and other systems). For each supplied id, frender will determine whether more reads can be demultiplexed using the forward (supplied) sequence or the reverse complement sequence. Barcodes will then be re-processed accordingly.
    - c: Number of cores to use for analysis, default = 1. Use 0 for all available, a number between 0 and 1 for a fraction of all available cores, or a number >=1 for a specified number of cores

OUTPUT/RETURNS:
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
"""


import os, argparse, re, csv, gzip
from itertools import islice, repeat, zip_longest
from multiprocessing import Pool
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


def analyze_barcode(idx1, idx2, all_idx1, all_idx2, all_ids, num_subs):
    """Determine which sample a combination of indexes belongs to.
    Inputs:
        - idx1: the first barcode (read 1 index/barcode)
        - idx2: the second barcode (read 2 index/barcode)
        - all_indexes: dict of lists "idx1", "idx2", and "id", which contain, in the same order,
            all index1s, index2s, and ids from the input csv.
        - num_subs: number of substitutions to allow when matching fastq barcodes to barcodes in the input csv
    Returns:
        A dict containing these entries:
            - "matched_idx1": the first index 1 in the input csv that matched
            - "matched_idx2": the first index 2 in the input csv that matched
            - "read_type": one of:
                - "index_hop": matches a supplied index1 and index2, but these matches are not associated with the same sample (e.g. maches index 1 for sample 'sample_x' and index 2 for sample 'sample_y')
                - "demuxable": maches one supplied index1 and index2 - an uniquely assignable read
                - "undetermined": no matches found for index1, index2, or both
                - "ambiguous": could be assigned to more than one sample, e.g. index 1 and index 2 both match 'sample x' and 'sample y'
            - "sample_name": if read_type is "demuxable", the associated sample id in the input csv
    """
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
    }


def analyze_barcodes_with_rc(
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
    temp = analyze_barcode(idx1, idx2, all_idx1, all_idx2, all_ids, num_subs)

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
        rc_all_idx2 = [reverse_complement(i) for i in all_idx2]
        rc_temp = analyze_barcode(idx1, idx2, all_idx1, rc_all_idx2, all_ids, num_subs)

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

        # in some cases, both forward index 2 and reverse complement index 2 could result in a valid demux call.
        # Test for this, and re-call as 'ambiguous' if this is the case.
        if (temp["read_type"] == "demuxable") & (rc_temp["read_type"] == "demuxable"):
            if (
                temp["sample_name"] == rc_temp["sample_name"]
            ):  # palindromic index 2, this is unusual but possible...
                pass
            else:  # This barcode is actually ambiguous:
                result.update(
                    {
                        "read_type": "ambiguous",
                        "sample_name": "",
                        "rc_read_type": "ambiguous",
                        "rc_sample_name": "",
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


def rc_mode_test(results_list):
    """Test if a list produced by read_scan_results includes reverse complement data"""
    return "rc_read_type" in results_list[0].keys()


def get_ids(results_list):
    """Get all sample ids from a list of dicts (same format as generated in frender_scan function)"""

    rc_mode = rc_mode_test(results_list)

    ids = []

    for entry in results_list:
        if rc_mode:
            if (entry["rc_sample_name"] != "") & (entry["rc_sample_name"] not in ids):
                ids += [entry["rc_sample_name"]]
        if (entry["sample_name"] != "") & (entry["sample_name"] not in ids):
            ids += [entry["sample_name"]]
    return ids


def call_rc_mode_per_id(results_list):
    """Given a list of dicts (same format as generated in frender_scan_function), for each sample id found, determine whether it should be demuxed with the forward or reverse complement index 2.
    Returns: a dictionary with each sample id and True (demux with rc index 2) or False (demux with forward index 2)
    """

    # Must have RC entries
    assert rc_mode_test(
        results_list
    ), "It looks like this frender result csv was not generated with the -rc flag. Either specify a different result csv, or run this command without setting the -rc flag."

    ids = {id: {"f": 0, "rc": 0, "demux_with_rc": ""} for id in get_ids(results_list)}

    for record in results_list:
        if record["sample_name"] != "":
            ids[record["sample_name"]]["f"] += int(record["reads"])
        if record["rc_sample_name"] != "":
            ids[record["rc_sample_name"]]["rc"] += int(record["reads"])

    for each in ids:
        if ids[each]["f"] >= ids[each]["rc"]:
            ids[each]["demux_with_rc"] = False
        else:
            ids[each]["demux_with_rc"] = True

    return {a: ids[a]["demux_with_rc"] for a in ids}


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
        cores = avail_cores
    elif 0 < cores < 1:
        cores = max(floor(cores * avail_cores), 1)
    else:  # cores >= 1
        cores = int(cores)

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
            print(f"Multiprocessing with {cores} cores")
            results = pool.starmap(
                analyze_barcodes_with_rc,
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
                analyze_barcodes_with_rc,
                barcode_counter,
                barcode_counter.values(),
                repeat(all_idx1),
                repeat(all_idx2),
                repeat(all_ids),
                repeat(num_subs),
                repeat(rc_mode),
            )
        )

    if not rc_mode:  # We're done, write out to csv.
        print(f"Analysis complete! Writing results to {out_csv_name}")
        keys = results[0].keys()
        with open(out_csv_name, "w", newline="") as output_file:
            dict_writer = csv.DictWriter(output_file, keys)
            dict_writer.writeheader()
            dict_writer.writerows(results)

    else:
        # Have preliminary results at this point. Now we can call whether to use forward or rc index 2 for each sample...
        rc_calls = call_rc_mode_per_id(results)
        rc_summary_file_name = out_csv_name.replace(
            "frender-scan-results_", "frender-index-2-calls_"
        )
        print("First round of analysis complete.")
        print(
            f"Based on the barcodes in the supplied fastq file, the following index 2 sequences will be used\n(also recorded in {rc_summary_file_name}):\n"
        )
        print(
            "Sample Name",
            "Supplied Index 2",
            "Final Index 2",
            "Reverse Complemented?",
            sep="\t",
        )
        for a in rc_calls:
            index = indexes["id"].index(a)
            print(
                a,
                indexes["idx2"][index],
                indexes["idx2"][index]
                if not rc_calls[a]
                else reverse_complement(indexes["idx2"][index]),
                "yes" if rc_calls[a] else "no",
                sep="\t",
            )
        with open(rc_summary_file_name, "w", newline="") as output_file:
            writer = csv.writer(output_file)
            writer.writerow(
                [
                    "sample_name",
                    "supplied_index_2",
                    "final_index_2",
                    "is_reverse_complemented",
                ]
            )
            for a in rc_calls:
                index = indexes["id"].index(a)
                writer.writerow(
                    [
                        a,
                        indexes["idx2"][index],
                        indexes["idx2"][index]
                        if not rc_calls[a]
                        else reverse_complement(indexes["idx2"][index]),
                        "TRUE" if rc_calls[a] else "FALSE",
                    ]
                )

        all_idx2 = [
            reverse_complement(indexes["idx2"][i])
            if rc_calls[id]
            else indexes["idx2"][i]
            for i, id in enumerate(indexes["id"])
        ]

        print(
            f"\nRe-analyzing {len(barcode_counter)} barcodes with corrected index 2 sequences...",
        )

        if cores > 1:
            with Pool(processes=cores) as pool:
                print(f"Multiprocessing with {cores} cores")
                results = pool.starmap(
                    analyze_barcodes_with_rc,
                    zip(
                        barcode_counter,
                        barcode_counter.values(),
                        repeat(all_idx1),
                        repeat(all_idx2),
                        repeat(all_ids),
                        repeat(num_subs),
                        repeat(False),
                    ),
                )

        else:
            results = list(
                map(
                    analyze_barcodes_with_rc,
                    barcode_counter,
                    barcode_counter.values(),
                    repeat(all_idx1),
                    repeat(all_idx2),
                    repeat(all_ids),
                    repeat(num_subs),
                    repeat(False),
                )
            )
        print(f"Analysis complete! Writing results to {out_csv_name}")
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
        help="In scan mode: scan for reverse complement of index 2 (to check for mistakes with e.g. HiSeq 4000 and other systems). In demultiplex mode, allow demultiplexing with reverse complement. ",
    )
    parser.add_argument(
        "-c",
        "--cores",
        help="Number of cores to use for analysis, default = 1. Use 0 for all available, a number between 0 and 1 for a fraction of all available cores, or a number >= for a specified number of cores",
        default=1,
        type=float,
        metavar="cores",
    )
    parser.add_argument(
        "-o",
        help="output csv file name",
        metavar="output.csv",
    )
    parser.add_argument(
        "-d",
        action="store_true",
        help="Demultiplex: write out reads to their assigned fastq files",
    )
    parser.add_argument(
        "-i",
        action="store_true",
        help="if -d is specified, also output fastq files with index-hopped reads",
    )
    parser.add_argument(
        "-a",
        action="store_true",
        help="if -d is specified, also output fastq files with ambiguous reads",
    )
    parser.add_argument(
        "-b",
        help="Barcode association table, csv format. Must include header with column names similar to 'index1', 'index2', and 'id'. Alternatively, you may provide a frender-formatted csv file as input when running in demultiplex mode.",
        required=True,
        metavar="barcodeAssociationTable.csv",
    )
    parser.add_argument(
        "-f",
        nargs=argparse.REMAINDER,
        help="Gzipped fastq file(s) to be scanned. If analyzing paired-end data, only one file will be used. If demultiplexing is performed (-d flag), both read 1 and read 2 files must be provided.",
        required=True,
        metavar="input.fastq.gz",
    )

    args = parser.parse_args()

    # Check proper specification of -d flag and csv file
    if check_frender_csv(args.b) & (not args.d):
        raise SystemExit(
            "It looks like the supplied csv was produced by frender. Did you forget to specify -d (demultiplex mode)?"
        )
    if args.d & (not check_frender_csv(args.b)):
        raise SystemExit(
            "-d (demultiplex mode) is specified but the csv supplied doesn't look like a frender result csv!"
        )

    if args.d:
        assert (
            len(args.f) == 2
        ), f"Exactly 2 input fastqs must be specified (found {len(args.f)})"

        print("Writing demultiplexed files...")
        frender_demux(args.f, read_scan_results(args.b), args.i, args.a, args.rc)

    # Provide a default output file name
    if not args.o:
        args.o = "frender-scan-results_" + os.path.basename(args.f[0]) + ".csv"
    if not args.o.endswith(".csv"):
        args.o = args.o + ".csv"

    rc_mode_text = (
        "both supplied and reverse-complement index 2 sequences"
        if args.reverse_complement
        else "index 2 sequences as supplied"
    )
    print(f"Scanning {args.f[0]} using {rc_mode_text}...")
    frender_scan(
        args.b,
        args.f[0],
        args.cores,
        args.n,
        args.reverse_complement,
        args.o,
    )
