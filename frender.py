import argparse, os, re, csv, gzip
from pathlib import Path
from math import floor
from itertools import islice, repeat, zip_longest
from multiprocessing import Pool
from datetime import datetime, timezone


def get_cores(cores):
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
    return cores


def find_barcode_file(dir):
    dir = Path(dir)
    assert Path.is_dir(dir), "The specified directory does not exist"
    files = []
    list_of_paths = list(dir.rglob("**/*"))
    for path in list_of_paths:
        if bool(re.search("barcode.*association", str(path), re.IGNORECASE)) | bool(
            re.search("sample.*sheet", str(path), re.IGNORECASE)
        ):
            files += [path]
    # If multiple files, pick the one with the shortest path

    filtered_files = []

    for path in files:
        if bool(re.search("\.csv$|\.txt$", str(path), re.IGNORECASE)):
            filtered_files += [path]
    filtered_files.sort(reverse=True)

    if not filtered_files:  # No barcode file identified
        raise SystemExit(
            "I couldn't find a barcode table in that directory. Please either specify one with the argment -b or specify a directory including a barcode table. File names matching '.*barcode.*association.*' or '.*sample.*sheet.*' (case insensitive) are accepted."
        )
    print(f"Found barcode association file {os.path.basename(filtered_files[0])}")
    return filtered_files[0]


def handle_illumina_csv(barcode_file):
    # Return number of lines to skip from Illumina-format Sample Sheet file
    with open(barcode_file, "r") as f:
        header = next(csv.reader(f))
        if bool(re.search("\[Header\]", header[0])):
            i = 1
            while not bool(re.search("\[Data\]", next(csv.reader(f))[0])):
                i += 1
            return i + 1
        else:
            return 0


def get_col(match_pattern, cols, discard_pattern=None):
    """Return the index of the first entry in cols matching pattern, case insensitive"""
    if discard_pattern:
        a = [
            (
                bool(re.search(match_pattern, string, flags=re.IGNORECASE))
                & (not bool(re.search(discard_pattern, string, flags=re.IGNORECASE)))
            )
            for string in cols
        ]
    else:
        a = [
            bool(re.search(match_pattern, string, flags=re.IGNORECASE))
            for string in cols
        ]
    result = [i for i, x in enumerate(a) if x]

    if result:
        return result[0]
    else:
        raise ValueError(
            f"""Couldn't find column matching "{match_pattern}"{' but not "'+ discard_pattern + '"' if discard_pattern != None else ''} in csv header {cols}"""
        )


def get_indexes(barcode_file):
    """Returns: dict of 3 lists (id, idx1, idx2) containing the values from their respective columns in the csv. Entries are in the same order."""
    # TODO: need some more robust error handling here.
    skip_lines = handle_illumina_csv(barcode_file)

    with open(barcode_file, "r") as f:
        for _ in range(skip_lines):
            next(csv.reader(f))

        header = next(csv.reader(f))

        try:
            id_col = get_col("id|name", header)
            idx1_col = get_col("index", header, "id|2")
            idx2_col = get_col("index.*2", header)
        except ValueError as e:
            print("Error finding columns in provided barcode file:")
            raise SystemExit(e)

        all_indexes = {"id": [], "idx1": [], "idx2": []}

        for row in csv.reader(f):
            all_indexes["id"] += [row[id_col]]
            all_indexes["idx1"] += [row[idx1_col]]
            all_indexes["idx2"] += [row[idx2_col]]

        return all_indexes


def parse_files(file_dict, just_r1):
    paths = []
    if list(file_dict.keys())[0] == "dir":
        print(
            f"Scanning {file_dict['dir']} for fastq files. {'Using read 1 files only for speed...' if just_r1 else ''}"
        )
        for each in Path(file_dict["dir"]).rglob("**/*"):
            if Path.is_file(each):
                paths += [each]
    elif list(file_dict.keys())[0] == "file":
        paths = (
            [Path(a) for a in file_dict["file"] if Path.is_file(Path(a))]
            if type(file_dict["file"]) == list
            else [file_dict["file"]]
        )

    # Sort out non-fastq files:
    filtered_paths = []

    for path in paths:
        if bool(re.search("\.f[ast]*q\.gz$", str(path), re.IGNORECASE)):
            filtered_paths += [path]
        else:
            print(f"Ignoring non-fastq file {str(os.path.basename(path))}")

    # If we're scanning a directory, pick only the read 1 files
    if (list(file_dict.keys())[0] == "dir") & (just_r1):
        filtered_paths = [
            path
            for path in filtered_paths
            if bool(re.search("R1", str(os.path.basename(path)), re.IGNORECASE))
        ]
    return filtered_paths


def scan_file(file, sample = None):
    file_barcodes = {}
    total_barcodes = {}
    name = str(os.path.basename(file))
    print(f"Tallying barcodes from {name}...", end="")
    with gzip.open(file, "rt") as read_file:
        actual_reads, new_barcodes = 0, 0
        for read_head in islice(read_file, 0, None, 4):

            if sample:
                if actual_reads >= sample:
                    break
            actual_reads += 1

            code = (
                read_head.rstrip("\n").split(" ")[1].split(":")[-1]
            )  # works for header format @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:AAAAAAAA+GGGGGGGG
            try:
                file_barcodes[code] += 1
                total_barcodes[code] += 1
            except KeyError:
                new_barcodes += 1
                file_barcodes[code] = 1
                total_barcodes[code] = 1
    print(
        f"found {new_barcodes} new barcode{'' if new_barcodes == 1 else 's'} in {actual_reads} reads."
    )
    return (name, file_barcodes, total_barcodes)

def tally_barcodes(cores, files, sample=None):
    print(f"Scanning {len(files)} files with {cores} core{'' if cores == 1 else 's'}...")
    if sample:
        assert sample >= 1, f"Number of reads to sample must be ≥ 1!"
        print(f"Sampling {sample} reads from the head of each file...")

    if cores > 1:
        with Pool(processes=cores) as pool:
            results = pool.starmap(
                scan_file, [(file, sample) for file in files]
            )
        print(type(results), len(results))
    else:
        results = [scan_file(file, sample = sample) for file in files]
        print(type(results), len(results))

    # combine all the 'total_barcodes' dictionaries into one dictionary, adding the values for any duplicate keys
    barcode_counter = {"total": {}}
    for d in [x[2] for x in results]:
        for k, v in d.items():
            barcode_counter["total"][k] = barcode_counter["total"].get(k, 0) + v
    for each in results:
        barcode_counter[each[0]] = each[1]

    return barcode_counter


def reverse_complement(string):
    return string.translate(str.maketrans("ATGCNatgcn", "TACGNtacgn"))[::-1]


def get_indexes_of_approx_matches(query, list_of_strings, hamming_dist):
    """Returns a list containing *indexes* of matches to query in list_of_strings within hamming_dist.
    Since all strings must be the same length, hamming_dist is equivalent to the number of substitutions/differences between strings.
    Case insensitive.
    **This function does the heavy lifting for barcode analysis**
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
    """Wrapper function to call analyze barcode. Handles special cases with rc_mode flag
    Inputs:
        - barcode: extracted barcode from fastq header line (matches [ACTG]+\+[ACTG]+)
        - num_reads: number of reads with this barcode in fastq file
        - all_idx1, all_idx2, all_ids: lists of index1/index2/sample ids from input csv file. Each list must be in the same order
        - num_subs: number of substitutions allowed when matching barcode to supplied indexes
        - rc_mode: if True, try to match using the reverse complement of index 2 as well as the sequence in the input csv

    """
    idx1, idx2 = barcode.split("+")[0:2]

    # analyze barcode using supplied idx2
    temp = analyze_barcode(idx1, idx2, all_idx1, all_idx2, all_ids, num_subs)

    temp["reads"] = num_reads
    result = temp

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


def call_rc_mode_per_id(results_list, ids):
    """Given a list of dicts (format generated in frender_scan_function), for each sample id found, determine whether it should be demuxed with the forward or reverse complement index 2.
    The 'forward' (supplied) index 2 sequence is preferred if it results in an equal or greater number of demuxable reads compared to the reverse compelement index 2 sequence.
    Also, if forward index 2 results in exactly 0 demuxable sequences, it is assumed that those reads have already been taken out of the 'undetermined' file;
    in this case, the forward index 2 sequence will be used.
    Returns: a dictionary with each sample id and True (demux with rc index 2) or False (demux with forward index 2)
    """

    # Must have RC entries
    assert (
        "rc_read_type" in results_list[0].keys()
    ), "It looks like this frender result csv was not generated with the -rc flag. Either specify a different result csv, or run this command without setting the -rc flag."

    ids = {id: {"f": 0, "rc": 0, "demux_with_rc": ""} for id in ids}

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

    return {
        a: {
            "call": ids[a]["demux_with_rc"],
            "reads_f": ids[a]["f"],
            "reads_rc": ids[a]["rc"],
        }
        for a in ids
    }


def process(cores, barcode_counter, indexes, num_subs, rc_mode):
    all_idx1 = indexes["idx1"]
    all_idx2 = indexes["idx2"]
    all_ids = indexes["id"]
    if cores > 1:

        with Pool(processes=cores) as pool:
            print(f"Multiprocessing with {cores} cores")
            temp = pool.starmap(
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
            results = dict(zip(barcode_counter.keys(), temp))

    else:
        results = {
            barcode: analyze_barcodes_with_rc(
                barcode,
                barcode_counter[barcode],
                all_idx1,
                all_idx2,
                all_ids,
                num_subs,
                rc_mode,
            )
            for barcode in barcode_counter
        }
    return results


def report_rc_call_info(rc_calls, indexes, out_csv_name):
    rc_summary_file_name = out_csv_name.replace(
        "frender-scan-results_", "frender-index-2-calls_"
    )
    print(
        f"Based on the barcodes in the supplied fastq file, the following index 2 sequences will be used\n(also recorded in {rc_summary_file_name}):\n"
    )
    print(
        "Sample Name",
        "Supplied Index 2",
        "Reads supporting (forward)",
        "Reverse complement Index 2",
        "Reads supporting (rev comp)",
        "Final call",
        sep="\t",
    )
    for a in rc_calls:
        index = indexes["id"].index(a)
        print(
            a,
            indexes["idx2"][index],
            rc_calls[a]["reads_f"],
            reverse_complement(indexes["idx2"][index]),
            rc_calls[a]["reads_rc"],
            "reverse complement" if rc_calls[a]["call"] else "forward",
            sep="\t",
        )
    with open(rc_summary_file_name, "w", newline="") as output_file:
        writer = csv.writer(output_file)
        writer.writerow(
            [
                "sample_name",
                "supplied_index_2",
                "reads_supplied_index_2",
                "rc_index_2",
                "reads_rc_index_2",
                "use_rc",
            ]
        )
        for a in rc_calls:
            index = indexes["id"].index(a)
            writer.writerow(
                [
                    a,
                    indexes["idx2"][index],
                    rc_calls[a]["reads_f"],
                    reverse_complement(indexes["idx2"][index]),
                    rc_calls[a]["reads_rc"],
                    "TRUE" if rc_calls[a]["call"] else "FALSE",
                ]
            )


def flatten_results(dict_of_dicts):
    """Given a dict of dicts in the frender format, 'flatten' it out into a list of dicts."""
    result = []

    for barcode in dict_of_dicts.keys():
        d = {"idx1": barcode.split("+")[0], "idx2": barcode.split("+")[1]}
        for key, val in dict_of_dicts[barcode].items():
            d[key] = val
        result.append(d)

    return result


def report_analysis(results, out_csv_name):
    print(f"Analysis complete! Writing results to {out_csv_name}")
    keys = results[0].keys()
    with open(out_csv_name, "w", newline="") as output_file:
        dict_writer = csv.DictWriter(output_file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(results)


def call_barcodes_correctly_distributed(barcode_counter, results, prefix):
    files = list(barcode_counter.keys())
    files.remove("total")

    mismatching_files = set()

    for barcode in list(results.keys()):
        read_type = results[barcode]["read_type"]
        a = []
        for file in files:
            # How many reads of this barcode in this file?
            try:
                reads = barcode_counter[file][barcode]
            except KeyError:
                reads = 0

            # Should this barcode be in this file?
            if read_type == "undetermined":
                match = bool(
                    re.search(re.compile("undetermined", re.I), file)
                )  # True if filename matches 'undetermined'; undetermined reads should only be in the undetermined file.
            elif read_type == "index_hop":
                match = bool(
                    re.search(
                        re.compile("undetermined|index-hop", re.I),
                        file,
                    )
                )  # True if filename matches 'undetermined' or 'index-hop'; these reads belong in one of these two files
            elif read_type == "ambiguous":
                match = bool(
                    re.search(
                        re.compile("undetermined|ambiguous", re.I),
                        file,
                    )
                )  # similar to above
            else:
                assert (
                    read_type == "demuxable"
                ), f"Strange read type ('{read_type}') found"
                match = bool(
                    re.search(
                        re.compile(
                            results[barcode]["sample_name"].removeprefix(prefix), re.I
                        ),
                        file,
                    )
                )

            a.append(
                (not (bool(reads))) | match
            )  # ok if there are zero reads OR it is a match
            good_demux = bool(len(a) == sum(a))
            results[barcode]["demux_ok"] = good_demux

        [
            mismatching_files.add(files[index])
            for index, match in enumerate(a)
            if not match
        ]

    return (results, mismatching_files)


def frender_scan(args):
    # Parse args: n=1, rc=True, c=1.0, s=None, o='test_name', p=None, b=None, files=['file1', 'file2', 'file3']
    num_subs = args.n
    rc_mode = args.rc
    cores = get_cores(args.c)
    sample = args.s
    user_infix = args.o if args.o else ""
    prefix = args.p if args.p else ""

    # barcode table must be present unless we can find one in the provided directory
    if args.b == None:
        if len(args.files) != 1:
            raise SystemExit(
                "You have not specified a barcode table. Please either specify one with the argment -b or specify a directory including a barcode table"
            )
        barcode_file = find_barcode_file(Path(args.files[0]))
    else:
        barcode_file = Path(args.b)
    indexes = get_indexes(barcode_file)

    if len(args.files) == 1:
        file = Path(args.files[0])
        if Path.is_dir(file):
            files = {"dir": file}
            out_csv_name = f"frender-scan-results_{num_subs}-mismatches_{user_infix}_{file.parts[-1]}.csv"
        elif Path.is_file(file):
            files = {"file": file}
            out_csv_name = f"frender-scan-results_{num_subs}-mismatches_{user_infix}_{file.name}.csv"
        else:
            raise SystemExit("Specified directory or file path doesn't seem to exist!")
    else:
        files = {"file": [Path(file) for file in args.files]}
        out_csv_name = f"frender-scan-results_{num_subs}-mismatches_{user_infix}_{datetime.strftime(datetime.now(timezone.utc), '%Y-%M-%d_%H%M_%Z')}.csv"

    out_csv_name = out_csv_name.replace("__", "_")

    # Filter out non-fastq files (and Read 2 files, if scanning dir...)
    files = parse_files(files, just_r1=True)

    barcode_counter = tally_barcodes(cores, files, sample)
    print(
        f"Scanning complete! Analyzing barcodes...",
    )
    results = process(cores, barcode_counter["total"], indexes, num_subs, rc_mode)

    if rc_mode:
        # Have preliminary results at this point. Now we can call whether to use forward or rc index 2 for each sample...
        rc_calls = call_rc_mode_per_id(flatten_results(results), indexes["id"])
        print("First round of analysis complete.")
        report_rc_call_info(rc_calls, indexes, out_csv_name)

        indexes["idx2"] = [
            reverse_complement(indexes["idx2"][i])
            if rc_calls[id]["call"]
            else indexes["idx2"][i]
            for i, id in enumerate(indexes["id"])
        ]

        print(
            f"\nRe-analyzing barcodes with corrected index 2 sequences...",
        )
        results = process(
            cores, barcode_counter["total"], indexes, num_subs, rc_mode=False
        )

    results, mismatching_files = call_barcodes_correctly_distributed(
        barcode_counter, results, prefix
    )

    if bool(mismatching_files):
        print("Incorrectly demultiplexed barcodes found! Affected files:")
        {print(a) for a in mismatching_files}
    else:
        print("It appears that all files are already correctly demultiplexed.")

    report_analysis(flatten_results(results), out_csv_name)


def parse_results_file(result_file):
    with open(result_file, newline="") as f:
        header = next(csv.reader(f))

        assert header[0:7] == [
            "idx1",
            "idx2",
            "reads",
            "matched_idx1",
            "matched_idx2",
            "read_type",
            "sample_name",
        ], f"${result_file} does not appear to be a valid frender result file!"
        results_dict = {}
        for line in csv.reader(f):
            results_dict[line[0] + "+" + line[1]] = {
                "read_type": line[5],
                "sample_id": line[6],
            }
    return results_dict


def open_files(name, out_dir):
    if not out_dir.endswith("/"):
        out_dir += "/"
    return {
        read: gzip.open(
            f"{out_dir}{name}_frender-demux_{args.o+'_' if args.o else ''}{read}.fq.gz",
            "wb",
        )
        for read in ["R1", "R2"]
    }


def close_files(list_of_file_dicts):
    for file_dict in list_of_file_dicts:
        if file_dict:
            [file_dict[file].close() for file in file_dict]


def is_read_mate(str1, str2):
    if len([0 for a, b in zip(str1, str2) if a != b]) != 1:
        return False
    r1 = int(re.search("_R[12]_", str1)[0].replace("_", "").replace("R", ""))
    r2 = int(re.search("_R[12]_", str2)[0].replace("_", "").replace("R", ""))
    if {r1, r2} == {1, 2}:
        return True
    else:
        return False


def get_paired_files(files_list):
    r1_files = [
        path for path in files_list if bool(re.search("_R1_", str(path), re.IGNORECASE))
    ]

    pairs = []
    for path in r1_files:
        r2_file = [
            a
            for a, b in enumerate(
                [is_read_mate(str(path), str(file)) for file in files_list]
            )
            if b
        ]
        if len(r2_file) > 1:
            raise SystemExit(f"Found more than one potential read 2 file for {path}")
        elif len(r2_file) == 0:
            raise SystemExit(f"Couldn't find a read 2 file for {path}")
        else:
            pairs += [(path, files_list[r2_file[0]])]
    return pairs


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def write_reads(records_tuple, files_dict):
    [
        [files_dict[file].write(str(line).encode("utf-8")) for line in record]
        for record, file in zip(records_tuple, files_dict)
    ]


def frender_demux(args):
    # no_index_hop=True, no_ambiguous=True, no_undeter=False, no_samples=False, o='test_name', d = 'outdir', r='result_file.csv', files=['file1', 'file2', 'file3']

    index_hop = not args.no_index_hop
    ambiguous = not args.no_ambiguous
    undeter = not args.no_undeter
    samples = not args.no_samples

    undeter_name = f"Undetermined{'-ambiguous' if ambiguous else ''}{'-index-hop' if index_hop else ''}"

    result_file = Path(args.r)
    if not Path.is_file(result_file):
        raise SystemExit(f"File {result_file} not found")

    results_dict = parse_results_file(result_file)

    ids = list(set([results_dict[a]["sample_id"] for a in results_dict.keys()]) - {""})
    if (not ids) & samples:
        print(
            "Warning: no demuxable sample ids found in the supplied frender result file!"
        )

    os.mkdir(args.d)
    sample_files = {id: open_files(id, args.d) for id in ids} if samples else None
    undeter_files = open_files(undeter_name, args.d) if undeter else None
    index_hop_files = open_files("Index-hop", args.d) if index_hop else undeter_files
    ambiguous_files = open_files("Ambiguous", args.d) if ambiguous else undeter_files

    if len(args.files) == 1:
        file = Path(args.files[0])
        if Path.is_dir(file):
            files = {"dir": file}
        elif Path.is_file(file):
            files = {"file": file}
        else:
            raise SystemExit("Specified directory or file path doesn't seem to exist!")
    else:
        files = {"file": [Path(file) for file in args.files]}

    input_files = get_paired_files(parse_files(files, just_r1=False))

    for read1_file, read2_file in input_files:
        print(f"Demultiplexing {read1_file.name}...")
        with gzip.open(read1_file, "rt") as read1, gzip.open(read2_file, "rt") as read2:
            for record_tuple in zip(grouper(read1, 4, ""), grouper(read2, 4, "")):
                code = record_tuple[1][0].split(":")[-1].rstrip("\n")
                try:  # should happen every time (code is in dict)
                    if (results_dict[code]["read_type"] == "demuxable") & bool(
                        sample_files
                    ):
                        write_reads(
                            (record_tuple),
                            sample_files[results_dict[code]["sample_id"]],
                        )

                    elif (results_dict[code]["read_type"] == "index_hop") & bool(
                        index_hop_files
                    ):
                        write_reads((record_tuple), index_hop_files)

                    elif (results_dict[code]["read_type"] == "ambiguous") & bool(
                        ambiguous_files
                    ):
                        write_reads((record_tuple), ambiguous_files)

                    elif (results_dict[code]["read_type"] == "undetermined") & bool(
                        undeter_files
                    ):
                        write_reads((record_tuple), undeter_files)
                    else:
                        raise SystemExit(
                            "Unrecognized read type found in supplied frender result file!"
                        )

                except KeyError:
                    raise SystemExit(
                        f"Couldn't find barcode {code} in supplied frender result file!"
                    )

    # Close files:
    close_files([sample_files[id] for id in sample_files])
    close_files([index_hop_files, ambiguous_files, undeter_files])


if __name__ == "__main__":

    # create the top-level parser
    parser = argparse.ArgumentParser(prog="frender.py")

    subparsers = parser.add_subparsers()

    # create the parser for the "a" command
    parser_scan = subparsers.add_parser(
        "scan", help="Scan file(s) or directory and compare to a supplied barcode table"
    )
    parser_scan.add_argument(
        "-n",
        metavar="[int]",
        type=int,
        required=True,
        help="REQUIRED: Number of mismatches allowed between supplied barcodes and fastq file(s)",
    )
    parser_scan.add_argument(
        "-rc",
        action="store_true",
        help="Scan/demultiplex using reverse complement of index 2 as well as forward sequence (to check for mistakes with e.g. HiSeq 4000 and other systems)",
    )
    parser_scan.add_argument(
        "-c",
        metavar="cores",
        type=float,
        default=1,
        help="Number of cores to use for analysis, default = 1. Use 0 for all available, a number between 0 and 1 for a fraction of all available cores, or a number >= 1 for a specified number of cores",
    )
    parser_scan.add_argument(
        "-s",
        metavar="sample",
        type=int,
        help="If set, sample an absolute number of reads from the head of each file (s >= 1)",
    )
    parser_scan.add_argument(
        "-o",
        metavar="output_name",
        help="name infix for output files",
    )
    parser_scan.add_argument(
        "-p",
        metavar="fix_prefix",
        help="When matching sample ids to filenames, remove this prefix from the sample id",
    )
    parser_scan.add_argument(
        "-b",
        metavar="barcode_table",
        help=".csv formatted file containing barcode associations with ids. REQUIRED unless you specify a directory already containing such a file.",
    )
    parser_scan.add_argument(
        "files",
        nargs="+",
        help="Fastq file, list of fastq files, or directory path containing fastq files (subdirectories will be searched as well)",
    )
    parser_scan.set_defaults(func=frender_scan)

    # create the parser for the "b" command
    parser_demux = subparsers.add_parser(
        "demux",
        help="Demultiplex reads into sample and undetermined files according to supplied frender scan results file",
    )
    parser_demux.add_argument(
        "-i",
        "--no-index-hop",
        action="store_true",
        help="don't split index hop reads into their own file (will be included in undetermined file unless -u is set)",
    )
    parser_demux.add_argument(
        "-a",
        "--no-ambiguous",
        action="store_true",
        help="don't split ambiguous reads into their own file (will be included in undetermined file unless -u is set)",
    )
    parser_demux.add_argument(
        "-u",
        "--no-undeter",
        action="store_true",
        help="do NOT produce undetermined files",
    )
    parser_demux.add_argument(
        "-s",
        "--no-samples",
        action="store_true",
        help="do NOT produce individual sample files",
    )
    parser_demux.add_argument(
        "-o",
        metavar="output_name",
        help="name infix for output files",
    )
    parser_demux.add_argument(
        "-d",
        metavar="output_dir",
        default=f"./frender-demux-output_{datetime.strftime(datetime.now(timezone.utc), '%Y-%M-%d_%H%M_%Z')}/",
        help="output directory (default: ./frender-demux-output_{date_time}/)",
    )
    parser_demux.add_argument(
        "-r",
        metavar="result_file",
        required=True,
        help="REQUIRED: frender scan result file (typically named 'frender-scan-result_n-mismatches_{output infix or file/directory name}.csv')",
    )
    parser_demux.add_argument(
        "files",
        nargs="+",
        help="Fastq file, list of fastq files, or directory path containing fastq files (subdirectories will be searched as well)",
    )
    parser_demux.set_defaults(func=frender_demux)

    args = parser.parse_args()

    args.func(args)