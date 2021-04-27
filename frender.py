"""
Fastq REad Name DemultiplexER

TODO: ADD ONE SENTENCE EXPLANATION

Creator:
    Nathan Spix

Date Created:
    September 2020

Requires:
    pandas 1.2.2
    regex 2020.11.13

Inputs:
    TODO: ADD DESCRIPTIONS OF INPUTS

OUTPUT/RETURNS:
    TODO: ADD DESCRIPTION OF OUTPUT FILES
"""
import argparse
from pathlib import Path
import pandas as pd
import gzip
from itertools import zip_longest, islice
import sys
from time import perf_counter

try:
    import regex
except ModuleNotFoundError:
    print("ERROR: Install regex module", file=sys.stderr)
    print(
        "\tpip install regex\tOR\tconda install -c conda-forge regex", file=sys.stderr
    )
    sys.exit(0)

counter_interval = 250000

# Print iterations progress
def printProgressBar(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=100,
    fill="█",
    printEnd="\r",
):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


def read_scan_results(csv_file):
    """Read a csv file produced by the frender_scan function and parse into an internally useful format.
    Inputs: csv file with (at least) the following columns: read type, idx1, idx2, sample_name
    Returns: tuple containing:
        - dict in the following format: {'CCTGTACT+ANACATCG': 'index_hop', 'CGCTACAG+ANACATCG': 'FT-SA81396'}
          Possible values for 'value' are: any sample name in barcodeAssociationTable, 'index_hop', 'ambiguous', 'undetermined'
        - list containing all sample names
    """
    info = pd.read_csv(csv_file)
    info = (
        info.assign(idx=info.idx1 + "+" + info.idx2)
        .loc[:, ("idx", "read_type", "sample_name")]
        .set_index("idx")
    )
    info["sample_name"] = info.apply(
        lambda x: x["read_type"] if pd.isnull(x["sample_name"]) else x["sample_name"],
        axis=1,
    )
    return (
        info.to_dict()["sample_name"],
        list(
            set(
                info[info.read_type == "demuxable"]
                .loc[:, "sample_name"]
                .values.flatten()
                .tolist()
            )
        ),
    )


def process_index(k):
    return tuple(k.split("+"))


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks.
       Ex. grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"

    Inputs -
        iterable  - string to group in blocks of length n
        n         - block size
        fillvalue - if block sizes of uneven length, fill missing values with
                    fillvalue
    Returns -
        iter object with fixed-length blocks
    """
    args = [iter(iterable)] * n

    return zip_longest(*args, fillvalue=fillvalue)


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


complement = str.maketrans("ATGCN", "TACGN")


def reverse_complement(string):
    return string.translate(complement)[::-1]


def write_read_to_file(checkvar, records, files):
    """Write records to files iff checkvar != "" """
    if checkvar != "":
        if isinstance(files, (tuple, list)):
            assert len(records) == len(files)
            for record, file in zip(records, files):
                for line in record:
                    file.write(str(line))
        else:
            for line in records:
                files.write(str(line))


def analyze_barcode(
    idx1, idx2, all_indexes, barcode_counts_df, num_subs, rc_flag=False
):

    orig_idx2 = reverse_complement(idx2) if rc_flag else idx2

    all_idx1 = all_indexes["Index1"].tolist()
    all_idx2 = all_indexes["Index2"].tolist()

    idx1_matches = fuz_match_list(idx1, all_idx1, num_subs)
    idx2_matches = fuz_match_list(idx2, all_idx2, num_subs)

    if bool(idx1_matches) and bool(idx2_matches):
        # Can find at least one barcode match for both indices
        match_isec = set(idx1_matches).intersection(idx2_matches)

        if len(match_isec) == 0:
            # this is an index hop
            matched_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
            matched_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()

            update_barcode_counts_df(
                barcode_counts_df,
                (idx1, orig_idx2),
                (matched_idx1, matched_idx2, "index_hop", "", rc_flag),
            )

        elif len(match_isec) == 1:
            # this is a good read
            matched_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
            matched_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
            sample_name = all_indexes.index[match_isec.pop()]

            update_barcode_counts_df(
                barcode_counts_df,
                (idx1, orig_idx2),
                (matched_idx1, matched_idx2, "demuxable", sample_name, rc_flag),
            )

        else:
            # this is an ambiguous read
            matched_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
            matched_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()

            update_barcode_counts_df(
                barcode_counts_df,
                (idx1, orig_idx2),
                (matched_idx1, matched_idx2, "ambiguous", "", rc_flag),
            )

    else:
        update_barcode_counts_df(
            barcode_counts_df,
            (idx1, orig_idx2),
            ("", "", "undetermined", "", rc_flag),
        )


def update_barcode_counts_df(df, indices, values):
    """Update the dataframe df at the location specified by tuple of indices (idx1, idx2) with the tuple of values
    (matched_idx1, matched_idx2, read_type, sample_name, idx2_is_reverse_complement)
    """
    if not values[4]:  # rc_flag is not set
        df.loc[
            indices,
            [
                "matched_idx1",
                "matched_idx2",
                "read_type",
                "sample_name",
            ],
        ] = values[0:4]
    else:
        df.loc[
            indices,
            [
                "matched_idx1",
                "matched_rc_idx2",
                "rc_read_type",
                "rc_sample_name",
            ],
        ] = values[0:4]


def frender(
    barcode,
    fastq_1,
    fastq_2,
    out_dir=".",
    preefix="",
    ihopped="",
    cnflict="",
    undeter="",
    num_subs=1,
):
    """Main module function.

    Inputs -
        barcode - CSV file containing barcodes
        fastq_1 - Read 1 FASTQ file
        fastq_2 - Read 2 FASTQ file
        out_dir - Output directory name [default: '.']
        preefix - Prefix to add to output files [default: '']
        ihopped - Name of FASTQ file to save index hopped reads to
                  [default: '']
        cnflict - Name of FASTQ file to save conflicting reads to
                  [default: '']
        undeter - Name of FASTQ file to save undetermined reads to
                  [default: '']
        num_subs - number of substitutions allowed in match
    Returns -
        Creates demultiplexed FASTQ files from fastq_1 and fastq_2 and a summary
        of the index hopped reads. Will also create index hopped FASTQ file,
        conflicting barcode FASTQ file, and undetermined FASTQ file if supplied
        during function call.
    """
    barcode_dict = None

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

    if test.values.flatten().tolist() == [
        "idx1",
        "idx2",
        "total_reads",
        "matched_idx1",
        "matched_idx2",
        "read_type",
        "sample_name",
        "idx2_is_reverse_complement",
    ]:
        print("frender scan csv file detected!", file=sys.stderr)
        barcode_dict, ids = read_scan_results(barcode)
        using_scan_results = True

    else:
        barcode_dict = {}
        using_scan_results = False

        # TODO: Handle single index
        # TODO: Check barcodes are all the same length, match [ATCGatcg]
        if test[0][0].startswith("Acc") or test[0][0].startswith("Client"):
            indexes = pd.read_csv(barcode, usecols=cols, index_col="AccessionID")
        else:
            indexes = pd.read_csv(
                barcode, header=None, names=cols, index_col="AccessionID"
            )

        ids = list(indexes.index)

        # these are only created if we are working from scratch (no scan results). They are not needed if we have scan results, because we'll never have to process a barcode
        all_idx1 = indexes["Index1"].tolist()
        all_idx2 = indexes["Index2"].tolist()
        # DataFrame to store index hopping information
        hops = pd.DataFrame()

    # Create output directory and prep prefix
    out_dir = out_dir.rstrip("/") + "/"  # Add trailing slash
    if out_dir != "./":
        Path(out_dir).mkdir(parents=True, exist_ok=True)

    if preefix != "" and not preefix.endswith("_"):
        preefix = preefix + "_"

    # Open output files
    r1_files = {}
    r2_files = {}

    for id in ids:
        r1_files[id] = open(f"{out_dir}{preefix}{id}_R1.fastq", "w")
        r2_files[id] = open(f"{out_dir}{preefix}{id}_R2.fastq", "w")

    if ihopped != "":
        r1_hop = open(f"{out_dir}{preefix}{ihopped}_R1.fastq", "w")
        r2_hop = open(f"{out_dir}{preefix}{ihopped}_R2.fastq", "w")
    if cnflict != "":
        r1_con = open(f"{out_dir}{preefix}{cnflict}_R1.fastq", "w")
        r2_con = open(f"{out_dir}{preefix}{cnflict}_R2.fastq", "w")
    if undeter != "":
        r1_und = open(f"{out_dir}{preefix}{undeter}_R1.fastq", "w")
        r2_und = open(f"{out_dir}{preefix}{undeter}_R2.fastq", "w")

    # Start processing data
    with gzip.open(fastq_1, "rt") as read1, gzip.open(fastq_2, "rt") as read2:
        reads_1 = grouper(read1, 4, "")
        reads_2 = grouper(read2, 4, "")

        record_count = 0
        for record_1, record_2 in zip(reads_1, reads_2):
            assert (len(record_1) == 4) and (len(record_2) == 4)
            if record_count % counter_interval == 1:
                print(f"Processed {record_count} reads", file=sys.stderr)
            record_count += 1

            # Parse record header line to compare name and extract barcode
            r1_head = record_1[0].rstrip("\n").split(" ")
            r2_head = record_2[0].rstrip("\n").split(" ")

            if r1_head[0] != r2_head[0]:
                print("WARNING: Skipping mismatched reads", end=" ", file=sys.stderr)
                print(f"{r1_head[0]} != {r2_head[0]} at line", end=" ", file=sys.stderr)
                print(f"{4*record_count} of {fastq_1} and {fastq_2}", file=sys.stderr)
                continue

            # Barcode indices
            code = r1_head[1].split(":")[-1]  # full index name

            if not using_scan_results:
                idx1 = code.split("+")[0]  # index 1 in full name
                idx2 = code.split("+")[1]  # index 2 in full name

            if code in barcode_dict:  # seen this combination before
                # could be an index hop:
                if barcode_dict[code] == "index_hop":

                    # write to files if requested
                    write_read_to_file(ihopped, (record_1, record_2), (r1_hop, r2_hop))

                    # update hop count table
                    if not using_scan_results:
                        idx1_matches = fuz_match_list(idx1, all_idx1, num_subs)
                        idx2_matches = fuz_match_list(idx2, all_idx2, num_subs)
                        i_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                        i_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                        hops.loc[i_idx1, i_idx2] += 1

                # could be a conflict:
                elif barcode_dict[code] == "ambiguous":
                    write_read_to_file(cnflict, (record_1, record_2), (r1_con, r2_con))

                # could be unkonwn:
                elif barcode_dict[code] == "undetermined":
                    write_read_to_file(undeter, (record_1, record_2), (r1_und, r2_und))

                # or it could be a valid demux:
                elif barcode_dict[code] not in [
                    "index_hop",
                    "ambiguous",
                    "undetermined",
                ]:
                    write_read_to_file(
                        True,
                        (record_1, record_2),
                        (r1_files[barcode_dict[code]], r2_files[barcode_dict[code]]),
                    )

            else:  # need to process this read; should never get here if using_scan_results
                assert not using_scan_results
                idx1_matches = fuz_match_list(idx1, all_idx1, num_subs)
                idx2_matches = fuz_match_list(idx2, all_idx2, num_subs)

                if bool(idx1_matches) and bool(idx2_matches):
                    # Can find at least one barcode match for both indices
                    match_isec = set(idx1_matches).intersection(idx2_matches)

                    if len(match_isec) == 0:
                        # No matches, so index hop
                        barcode_dict[code] = "index_hop"  # add to known barcodes

                        # This is currently imprecise. A 'hop' could match more than one barcode in either index,
                        # idx1_matches = [4] and idx2_matches = [2, 3]. In this case, the matched indices could
                        # be identical, or they could both be 1 substitution away from the target (bad index design, but that's not my problem here :-)
                        # This doesn't account for that second case.
                        i_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                        i_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                        try:
                            hops.loc[i_idx1, i_idx2] += 1
                        except KeyError:
                            # Combination of indices hasn't been initialized
                            hops.loc[i_idx1, i_idx2] = 1

                        # Write to output file (if applicable)
                        if ihopped != "":
                            for line in record_1:
                                r1_hop.write(str(line))
                            for line in record_2:
                                r2_hop.write(str(line))
                    elif len(match_isec) == 1:
                        # good read; idx1 and idx2 line up in exactly one spot
                        demux_id = indexes.index[match_isec.pop()]

                        barcode_dict[code] = demux_id
                        for line in record_1:
                            r1_files[barcode_dict[code]].write(str(line))
                        for line in record_2:
                            r2_files[barcode_dict[code]].write(str(line))
                    else:
                        # Read matches to more than one possible output file
                        barcode_dict[code] = "ambiguous"
                        if cnflict != "":
                            for line in record_1:
                                r1_con.write(str(line))
                            for line in record_2:
                                r2_con.write(str(line))
                else:
                    # Can't find a match for at least one barcode
                    barcode_dict[code] = "undetermined"
                    write_read_to_file(undeter, (record_1, record_2), (r1_und, r2_und))

    # Close output files
    for key in r1_files:
        r1_files[key].close()
    for key in r2_files:
        r2_files[key].close()

    if ihopped != "":
        r1_hop.close()
        r2_hop.close()
    if cnflict != "":
        r1_con.close()
        r2_con.close()
    if undeter != "":
        r1_und.close()
        r2_und.close()

    # Write index hopping report
    if not using_scan_results:
        hops = hops.reset_index().melt(
            id_vars="index", var_name="idx2", value_name="num_hops_observed"
        )
        hops = hops[hops["num_hops_observed"] > 0]
        hops.columns = ["idx1", "idx2", "num_hops_observed"]
        hops.to_csv(f"{out_dir}{preefix}barcode_hops.csv", index=False)


def frender_se(
    barcode,
    fastq,
    out_dir=".",
    preefix="",
    ihopped="",
    cnflict="",
    undeter="",
    num_subs=1,
):
    """Demultiplexes a single-end read FASTQ file based on barcodes in read name.

    Inputs -
        barcode - CSV file containing barcodes
        fastq - FASTQ file
        out_dir - Output directory name [default: '.']
        preefix - Prefix to add to output files [default: '']
        ihopped - Name of FASTQ file to save index hopped reads to
                  [default: '']
        cnflict - Name of FASTQ file to save conflicting reads to
                  [default: '']
        undeter - Name of FASTQ file to save undetermined reads to
                  [default: '']
        num_subs - Number of substitutions allowed
    Returns -
        Creates demultiplexed FASTQ files and a summary
        of the index hopped reads. Will also create index hopped FASTQ file,
        conflicting barcode FASTQ file, and undetermined FASTQ file if supplied
        during function call.
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

    if test.values.flatten().tolist() == [
        "idx1",
        "idx2",
        "total_reads",
        "matched_idx1",
        "matched_idx2",
        "read_type",
        "sample_name",
        "idx2_is_reverse_complement",
    ]:
        print("frender scan csv file detected!", file=sys.stderr)
        barcode_dict, ids = read_scan_results(barcode)
        using_scan_results = True

    else:
        barcode_dict = {}
        using_scan_results = False

        # TODO: Handle single index
        # TODO: Check barcodes are all the same length, match [ATCGatcg]
        if test[0][0].startswith("Acc") or test[0][0].startswith("Client"):
            indexes = pd.read_csv(barcode, usecols=cols, index_col="AccessionID")
        else:
            indexes = pd.read_csv(
                barcode, header=None, names=cols, index_col="AccessionID"
            )

        ids = list(indexes.index)

        # these are only created if we are working from scratch (no scan results). They are not needed if we have scan results, because we'll never have to process a barcode
        all_idx1 = indexes["Index1"].tolist()
        all_idx2 = indexes["Index2"].tolist()
        # DataFrame to store index hopping information
        hops = pd.DataFrame()

    # Create output directory and prep prefix
    out_dir = out_dir.rstrip("/") + "/"  # Add trailing slash
    if out_dir != "./":
        Path(out_dir).mkdir(parents=True, exist_ok=True)

    if preefix != "" and not preefix.endswith("_"):
        preefix = preefix + "_"

    # Open output files
    r1_files = {}

    for id in ids:
        r1_files[id] = open(f"{out_dir}{preefix}{id}.fastq", "w")

    if ihopped != "":
        r1_hop = open(f"{out_dir}{preefix}{ihopped}.fastq", "w")
    if cnflict != "":
        r1_con = open(f"{out_dir}{preefix}{cnflict}.fastq", "w")
    if undeter != "":
        r1_und = open(f"{out_dir}{preefix}{undeter}.fastq", "w")

    # Start processing data
    with gzip.open(fastq, "rt") as read1:
        reads = grouper(read1, 4, "")

        record_count = 0
        for record in reads:
            assert len(record) == 4
            if record_count % counter_interval == 1:
                print(f"Processed {record_count} reads", file=sys.stderr)
            record_count += 1

            # Parse record header line to compare name and extract barcode
            r1_head = record[0].rstrip("\n").split(" ")

            # Barcode indices
            code = r1_head[1].split(":")[-1]  # full index name

            if not using_scan_results:
                idx1 = code.split("+")[0]  # index 1 in full name
                idx2 = code.split("+")[1]  # index 2 in full name

            if code in barcode_dict:  # seen this combination before
                # could be an index hop:
                if barcode_dict[code] == "index_hop":

                    # write to files if requested
                    write_read_to_file(ihopped, record, r1_hop)

                    # update hop count table
                    if not using_scan_results:
                        idx1_matches = fuz_match_list(idx1, all_idx1, num_subs)
                        idx2_matches = fuz_match_list(idx2, all_idx2, num_subs)
                        i_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                        i_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                        hops.loc[i_idx1, i_idx2] += 1

                # could be a conflict:
                elif barcode_dict[code] == "ambiguous":
                    write_read_to_file(cnflict, record, r1_con)

                # could be unkonwn:
                elif barcode_dict[code] == "undetermined":
                    write_read_to_file(undeter, record, r1_und)

                # or it could be a valid demux:
                elif barcode_dict[code] not in [
                    "index_hop",
                    "ambiguous",
                    "undetermined",
                ]:
                    write_read_to_file(True, record, r1_files[barcode_dict[code]])

            else:  # need to process this read
                assert not using_scan_results
                idx1_matches = fuz_match_list(idx1, all_idx1, num_subs)
                idx2_matches = fuz_match_list(idx2, all_idx2, num_subs)

                if bool(idx1_matches) and bool(idx2_matches):
                    # Can find at least one barcode match for both indices
                    match_isec = set(idx1_matches).intersection(idx2_matches)

                    if len(match_isec) == 0:
                        # No matches, so index hop
                        barcode_dict[code] = "index_hop"  # add to known barcodes

                        # This is currently imprecise. A 'hop' could match more than one barcode in either index,
                        # idx1_matches = [4] and idx2_matches = [2, 3]. In this case, the matched indices could
                        # be identical, or they could both be 1 substitution away from the target (bad index design, but that's not my problem here :-)
                        # This doesn't account for that second case.
                        i_idx1 = set([all_idx1[i] for i in idx1_matches]).pop()
                        i_idx2 = set([all_idx2[i] for i in idx2_matches]).pop()
                        try:
                            hops.loc[i_idx1, i_idx2] += 1
                        except KeyError:
                            # Combination of indices hasn't been initialized
                            hops.loc[i_idx1, i_idx2] = 1

                        # Write to output file (if applicable)
                        write_read_to_file(ihopped, record, r1_hop)

                    elif len(match_isec) == 1:
                        # good read; idx1 and idx2 line up in exactly one spot
                        barcode_dict[code] = indexes.index[match_isec.pop()]
                        write_read_to_file(True, record, r1_files[barcode_dict[code]])

                    else:
                        # Read matches to more than one possible output file
                        barcode_dict[code] = "ambiguous"
                        write_read_to_file(cnflict, record, r1_con)

                else:
                    # Can't find a match for at least one barcode
                    barcode_dict[code] = "undetermined"
                    write_read_to_file(undeter, record, r1_und)

    # Close output files
    for key in r1_files:
        r1_files[key].close()

    if ihopped != "":
        r1_hop.close()
    if cnflict != "":
        r1_con.close()
    if undeter != "":
        r1_und.close()

    # Write index hopping report
    if not using_scan_results:
        hops = hops.reset_index().melt(
            id_vars="index", var_name="idx2", value_name="num_hops_observed"
        )
        hops = hops[hops["num_hops_observed"] > 0]
        hops.columns = ["idx1", "idx2", "num_hops_observed"]
        hops.to_csv(f"{out_dir}{preefix}barcode_hops.csv", index=False)


def frender_scan(rc_mode, barcode, fastq_1, out_dir=".", preefix="", num_subs=1):
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

    # Create output directory and prep prefix
    out_dir = out_dir.rstrip("/") + "/"  # Add trailing slash
    if out_dir != "./":
        Path(out_dir).mkdir(parents=True, exist_ok=True)

    if preefix != "" and not preefix.endswith("_"):
        preefix = preefix + "_"

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
        f"Scanning complete! Analyzing {len(barcode_counter)} barcodes...",
        file=sys.stderr,
    )
    # Turn the barcode_counter dict into a pd dataframe
    barcode_counts = pd.DataFrame.from_dict(
        barcode_counter, orient="index", columns=["total_reads"]
    )
    # Turn the single index into a pd MultiIndex (split barcodes on "+")
    barcode_counts["idx1"], barcode_counts["idx2"] = zip(
        *map(process_index, barcode_counts.index)
    )
    barcode_counts.set_index(["idx1", "idx2"], inplace=True)

    all_found_indexes = barcode_counts.index.values
    barcode_count = 0

    for i in all_found_indexes:
        if barcode_count % (counter_interval / 2500) == 1:
            printProgressBar(
                barcode_count,
                len(barcode_counter),
                prefix="Progress:",
                suffix="complete",
                length=50,
            )
        barcode_count += 1

        analyze_barcode(
            i[0], i[1], indexes, barcode_counts, num_subs=num_subs, rc_flag=False
        )

        if rc_mode:
            analyze_barcode(
                i[0],
                reverse_complement(i[1]),
                indexes,
                barcode_counts,
                num_subs=num_subs,
                rc_flag=True,
            )

        # Write report
        print(barcode_counts.to_csv())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Demultiplex Undetermined FastQ files with given barcodes."
    )
    parser.add_argument(
        "-o",
        metavar="output_dir",
        default=".",
        help="output files directory [default: '.']",
    )
    parser.add_argument(
        "-p",
        metavar="file_prefix",
        default="",
        help="prefix added to output files [default: " "]",
    )
    parser.add_argument(
        "-i",
        metavar="index_hopped.fq",
        default="",
        help="output fastq for index hopped reads [default: " "]",
    )
    parser.add_argument(
        "-c",
        metavar="conflicting.fq",
        default="",
        help="output fastq for reads with conflicting barcodes [default: " "]",
    )
    parser.add_argument(
        "-u",
        metavar="undetermined.fq",
        default="",
        help="output fastq for remaining undetermined reads [default: " "]",
    )
    parser.add_argument(
        "-b",
        metavar="barcode.txt",
        help="Barcode association table, csv format",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--numsubs",
        metavar="num_subs",
        help="Number of substitutions allowed in barcode when matching",
        type=int,
        required=True,
    )
    parser.add_argument(
        "-s",
        "--scan",
        action="store_true",
        help="""Scan a single fastq file, counting exact and inexact barcode matches, conflicting barcodes, index hops, and undetermined reads. 
                No demultiplexing is performed. 
                ****** 
                Streams CSV formatted output (idx1, idx2, total_reads, matched_idx1, matched_idx2,read_type, sample_name, idx2_is_reverse_complement)
                to stdout (make sure to redirect to a file!) 
                ******
             """,
    )
    parser.add_argument(
        "-rc",
        action="store_true",
        help="When used with --scan, also scan for reverse complement of index 2 (to check for mistakes with e.g. HiSeq 4000 and other systems)",
    )
    parser.add_argument(
        "fastqs",
        nargs=argparse.REMAINDER,
        help="""Gzipped fastq files to be scanned or demultiplexed. If -s is specified, only one fastq is used. 
                  If -s is NOT specified, and one fastq is supplied, one set of fastqs will be produced (single-end mode)
                  If -s is NOT specified, and two fastqs are supplied, two sets of fastqs will be produced (R1 and R2; paired-end mode)
               """,
    )

    args = parser.parse_args()

    if args.scan:
        rc_mode_text = (
            "both supplied and reverse-complement index 2 sequences"
            if args.rc
            else "index 2 sequences as supplied"
        )
        print(f"Scanning {args.fastqs[0]} using {rc_mode_text}...", file=sys.stderr)
        frender_scan(args.rc, args.b, args.fastqs[0], args.o, args.p, args.numsubs)

    else:
        # single end case
        if len(args.fastqs) == 1:
            print("Only one input fastq detected. Running in single-end mode...")
            frender_se(
                args.b,
                args.fastqs[0],
                args.o,
                args.p,
                args.i,
                args.c,
                args.u,
                args.numsubs,
            )

        # paired-end case
        elif len(args.fastqs) == 2:
            print("Two input fastqs detected. Running in paired-end mode...")
            frender(
                args.b,
                args.fastqs[0],
                args.fastqs[1],
                args.o,
                args.p,
                args.i,
                args.c,
                args.u,
                args.numsubs,
            )
        else:
            raise TypeError(
                "Wrong number of fastq files provided. Please specify only one or two fastq files to use as input!"
            )

# TODO:
# Have an option to proceed with standard frender() if any un-demultiplexed reads are found.
