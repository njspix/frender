import argparse
import os
import pandas as pd
import gzip
from itertools import zip_longest
import regex

parser = argparse.ArgumentParser()
parser.add_argument("-b", help="Barcode association table, csv format", required=True)
parser.add_argument("-o", help="output files directory", default=".")
parser.add_argument(
    "-p", help="prefix to be added to demuxed output files and index hop summary"
)
parser.add_argument(
    "-i", help="(optional) output fastq for index hopped reads", default=""
)
parser.add_argument(
    "-u", help="(optional) output fastq for remaining undetermined reads", default=""
)
parser.add_argument(
    "-c", help="(optional) output fastq for reads with conflicting barcodes", default=""
)
parser.add_argument("fastqs", nargs=2)

args = parser.parse_args()

# read barcode association file:
indexes = pd.read_csv(args.b, names=["id", "idx1", "idx2"]).set_index("id")
all_idx1 = indexes["idx1"].tolist()
all_idx2 = indexes["idx2"].tolist()
ids = list(indexes.index)
# todo: handle single index
# todo: check that barcodes are all same length, match [ACTGactg]

# create output dir
args.o = args.o.rstrip("/") + "/"  # add trailing slash just to make sure
# todo: what if this doesn't start with './'?
try:
    os.makedirs(args.o)
except OSError as e:
    if e.errno != 17:  # 17 = file exists
        raise

if args.p != "":
    args.p = args.p + "_"

r1_files = []
r2_files = []

# open output files:
for i in range(len(ids)):
    r1_files.append(gzip.open(f"{args.o}{args.p}{ids[i]}_R1.fq.gz", "wb"))
    r2_files.append(gzip.open(f"{args.o}{args.p}{ids[i]}_R2.fq.gz", "wb"))

if args.u != "":
    r1_undeter = gzip.open(f"{args.o}{args.u}_R1.fq.gz", "wb")
    r2_undeter = gzip.open(f"{args.o}{args.u}_R2.fq.gz", "wb")

if args.i != "":
    r1_hop = gzip.open(f"{args.o}{args.i}_R1.fq.gz", "wb")
    r2_hop = gzip.open(f"{args.o}{args.i}_R2.fq.gz", "wb")

if args.c != "":
    r1_conflict = gzip.open(f"{args.o}{args.c}_R1.fq.gz", "wb")
    r2_conflict = gzip.open(f"{args.o}{args.c}_R2.fq.gz", "wb")


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def fuz_match_list(pattern, set_of_strings):
    "Given a query string and a list of strings, return a list of indices where a match (+/- 1 substitution) is found"
    pattern = regex.compile("(?:" + pattern + "){s<=1}", regex.IGNORECASE)
    matches = [
        bool(regex.match(pattern, set_of_strings[i]))
        for i in range(len(set_of_strings))
    ]
    matches = [i for i, val in enumerate(matches) if val]
    return matches


hops = pd.DataFrame()

with gzip.open(args.fastqs[0], "rt") as read1, gzip.open(args.fastqs[1], "rt") as read2:

    reads_1 = grouper(read1, 4, "")
    reads_2 = grouper(read2, 4, "")

    barcode_dict = {}

    for record_1, record_2 in zip(reads_1, reads_2):
        assert len(record_1) == 4
        # todo: check that read 1 and read 2 are identical (add switch for this?)
        idx1 = record_1[0].split(":")[-1].split("+")[0].rstrip("\n")
        idx2 = record_1[0].split(":")[-1].split("+")[1].rstrip("\n")

        if idx1 + "+" + idx2 in barcode_dict:  # we have seen this combination before...
            if (barcode_dict[idx1 + "+" + idx2] == "undetermined") & (args.u != ""):
                for line in record_1:
                    r1_undeter.write(str(line).encode("utf-8"))
                for line in record_2:
                    r2_undeter.write(str(line).encode("utf-8"))

            elif (barcode_dict[idx1 + "+" + idx2] == "hop") & (args.i != ""):
                for line in record_1:
                    r1_hop.write(str(line).encode("utf-8"))
                for line in record_2:
                    r2_hop.write(str(line).encode("utf-8"))

            elif (barcode_dict[idx1 + "+" + idx2] == "conflict") & (args.c != ""):
                for line in record_1:
                    r1_conflict.write(str(line).encode("utf-8"))
                for line in record_2:
                    r2_conflict.write(str(line).encode("utf-8"))

            elif barcode_dict[idx1 + "+" + idx2] not in [
                "hop",
                "undetermined",
                "conflict",
            ]:
                for line in record_1:
                    r1_files[barcode_dict[idx1 + "+" + idx2]].write(
                        str(line).encode("utf-8")
                    )
                for line in record_2:
                    r2_files[barcode_dict[idx1 + "+" + idx2]].write(
                        str(line).encode("utf-8")
                    )

        else:  # need to process this read
            idx1_matches = fuz_match_list(idx1, all_idx1)
            idx2_matches = fuz_match_list(idx2, all_idx2)

            if bool(idx1_matches) & bool(
                idx2_matches
            ):  # can find at least one match for both idx1 and idx2 somewhere

                match_intersection = set(idx1_matches).intersection(idx2_matches)

                if len(match_intersection) == 0:  # this is an index hop

                    barcode_dict[idx1 + "+" + idx2] = "hop"  # add to known barcodes

                    # record in index hop report:
                    try:
                        hops.loc[
                            set([all_idx1[i] for i in idx1_matches]).pop(),
                            set([all_idx2[i] for i in idx2_matches]).pop(),
                        ] += 1
                        # this is a bit of a hack. Could fail if there is a 'perfectly bad' index read that matches > 1 index because it is 1 sub away from both of them.
                    except KeyError:  # this combination of indices hasn't been initialized
                        hops.loc[
                            set([all_idx1[i] for i in idx1_matches]).pop(),
                            set([all_idx2[i] for i in idx2_matches]).pop(),
                        ] = 1

                    # optional write to output file
                    if args.i != "":
                        for line in record_1:
                            r1_hop.write(str(line).encode("utf-8"))
                        for line in record_2:
                            r2_hop.write(str(line).encode("utf-8"))

                elif (
                    len(match_intersection) == 1
                ):  # good read; idx1 and idx2 line up in exactly one spot

                    demux_id = indexes.index[match_intersection.pop()]

                    barcode_dict[idx1 + "+" + idx2] = indexes.index.get_loc(demux_id)

                    for line in record_1:
                        r1_files[indexes.index.get_loc(demux_id)].write(
                            str(line).encode("utf-8")
                        )
                    for line in record_2:
                        r2_files[indexes.index.get_loc(demux_id)].write(
                            bytes(str(line).encode("utf-8"))
                        )

                else:  # something is weird, we have more than one possible output file this read could belong to

                    barcode_dict[idx1 + "+" + idx2] = "conflict"

                    if args.c != "":
                        for line in record_1:
                            r1_conflict.write(str(line).encode("utf-8"))
                        for line in record_2:
                            r2_conflict.write(str(line).encode("utf-8"))

            else:  # can't match both barcodes

                barcode_dict[idx1 + "+" + idx2] = "undetermined"

                if args.u != "":
                    for line in record_1:
                        r1_undeter.write(str(line).encode("utf-8"))
                    for line in record_2:
                        r2_undeter.write(str(line).encode("utf-8"))

# close output files:
for i in range(len(r1_files)):
    r1_files[i].close()
for i in range(len(r2_files)):
    r2_files[i].close()

if args.u != "":
    r1_undeter.close()
    r2_undeter.close()

if args.i != "":
    r1_hop.close()
    r2_hop.close()

if args.c != "":
    r1_conflict.close()
    r2_conflict.close()