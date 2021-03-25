def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # last is empty: the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:  # if last is empty
            break
        name, seqs, last = (last[1:].partition(" ")[0], [], None)  # get sequence name, create empty list for sequence, and reset 'last' - we've read all the information out of this (header) line
        for l in fp:  # get the next line, twice: once for the sequence, and once for the separator
            if l[0] in "@+>": # are we reading the sequence line or a separator/header?
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last:  # if last is empty
            break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []  #join multiple sequences together 
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, "".join(seqs)
                    # yield a fastq record
                    break
            if last:  # last is not empty: reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def read_fastq_barcodes(file):
    last = None
    while True:
        if not last:
            for line in file:
                if line[0] in ">@":
                    last = line[:-1]
                    break
        yield last.split(" ")[-1].split(":")[-1]
        last = None


# for name, seq, qual in readfq(open("test.fastq")):
#     print(seq)

for barcode in read_fastq_barcodes(open("test.fastq")):
    print(barcode)
