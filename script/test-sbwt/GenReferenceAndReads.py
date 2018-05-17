from sys import argv
import os
import random
from GenRandomDnas import GenRandomDnas


def IsDna(c):
    return c == 'A' or c == 'C' or c == 'G' or c == 'T'


def Mutation(reads):
    d = {'A':"CGT", 'C':"AGT", "G":"ACT", "T":"ACG"}
    idxs = range(0, len(reads))
    random.shuffle(idxs)
	# max 7 mismatches
    num = random.randrange(0, 8)
    #num = random.randrange(1, 4)
    for i in xrange(0, num):
        i = idxs[i]
        c = reads[i]
        reads = reads[:i:] + random.choice(d[c]) + reads[i+1::]
    return reads, num


def generate_reads_kmer(file_name, kmer, size):
    """
	Assuming there exists only one read in fasta file.
    """

    size0 = size
    ref = ""
    with open(file_name) as mf:
        for ml in mf:
            if ml.startswith('>'):
                continue
            ref = ref + ml.rstrip(os.linesep).replace("N", "A")
    if len(ref) < kmer or size < kmer:
        exit(1)
    for i in xrange(0, size - kmer):
        line, num = Mutation(ref[i:i+kmer])
        print(">r"+str(i)+"-"+str(num)+"-"+str(i))
        print(line)


def CleanN_Fasta(file_name):
    line_width = 70
    line = ""

    with open(file_name) as mf:
        for ml in mf:
            if ml[0] == '>':
                if len(line) > 0:
                    print line
                    line = ""
                    print ml,
                else:
                    for c in ml:
                        if IsDna(c):
                            if len(line) == line_width:
                                print line
                                line = ""+c
                            else:
                                line += c
        if len(line) > 0:
            print line


def generate_reads_kmer_main():
    if len(argv) != 4:
        print ("usage: " + argv[0] + " [.fa] [length (150)] [size 1024*1024]" + os.linesep)
        exit()
    kmer = int(argv[2])
    size = int(argv[3])
    generate_reads_kmer(argv[1], kmer, size)


def print_hamming_weight(r, f):
    num_diff = 0
    line0 = ""

    c = 0
    for i in xrange(len(r)):
        if i % 32 == 0:
            line0 += str(c)
            c += 1
        else:
            line0 += " "

    line = ""

    for i, j in zip(f, r):
        if i != j:
            num_diff += 1
            line += "|"
        else:
            line += " "
    print(str(num_diff))
    print(line0)
    print(r)
    print(line)
    print(f)


def foo():
    f = ""
    with open("f1.fa") as mf:
        for ml in mf:
            if ml[0] == ">":
                continue
            else:
                f = ml[:-1:]
    rs = []
    with open("r1.fa") as mf:
        for ml in mf:
            if ml[0] == ">":
                continue
            else:
                rs.append(ml[:-1:])

    for i in xrange(len(rs)):
        print_hamming_weight(rs[i], f[i:150+i:])
        print


if __name__ == '__main__':
    generate_reads_kmer_main()

