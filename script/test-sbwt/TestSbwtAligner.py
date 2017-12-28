from sys import argv


def IsDna(c):
        return c == 'A' or c == 'C' or c == 'G' or c == 'T'


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


if __name__ == '__main__':
    if len(argv) < 2:
        print "usage: ", argv[0], " [.fa] [length (150)]"
        exit()
        kmer = 150
        if len(argv) > 2:
            kmer = int(argv[2])
        CleanN_Fasta(argv[1])
