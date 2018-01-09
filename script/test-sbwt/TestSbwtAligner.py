from sys import argv
import random



def IsDna(c):
        return c == 'A' or c == 'C' or c == 'G' or c == 'T'


def Mutation(reads):
	d = {'A':"CGT", 'C':"AGT", "G":"ACT", "T":"ACG"}
	idxs = range(0, len(reads))
	num = random.randrange(0, len(reads)+1)
	for i in xrange(0, num):
		c = reads[i]
		reads = reads[:i:] + random.choice(d[c]) + reads[i+1::]
	return reads, num


def GenerateReadsKmer(file_name, kmer, size):
	"""
	Assuming there exists only one read in fasta file.
	"""

	size0 = size

	flag = False
	with open(file_name) as mf:
		for ml in mf:
			size = min(len(ml), size0+kmer)

			if ml[0] == '>':
				if flag:
					break
				flag = True
				continue
			else:
				ml = ml[:-1:].replace("N", "")
				if len(ml) < kmer or size < kmer:
					exit()

				for i in xrange(0, size - kmer):
					line, num = Mutation(ml[i:i+kmer])
					print "r"+str(i)+"-"+str(num)
					print line

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
	if len(argv) != 4:
		print "usage: ", argv[0], " [.fa] [length (150)] [size 1024*1024] \n\r(Asumming there exits only one reads in .fa)"
		exit()
	kmer = int(argv[2])
	size = int(argv[3])
	GenerateReadsKmer(argv[1], kmer, size)
