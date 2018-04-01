from random import randint
from random import shuffle


if __name__ == '__main__':

        size_seed = 50
        period = 3
        size_ref = 150*10000
        seed = "CTGAGTGGCCATTAAATTTGTATTTCAAGAACTGGGTTGACGAAGCCCTA"
        r = ["A","C","G","T"]*(size_ref/4)
        shuffle(r)

        i = 0
        while (i < size_ref/10):
                for j in range(0, size_seed):
                        r[i+j*period+2] = seed[j]
                i += 150

        
        r = ''.join(r)
        wref = open("rep.norm.fa", "w")
        wref.write(">chr1\n")
        wref.write(r)
        wref.write("\n")

        wread = open("rep.norm.read.fa", 'w')
        wread.write(">1\n")
        wread.write(r[150*45:150*46]+"\n")


