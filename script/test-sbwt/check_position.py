from sys import argv
import os

def reverse_complement(s):
        mp = {'A': 'T', 'C':'G', 'T':'A', 'G':'C', 'N':'A'}
        st = list(s)
        size = len(s)
        t = ['A'] * size
        for i in range(0, size):
                t[size-1-i] = mp[s[i]]
        return ''.join(t)

if __name__ == '__main__':
        if len(argv) != 3:
                print "usage: [alignment files] [reference files]"
                exit(1)
        ref = ""
        with open(argv[2])  as mf:
                for ml in mf:
                        if ml.startswith('>'):
                                continue

                        ml = ml.rstrip(os.linesep)
                        ref += ml
        with open (argv[1]) as mf:
                for ml in mf:
                        ml.rstrip(os.linesep)
                        read = ml.split('\t')[4].rstrip(os.linesep)
                        pos  = int(ml.split('\t')[3])
                        flag = ml.split('\t')[1]
                        s = ref[pos:pos+len(read)]
                        if flag == '-':
                                read = reverse_complement(read)
                        count = 0
                        for j in range(0, len(read)):
                                if read[j] != s[j]:
                                        count += 1
                        print count

