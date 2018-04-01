from sys import argv
import os

if __name__ == '__main__':
        if len(argv) != 3:
                print "usage: [alignment files] [reference files]"
                exit(1)
        ref = ""
        with open(argv[2])  as mf:
                for ml in mf:
                        if ml.startswith('>'):
                                continue

                        ml.rstrip(os.linesep)
                        ref += ml
        with open (argv[1]) as mf:
                for ml in mf:
                        ml.rstrip(os.linesep)
                        read = ml.split('\t')[4]
                        pos  = int(ml.split('\t')[3])
                        s = ref[pos:pos+len(read)]
                        count = 0
                        for j in range(0, len(read)):
                                if read[j] != s[j]:
                                        count += 1
                        print count

