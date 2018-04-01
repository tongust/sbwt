from sys import argv
from random import randint
import os


if __name__ == '__main__':
        
        if len(argv) != 2:
                print("usage: xx.log")
                exit()

        mymap = {'A': 'CGT',
                'C': 'AGT',
                'G': 'ACT',
                'T': 'ACG'}
        with open(argv[1]) as mf:
                for ml in mf:
                        if len(ml)>=30 and (ml.startswith('A') or
                                           ml.startswith('B') or
                                           ml.startswith('C') or
                                           ml.startswith('G')):
                                ml = ml.replace(' ', '')
                                ml = ml.rstrip(os.linesep)
                                ml_bcakup = ml
                                for i in range(0, 1000):
                                        print '>'+ml_bcakup
                                        ml = ml_bcakup
                                        ml = list(ml)
                                        i0 = randint(0, len(ml)-1)
                                        i1 = randint(0, len(ml)-1)
                                        while i0 == i1:
                                                i1 = randint(0, len(ml)-1)
                                        ml[i0] = mymap[ml[i0]][randint(0, 2)]
                                        ml[i1] = mymap[ml[i1]][randint(0, 2)]
                                        i0 = randint(0, len(ml)-1)
                                        ml[i0] = mymap[ml[i0]][randint(0, 2)]
                                        print ''.join(ml)
