import numpy as np
from random import randint
from sys import argv
if __name__ == '__main__':
    if len(argv) != 2:
        NN = 1000
    else:
        NN = int(argv[1])
    BASE = {0:'A', 1:'C', 2:'G', 3:'T'}
    LL = randint(4, NN)
    LL = NN
    print ''.join([BASE[randint(0,3)] for i in xrange(LL)])
