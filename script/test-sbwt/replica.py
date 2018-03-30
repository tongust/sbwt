from sys import argv
from random import randint
import random
import os


if __name__ == '__main__':

        size_rep = 10
        num_rep = 500
        size_read = 10
        period = 3

        size_ref = 102400
        size1 = size_ref / 4
        a = "A"*size1 + "C"*size1 + "G"*size1 + "T"*size1
        s = list(a)
        random.shuffle(s)

        for i in range(0, size_rep):
                size1 = 12
                a = "A"*size1 + "C"*size1 + "G"*size1 + "T"*(size1+2)
                a = list(a)
                random.shuffle(a)
                a = ''.join(a)

                for j in range(0, num_rep):
                        pos = randint(0, size_ref-size_read*4)
                        for k in range(0, size_read):
                                s[pos + k*3] = a[k]

        s = ''.join(s)
        print s
        print

