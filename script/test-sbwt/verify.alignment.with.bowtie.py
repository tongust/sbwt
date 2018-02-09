import os
from sys import argv


if __name__ == '__main__':
    
    if len(argv) != 3:
        print("usage: [sbwt.log] [others.log]")
        exit(1)
    
    d = set([])
    
    with open(argv[1]) as mf:
        for ml in mf:
            if ml.startswith('>'):
                d.add(ml.rstrip(os.linesep).split('\t')[0][1::])
            else:
                d.add(ml.rstrip(os.linesep).split('\t')[0])

    with open(argv[2]) as mf:
        for ml in mf:
            if ml.startswith('>'):
                key = ml.rstrip(os.linesep).split('\t')[0][1::]
            else:
                key = ml.rstrip(os.linesep).split('\t')[0]
            if key not in d:
                print(ml.rstrip(os.linesep))
