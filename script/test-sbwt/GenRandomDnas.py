import random

from sys import argv


def GenRandomDnas(size):
    if (size/4) < 1:
        print("size is too small")
        return ""
    size = int(size / 4)
    a = "A"*size + "C"*size + "G"*size + "T"*size
    a = list(a)
    random.shuffle(a)

    return "".join(a)


if __name__ == "__main__":

    if len(argv) != 2:
        print("usage: <size>")
        exit(1)

    size = int(argv[1])

    print(">1")
    print(GenRandomDnas(size))
    print
