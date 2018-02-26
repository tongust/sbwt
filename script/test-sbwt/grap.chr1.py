from sys import argv


if __name__ == '__main__':
    if len(argv) != 2:
        print("usage: [fa]")
        exit()

    count_chr = 0
    with open(argv[1]) as mf:
        for ml in mf:
            if ml[0] == '>':
                if count_chr == 0:
                    count_chr += 1
                else:
                    break
                print ml,
            elif count_chr == 1:
                print ml,
            else:
                break
