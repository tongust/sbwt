import os


def prepare_data(num_test):
    print("Generating reference fasta...")
    os.system("python ./GenRandomDnas.py " +
              str(num_test+10000) +
              "> f.fa")

    print("Generating faked reads fasta...")
    os.system("python ./GenReferenceAndReads.py f.fa 150 " +
              str(num_test) +
              " > r.fa")


def load_data(file_r, file_f):
    os.system("./sbwt " + file_r + " " + file_f + " > run.check.fa")

    ref_data = []
    run_data = []

    c = False
    with open(file_r, "r") as mf:
        for ml in mf:
            if c:
                c = False
                continue
            else:
                c = True
                ref_data.append(
                    int(
                        ml.split("-")[1].rstrip('\n')
                    )
                )
    with open("run.check.fa", "r") as mf:
        for ml in mf:
            run_data.append(
                int(
                    ml.rstrip('\n')
                )
            )
    # clean files
    os.system("rm -rf run.check.fa")
    return run_data, ref_data


def check_nums(r, f):
    c = 0
    for i, j in zip(r, f):
        if i != j:
            c += 1
    return c


if __name__ == '__main__':
    num_test = 10000000
    prepare_data(num_test)
    print("Test "+str(num_test) + " reads")

    exit()

    r, f = load_data("r.fa", "f.fa")
    c = check_nums(r, f)

    print("Number of mismatches: " + str(c))
