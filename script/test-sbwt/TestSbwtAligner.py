import os
import random


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


def test_bitset():
    num_test = 10000000
    prepare_data(num_test)
    print("Test "+str(num_test) + " reads")

    exit()

    r, f = load_data("r.fa", "f.fa")
    c = check_nums(r, f)

    print("Number of mismatches: " + str(c))


def gen_reads(period, ref, reads):
    num_reads = 100
    size_read = random.randrange(7, 128)
    fopen = open(reads, 'w')
    with open(ref) as mf:
        for ml in mf:
            if ml.startswith(">"):
                continue
            else:
                ml = ml.rstrip("\n")
                size = len(ml)
                for i in range(0, size - period*size_read):
                    line = ""
                    first = True
                    for j in range(0, size_read):
                        # with error
                        if random.randrange(0, 10) > 5 and first:
                            line += "A"
                            first = False
                        else:
                            line += ml[i+j*period]
                    fopen.write(">\n")
                    fopen.write(line+"\n")

                    num_reads -= 1
                    if num_reads < 0:
                        pass
                        break

    fopen.close()


def spaced_exact_match(period, ref, reads, res):
    fw = open(res, "w")

    q = ""
    with open(ref) as mf:
        for ml in mf:
            if ml.startswith(">"):
                continue
            else:
                q = ml.rstrip("\n")

    size = len(q)
    with open(reads) as mf:
        for ml in mf:
            p = ""
            if ml.startswith(">"):
                continue
            else:
                p = ml.rstrip('\n')
            end = size - len(p)*period

            match_res = []

            for i in range(0, end):
                match = True
                for j in range(0, len(p)):
                    if p[j] != q[i+j*period]:
                        match = False
                        break
                if match:
                    match_res.append(i)
            if len(match_res) == 0:
                match_res.append(-1)
            fw.write(','.join(str(u) for u in match_res) + "\n")
    fw.close()
    return


def verify_match(py_log, cpp_log):

    with open(py_log) as p, open(cpp_log) as q:
        for a, b in zip(p, q):
            a = a.rstrip(os.linesep).split(",")
            b = b.rstrip(os.linesep).split(",")
            a = [int (i) for i in a]
            b = [int (i) for i in b]
            a.sort()
            b.sort()

            if a != b:
                if len(b) < len(a):
                    return False
                else:
                    b = b[0:len(a):]
                    if a != b:
                        return False
    return True


def test_exact_match():
    flog = open("test.exact.match.log", "w")

    for t in range(0,10000):
        print("---------------------------- " + str(t))
        os.system("rm -rf genome.fa* reads.fa *.log")

        num_ref = random.randrange(128*21, 100*1024)
        period = random.randrange(1, 20)
        reads = "reads.fa"
        ref = "genome.fa"

        os.system("python ./GenRandomDnas.py " + str(num_ref) + " > genome.fa")
        os.system("./build_index " + str(period) + " " + ref)
        gen_reads(period, ref, reads)
        print("sbwt")
        os.system("./sbwt_test reads.fa genome.fa > res.cpp.log")
        print("...")
        spaced_exact_match(period, ref, reads, "res.py.log")
        print("naive method")
        if verify_match("res.py.log", "res.cpp.log"):
            print("Successful")
        else:
            print("Failed!")
            break
        flog.close()


if __name__ == '__main__':
    test_exact_match()

