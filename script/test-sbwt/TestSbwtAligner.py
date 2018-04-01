import os
import sys
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


def verify_sbwt_from_log(res_log, reads, period):

    # check position
    check1 = True
    log_set = set()

    with open(res_log) as mf:
        for ml in mf:
            if len(ml) > 0:
                ml.rstrip(os.linesep)
                ml = ml.split("\t")
                posl = ml[0].split("-")[-1]
                posr = ml[3]
                log_set.add(ml[0].split("-")[0])
                if posl != posr:
                    print "Failure in phase 1!"
                    check1 = False

    check2 = True

    with open(reads) as mf:
        for ml in mf:
            if ml.startswith(">"):
                ml.rstrip(os.linesep)
                a = ml[1::].split('-')[0]
                b = int(ml.split('-')[1])
                if b < period:
                    if a not in log_set:
                        check2 = False
                        print period,a
                        print "Failure in phase 2!"

    return check1 and check2

def test_exact_match():

    for t in range(0,100):
        print("---------------------------- " + str(t))
        os.system("rm -rf genome-test.fa* reads-test.fa *.log")
        num_ref = 150000
        size_read = (num_ref / 150) - 10;
        period = random.randrange(2, 5)
        period = 3
        reads = "reads-test.fa"
        ref = "genome-test.fa"

        os.system("python ./GenRandomDnas.py " + str(num_ref) + " > genome-test.fa")
        print("Generating faked reads fasta...")
        os.system("python ./GenReferenceAndReads.py genome-test.fa 150 " + str(size_read) + " > reads-test.fa")

        os.system("./build_index " + ref + " " + str(period) + " 50" + " > /dev/null 2>&1")

        os.system("./sbwt reads-test.fa genome-test.fa > res.sbwt.log")

        if  not verify_sbwt_from_log("res.sbwt.log", "reads-test.fa", period):
            print "ERROR"
            break


    print "================================================="

if __name__ == '__main__':
    test_exact_match()

