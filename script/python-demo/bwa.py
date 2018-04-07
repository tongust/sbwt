# Burrows-Wheeler Transform
import numpy as np
import mbwa
import time
from sys import argv
import os
from random import randint
MCCC = 0;
def rotation(t):
    """ Return list of rotations of input stringt """
    tt = t*2
    return [ tt[i:i+len(t)] for i in xrange(0, len(t))]


def bwm(t):
    """ Return lexicographically sorted list of t's rotations """
    return sorted(rotation(t))


def bwtViaBwm(t):
    """ Given T, returns BWT(T) by the way of the BWM """
    return ''.join(map(lambda x: x[-1], bwm(t)))


def calCO(B):
    """ Let C[a] be the the number of symbols in X[0, n-2] taht are lexicographically smaller than a
        and O[a,i] be the # of occurrences of a in B[0, i].
    """
    C = []
    N = len(B) - 1
    O = np.zeros([4,N+2]).tolist()# record A C G expect T
    mc = 0
    ac = cc = gc = 0
    for i in B:
        mc += 1
        if i == 'A':
            ac += 1
            O[0][mc] = O[0][mc-1] + 1# A changing
            O[1][mc] = O[1][mc-1]# C
            O[2][mc] = O[2][mc-1]# G
            O[3][mc] = O[3][mc-1]# T
            continue
        elif i == 'C':
            cc += 1
            O[0][mc] = O[0][mc-1]# A
            O[1][mc] = O[1][mc-1] + 1# C changing
            O[2][mc] = O[2][mc-1]# G
            O[3][mc] = O[3][mc-1]# T
            continue
        elif i == 'G':
            gc += 1
            O[0][mc] = O[0][mc-1]# A
            O[1][mc] = O[1][mc-1]# C
            O[2][mc] = O[2][mc-1] + 1# G changing
            O[3][mc] = O[3][mc-1]# T
            continue
        elif i == 'T': # T changing
            O[0][mc] = O[0][mc-1]# A
            O[1][mc] = O[1][mc-1]# C
            O[2][mc] = O[2][mc-1]# G
            O[3][mc] = O[3][mc-1] + 1# T Changing
            continue
        else:
            O[0][mc] = O[0][mc-1]# A
            O[1][mc] = O[1][mc-1]# C
            O[2][mc] = O[2][mc-1]# G
            O[3][mc] = O[3][mc-1]# T
    for i in [0,1,2,3]:
        O[i][::] = O[i][1::]
    for i in xrange(0, len(O)):
        for j in xrange(0, len(O[0])):
            O[i][j] = int(O[i][j])
    C = [0,0,0,0]
    C[0] = 0
    C[1] = C[0] + ac
    C[2] = C[1] + cc
    C[3] = C[2] + gc
    for i in range(0,0):###############################################
        O[i][0] = 0;
    return C, O


def CalculateD(W, X, C, O_, BASE):
    D = [0]*len(W)
    k = 0
    N_X = len(X) - 1
    l = N_X
    z = 0
    for i in range(0, len(W)):
        if (k == 0) :
            k = C[BASE[W[i]]] + 1
        else:
            k = C[BASE[W[i]]] + O_[BASE[W[i]]][k-1]+1
        if l == -1 :
            l = C[BASE[W[i]]]
        else:
            l = C[BASE[W[i]]] + O_[BASE[W[i]]][l]
        k, l = int(k), int(l)
        if k > l:
            k = 0
            l = N_X
            z += 1
        D[i] = z
    return D


def InexRecur1(W, i, z, k, l, D, C, O, BASE):
    global MCCC
    MCCC += 1;
    # condition: i < 0
    if i < 0:
        if z >= 0:
            return [k,l]
        else:
            return []
    if z < D[i]:
        return []

    I = []
#    I.extend(InexRecur1(W, i-1, z-1, k, l, D, C, O, BASE)) # insertion
    N_O = len(O)-1
    k_ = k
    l_ = l ################
    for ib in ['A', 'C', 'G', 'T']:
        if k_ == 0:
            k = C[BASE[ib]] + 1

        else:
            k = C[BASE[ib]] + O[BASE[ib]][k_-1]+1
        if l_ == -1:
            l = C[BASE[ib]]
        else:
            l = C[BASE[ib]] + O[BASE[ib]][l_]
        k,l = int(k), int(l)
        if k <= l:
#            I.extend(InexRecur1(W, i, z-1, k, l, D, C, O, BASE))  # Deletion of W
            if ib == W[i]:
                I.extend(InexRecur1(W, i-1, z, k, l, D, C, O, BASE)) # match
                if len(I) >= 2:
                    return I
            else:
                I.extend(InexRecur1(W, i-1, z-1, k, l, D, C, O, BASE)) # Mismatch
                if len(I) >= 2:
                    return I
    return I


def InexRecur1_gap(W, i, z, k, l, D, C, O, BASE):
    #global MCCC
    #MCCC += 1;
    # condition: i < 0
    if i < 0:
        if z >= 0:
            return [k,l]
        else:
            return []
    if z < D[i]:
        return []

    I = []
    I.extend(InexRecur1_gap(W, i-1, z-1, k, l, D, C, O, BASE)) # insertion
    N_O = len(O)-1
    k_ = k
    l_ = l
    for ib in ['A', 'C', 'G', 'T']:
        if k_ == 0:
            k = C[BASE[ib]] + 1

        else:
            k = C[BASE[ib]] + O[BASE[ib]][k_-1]+1
        if l_ == -1:
            l = C[BASE[ib]]
        else:
            l = C[BASE[ib]] + O[BASE[ib]][l_]
        k,l = int(k), int(l)
        if k <= l:
            I.extend(InexRecur1_gap(W, i, z-1, k, l, D, C, O, BASE))  # Deletion of W
            if ib == W[i]:
                I.extend(InexRecur1_gap(W, i-1, z, k, l, D, C, O, BASE)) # match
            else:
                I.extend(InexRecur1_gap(W, i-1, z-1, k, l, D, C, O, BASE)) # Mismatch
    return I


def ExMatch(W, X, C, O, BASE):
    k = 0
    N_X = len(X) - 1
    l = N_X
    z = 0
    for i in range(0, len(W)):
        if (k == 0) :
            k = C[BASE[W[i]]] + 1
        else:
            k = C[BASE[W[i]]] + O[BASE[W[i]]][k-1]+1

        if l == -1 :
            l = C[BASE[W[i]]]

        else:
            l = C[BASE[W[i]]] + O[BASE[W[i]]][l]
        k, l = int(k), int(l)
        if k > l:
            return [k,l]
    return [k,l]


def suffixArray(s):
    satups = sorted((s[i:], i) for i in xrange(0, len(s)+1))
    return map(lambda x: x[1], satups)


def InexSearch_gap(W, X, z):
    #z = 2 # the upper bound on # of the differences (mismatches / gaps)
    X += "$"
    BASE = {'A':0, 'C':1, 'G':2, 'T':3}
    ##### Precalculation:
    B = bwtViaBwm(X)                         # BWT of X
    C, O = calCO(B)
    X_ = X[len(X)-2::-1] # the reverse string X
    X_ += '$'
    B_ = bwtViaBwm(X_)
    C, O_ = calCO(B_)
    #print X,'\n',B,'\n',X_,'\n',B_
    #print C
    #for i in O:
    #    print i

    #for i in O_:
    #    print i
    D = CalculateD(W, X, C, O_, BASE)
    #D = D - [0,1,1,1,2]
    #D = [0,0,0,0,1]
    print D
    Xs = rotation(X)

    MC = -1
    #for i in Xs:
    #    MC += 1
    #    print i, MC
    Xs = bwm(X)
    MC = -1
    #for i in Xs:
    #    MC += 1
    #    print i, MC
    return InexRecur1_gap(W, len(W)-1, z, 0, len(X)-1, D, C, O, BASE)


def InexSearch_nongap(Ws, X, z):
    X += "$"
    BASE = {'A':0, 'C':1, 'G':2, 'T':3}
    ##### Precalculation:
    B = bwtViaBwm(X)                         # BWT of X
    C, O = calCO(B)
    X_ = X[len(X)-2::-1] # the reverse string X
    X_ += '$'
    B_ = bwtViaBwm(X_)
    C, O_ = calCO(B_)
    # cal SA
    mSA = suffixArray(X)
    mres = [-1]*len(Ws)
    t00 = time.clock()
    cic = -1
    for W in Ws:
        cic += 1
        D = CalculateD(W, X, C, O_, BASE)
        resIntervals = InexRecur1(W, len(W)-1, z, 0, len(X)-1, D, C, O, BASE)
        nRes = int(len(resIntervals)/2)

#        for i in xrange(0,nRes):
#            for j in xrange(resIntervals[i], resIntervals[i+1]+1):
#                mres.extend([mSA[j+1]])# the array in python start at zero.
###standard
        if nRes > 0 and resIntervals[0] <= resIntervals[1]:
            mres[cic] = mSA[resIntervals[0]+1]
    d01 = time.clock() - t00
    #print d01
    return mres, d01


def bwaTEST():
    #X = "ACGGGGCGGGCTGGTCGTAGAGCGTAGGCGTACGTGCT"
    X = "ACGAAAACCGGCGCGCACGTACGAGGTCGTAGGTCGTAGGTCGTA"
    W =           "GCGCGTACGTACGA"
    Ws = []
    for i in xrange(0,100000):
        Ws.extend([W])
    res = InexSearch_nongap(Ws, X, 1)#mbwa.mBWT(X)
    #print res
    #InexSearch(W,X,1)


def TESTall():
    W =             "CGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACC"#0
    #X = "ACCGGCGTTCGATCCT"
    #W =      "CTTTCGATCC"
    #              0 0 0 0 0 0 0 0
    print len(X), len(W)
    z = 1
    Ws = []
    for i in xrange(0,10000):
        Ws.extend([W])
#    Ws.extend([W])
#    Ws.extend([X[100:149]+"T"])
#    print Ws
    print len(Ws)
    # bwa InexSearch
    t00 = time.clock()
    res0 = InexSearch_nongap(Ws, X, z)
    t01 = time.clock()
    d1 = t01 - t00
    print "BWA costs: ", (t01 - t00), " s"
    print len(res0)
    # supsample Burrows-Wheeler transform
    t00 = time.clock()
    res1 = mbwa.sBwtInexactMatch_nongap(Ws, X, z)
    d2 = time.clock() - t00
    print "sBWA costs: ", (time.clock() - t00), " s"
    print len(res1)
    print "speed up: ", (d1/d2)
    print len(res1),len(res1),min(res1),min(res1),max(res1),max(res1)

    
def TESTmis_1():
    if len(argv) != 5:
        print "usage: xx.fa begin end output (only line)"
        exit(1)

    beg = int(argv[2])
    end = int(argv[3])
    with open(argv[1]) as mf:
        for ml in mf:
            if not ml.startswith('>'):
                ml = ml.rstrip(os.linesep)
                X = ml
                break

    len_ref = len(X)
    BASE = {'A':0, 'C':1, 'G':2, 'T':3}
    z = 1
    num_tst = 100000
    mTable = ["CGT","AGT", "ACT", "ACG"]
    mRecord = open(argv[4],'w')
    mRecord.write("len,bwa,sbwa,speedup,tests,errors\n")
    print "len speedup  bwa sbwa postivie true errors"
    for t1 in xrange(beg, end):#500
        Ws = []
        sts_ref = []
        #len_pattern = t1*8 + 10
        len_pattern = t1
        if t1 >= end:
            break

        #len_pattern = 100
        #if len_pattern > 600:
            #break
        region_slect = len_ref - len_pattern + 1# [0, region_slect)
        mcc1 = 0
        start1 = -1
        while True:
            mcc1 += 1
            start1 += 1
            start1 = start1%region_slect
            #print start1, len_pattern
            pos_mutation = randint(0, len_pattern-1)
            W = ""+X[start1:start1+pos_mutation:]
            W += mTable[BASE[X[start1+pos_mutation]]][randint(0,2)] # MUTATION
            W += X[start1+pos_mutation+1:start1+len_pattern:]
            Ws.extend([W])
            sts_ref.extend([start1])
            if mcc1 >= num_tst :
                break
        # Ws
        # bounded bwa backward
        res0,d0 = InexSearch_nongap(Ws, X, z)
        # sbwa
        res1,d1 = mbwa.sBwtInexactMatch_nongap(Ws, X, z)
        # verify
        dis_sbwa = 0
        for iver in xrange(0, num_tst):
            mc_sbwa = 0
            t363 = res1[iver]
            if t363 >= 0:
                if sum(a!=b for a,b in zip(X[t363:len_pattern+t363],Ws[iver])) > 5:
                    dis_sbwa += 1
            else:
                dis_sbwa += 1
        print len_pattern, d0/(d1+10**-14), d0, d1, dis_sbwa
        mRecord.write(str(len_pattern)+","+str(d0)+","+str(d1)+","+str(d0/(d1+10**-14))+","+str(num_tst)+","+str(dis_sbwa)+"\n")
    mRecord.close()
    return


def nonfun():
    z = 2 # the upper bound on # of the differences (mismatches / gaps)
    X += "$"
    BASE = {'A':0, 'C':1, 'G':2, 'T':3}
    ##### Precalculation:
    B = bwtViaBwm(X)                         # BWT of X
    C, O = calCO(B)
    X_ = X[len(X)-2::-1] # the reverse string X
    X_ += '$'
    B_ = bwtViaBwm(X_)
    C, O_ = calCO(B_)
    print X,'\n',B,'\n',X_,'\n',B_
    print C
    for i in O:
        print i
    print
    for i in O_:
        print i
    D = CalculateD(W, X, C, O_, BASE)
    #D = D - [0,1,1,1,2]
    #D = [0,0,0,0,1]
    print D
    Xs = rotation(X)

    MC = -1
    for i in Xs:
        MC += 1
        print i, MC
    Xs = bwm(X)
    MC = -1
    for i in Xs:
        MC += 1
        print i, MC

    print InexRecur1(W, len(W)-1, z, 0, len(X)-1, D, C, O, BASE)
    #print ExMatch("GCA", X, C, O, BASE)
    print "TOTAL: " + str(MCCC)


if __name__ == "__main__":
    TESTmis_1()
