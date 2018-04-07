import numpy as np
from bwa import rotation
from bwa import calCO
import time 
def mSA(s):# interval: odd and even. 
    if s[-1] != "$":
        s += "$$"
        #s = s[0:-1]# cut off the last "$"
#    ml = len(s)

    satups = sorted((s[i::2], i) for i in xrange(0, len(s)))
    #print satups 
    sa =  map(lambda x: x[1], satups)
    return sa
    #print sa
def mBWT(s):# construct sbwt from sa
    if s[-1] != '$':
        s += "$$"
    ml = len(s)
    sbwt = ''
    sa =  mSA(s)
    for i in xrange(0, ml):
        if sa[i] == 0 or sa[i] == 1:
            sbwt += '$'
        else:
            sbwt += s[sa[i]-2]
#    print sa
#    arr = rotation(s)
#    tc = -1
#    print "sbwt: \n--------------"
#    for i in sa:
#        tc += 1
#        print arr[i], tc
#    print "--------------\nend of sbwt\n"
    return sbwt
def calCO_sbwt(B):
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
    C[0] = 1
    C[1] = C[0] + ac
    C[2] = C[1] + cc
    C[3] = C[2] + gc
    return C, O

def sBwaExactMatch_backup(W, X, C, O, BASE):
    k = 0
    l = len(X)+1
    z = 0
    mflag = 0
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
#        print "sbwt exaxt match", i, " ", [k,l]
        if k > l:
            return [k,l]
        if k == l:
            mflag += 1
            if mflag == 2:
                return [k,l]

    return [k,l]


def sBwaExactMatch(W, X, C, O, BASE):
    k = 0
    l = len(X)+1
    z = 0
    mflag = 0
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
#        print "sbwt exaxt match", i, " ", [k,l]
        if k > l:
            return [k,l,i]
        if k == l:
            mflag += 1
            if mflag == 2:
                return [k,l,i]

    return [k,l,len(W)-1]
def sBwt_nongap(w,x,z,BASE, msSA, C, O):
    mtc = 0
    mres = []# the resault of return
    len_w = len(w)
    len_half_w = int(len_w/2) - 1
    w0 = w[-1::-2]#111
    w1 = w[-2::-2]#000
    
    #print msSA
    # even
    kl0 = sBwaExactMatch(w0, x, C, O, BASE)
    if len(kl0) == 2:
        kl0.extend([len_half_w])
    for i in xrange(kl0[0], kl0[1]+1):
        mtc = 0
        pos_start_x = msSA[i]-1-(len_half_w-kl0[2])*2
        pos_end_x = pos_start_x + len_w 

        if (sum(a==b for a,b in zip(x[pos_start_x  :  pos_end_x  :2], w[0::2])) >= len_half_w):
            #mres.extend([pos_start_x])# at 1
            return pos_start_x
    # Greedly
    #if len(mres) != 0:
    #    return mres
    # odd
    kl1 = sBwaExactMatch(w1, x, C, O, BASE)
    if len(kl1) == 2:
        kl1.extend([len_w])
    for i in xrange(kl1[0], kl1[1]+1):
        mtc = 0
        pos_start_x = msSA[i]+1-(len_half_w-kl1[2])*2
        pos_end_x = pos_start_x+len_w-1
        if sum(a==b for a,b in zip(x[pos_start_x  : pos_end_x :2], w[1::2])) >= len_half_w :
            #mres.extend([pos_start_x-1])# at 0
            return pos_start_x-1

    return -1
def sBwt_nongap_backup(w,x,z,BASE, msSA, C, O):
    mtc = 0
    mres = []# the resault of return
    w0 = w[-1::-2]#111
    w1 = w[-2::-2]#000
    
    #print msSA
    # even
    kl0 = sBwaExactMatch(w0, x, C, O, BASE)
    
#    print w0[::-1]
#    print w0
    for i in xrange(kl0[0], kl0[1]+1):
        #print x[msSA[i]:msSA[i]+14:2]
        mtc = 0
        for vf,vq in zip(x[msSA[i]-1:msSA[i]+13:2], w[0::2]):
            if vf != vq:
                mtc += 1
            if mtc > z:
                break
        if mtc <= z:
            mres.extend([msSA[i]-1])# at 1
    # Greedly
    if len(mres) != 0:
        return mres
    # odd
    kl1 = sBwaExactMatch(w1, x, C, O, BASE)
#    print kl1
    for i in xrange(kl1[0], kl1[1]+1):
        #print x[msSA[i]+1:msSA[i]+14:2]
        mtc = 0
        for vf,vq in zip(x[msSA[i]+1:msSA[i]+14:2], w[1::2]):
            if vf != vq:
                mtc += 1
            if mtc > z:
                break
        if mtc <= z:
            mres.extend([msSA[i]])# at 0
    return mres
def sBwtInexactMatch_nongap(ws, x, z):
    BASE = {'A':0, 'C':1, 'G':2, 'T':3}
    sB = mBWT(x)
    C, O = calCO_sbwt(sB)
    msSA = mSA(x)
    res = [-1]*len(ws)
    t00 = time.clock()
    cic = -1
    for w in ws:
        cic += 1
        res[cic] = sBwt_nongap(w,x,z,BASE, msSA, C, O)
        #res += (ta, )
    d01 = time.clock() - t00
    #print d01
    return res, d01
    
def sbwaTEST():
    BASE = {'A':0, 'C':1, 'G':2, 'T':3}
    X = "ACGAAAACCGGCGCGCACGTACGAGGTCGTAGGTCGTAGGTCGTA"
    W =           "GCGCGTACGTACGA"
    Ws = []
    for i in xrange(0,100000):
        Ws.extend([W])
    t00 = time.clock()
    mres = sBwtInexactMatch_nongap(Ws,X, 1)
    print "sbwt: " , (time.clock() - t00)
    
