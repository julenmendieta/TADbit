"""
20 juin 2014

Implementation of iterative correction Imakaev 2012

Schematic flow chart for five iterations (it = 0->5) starting from the symmetric
matrix W of size N:

   +---------->  Wij
   |    
   |              |
   |              v
   |             __
   |             \
   |        Si = /_    Wij
   |            j=0->N
   |    
   |              |
   |              v
   |    
   |                   _
   |        DBi = Si / S
   |    
   |              |
   |              v
   |    
   |        Bi = Bi x DBi             ---> keep track, used as expected value
   |    
   |              |
   |              v
   |    
   |                 Wij
   |       Wij = -----------
   |              DBi x DBj  
   |    
   |              |
   |        it<5 / \ it=5
   |____________/   \_________   TADbit          _           
    it++                     \`----------> Wij / S    meaning that: Si = O(1)
                             |                        ('Si' tends towards one
                             |                         when 'it' -> infinite)
                             |Strict Imakaev
                             |
                             v
                             
                            Wij
                          -------  meaning that: Si = 1
                           ___
                           \
                           /__ Wi

"""

def _update_S(W):
    S = {}
    meanS = 0.0
    for bin1 in W:
        S[bin1] = sum(W[bin1].values())
        meanS += S[bin1]
    meanS /= len(W)
    return S, meanS

def _updateDB(S, meanS, B):
    DB = {}
    for bin1 in S:
        DB[bin1] = float(S[bin1]) / meanS
        B[bin1] *= DB[bin1]
    return DB

def _update_W(W, DB):
    for bin1 in W:
        DBbin1 = DB[bin1]
        W1 = W[bin1]
        for bin2 in W1:
            try:
                W1[bin2] /= DBbin1 * DB[bin2]
            except ZeroDivisionError: # whole row is empty
                continue

def iterative(hic_data, remove=None, iterations=0):
    """
    """
    size = len(hic_data)
    remove = remove or tuple([int(hic_data[i+i*size]==0) for i in xrange(size)])
    W = {}
    for i in xrange(size):
        if remove[i]:
            continue
        W[i] = {}
        for j in xrange(size):
            if remove[j]:
                continue
            W[i][j] = hic_data[i, j]
    B = dict([(b, 1.) for b in W])
    for _ in xrange(iterations + 1):
        S, meanS = _update_S(W)
        DB = _updateDB(S, meanS, B)
        _update_W(W, DB)
    for i in xrange(size):
        try:
            if B[i]:
                B[i] *= meanS**.5
            else:
                B[i] = 1.
        except KeyError:
            B[i] = 1.
    return B