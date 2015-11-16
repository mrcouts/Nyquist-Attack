# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 09:13:11 2015

@author: couto
"""

from sympy import *
#import math
#from operator import concat
init_printing(use_unicode=True)

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

class SMatrix(object):
    """ Matrizes acessadas por simbolos """
    def __init__(self, M_, rowl_, coll_=['1']):
        self.rowl_ = rowl_
        self.coll_ = coll_
        self.dic_rowl = {rowl_[i]:i for i in xrange(len(rowl_))}
        self.dic_coll = {coll_[i]:i for i in xrange(len(coll_))}
        if M_ == 0:
            self.M_ = zeros(len(rowl_),len(coll_))
        elif M_ == 1:
            self.M_ = Matrix([[1 if i==j else 0 for j in coll_] for i in rowl_] )
        else:
            self.M_ = M_
        
        
    def S(self,rowl,coll):
        if rowl in self.dic_rowl and coll in self.dic_coll:
            return self.M_[self.dic_rowl[rowl],self.dic_coll[coll]]
        else:
            return 0
    
    def S_(self,rowl_,coll_):
        nrows = len(rowl_)
        ncols = len(coll_)
        M = zeros(nrows,ncols)
        for i in xrange(nrows):
            for j in xrange(ncols):
                M[i,j] = self.S(rowl_[i],coll_[j])
        return M
        
    def __add__(self, other):
        rowl_ = union(self.rowl_,other.rowl_)
        coll_ = union(self.coll_,other.coll_)
        M_ = self.S_(rowl_,coll_) + other.S_(rowl_,coll_)
        return SMatrix(M_,rowl_,coll_)
        
    def __sub__(self, other):
        rowl_ = union(self.rowl_,other.rowl_)
        coll_ = union(self.coll_,other.coll_)
        M_ = self.S_(rowl_,coll_) - other.S_(rowl_,coll_)
        return SMatrix(M_,rowl_,coll_)
        
    def __mul__(self, other):
        col1row2_ = union(self.coll_,other.rowl_)
        M_ = self.S_(self.rowl_,col1row2_) * other.S_(col1row2_,other.coll_)
        return SMatrix(M_,self.rowl_,other.coll_)
        
    def T(self):
        return SMatrix(self.M_.T,self.coll_,self.rowl_)
        
    def inv(self):
        return SMatrix(self.M_**-1,self.coll_,self.rowl_)
        
    def ns(self):
        ns_ = self.M_.nullspace()
        M_ns_ = Matrix([ns_i.T for ns_i in ns_]).T
        rowl_ = self.coll_
        coll_ = [symbols('pi_' + str(i+1)) for i in xrange(M_ns_.cols)]
        return SMatrix(M_ns_,rowl_,coll_)
    
M1_ = Matrix([[1000,2000],[3000,4000]])
rowl1_ = ['a','b']
coll1_ = ['c','d']
M2_ = Matrix([[5,6],[7,8]])
rowl2_ = ['b','e']
coll2_ = ['d','f']
sM1_ = SMatrix(M1_,rowl1_,coll1_)
sM2_ = SMatrix(M2_,rowl2_,coll2_)
sM3_ = sM1_ + sM2_
sM4_ = sM1_.T() * sM2_
pprint(sM1_.S_(['a','b','e'],['c','d','f']))
pprint(sM2_.S_(['a','b','e'],['c','d','f']))
pprint(sM3_.M_)
pprint(sM4_.M_)

M5_ = Matrix([symbols('a'),symbols('b')])
rowl5_ = ['a','b']
sM5_ = SMatrix(M5_,rowl5_)

M6_ = Matrix([symbols('b'),symbols('c')])
rowl6_ = ['b','c']
sM6_ = SMatrix(M6_,rowl6_)

sM7_ = sM5_ - sM6_
pprint(sM7_.M_)
print(sM7_.rowl_)

MA_ = Matrix([
    [symbols('A_aw'), symbols('A_ax'), symbols('A_ay')],
    [symbols('A_bw'), symbols('A_bx'), symbols('A_by')],
    [symbols('A_cw'), symbols('A_cx'), symbols('A_cy')],
    [symbols('A_dw'), symbols('A_dx'), symbols('A_dy')]
    ])
rowlA_=['a','b','c','d']
collA_=['w','x','y']
sMA_ = SMatrix(MA_, rowlA_, collA_)

MB_ = Matrix([
    [symbols('B_aq'), symbols('B_aw'), symbols('B_ay'), symbols('B_az')],
    [symbols('B_cq'), symbols('B_cw'), symbols('B_cy'), symbols('B_cz')],
    [symbols('B_rq'), symbols('B_rw'), symbols('B_ry'), symbols('B_rz')],
    [symbols('B_sq'), symbols('B_sw'), symbols('B_sy'), symbols('B_sz')],
    ])
rowlB_=['a','c','r','s']
collB_=['q','w','y','z']
sMB_ = SMatrix(MB_, rowlB_, collB_)

# pprint(sMA_.M_)
# pprint(sMB_.M_)

sMC_ = sMA_ + sMB_
pprint(sMC_.M_)
pprint(sMC_.rowl_)
pprint(sMC_.coll_)

sMD_ = sMB_ * (sMA_.T())
pprint(sMD_.M_)
pprint(sMD_.rowl_)
pprint(sMD_.coll_)

sMx_ = SMatrix(Matrix([1,2,3]).T,['1'],['a','b','c'])
sMx_ns_ = sMx_.ns()
pprint(sMx_ns_.M_)

One_ = SMatrix(1,['a','b','c','d','e'],['c','b','d','a'])
pprint(One_.M_)