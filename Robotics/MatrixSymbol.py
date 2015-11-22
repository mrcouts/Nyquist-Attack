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
    
def lsolve(self,other):
    """ LDLsolve for SMatrix """
    row1row2_ = union(self.rowl_,other.rowl_)
    M_ = self.S_(row1row2_,self.coll_).LDLsolve( other.S_(row1row2_,other.coll_) )
    return SMatrix(M_,self.coll_,other.coll_)    

class SMatrix(object):
    """ Matrizes acessadas por simbolos """
    def __init__(self, M_, rowl_=['svector'], coll_=['vector']):
        if rowl_ == ['svector'] and M_.cols ==1:
            rowl_ = M_
        rowl_ = list(rowl_)
        coll_ = list(coll_)
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
        
    def __repr__(self):
        if self.rowl_ == list(self.M_):
            return pretty(self.M_)
        elif self.coll_ == ['vector']:
            return pretty(Matrix([ Matrix(self.rowl_).T, self.M_.T ]).T)
        else:
            return pretty(Matrix([ Matrix(['_']+self.rowl_).T, Matrix([ Matrix(self.coll_).T, self.M_ ]).T  ]).T)
        
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
        if other == 0:
            return self
        else:
            rowl_ = union(self.rowl_,other.rowl_)
            coll_ = union(self.coll_,other.coll_)
            M_ = self.S_(rowl_,coll_) + other.S_(rowl_,coll_)
            return SMatrix(M_,rowl_,coll_)
            
    def __radd__(self, other):
        if other == 0:
            return self
        else:
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
        if not isinstance(other, SMatrix):
            return SMatrix(other*self.M_,self.rowl_,self.coll_)
        else:
            col1row2_ = union(self.coll_,other.rowl_)
            M_ = self.S_(self.rowl_,col1row2_) * other.S_(col1row2_,other.coll_)
            return SMatrix(M_,self.rowl_,other.coll_)
            
    def __rmul__(self, other):
        if not isinstance(other, SMatrix):
            return SMatrix(other*self.M_,self.rowl_,self.coll_)
        else:
            col1row2_ = union(self.coll_,other.rowl_)
            M_ = self.S_(self.rowl_,col1row2_) * other.S_(col1row2_,other.coll_)
            return SMatrix(M_,self.rowl_,other.coll_)
            
    def __mod__(self, other):
        rowl_ = union(self.rowl_,other.rowl_)
        coll_ = union(self.coll_,other.coll_)
        M_ = self.S_(rowl_,coll_).multiply_elementwise(other.S_(rowl_,coll_))
        return SMatrix(M_,rowl_,coll_)
        
    def __div__(self,other):
        return SMatrix(self.M_ / other, self.rowl_, self.coll_)
        
    def __truediv__(self,other):
        return SMatrix(self.M_ / other, self.rowl_, self.coll_)
        
    def T(self):
        return SMatrix(self.M_.T,self.coll_,self.rowl_)
        
    def inv(self):
        return SMatrix(self.M_.inv(),self.coll_,self.rowl_)
        
    def pinv(self):
        return SMatrix(self.M_.pinv(),self.coll_,self.rowl_)
        
    def nullspace(self, symplify=False):
        ns_ = self.M_.nullspace(symplify)
        M_ns_ = Matrix([ns_i.T for ns_i in ns_]).T
        rowl_ = self.coll_
        coll_ = [symbols('pi_' + str(i+1)) for i in xrange(M_ns_.cols)]
        return SMatrix(M_ns_,rowl_,coll_)
        
    def LDLsolve(self,other):
        row1row2_ = union(self.rowl_,other.rowl_)
        M_ = self.S_(row1row2_,self.coll_).LDLsolve( other.S_(row1row2_,other.coll_) )
        return SMatrix(M_,self.coll_,other.coll_)   
        
    def diff(self, *symbols, **kwargs):
        return SMatrix(self.M_.diff(*symbols, **kwargs), self.rowl_, self.coll_)
        
    def subs(self, *args, **kwargs):
        return SMatrix(self.M_.subs(*args, **kwargs), self.rowl_, self.coll_)
        
    def cols(self):
        return self.M_.cols
        
    def rows(self):
        return self.M_.rows
        
    def extract(self, rowl_, coll_):
        return SMatrix(self.S_(rowl_,coll_), rowl_, coll_)