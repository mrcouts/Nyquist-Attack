from sympy import *
import math
from operator import concat
init_printing(use_unicode=True)

###############################################################################

t = symbols('t')

Rotx = lambda T: Matrix([[1, 0, 0],[0, cos(T),-sin(T)],[0, sin(T), cos(T)]])
Roty = lambda T: Matrix([[cos(T), 0, sin(T)],[0, 1, 0],[-sin(T), 0, cos(T)]])
Rotz = lambda T: Matrix([[cos(T), -sin(T), 0],[sin(T), cos(T), 0], [0, 0, 1]])

Hx = lambda T,dx,dy,dz: Rotx(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
Hy = lambda T,dx,dy,dz: Roty(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
Hz = lambda T,dx,dy,dz: Rotz(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def H(axis, T, dx, dy, dz):
    if axis == 'x' :
        return Hx(T,dx,dy,dz)
    elif axis == 'y' :
        return Hy(T,dx,dy,dz)
    elif axis == 'z' :
        return Hz(T,dx,dy,dz)
    elif axis == '0':
        return Hz(0,dx,dy,dz)
    else:
        raise ValueError("Soh pode ser string x, y, z ou 0!")
        
###############################################################################

def H_d(mat):
    _H_d = lambda a, alpha, d, theta: H('z',theta,0,0,0)*H('0',0,0,0,d)*H('0',0,a,0,0)*H('x',alpha,0,0,0)
    return [_H_d(*mat[i,0:4]) for i in range(mat.rows)]

class Serial(object):
    """Serial robots dynamics."""
    def __init__(self, name, ID, DH_):
        self.name = name
        self.ID = ID
        self.dof = DH_.rows
        self.DH_ = DH_

        Id = '' if ID == '' else '_' + str(ID)
        
        self.Hr__ = H_d(DH_[:,0:4])
        
        self.H__ = []
        self.H__.append(self.Hr__[0])
        for i in range(1,self.dof):
            self.H__.append(simplify(self.H__[i-1]*self.Hr__[i]))
        
        self.z__ = [Matrix([0,0,1])] + [self.H__[i][0:3,2] for i in range(self.dof)]
        self.o__ = [Matrix([0,0,0])] + [self.H__[i][0:3,3] for i in range(self.dof)]
        self.og__ = [ simplify(Matrix([(self.H__[i]*Matrix([DH_[i,4:7].T,[1]]))[0:3] ]).T) for i in range(self.dof)]
        
        self.Jv__ = [simplify( Matrix([self.z__[i].cross(self.og__[j]- self.o__[i]).T if str(DH_[i,7]) == 'R' and i <= j else ( self.z__[i].T if i <= j else zeros(1,3) ) for i in range(self.dof)]).T ) for j in range(self.dof)]
        self.Jw__ = [simplify( self.H__[j][0:3,0:3].T*Matrix([self.z__[i].T if str(DH_[i,7]) == 'R' and i <= j else zeros(1,3) for i in range(self.dof)]).T ) for j in range(self.dof)]
        self.Jg_ = Matrix( [Matrix([self.Jv__[i],self.Jw__[i]]) for i in range(self.dof)] )
        
        self.Jv_ = simplify(Matrix( [self.z__[i].cross(self.o__[self.dof] - self.o__[i]).T if str(DH_[i,7]) == 'R' else self.z__[i].T for i in range(self.dof) ]).T )
        self.Jw_ = simplify( self.H__[self.dof-1][0:3,0:3].T*Matrix( [self.z__[i].T if str(DH_[i,7]) == 'R' else zeros(1,3) for i in range(self.dof) ]).T )
        self.J_ = Matrix([self.Jv_,self.Jw_])
        
        #Dinamica
        self.I__ = [Matrix([[symbols('Jx'+str(i+1)), symbols('Jxy'+str(i+1)),symbols('Jxz'+str(i+1))],
                            [symbols('Jxy'+str(i+1)),symbols('Jy'+str(i+1)), symbols('Jyz'+str(i+1))],
                            [symbols('Jxz'+str(i+1)),symbols('Jyz'+str(i+1)),symbols('Jz'+str(i+1))]]) for i in range(self.dof)]
                            
        self.w__ = [Matrix([symbols('wy'+str(i+1)),symbols('wy'+str(i+1)),symbols('wz'+str(i+1))]) for i in range(self.dof)]        
        
        self.M__ = [diag(eye(3)*symbols('m'+str(i+1)), self.I__[i]) for i in range(self.dof)]
        self.v__ = [Matrix([zeros(3,1),self.w__[i].cross(self.I__[i]*self.w__[i])]) for i in range(self.dof)]
        
        self.M_ = diag(zeros(self.dof),*self.M__)
        self.v_ = Matrix([zeros(self.dof,1), Matrix([R.v__[i] for i in range(self.dof)]) ])

        self.C_ = Matrix([eye(self.dof),self.Jg_])              
        

#RRP        
DH_ = Matrix([
[0, +pi/2, symbols('l_1')               , symbols('theta_1')+pi/2, 0, -symbols('l_1')+symbols('lg_1'), 0                              , 'R'],
[0, +pi/2, 0                            , symbols('theta_2')+pi/2, 0, 0                              , symbols('lg_2')                , 'R'], 
[0,  0   , symbols('l_2')+symbols('d_3'), 0                      , 0, 0                              , -symbols('l_3')+symbols('lg_3'), 'P']
])
        
DH2_ = Matrix([
[0,  pi/2, 0                            , symbols('theta_1')+pi/2, 0, 0                             , symbols('lg_1'), 'R'],
[0, -pi/2, symbols('l_1')+symbols('l_2'), symbols('theta_2')     , 0, symbols('l_2')-symbols('lg_2'), 0              , 'R'], 
[0,  pi/2, 0                            , symbols('theta_3')     , 0, 0                             , symbols('lg_3'), 'R'],
[0, -pi/2, symbols('l_3')+symbols('l_4'), symbols('theta_4')     , 0, symbols('l_4')-symbols('lg_4'), 0              , 'R'] 
])

R = Serial('RRP',0,DH_)