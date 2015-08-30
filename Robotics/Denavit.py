from sympy import *
import math
from operator import concat
init_printing(use_unicode=True)

###############################################################################

t = symbols('t')

def Rotx(T):
    return Matrix([[1, 0, 0],[0, cos(T),-sin(T)],[0, sin(T), cos(T)]])

def Roty(T):
    return Matrix([[cos(T), 0, sin(T)],[0, 1, 0],[-sin(T), 0, cos(T)]])

def Rotz(T):
    return Matrix([[cos(T), -sin(T), 0],[sin(T), cos(T), 0], [0, 0, 1]])

def Hx(T,dx,dy,dz):
    return Rotx(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def Hy(T,dx,dy,dz):
    return Roty(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def Hz(T,dx,dy,dz):
    return Rotz(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
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
        
H_d = lambda a, alpha, d, theta: H('z',theta,0,0,0)*H('0',0,0,0,d)*H('0',0,a,0,0)*H('x',alpha,0,0,0)
H_d2 = lambda lista:H_d(lista[0],lista[1],lista[2],lista[3])
H_d3 = lambda mat:[H_d2(mat[i,:]) for i in range(mat.rows)]
    

#a = Matrix([0,0,0,0,0,0])
#alpha = Matrix([pi/2,-pi/2,pi/2,-pi/2,pi/2,0])
#d = Matrix([0,symbols('l_1')+symbols('l_2'),0,symbols('l_3')+symbols('l_4'),0,symbols('l_5')+symbols('l_6')])
#theta = Matrix([symbols('theta_1')+pi/2,symbols('theta_2'),symbols('theta_3'),symbols('theta_4'),symbols('theta_5'),symbols('theta_6')])

#a = Matrix([0,0,0,0])
#alpha = Matrix([pi/2,-pi/2,pi/2,-pi/2])
#d = Matrix([0,symbols('l_1')+symbols('l_2'),0,symbols('l_3')+symbols('l_4')])
#theta = Matrix([symbols('theta_1')+pi/2,symbols('theta_2'),symbols('theta_3'),symbols('theta_4')])
#cx = Matrix([0,0,0,0])
#cy = Matrix([0,symbols('l_2')-symbols('lg_2'),0,symbols('l_4')-symbols('lg_4')])
#cz = Matrix([symbols('lg_1'),0,symbols('lg_3'),0])
#Denavit_matrix = Matrix([a.T,alpha.T,d.T,theta.T, cx.T, cy.T, cz.T]).T

Denavit_matrix = Matrix([
[0,  pi/2, 0                            , symbols('theta_1')+pi/2, 0, 0                             , symbols('lg_1')],
[0, -pi/2, symbols('l_1')+symbols('l_2'), symbols('theta_2')     , 0, symbols('l_2')-symbols('lg_2'), 0              ], 
[0,  pi/2, 0                            , symbols('theta_3')     , 0, 0                             , symbols('lg_3')],
[0, -pi/2, symbols('l_3')+symbols('l_4'), symbols('theta_4')     , 0, symbols('l_4')-symbols('lg_4'), 0              ] 
])

#Matrizes de transformacao homogenea relativas  
vH_ = H_d3(Denavit_matrix[:,0:4])

#Matrizes de transformacao homogenea absolutas
vI_ = []
vI_.append(vH_[0])
for i in range(1,Denavit_matrix.rows):
    vI_.append(simplify(vI_[i-1]*vH_[i]))
    
z_ = [Matrix([0,0,1])] + [vI_[i][0:3,2] for i in range(Denavit_matrix.rows)]
O_ = [Matrix([0,0,0])] + [vI_[i][0:3,3] for i in range(Denavit_matrix.rows)]
Og_ = [ simplify(Matrix([(vI_[i]*Matrix([Denavit_matrix[i,4:7].T,[1]]))[0:3] ]).T) for i in range(Denavit_matrix.rows)]

Jv__ = [  simplify( Matrix([z_[i].cross(Og_[j]- O_[i]).T if i <= j else zeros(1,3)  for i in range(Denavit_matrix.rows)]).T ) for j in range(Denavit_matrix.rows)]
Jw__ = [  simplify( vI_[j][0:3,0:3].T*Matrix([z_[i].T if i <= j else zeros(1,3)  for i in range(Denavit_matrix.rows)]).T ) for j in range(Denavit_matrix.rows)]
J__ = Matrix( [Matrix([Jv__[i],Jw__[i]]) for i in range(Denavit_matrix.rows)] )

dO_ = [simplify(O_[Denavit_matrix.rows] - O_[i]) for i in range(Denavit_matrix.rows)]

Jv_ = simplify(Matrix( [z_[i].cross(dO_[i]).T for i in range(Denavit_matrix.rows) ]).T )
Jw_ = simplify( vI_[Denavit_matrix.rows-1][0:3,0:3].T*Matrix( [z_[i].T for i in range(Denavit_matrix.rows) ]).T )

J_ = Matrix([Jv_,Jw_])