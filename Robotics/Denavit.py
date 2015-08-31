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
#RRP        
#DH = Matrix([
#[0, +pi/2, symbols('l_1')               , symbols('theta_1')+pi/2, 0, -symbols('l_1')+symbols('lg_1'), 0                              , 'R'],
#[0, +pi/2, 0                            , symbols('theta_2')+pi/2, 0, 0                              , symbols('lg_2')                , 'R'], 
#[0,  0   , symbols('l_2')+symbols('d_3'), 0                      , 0, 0                              , -symbols('l_3')+symbols('lg_3'), 'P']
#])
        
DH = Matrix([
[0,  pi/2, 0                            , symbols('theta_1')+pi/2, 0, 0                             , symbols('lg_1'), 'R'],
[0, -pi/2, symbols('l_1')+symbols('l_2'), symbols('theta_2')     , 0, symbols('l_2')-symbols('lg_2'), 0              , 'R'], 
[0,  pi/2, 0                            , symbols('theta_3')     , 0, 0                             , symbols('lg_3'), 'R'],
[0, -pi/2, symbols('l_3')+symbols('l_4'), symbols('theta_4')     , 0, symbols('l_4')-symbols('lg_4'), 0              , 'R'] 
])
                
def H_d(mat):
    _H_d = lambda a, alpha, d, theta: H('z',theta,0,0,0)*H('0',0,0,0,d)*H('0',0,a,0,0)*H('x',alpha,0,0,0)
    return [_H_d(*mat[i,0:4]) for i in range(mat.rows)]

#Matrizes de transformacao homogenea relativas  
vH_ = H_d(DH[:,0:4])

#Matrizes de transformacao homogenea absolutas
vI_ = []
vI_.append(vH_[0])
for i in range(1,DH.rows):
    vI_.append(simplify(vI_[i-1]*vH_[i]))
    
z_ = [Matrix([0,0,1])] + [vI_[i][0:3,2] for i in range(DH.rows)]
O_ = [Matrix([0,0,0])] + [vI_[i][0:3,3] for i in range(DH.rows)]
Og_ = [ simplify(Matrix([(vI_[i]*Matrix([DH[i,4:7].T,[1]]))[0:3] ]).T) for i in range(DH.rows)]

Jv__ = [simplify( Matrix([z_[i].cross(Og_[j]- O_[i]).T if str(DH[i,7]) == 'R' and i <= j else ( z_[i].T if i <= j else zeros(1,3) ) for i in range(DH.rows)]).T ) for j in range(DH.rows)]
Jw__ = [simplify( vI_[j][0:3,0:3].T*Matrix([z_[i].T if str(DH[i,7]) == 'R' and i <= j else zeros(1,3) for i in range(DH.rows)]).T ) for j in range(DH.rows)]
J__ = Matrix( [Matrix([Jv__[i],Jw__[i]]) for i in range(DH.rows)] )

Jv_ = simplify(Matrix( [z_[i].cross(O_[DH.rows] - O_[i]).T for i in range(DH.rows) ]).T )
Jw_ = simplify( vI_[DH.rows-1][0:3,0:3].T*Matrix( [z_[i].T for i in range(DH.rows) ]).T )
J_ = Matrix([Jv_,Jw_])