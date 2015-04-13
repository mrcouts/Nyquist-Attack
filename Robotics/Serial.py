from sympy import *
import math
init_printing(use_unicode=True)

def Rotx(T):
    return Matrix([[1, 0, 0],[0, cos(T),-sin(T)],[0, sin(T), cos(T)]])

def Roty(T):
    return Matrix([[ cos(T), 0, sin(T)],[0, 1, 0],[-sin(T), 0, cos(T)]])

def Rotz(T):
    return Matrix([[cos(T),-sin(T),0],[sin(T), cos(T),0], [0, 0, 1]])

def Hx(T,dx,dy,dz):
    return Rotx(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def Hy(T,dx,dy,dz):
    return Roty(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def Hz(T,dx,dy,dz):
    return Rotz(T).col_insert(3,Matrix([[dx,dy,dz]]).T).row_insert(3,Matrix([[0,0,0,1]]))
    
def vec2h(v):
    return v.row_insert(3,Matrix([[1]]))
    
def h2vec(hv):
    return hv.row(range(3))
    
def h2rot(M):
    return M.row(range(3)).col(range(3))
    
def nrows(M):
    return M.shape[0]
    
def ncols(M):
    return M.shape[1]
    
T = symbols('T', cls=Function)
t = symbols('t')
theta = symbols('theta', cls=Function)
x = symbols('x', cls=Function)
y = symbols('y', cls=Function)
z = symbols('z', cls=Function)