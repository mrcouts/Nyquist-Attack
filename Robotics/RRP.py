from Denavit import *

dof = 3

#Matriz de parametros de Denavit-Hartemberg
fDH_ = lambda q_,l_,lg_: Matrix([
[0, +pi/2, l_[0], q_[0]+pi/2, 0, -l_[0] + lg_[0], 0            , 'R'],
[0, +pi/2, 0    , q_[1]+pi/2, 0, 0              , lg_[1]       , 'R'], 
[0,  0   , l_[1] + q_[2], 0 , 0, 0              , -l_[2]+lg_[2], 'P']
])

R = Serial('RRP', '', dof, fDH_, Matrix([0,0,-1]))