from Denavit import *

ID = ''
Id = '' if ID == '' else '_' + str(ID)

dof = 2

#Parametros geometricos
l_ = [symbols('l_' +str(i+1)) for i in range(dof)]
lg_= [symbols('lg_'+str(i+1)) for i in range(dof)]

#Parametros de inercia
m_ = [symbols('m_'+str(i+1)) for i in range(dof)]
I__ = [diag(symbols('Jx_'+str(i+1)),
            symbols('Jy_'+str(i+1)),
            symbols('Jz_'+str(i+1)) )
            for i in range(dof)]

#Matriz de parametros de Denavit-Hartemberg
DH_ = Matrix([
[l_[0], 0, 0, Function('theta_1'+Id)(t), -l_[0]+lg_[0], 0, 0, 'R'],
[l_[1], 0, 0, Function('theta_2'+Id)(t), -l_[1]+lg_[1], 0, 0, 'R']
])

R = Serial('RR',ID,DH_, m_, I__, Matrix([0,-1,0]))

SUBS = [
    (sin(R.q_[0]), symbols('s_1')),
    (cos(R.q_[0]), symbols('c_1')),
    (sin(R.q_[1]), symbols('s_2')),
    (cos(R.q_[1]), symbols('c_2')),
    (sin(R.q_[0]+R.q_[1]), symbols('s_12')),
    (cos(R.q_[0]+R.q_[1]), symbols('c_12')),
    (R.q_[0], symbols('theta_1')),
    (R.q_[1], symbols('theta_2'))]
