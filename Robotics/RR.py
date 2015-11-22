from Denavit import *

ID = ''
Id = '' if ID == '' else '_' + str(ID)

dof = 2

#Parametros geometricos
l_ = [symbols('l_' +str(i+1)) for i in xrange(dof)]
lg_= [symbols('lg_'+str(i+1)) for i in xrange(dof)]

#Parametros de inercia
m_ = [symbols('m_'+str(i+1)) for i in xrange(dof)]
I__ = [diag(symbols('Jx_'+str(i+1)),
            symbols('Jy_'+str(i+1)),
            symbols('Jz_'+str(i+1)) )
            for i in xrange(dof)]

#Parametros de atrito
c_     = [symbols('b_'    +str(i+1)) for i in xrange(dof)]
gamma_ = [symbols('gamma_'+str(i+1)) for i in xrange(dof)]
n_     = [symbols('n_'    +str(i+1)) for i in xrange(dof)]

#Matriz de parametros de Denavit-Hartemberg
fDH_ = lambda q_,l_,lg_: Matrix([
[l_[0], 0, 0, q_[0], -l_[0]+lg_[0], 0, 0, 'R'],
[l_[1], 0, 0, q_[1], -l_[1]+lg_[1], 0, 0, 'R']
])

R = Serial('RR', ID, dof, fDH_, Matrix([0,-1,0]), l_, lg_, m_, I__, c_, gamma_, n_)

"""SUBS = [
    (sin(R.q_[0]), symbols('s_1')),
    (cos(R.q_[0]), symbols('c_1')),
    (sin(R.q_[1]), symbols('s_2')),
    (cos(R.q_[1]), symbols('c_2')),
    (sin(R.q_[0]+R.q_[1]), symbols('s_12')),
    (cos(R.q_[0]+R.q_[1]), symbols('c_12')),
    (R.q_[0], symbols('theta_1')),
    (R.q_[1], symbols('theta_2'))]"""
