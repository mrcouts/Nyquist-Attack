from Denavit import *

ID = ''
Id = '' if ID == '' else '_' + str(ID)

dof = 3

#Parametros geometricos
l_ = [symbols('l_' +str(i+1)) for i in range(dof)]
lg_= [symbols('lg_'+str(i+1)) for i in range(dof)]

#Parametros de inercia
m_ = [symbols('m_'+str(i+1)) for i in range(dof)]
I__ = [diag(symbols('Jx_'+str(i+1)),
            symbols('Jy_'+str(i+1)),
            symbols('Jz_'+str(i+1)) )
            for i in range(dof)]

#RRP        
DH_ = Matrix([
[0, +pi/2, l_[0]                        , Function('theta_1'+Id)(t)+pi/2, 0, -l_[0] + lg_[0], 0            , 'R'],
[0, +pi/2, 0                            , Function('theta_2'+Id)(t)+pi/2, 0, 0              , lg_[1]       , 'R'], 
[0,  0   , l_[1] + Function('d_3'+Id)(t), 0                             , 0, 0              , -l_[2]+lg_[2], 'P']
])

R = Serial('RRP',ID,DH_, m_, I__, Matrix([0,0,-1]))