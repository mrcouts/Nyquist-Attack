from Denavit import *

ID = 0
Id = '' if ID == '' else '_' + str(ID)

dof = 4

#Parametros geometricos
l_ = [symbols('l_' +str(i+1)) for i in range(dof)]
lg_= [symbols('lg_'+str(i+1)) for i in range(dof)]

#Parametros de inercia
m_ = [symbols('m_'+str(i+1)) for i in range(dof)]
I__ = [diag(symbols('Jx_'+str(i+1)),
            symbols('Jy_'+str(i+1)),
            symbols('Jz_'+str(i+1)) )
            for i in range(dof)]
        
DH_ = Matrix([
[0,  pi/2, 0          , symbols('theta_1'+Id)+pi/2, 0, 0           , l_[0], 'R'],
[0, -pi/2, l_[0]+l_[1], symbols('theta_2'+Id)     , 0, l_[1]-lg_[1], 0    , 'R'], 
[0,  pi/2, 0          , symbols('theta_3'+Id)     , 0, 0           , l_[2], 'R'],
[0, -pi/2, l_[2]+l_[3], symbols('theta_4'+Id)     , 0, l_[3]-lg_[3], 0    , 'R'] 
])


R = Serial('RR',ID,DH_, m_, I__, Matrix([0,-1,0]))