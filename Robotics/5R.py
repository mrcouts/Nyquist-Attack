from Denavit import *

dof = 2

class PontualBody2D(object):
    """Serial robots dynamics."""
    def __init__(self):
        self.q_ = SMatrix( Matrix([Function('x')(t),Function('y')(t)]) )
        self.dq_ = SMatrix( self.q_.M_.diff(t) )
        self.C_ = SMatrix(1, self.dq_.rowl_, self.dq_.rowl_)
        self.M_ = SMatrix( zeros(2), self.dq_.M_, self.dq_.M_ )
        self.v_ = SMatrix( zeros(2,1), self.dq_.M_)
        self.g_ = SMatrix( zeros(2,1), self.dq_.M_)

#Matriz de parametros de Denavit-Hartemberg
fDH_ = lambda q_,l_,lg_: Matrix([
[l_[0], 0, 0, q_[0], -l_[0]+lg_[0], 0, 0, 'R'],
[l_[1], 0, 0, q_[1], -l_[1]+lg_[1], 0, 0, 'R']
])

RR1 = Serial('RR1', '1', dof, fDH_, Matrix([0,-1,0]))
RR2 = Serial('RR2', '2', dof, fDH_, Matrix([0,-1,0]))
P   = PontualBody2D()

qh_ = P.q_
qo_ = RR1.q_ + RR2.q_
q_ = qh_ + qo_ 

ph_ = P.dq_
po_ = RR1.p_ + RR2.p_
p_ = ph_ + po_

rhoh_ = P.dq_
rhoo_ = RR1.dq_ + RR2.dq_
rho_ = rhoh_ + rhoo_

H1 = H('y',0,symbols('l_0'),0,0)
H2 = H('y',pi,-symbols('l_0'),0,0)

#Vinculos de posicao:
_q_ = SMatrix( Matrix([ qh_.M_ - (H1*Matrix([RR1.o__[2], Matrix([1]) ]))[0:2,0], qh_.M_ - (H2*Matrix([RR2.o__[2], Matrix([1]) ]))[0:2,0] ]) , range(4) )

Jh_ = SMatrix( _q_.M_.jacobian(qh_.M_), _q_.rowl_, rhoh_.rowl_ )
Jo_ = SMatrix( _q_.M_.jacobian(qo_.M_), _q_.rowl_, rhoo_.rowl_ )

Cch_ = SMatrix(1, ph_.rowl_, ph_.rowl_ ) + (-1*Jo_.inv()*Jh_).simplify()
C_ = ( (RR1.C_ + RR2.C_ + P.C_)*Cch_ ).simplify()

M_ = RR1.M_ + RR2.M_ + P.M_
v_ = RR1.v_ + RR2.v_ + P.v_
g_ = RR1.g_ + RR2.g_ + P.g_