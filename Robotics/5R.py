from Denavit import *
from math import pi as Pi

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
dC_ = SMatrix( C_.M_.diff(symbols('t')), C_.rowl_, C_.coll_ )

M_ = RR1.M_ + RR2.M_ + P.M_
v_ = RR1.v_ + RR2.v_ + P.v_
g_ = RR1.g_ + RR2.g_ + P.g_

uh_ = SMatrix( Matrix([symbols('tau_1'),symbols('tau_2')]), [RR1.dq_.rowl_[0],RR2.dq_.rowl_[0]] )
u_ = SMatrix( 0, rho_.rowl_ ) + uh_

Cch_T_u_ = Cch_.T()*u_
Z_ = SMatrix( Cch_T_u_.M_.jacobian(uh_.M_), Cch_T_u_.rowl_, uh_.rowl_ )

l0,l1,l2=0.05,0.12,0.15

E1=lambda x,y: -2*x*l1+2*l0*l1
F1=lambda x,y: -2*y*l1
G1=lambda x,y: x**2+l0**2+l1**2+y**2-2*x*l0-l2**2

E2=lambda x,y: 2*x*l1+2*l0*l1
F2=lambda x,y: -2*y*l1
G2=lambda x,y: x**2+l0**2+l1**2+y**2+2*x*l0-l2**2

theta1=lambda x,y: 2*atan((-F1(x,y)-sqrt(E1(x,y)**2+F1(x,y)**2-G1(x,y)**2))/(G1(x,y)-E1(x,y))) if G1(x,y)!=E1(x,y) else 2*atan(-E1(x,y)/F1(x,y)) 
theta2=lambda x,y: 2*atan((-F2(x,y)-sqrt(E2(x,y)**2+F2(x,y)**2-G2(x,y)**2))/(G2(x,y)-E2(x,y))) if G2(x,y)!=E2(x,y) else 2*atan(-E2(x,y)/F2(x,y))

def thetai(x,y):
    theta_11 = theta1(x,y)
    theta_21 = theta2(x,y)
    x1 = l0 + l1*cos(theta_11)
    y1 = l1*sin(theta_11)
    x2 = -l0 - l1*cos(theta_21)
    y2 = l1*sin(theta_21)
    theta_12 =  (atan2(y - y1, x - x1) - theta_11).evalf()
    theta_22 = (atan2(y - y2, x - x2) - theta_21).evalf()
    return [x, y, theta_11, theta_12, theta_21, theta_22]
    
PARAMETROS = [(symbols('l_0'), l0),
              (symbols('l_1'), l1),
              (symbols('l_2'), l2),
              (symbols('lg_1'), 0.5*l1),
              (symbols('lg_2'), 0.5*l2),
              (symbols('m_1'), 0.143),
              (symbols('m_2'), 0.171),
              (symbols('Jz_1'), (0.143*l1**2)/12.0),
              (symbols('Jz_2'), (0.171*l2**2)/12.0),
              (symbols('g'), 9.8)         
]
                  
rep = lambda vec1,vec2: [(vec1[i],vec2[i]) for i in xrange(len(vec1))]
    
C_n = lambda q_n: SMatrix( (C_.M_.subs(rep( [q_.M_[2],q_.M_[5],q_.M_[1],q_.M_[0],q_.M_[4],q_.M_[3]] ,q_n)).subs(PARAMETROS)).evalf(), C_.rowl_, C_.coll_)
dC_n = lambda q_n, p_n: SMatrix( (dC_.M_.subs( rep(list(p_.M_), p_n) ).subs(rep( [q_.M_[2],q_.M_[5],q_.M_[1],q_.M_[0],q_.M_[4],q_.M_[3]] ,q_n)).subs(PARAMETROS)).evalf(), dC_.rowl_, dC_.coll_)  
Z_n = lambda q_n: SMatrix( (Z_.M_.subs(rep( [q_.M_[2],q_.M_[5],q_.M_[1],q_.M_[0],q_.M_[4],q_.M_[3]] ,q_n)).subs(PARAMETROS)).evalf(), Z_.rowl_, Z_.coll_)  
    
#Jacobianos

def Jx(x,y):
    t1=theta1(x,y)
    t2=theta2(x,y)
    return Matrix([[x-l0-cos(t1)*l1, y-sin(t1)*l1],[x+l0+cos(t2)*l1, y-sin(t2)*l1]])
    
def Jt(x,y):
    t1=theta1(x,y)
    t2=theta2(x,y)
    return Matrix([[l1*(-y*cos(t1)+sin(t1)*(x-l0)), 0],[0,-(y*cos(t2)+sin(t2)*(x+l0))*l1]])

lx=0.25
ly=0.28
dl=0.0005
nx=int(lx/dl)
ny=int(ly/dl)
M=zeros(ny,nx)
epsx=0.0004
epst=0.00007

r = 0.085
x0 = 0.0
y0 = 0.16

v=1.0
w=v/r
p=0.6
dt=0.01
nt=int(p/dt)+1
t=[i*dt for i in xrange(nt)]

xt=[x0+r*cos(w*T) for T in t]
yt=[y0+r*sin(w*T) for T in t]
Xt=[[xt[i],yt[i]] for i in xrange(nt)]

dxt=[-w*r*sin(w*T) for T in t]
dyt=[w*r*cos(w*T) for T in t]
dXt=[[dxt[i],dyt[i]] for i in xrange(nt)]

d2xt=[-(w**2)*r*cos(w*T) for T in t]
d2yt=[-(w**2)*r*sin(w*T) for T in t]
d2Xt=[[d2xt[i],d2yt[i]] for i in xrange(nt)]

t1t=[theta1(*XT) for XT in Xt]
t2t=[theta2(*XT) for XT in Xt]

Jxt=[Jx(*XT) for XT in Xt]
Jtt=[Jt(*XT) for XT in Xt]

#qt=[thetai(*XT) for XT in Xt]
#C_t=[C_n(*XT) for XT in Xt]

q_t = []
C_t = []
dC_t = []
Z_t = []
p_t = []
dp_t = []
tau_t = []

M_n = SMatrix( M_.M_.subs(PARAMETROS), M_.rowl_, M_.coll_)
v_n = SMatrix( v_.M_.subs(PARAMETROS), v_.rowl_)
g_n = SMatrix( g_.M_.subs(PARAMETROS), g_.rowl_)

for i in xrange(nt):
    q_t.append(thetai(*Xt[i]))
    C_t.append(C_n(q_t[i]))
    Z_t.append(Z_n(q_t[i]))
    p_t.append(C_t[i]*SMatrix( Matrix(dXt[i]),ph_.rowl_))
    dC_t.append(dC_n(q_t[i],list(p_t[i].M_)))
    dp_t.append( dC_t[i]*SMatrix( Matrix(dXt[i]),ph_.rowl_) + C_t[i]*SMatrix( Matrix(d2Xt[i]),ph_.rowl_) )
    tau_t.append( Z_t[i].inv()*C_t[i].T()*( M_n*dp_t[i] + v_n + g_n ) )
    

dtheta=[Jt(*Xt[i]).inv()*Jx(*Xt[i])*Matrix(dXt[i]) for i in xrange(nt)]

dtheta_max = max([max([abs(dtheta[i][0]) for i in range(nt) ]), max([abs(dtheta[i][1]) for i in range(nt) ]) ])
dtheta_max_rpm = 60*dtheta_max/(2*Pi)

import matplotlib.pyplot as plt
import numpy as np

t_np = np.linspace(t[0], t[-1], len(t))
t1t_np = t_np.copy()
t2t_np = t_np.copy()
dtheta1_np = t_np.copy()
dtheta2_np = t_np.copy()
tau1_np = t_np.copy()
tau2_np = t_np.copy()

for i in np.arange(np.size(t_np)):
    t1t_np[i] = t1t[i]
    t2t_np[i] = t2t[i]
    dtheta1_np[i] = 60*dtheta[i][0]/(2*Pi)
    dtheta2_np[i] = 60*dtheta[i][1]/(2*Pi)
    tau1_np[i] = tau_t[i].M_[0]
    tau2_np[i] = tau_t[i].M_[1]
    

plt.figure()
plt.plot(t_np, t1t_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, t2t_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, dtheta1_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, tau1_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

plt.figure()
plt.plot(t_np, tau2_np, 'r')
#plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()

"""
for i in range(ny):
    for j in range(nx):
        x=j*dl+0.5*dl
        y=i*dl+0.5*dl
        
        if im(theta1(x,y))==0 and im(theta2(x,y))==0:
            if abs(det(Jx(x,y)))<epsx or abs(det(Jt(x,y)))<epst:
                M[i,j]=2
            else:
                M[i,j]=1
                
        if (x - x0)**2 + (y - y0)**2 < r**2:
            M[i,j]=3

M2=M.extract([ny-1-i for i in range(ny)],range(nx))
pprint(M2)

import matplotlib.pyplot as plt
import numpy as np

M3 = np.zeros((ny,nx))
for i in range(ny):
    for j in range(nx):
        M3[i,j] = int(M2[i,j])
        
plt.matshow(M3)
"""