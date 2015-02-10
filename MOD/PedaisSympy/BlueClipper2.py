from sympy import *
import math
init_printing(use_unicode=True)

s = symbols('s') #Laplace variable
z = symbols('z') #Z-transform variable
iz = symbols('iz') #iz = z**-1
w = symbols('w', real=True) #Angular frequency
t = symbols('t') #Time variable
k = symbols('k') #Discrete time variable
c = symbols('c')
T = symbols('T') #Sampling period

Bilinear = [(s, c*(1-iz)/(1+iz))]

fs = 48000.0

def BilinearConstant(f0,fs):
    if f0 == 0:
        return 2*fs
    else:
        return (2*math.pi*f0)/(math.tan(math.pi*f0/fs))


#Stage 1 ---------------------------------------------------------------------

C1, R1 = symbols('C_1 R_1') #Component parameters
Z1, Z2 = symbols('Z_1 Z_2') #Component impedances
I1, I2 = symbols('I_1 I_2') #Currents through the components
Vi, V1 = symbols('V_i V_1') #Circuit node tensions

#Current Laws
CURRENTS_1 = Matrix(
    [I1 - I2]
)

#Ohm Laws
OHM_1 = [ 
    (I1,(Vi-V1)/Z1),
    (I2,(V1-0)/Z2)
]

#Impedances of the components
IMPEDANCES_1 = [
    (Z1, 1/(C1*s)),
    (Z2, R1)
]

#Component Values
SDATA_1 = [
    (C1, 47e-9),
    (R1, 0.5*470e3),
    (c, BilinearConstant(14.444,fs))
]

#Transfer Function
G1 = V1.subs(
    solve
    (
        CURRENTS_1.subs(OHM_1).subs(IMPEDANCES_1),
        V1
    )
).diff(Vi)

G1_ = together(G1.subs(Bilinear).subs(SDATA_1) )
Num1 = numer( G1_ )
Den1 = denom( G1_ )

Num1_ = Poly(Num1,iz).coeffs()
Den1_ = Poly(Den1,iz).coeffs()

a0 = Den1_[-1]

for i in range(len(Num1_)):
    Num1_[i] = Num1_[i]/a0
    
for i in range(len(Den1_)):
    Den1_[i] = Den1_[i]/a0

pprint(  Num1_ )
pprint(  Den1_ )
print(' ')

#Stage 2 ----------------------------------------------------------------------

C2, R3, R4, tau, A0, Vsat = symbols('C_2 R_3 R_4 tau A_0 Vsat') #Component parameters
A = symbols('A', cls=Function)
I3, I4, I5 = symbols('I_3 I_4 I_5') #Currents through the components
V1, V2, V3, V4 = symbols('V_1 V_2 V_3 V_4', cls=Function) #Circuit node tensions
dV3, dV4 = symbols('dV_3 dV_4', cls=Function) #Derivate of the circuit node tensions

def A(x):
    return Piecewise( (A0*( 1 - ((x)/(Vsat))**2 ) , Abs(x) < Vsat ) , (0, True) ) 

#Current Laws
CURRENTS_2 = Matrix([
    I3 - I4,
    I4 - I5
])

#AmpOp Law
AMPOPLAW_2 = Matrix([
    tau*dV3(t) + V3(t) - A(V3(t))*(V1(t)-V2(t))
])

#Ohm Laws
OHM_2 = [ 
    (I3, (V3(t)-V2(t))/R4),
    (I4, (V2(t)-V4(t))/R3),
    (I5, C2*(dV4(t)-0))  
]

#Component Values
SDATA_2 = [
    (C2, 4.7e-6),
    (R3, 2.4e3),
    (R4, 150e3),
    (tau,  1/(2*math.pi*100)),
    (A0, 10.0**(100/20)),
    (Vsat, 4.5),
    (T, 1/(8*48000.0))
]

# 0: Algebrica linear -> V2
# 1: Diferencial
# 2: Diferencial

#Eliminando a equacao algebrica linear:
V2Sol = solve(
        Matrix([CURRENTS_2,AMPOPLAW_2]).subs(OHM_2)[0,:],
        V2(t)
)

DiffEq_2 = simplify(Matrix([CURRENTS_2,AMPOPLAW_2]).subs(OHM_2)[1:3,:].subs(V2Sol)) #Non-Linear Algebric-Differencial Equation System

pprint(  Matrix([CURRENTS_2,AMPOPLAW_2]).subs(OHM_2) )
pprint(V2Sol)
pprint(DiffEq_2)

#Separando as derivadas no lado direito
X_2 = Matrix([V3(t),V4(t)]) #Nos que "tem derivadas"
dX_2 = Matrix([dV3(t),dV4(t)])
J_2 = DiffEq_2.jacobian(dX_2)
DiffEqR_2 = -J_2*X_2 #Na verdade essa eh a integral do lado direito da equacao
DiffEqL_2 = simplify(expand(DiffEq_2 - J_2*dX_2))


pprint(DiffEqR_2)
pprint(DiffEqL_2)

DifferenceEqR_2 = simplify(  DiffEqR_2.subs(t,k+1) - (T/2)*(DiffEqL_2.subs(t,k+1) ) )
DifferenceEqL_2 = simplify(  DiffEqR_2.subs(t,k) +   (T/2)*(DiffEqL_2.subs(t,k) ) )

pprint(DifferenceEqR_2)
pprint(DifferenceEqL_2)

Y_2 = Matrix([V3(t),V4(t)]) #States
U_2 = Matrix([V1(t)]) #Input

f_2 = DifferenceEqR_2
g_2 = DifferenceEqL_2

fk_2 = f_2.subs(k+1,k)
fkLin_2 = fk_2.subs(A(V3(t)),0)
fkNLin_2 = fk_2 - fkLin_2

gLin_2 = g_2.subs(A(V3(t)),0)
gNLin_2 = g_2 - gLin_2

dfdY_2 = f_2.jacobian(Y_2.subs(t,k+1)).subs(k+1,k)
dfdU_2 = f_2.jacobian(U_2.subs(t,k+1)).subs(k+1,k)
dgdY_2 = gLin_2.jacobian(Y_2.subs(t,k))
dgdU_2 = gLin_2.jacobian(U_2.subs(t,k))

dfkLindY_2 = fkLin_2.jacobian(Y_2.subs(t,k))
dfkLindU_2 = fkLin_2.jacobian(U_2.subs(t,k))

A_2 = simplify(dfdY_2)
B_2 = simplify( dgdY_2+dfdY_2-dfkLindY_2)
C_2 = simplify(-dfdU_2)
D_2 = simplify( dgdU_2+dfdU_2-dfkLindU_2)
E_2 = simplify(-fkNLin_2 + gNLin_2)

pprint(A_2)
pprint(B_2)
pprint(C_2)
pprint(D_2) #D=C
pprint(E_2) #E=0

def NA_2(v1,v3,v4):
    return (A_2.subs(SDATA_2).subs([(V1(k),v1),(V3(k),v3),(V4(k),v4)])).evalf()
    
def NB_2(v1,v3,v4):
    return (B_2.subs(SDATA_2).subs([(V1(k),v1),(V3(k),v3),(V4(k),v4)])).evalf()
    
def NC_2(v1,v3,v4):
    return (C_2.subs(SDATA_2).subs([(V1(k),v1),(V3(k),v3),(V4(k),v4)])).evalf()

    
def NAB_2(v1,v3,v4):
    return NA_2(v1,v3,v4)**-1 * NB_2(v1,v3,v4)
    
def NAC_2(v1,v3,v4):
    return NA_2(v1,v3,v4)**-1 * NC_2(v1,v3,v4)
    
    
pprint( NA_2(1,2,3) )
pprint( NB_2(1,2,3) )
pprint( NC_2(1,2,3) )
pprint( NAB_2(1,2,3) )
pprint( NAC_2(1,2,3) )



#Stage 3 ---------------------------------------------------------------------

C3, C4, R5, R6, R7, Is, Vt = symbols('C_3 C_4 R_5 R_6 R_7 I_s V_t') #Component values
I5, I6, I7, I8, I9, I10, I11, I12 = symbols('I_5 I_6 I_7 I_8 I_9 I_10 I_11 I_12') #Currents through the components
V3, V4, V5, V6, Vo = symbols('V_3 V_4 V_5 V_6 V_o', cls=Function) #Circuit node tensions
dV3, dV4, dV6 = symbols('dV_3 dV_4 dV_6', cls=Function) #Derivate of the circuit node tensions
vol = symbols('vol')

#Current Laws
CURRENTS_3 = Matrix([
    I5 - I6,
    I6 + I7 - I8 - I9,
    I9 - I10 - I11,
    I11 - I12
])

#Ohm Laws
OHM_3 = [ 
    (I5, C3*(dV3(t)-dV4(t))),
    (I6, (V4(t)-V5(t))/R5),
    (I7, Is*(exp((0-V5(t))/Vt) -1) ),
    (I8, Is*(exp((V5(t)-0)/Vt) -1) ),
    (I9, (V5(t)-V6(t))/R6),
    (I10, C4*(dV6(t)-0)),
    (I11, (V6(t)-Vo(t))/( (1-vol)*R7 )),
    (I12, (Vo(t)-0)/(vol*R7))    
]

#Component Values
SDATA_3 = [
    (C3, 4.7e-6),
    (C4, 33e-6),
    (R5, 8.2e3),
    (R6, 10e3),
    (R7, 50e3),
    (T, 1/(8*fs)),
    (Is, 2.52e-9),
    (Vt, 45.3e-3)
]
    
# 0: Diferencial
# 1: Algebrica nao linear
# 2: Diferencial
# 3: Algebrica linear -> Vo, V6    
    
#Eliminando a equacao algebrica linear:
VoSol = solve(
        (CURRENTS_3.subs(OHM_3))[3,:],
        Vo(t)
)

pprint(VoSol)
    
DiffEq = simplify((CURRENTS_3.subs(OHM_3))[0:3,:].subs(VoSol)) #Non-Linear Algebric-Differencial Equation System

#Separando as derivadas no lado direito
X = Matrix([V3(t),V4(t),V6(t)]) #Nos que "encostam" em capacitores
dX = Matrix([dV3(t),dV4(t),dV6(t)])
J = DiffEq.jacobian(dX)
DiffEqR = -J*X #Na verdade essa eh a integral do lado direito da equacao
DiffEqL = simplify(DiffEq - J*dX)

DefferentialEqR = DiffEqR[(0,2),:]
DefferentialEqL = DiffEqL[(0,2),:]

AlgebricEqR = DiffEqL[1,:]
AlgebricEqL = DiffEqR[1,:]

DifferenceEq1R = simplify(  DefferentialEqR.subs(t,k+1) - (T/2)*(DefferentialEqL.subs(t,k+1) ) )
DifferenceEq1L = simplify(  DefferentialEqR.subs(t,k) +   (T/2)*(DefferentialEqL.subs(t,k) ) )
DifferenceEq2R = simplify( AlgebricEqR.subs(t,k+1) )
DifferenceEq2L = simplify( AlgebricEqL)

DifferenceEqR = Matrix([DifferenceEq1R,DifferenceEq2R])
DifferenceEqL = Matrix([DifferenceEq1L,DifferenceEq2L])

Y = Matrix([V4(t),V5(t),V6(t)]) #States
U = Matrix([V3(t)]) #Input

f = DifferenceEqR
g = DifferenceEqL

# f[k+1] = g[k]
# f[k] + (df[k+1]/dX[k+1])[k] dX + (df[k+1]/dU[k+1])[k] dU = gNLin[k] + gLin[k]
#(df[k+1]/dX[k+1])[k] X[k+1] = -fNLin[k] -fLin[k] + gNLin[k] + gLin[k] + (df[k+1]/dX[k+1])[k] X[k] - (df[k+1]/dU[k+1])[k] dU
#gLin[k] = (dgLin[k]/dX[k])[k] X[k] + (dgLin[k]/dU[k])[k] U[k]
#fLin[k] = (dfLin[k]/dX[k])[k] X[k] + (dfLin[k]/dU[k])[k] U[k]



fk = f.subs(k+1,k)
fkLin = fk.subs(Is,0)
fkNLin = fk - fkLin

gLin = g.subs(Is,0)
gNLin = g - gLin

dfdY = f.jacobian(Y.subs(t,k+1)).subs(k+1,k)
dfdU = f.jacobian(U.subs(t,k+1)).subs(k+1,k)
dgdY = gLin.jacobian(Y.subs(t,k))
dgdU = gLin.jacobian(U.subs(t,k))

dfkLindY = fkLin.jacobian(Y.subs(t,k))
dfkLindU = fkLin.jacobian(U.subs(t,k))

#A Y(k+1) = B Y(k) + C U(k+1) + D U(k) + E

A = simplify(dfdY)
B = simplify( dgdY+dfdY-dfkLindY)
C = simplify(-dfdU)
D = simplify( dgdU+dfdU-dfkLindU)
E = simplify(-fkNLin + gNLin)

pprint(A)
pprint(B)
pprint(C)
pprint(D)
pprint(E)

#------------------------------------------------------------------------------

def NA(x):
    return (A.subs(SDATA_3).subs(V5(k),x)).evalf()
    
def NB(x):
    return (B.subs(SDATA_3).subs(V5(k),x)).evalf()
    
def NC(x):
    return (C.subs(SDATA_3).subs(V5(k),x)).evalf()
    
def ND(x):
    return (D.subs(SDATA_3).subs(V5(k),x)).evalf()
    
def NE(x):
    return (E.subs(SDATA_3).subs(V5(k),x)).evalf()

def NAB(x):
    return NA(x)**-1 * NB(x)
    
def NAC(x):
    return NA(x)**-1 * NC(x)
    
def NAD(x):
    return NA(x)**-1 * ND(x)
    
def NAE(x):
    return NA(x)**-1 * NE(x)
    
pprint(NAB(-1))
pprint(NAC(-1))
pprint(NAD(-1))
pprint(NAE(-1))

#-----------------------------Table generation---------------------------------

import sys
from mpmath import mp
mp.dps = 50

N = 40
inicio = -4.5
fim = 4.5
dx = mp.mpf((fim-inicio)/(N-1))
Idx = mp.mpf(1/dx)
print(Idx)

with open('MatrixAB_2.h', 'w') as f:

    f.write("\n".join([
        '#include <armadillo>',
        '',
        'using namespace arma;',
        '',
        '#define MatrixAB2_N {N}',
        '#define MatrixAB2_Idx {Idx}',
        '#define MatrixAB2_inicio {inicio}',
        '#define MatrixAB2_fim {fim}',
        '',
        '',
    ]).format(**locals()))

    f.write('const mat MatrixAB2[MatrixAB2_N][MatrixAB2_N][MatrixAB2_N] = \n{ \n')

    for i in range(N):
        f.write('    { \n')
        for j in range(N):
            f.write('        { \n')            
            for l in range(N):
                Matrix = NAB_2(i*dx + inicio,j*dx + inicio,l*dx + inicio)
                f.write(
                '            mat(" ' + str(Matrix[0,0]) + ' ' + str(Matrix[0,1]) + '; '
                                     + str(Matrix[1,0]) + ' ' + str(Matrix[1,1]) + ' ")')
                if l != N-1:
                    f.write(',\n')
                else:
                    f.write('\n')
            if j != N-1:            
                f.write('        },\n')
            else:
                f.write('        }\n')
        if i != N-1:            
            f.write('    },\n')
        else:
            f.write('    }\n')            
    f.write('};')
    
with open('MatrixAC_2.h', 'w') as f:

    f.write("\n".join([
        '#include <armadillo>',
        '',
        'using namespace arma;',
        '',
        '#define MatrixAC2_N {N}',
        '#define MatrixAC2_Idx {Idx}',
        '#define MatrixAC2_inicio {inicio}',
        '#define MatrixAC2_fim {fim}',
        '',
        '',
    ]).format(**locals()))

    f.write('const mat MatrixAC2[MatrixAC2_N][MatrixAC2_N][MatrixAC2_N] = \n{ \n')

    for i in range(N):
        f.write('    { \n')
        for j in range(N):
            f.write('        { \n')            
            for l in range(N):
                Matrix = NAC_2(i*dx + inicio,j*dx + inicio,l*dx + inicio)
                f.write(
                '            mat(" ' + str(Matrix[0,0]) + '; '
                                     + str(Matrix[1,0]) + ' ")')
                if l != N-1:
                    f.write(',\n')
                else:
                    f.write('\n')
            if j != N-1:            
                f.write('        },\n')
            else:
                f.write('        }\n')
        if i != N-1:            
            f.write('    },\n')
        else:
            f.write('    }\n')            
    f.write('};')

with open('MatrixAB.h', 'w') as f:

    f.write("\n".join([
        '#include <armadillo>',
        '',
        'using namespace arma;',
        '',
        '#define MatrixAB_N {N}',
        '#define MatrixAB_Idx {Idx}',
        '#define MatrixAB_inicio {inicio}',
        '#define MatrixAB_fim {fim}',
        '',
        '',
    ]).format(**locals()))

    f.write('const mat MatrixAB[MatrixAB_N] = \n{ \n')

    for i in range(N):
        Matrix = NAB(i*dx + inicio)
        f.write(
        '    mat(" ' + str(Matrix[0,0]) + ' ' + str(Matrix[0,1]) + ' ' + str(Matrix[0,2]) + '; '
                     + str(Matrix[1,0]) + ' ' + str(Matrix[1,1]) + ' ' + str(Matrix[1,2]) + '; '
                     + str(Matrix[2,0]) + ' ' + str(Matrix[2,1]) + ' ' + str(Matrix[2,2]) + '")')
        if i != N-1: f.write(',\n')
    f.write('\n};')
    
with open('MatrixAC.h', 'w') as f:

    f.write("\n".join([
        '#include <armadillo>',
        '',
        'using namespace arma;',
        '',
        '#define MatrixAC_N {N}',
        '#define MatrixAC_Idx {Idx}',
        '#define MatrixAC_inicio {inicio}',
        '#define MatrixAC_fim {fim}',
        '',
        '',
    ]).format(**locals()))

    f.write('const mat MatrixAC[MatrixAC_N] = \n{ \n')

    for i in range(N):
        Matrix = NAC(i*dx + inicio)
        f.write(
        '    mat(" ' + str(Matrix[0,0]) + '; '
                     + str(Matrix[1,0]) + '; '
                     + str(Matrix[2,0]) + '")')
        if i != N-1: f.write(',\n')
    f.write('\n};')
    
with open('MatrixAD.h', 'w') as f:

    f.write("\n".join([
        '#include <armadillo>',
        '',
        'using namespace arma;',
        '',
        '#define MatrixAD_N {N}',
        '#define MatrixAD_Idx {Idx}',
        '#define MatrixAD_inicio {inicio}',
        '#define MatrixAD_fim {fim}',
        '',
        '',
    ]).format(**locals()))

    f.write('const mat MatrixAD[MatrixAD_N] = \n{ \n')

    for i in range(N):
        Matrix = NAD(i*dx + inicio)
        f.write(
        '    mat(" ' + str(Matrix[0,0]) + '; '
                     + str(Matrix[1,0]) + '; '
                     + str(Matrix[2,0]) + '")')
        if i != N-1: f.write(',\n')
    f.write('\n};')
    
with open('MatrixAE.h', 'w') as f:

    f.write("\n".join([
        '#include <armadillo>',
        '',
        'using namespace arma;',
        '',
        '#define MatrixAE_N {N}',
        '#define MatrixAE_Idx {Idx}',
        '#define MatrixAE_inicio {inicio}',
        '#define MatrixAE_fim {fim}',
        '',
        '',
    ]).format(**locals()))

    f.write('const mat MatrixAE[MatrixAE_N] = \n{ \n')

    for i in range(N):
        Matrix = NAE(i*dx + inicio)
        f.write(
        '    mat(" ' + str(Matrix[0,0]) + '; '
                     + str(Matrix[1,0]) + '; '
                     + str(Matrix[2,0]) + '")')
        if i != N-1: f.write(',\n')
    f.write('\n};')
    
#------------------------------------------------------------------------------

def G1abs(f):
    return 20*log(Abs((G1.subs(SDATA_1)).subs(s,2*math.pi*f*I)),10).evalf()
    
print(G1abs(14.0873))


import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(1, 1000, 10)
y = x.copy()

for i in np.arange(np.size(x)):
    y[i] = G1abs(x[i])

plt.figure()
plt.plot(x, y, 'r')
plt.xscale('log')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()