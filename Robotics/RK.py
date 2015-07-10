from sympy import *
import math
from numpy import dot
from numpy.linalg.linalg import norm
from numpy.linalg.linalg import inv

init_printing(use_unicode=True)

class RK(object):
    """Runge-Kutta Methods."""
    def __init__(self, method='RK4'):
        if method == 'Euler':
            self.N = 1
            self.c = [1]
        elif method == 'Heun':
            self.N = 2
            self.a = [[1.0]]
            self.c = [0.5, 0.5]
        elif method == 'Mid':
            self.N = 2
            self.a = [[0.5]]
            self.c = [0, 1.0]
        elif method == 'Ralston':
            self.N = 2
            self.a = [[(2.0/3)]]
            self.c = [0.25, 0.75]
        elif method == 'RK3':
            self.N = 3
            self.a = [[0.5],[-1.0, 2.0]]
            self.c = [1.0/6, 2.0/3, 1.0/6]
        elif method == 'RK4':
            self.N = 4
            self.a = [[0.5],[0, 0.5],[0, 0, 1.0]]
            self.c = [1.0/6, 1.0/3, 1.0/3, 1.0/6]
        elif method == 'RK4_3/8':
            self.N = 4
            self.a = [[1.0/3],[-1.0/3, 1.0],[1.0, -1.0, 1.0]]
            self.c = [1.0/8, 3.0/8, 3.0/8, 1.0/8]
        elif method == 'RK5':
            self.N = 6
            self.a = [[0.5],[3.0/16, 1.0/16],[0, 0, 0.5],[0, -3.0/16, 6.0/16, 9.0/16],[1.0/7, 4.0/7, 6.0/7, -12.0/7, 8.0/7 ]]
            self.c = [7.0/90, 0, 32.0/90, 12.0/90, 32.0/90, 7.0/90]
        else:
            self.N = 4
            self.a = [[0.5],[0, 0.5],[0, 0, 1.0]]
            self.c = [1.0/6, 1.0/3, 1.0/3, 1.0/6]
        
        if self.N > 1:
            self.b = [sum(self.a[i]) for i in range(self.N-1)]
            
    def RKX(self, f, t0, Y0, n, tf):
        h = (tf-t0)/(1.0*n)
        Y = zeros(len(Y0),n+1)
        Y[:,0] = Y0
        t = t0
        k = [zeros(len(Y0),1) for i in range(self.N)]
        for i in range(n):
            k[0] = f(t, Y[:,i])
            for j in range(self.N-1):
                k[j+1] = f(t + h*self.b[j], Y[:,i] + h*dot(self.a[j],k[0:j+1]) )
            Y[:,i+1] = Y[:,i] + h*dot(self.c,k) 
            t += h
        return Y
        
gr=(math.sqrt(5)-1)/2
def gss(f,a,b,tol=1e-5):
    '''golden section search
to find the minimum of f on [a,b]
f: a strictly unimodal function on [a,b]

example:
>>> f=lambda x:(x-2)**2
>>> x=gss(f,1,5)
>>> x
2.000009644875678
'''
    c=b-gr*(b-a)
    d=a+gr*(b-a)
    n = 0
    while abs(c-d)>tol and n < 100: 
        n = n+1
        fc=f(c);fd=f(d)
        if fc<fd:
            b=d
            d=c  #fd=fc;fc=f(c)
            c=b-gr*(b-a)
        else:
            a=c
            c=d  #fc=fd;fd=f(d)
            d=a+gr*(b-a)
    return [(b+a)/2, n]
    
def gnewton(f,x0,tol=1e-5):
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0 or f(x0) == 0:
        return [x0,f(x0), 0, False, False]
    f0 = f(x0)
    k1 = F(x0,x0)
    x = x0 + 1.0*k1
    fx = f(x)
    n = 0
    fg = False
    fn = False
    while abs(x-x0)>tol and fx != 0 and n < 100 and fx < f0:
        n = n+1
        x0 = x
        f0 = fx
        k1 = F(x0,x0)
        k2 = F(x0 + 0.5*k1,x0)
        k3 = F(x0 + 0.5*k2,x0)
        k4 = F(x0 + 1.0*k3,x0)
        x = x0 + 1.0*(k1+2*k2+2*k3+k4)/6
        fx = f(x)
        fn = True 
    if fx >= f0 and abs(x-x0)>tol and fx != 0:
        if df(x0)*(x-x0) < 0:
            sol = gss(f,min(x0,x),max(x0,x),tol)
            x = sol[0]
            n = n + sol[1]
            fg = True
        else:
            d = x - x0
            k = 0
            x = x0 - (2**k)*d
            fx = f(x)
            while df(x)*df(x0) > 0 and k < 100:
                k = k+1
                x = x0 - (2**k)*d
            if df(x)*df(x0) > 0:
                return "Fudeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeu"
            else:
                sol = gss(f,min(x0,x),max(x0,x),tol)
                x = sol[0]
                n = n + sol[1]
                fg = True
    return [x, f(x), n, fn, fg]
        
def rk4newton(f,x0,tol=1e-5):
    X = Matrix([ symbols('x_'+str(i+1)) for i in range(len(x0)) ])
    J = lambda x: f(X).jacobian(X).subs([(X[i],x[i]) for i in range(len(x0))]).evalf()
    F = lambda x,x0: -inv(J(x))*f(x0)
    if norm(f(x0),1) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,11):
            k1 = F(x0,x0)
            k2 = F(x0 + 0.5*k1,x0)
            k3 = F(x0 + 0.5*k2,x0)
            k4 = F(x0 + 1.0*k3,x0)
            x = x0 + 1.0*(1.0*k1 + 2.0*k2+ 2.0*k3 + 1.0*k4)/6.0
            if norm(f(x),1) > norm(f(x0),1):
                s = x - x0
                f2 = lambda Y: (f(Y).T*f(Y))[0]
                f2_= lambda alpha:f2(x0 + alpha*s)
                alpha = gnewton(f2_,0)[0]
                x = x0 + alpha*s
            if norm(x-x0, 1)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
        
class TR(object):
    x = symbols('x')
    """Trapezoidal Rule Method."""
    def __init__(self, method='RK5'):
        self.RK = RK(method)
        
    def TRX(self, f, t0, Y0, n, tf):
        h = (tf-t0)/(1.0*n)
        Yrk = zeros(len(Y0),n+1)
        Y = zeros(len(Y0),n+1)
        Yrk[:,0] = Y0
        Y[:,0] = Y0
        t = t0
        for i in range(n):
            Yrk[:,i+1] = self.RK.RKX(f, t, Yrk[:,i],1,t+h)[:,1]
            F0 = Y[:,i] + 0.5*h*f(t,Y[:,i])
            F = lambda Y: Y - 0.5*h*f(t+h, Y) - F0
            #Y[:,i+1] = rk4newton(F, Yrk[:,i+1])[0]
            if( norm(F(Yrk[:,i+1]),1) < norm(F(Y[:,i]),1) ):
                sol = rk4newton(F, Yrk[:,i+1])
            else:
                sol = rk4newton(F, Y[:,i])
            Y[:,i+1] = sol[0]
            print(sol[2])
            Yrk[:,i+1] = Y[:,i+1]
            t += h
        return Y
    
f = lambda t,Y: Matrix([
Y[1],
-20*Y[1] - 80.0*tanh(100.0*Y[1]) - 100*Y[0]+100*sin(t)
])

t0 = 0
tf = 20
n = 400  
Y = TR('RK5').TRX(f, t0, Matrix([1,1]), n, tf )

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(t0, tf, n)
y = x.copy()

for i in np.arange(np.size(x)):
    y[i] = Y[0,i]

plt.figure()
plt.plot(x, y, 'r')
plt.xlabel('x')
plt.ylabel('y')
plt.title('title')
plt.show()