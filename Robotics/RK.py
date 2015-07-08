from sympy import *
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
    
f = lambda t,Y: Matrix([
Y[1],
-20*Y[1] - 40.0*tanh(100.0*Y[1]) - 100*Y[0]+100*sin(t)
])

t0 = 0
tf = 20
n = 2000  
Y = RK().RKX(f, t0, Matrix([1,1]), n, tf )
pprint(Y)

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