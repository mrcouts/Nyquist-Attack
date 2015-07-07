# python program for golden section search

from sympy import *

gr=(math.sqrt(5)-1)/2
X = symbols('x')
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
    
def newton(f,x0,tol=1e-5):
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,101):
            k1 = F(x0,x0)
            x = x0 + 1.0*k1
            if abs(x-x0)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
    
def heunnewton(f,x0,tol=1e-5):
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,101):
            k1 = F(x0,x0)
            k2 = F(x0 + 1.0*k1,x0)
            x = x0 + 1.0*(1.0*k1 +1.0*k2)/2.0
            if abs(x-x0)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
    
def midnewton(f,x0,tol=1e-5):
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,101):
            k1 = F(x0,x0)
            k2 = F(x0 + 1.0*k1,x0)
            x = x0 + 1.0*k2
            if abs(x-x0)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
    
def ralstonnewton(f,x0,tol=1e-5):
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,101):
            k1 = F(x0,x0)
            k2 = F(x0 + (2.0/3.0)*k1,x0)
            x = x0 + 1.0*(1.0*k1 +3.0*k2)/4.0
            if abs(x-x0)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
    
def rk3newton(f,x0,tol=1e-5):
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,101):
            k1 = F(x0,x0)
            k2 = F(x0 + 0.5*k1,x0)
            k3 = F(x0 - 1.0*k1 + 2.0*k2,x0)
            x = x0 + 1.0*(1.0*k1 + 4.0*k2+ 1.0*k3)/6.0
            if abs(x-x0)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
    
def rk4newton(f,x0,tol=1e-5):
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,101):
            k1 = F(x0,x0)
            k2 = F(x0 + 0.5*k1,x0)
            k3 = F(x0 + 0.5*k2,x0)
            k4 = F(x0 + 1.0*k3,x0)
            x = x0 + 1.0*(1.0*k1 + 2.0*k2+ 2.0*k3 + 1.0*k4)/6.0
            if abs(x-x0)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
    
def rk438newton(f,x0,tol=1e-5):
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,101):
            k1 = F(x0,x0)
            k2 = F(x0 + (1.0/3.0)*k1,x0)
            k3 = F(x0 - (1.0/3.0)*k1 + 1.0*k2,x0)
            k4 = F(x0 + 1.0*k1 - 1.0*k2 + 1.0*k3,x0)
            x = x0 + 1.0*(1.0*k1 + 3.0*k2+ 3.0*k3 + 1.0*k4)/8.0
            if abs(x-x0)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
    
def rk4newton2(f,x0, N=2,tol=1e-5):
    h = 1.0/N
    X = symbols('x')
    df = lambda x: f(X).diff(X).subs(X,x).evalf()
    d2f= lambda x: df(X).diff(X).subs(X,x).evalf()
    F = lambda x,x0: -df(x0)/d2f(x)
    if df(x0) == 0:
        return [x0,f(x0), 0]
    else:
        for n in range(1,101):
            k1 = F(x0,x0)
            k2 = F(x0 + 0.5*h*k1,x0)
            k3 = F(x0 + 0.5*h*k2,x0)
            k4 = F(x0 + 1.0*h*k3,x0)
            x = x0 + h*(k1+2*k2+2*k3+k4)/6
            for j in range(N-1):
                k1 = F(x,x0)
                k2 = F(x + 0.5*h*k1,x0)
                k3 = F(x + 0.5*h*k2,x0)
                k4 = F(x + 1.0*h*k3,x0)
                x = x + h*(k1+2*k2+2*k3+k4)/6
            if abs(x-x0)<tol:
                break
            else:
                x0 = x
    return [x,f(x), n]
    
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

#ds/dt = - k*sign(s)
#Trapezios: s_(k+1) + 0.5*T*k*sign(s_(k+1)) = s_(k) - 0.5*T*k*sign(s_(k))
#Cade passo de integracao pede a resolucao de um sistema nao linear
#Transformamos o problema em um equivalente de minimizacao,
#substituinfo sign(s) por tanh(n*s)
  

k = 1000*rand()
s0 = rand() - 0.5
n = 10000*rand()
T = 0.1*rand()
f0 = s0 - 0.5*T*k*tanh(n*s0)
F = 0.5*(X+0.5*T*k*tanh(n*X)-f0)**2

f=lambda x:(x-2)**8
g=lambda x:F.subs(X,x).evalf()
    
x=midnewton(f,4)
print(x)
