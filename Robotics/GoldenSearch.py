# python program for golden section search

from sympy import *

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
    
def newton(df,d2f,x0,tol=1e-5):
    x = x0 - df(x0)/d2f(x0)
    n = 0
    while abs(x-x0)>tol and n < 100: 
        n = n+1
        x0 = x
        x = x0 - df(x0)/d2f(x0)
    return [x, n]
    
def gnewton(f,df,d2f,x0,tol=1e-5):
    f0 = f(x0)
    x = x0 - df(x0)/d2f(x0)
    fx = f(x)
    n = 0
    fg = False
    fn = False
    while abs(x-x0)>tol and n < 100 and fx < f0: 
        n = n+1
        x0 = x
        f0 = fx
        x = x0 - df(x0)/d2f(x0)
        fx = f(x)
        fn = True 
    if fx >= f0:
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
            while fx < f0 and k < 10:
                k = k+1
                x = x0 - (2**k)*d
                fx = f(x)
            if fx < f0:
                return "Fudeu"
            else:
                sol = gss(f,min(x0,x),max(x0,x),tol)
                x = sol[0]
                n = n + sol[1]
                fg = True
    return [x,n, fn, fg]
    
X = symbols('x')
k = 100
s0 = 0.1
n = 100
T = 0.01
f0 = s0 - 0.5*T*k*tanh(n*s0)
F = 0.5*(X+0.5*T*k*tanh(n*X)-f0)**2
dF = F.diff(X)
d2F = dF.diff(X)

f=lambda x:(x-2)**2
df=lambda x:2*(x-2)
d2f=lambda x:2

g  =lambda x:F.subs(X,x).evalf()
dg =lambda x:dF.subs(X,x).evalf()
d2g=lambda x:d2F.subs(X,x).evalf()
    
x=gnewton(g,dg,d2g,s0)
print(x)
print(g(x[0]))
print(sqrt(2*g(x[0])))
