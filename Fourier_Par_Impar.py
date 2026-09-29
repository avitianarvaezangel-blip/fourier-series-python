import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# variables y funciones simbolicas
x = sp.Symbol('x')

n = sp.Symbol('n', integer=True , positive=True)

f = x**2

# periodo

T = 2 * sp.pi

# Omega 0

w0 = (2 * sp.pi)/T

# f(-x)

f_x = f.subs(x, -x)

print("f(x) =", f)

print("f(-x) =", f_x)

if f == f_x:
    print("f(x) es par")
    
    a0 = (4/T)*sp.integrate(f,(x,0, T/2))
    
    an = (4/T)*sp.integrate(f * sp.cos(n * w0 * x),(x,0, T/2))
    
    bn = 0
    
    a0 = sp.simplify(a0)
    
    an = sp.simplify(an)
    
    print("a0 =", a0)
    
    print("an =", an)
    
    print("bn =", bn)

elif f == - f_x:
    print("f(x) es impar")
    
    a0 = 0
    
    an = 0
    
    bn = (4/T)*sp.integrate(f * sp.sin(n * w0 * x),(x,0, T/2))
    
    a0 = sp.simplify(a0)
    
    an = sp.simplify(an)
    
    print("a0 =", a0)
    
    print("an =", an)
    
    print("bn =", bn)
    
else:
    print("No es ninguna de las dos")
    
    
N = 5


S = a0/2 

for i in range(1, N + 1):
    
    if an != 0:
        ai = an.subs(n, i)
        
    else:
        
        ai = 0
    
    if bn != 0:
       bi = bn.subs(n, i)
       
    else:
        
        bi = 0
    
    S = S + ai * sp.cos(i * w0 * x) + bi * sp.sin(i * w0 * x)
    
S = sp.simplify(S)

print("S = ",S)


fn = sp.lambdify(x, f, 'numpy')

Sn = sp.lambdify(x, S, 'numpy')

xi = np.linspace(- np.pi, np.pi, 100)

plt.plot(xi, fn(xi), label=r'$f(x) = x^{2}$', color='b')

plt.plot(xi, Sn(xi), 
         label=r'$\frac{\pi^{2}}{3} + \sum_{n = 1}^{\infty} \frac{4(-1)^{n}}{n^2} \cos{nx}$',
         color='r')

plt.legend(fontsize=12)

plt.title("Serie de Fourier", fontsize=20)
    
plt.grid()

plt.show()
