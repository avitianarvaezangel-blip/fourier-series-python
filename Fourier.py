import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


# Variable x simbolica
x = sp.Symbol('x')

# Variable n simbolica y entera
n = sp.Symbol('n', integer=True)

# funcion simbolica
f = sp.exp(x)

# Periodo
T = 2 * sp.pi

# Omega 0
w0 = (2 * sp.pi)/T

# Coeficientes de Fourier

a0 = (2/T)*sp.integrate(f,(x, -sp.pi, sp.pi))

an = (2/T)*sp.integrate(f * sp.cos(n * w0 * x),(x, -T/2, T/2))

bn = (2/T)*sp.integrate(f * sp.sin(n * w0 * x),(x, -T/2, T/2))

# Coeficientes Simplificados

a0 = sp.simplify(a0)
an = sp.simplify(an)
bn = sp.simplify(bn)

# Moatrar coeficientes 
print("a0 =", a0)
print("an =", an)
print("bn =", bn)

# Serie de Fourier

# numero de terminos
N = 20

# Inicio con el termino a0/2

S = a0/2 

for i in range(1, N + 1):
    
    ai = an.subs(n, i)
    
    bi = bn.subs(n, i)
    
    S = S + ai * sp.cos(i * w0 * x) + bi * sp.sin(i * w0 * x)
    
S = sp.simplify(S)

print("S = ",S)

fn = sp.lambdify(x, f, 'numpy')

Sn = sp.lambdify(x, S, 'numpy')

xi = np.linspace(- np.pi, np.pi, 100)

plt.plot(xi, fn(xi), label=r'$f(x) = e^{x}$', color='b')

plt.plot(xi, Sn(xi), 
         label=r'$ \frac{2 \sinh{\pi}}{\pi} \left\{ \frac{1}{2} + \sum_{1}^{\infty} \frac{(-1)^n}{n^2 + 1} [\cos(n x) - n \sin(nx)] \right\}$',
         color='r')

plt.legend(fontsize=12)

plt.title("Serie de Fourier", fontsize=20)
    
plt.grid()

plt.show()
    




