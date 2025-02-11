# Créé par christophe.frings, le 21/12/2023 en Python 3.7
from math import exp

def f(x):
    return 4/(1+exp(12-3*x))

def g(x):
    return 12*exp(12-3*x)/(1+2*exp(12-3*x)+exp(24-6*x))

def h(x):
    return 36*(exp(12-3*x)-1)/(3+3*exp(12-3*x)+exp(24-6*x)+exp(3*x-12))

f_out = []
g_out = []
h_out = []

for i in range(101):
    f_out.append((8*i/100, f(8*i/100)))
    g_out.append((8*i/100, g(8*i/100)))
    h_out.append((8*i/100, h(8*i/100)))


print("--------------------------")
print(" ".join([str(p) for p in f_out]))
print("--------------------------")
print(" ".join([str(p) for p in g_out]))
print("--------------------------")
print(" ".join([str(p) for p in h_out]))
print("--------------------------")
