import matplotlib.pyplot as plt
import numpy as np
import sys

# Plots the Fourier series approximation to the function, with the coefficients
# calculated by fourier_series.f90.
# Usage: python fourier_series.py datafile.txt

data = open(sys.argv[1],'r').read()
blocks = data.split("---")

coeffs = np.loadtxt(blocks[0].split("\n"), dtype=float)
N = 2*coeffs[1:,0].shape[0]+1 # N should be odd
a = coeffs[:,0]
b = coeffs[:,1]

x = np.linspace(-np.pi,np.pi,N)
#       the period is assumed to be 2*pi, try other values of N

f = np.zeros_like(x,dtype=float)

f += a[0]/2

for k in range(1,N//2):
    f += a[k]*np.cos(k*x)
    f += b[k]*np.sin(k*x)

sampled_function = np.loadtxt(blocks[1].split(), dtype=float)

plt.plot(x,sampled_function,label="original data")

plt.plot(x,f,label="series approx.")
# print(x)
plt.legend()
plt.show()