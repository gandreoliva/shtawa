import matplotlib.pyplot as plt
import numpy as np
import sys

# Plots the Fourier series approximation to the function, with the coefficients
# calculated by fourier_series.f90.
# Usage: python fourier_series.py datafile.txt

data = np.loadtxt(sys.argv[1], dtype=float)
N = data.shape[0]
fourier = data[:,0] + 1j*data[:,1]

x = np.linspace(0,2*np.pi,N)

fig, axs = plt.subplots(ncols=2,nrows=1, figsize=(8,4))
axs[0].plot(x,np.real(fourier),label=r"Re[F(f(x))]")
axs[0].plot(x,np.imag(fourier),label=r"Im[F(f(x))]")
axs[0].legend()


sampled_function = data[:,2]
inverse_fourier = data[:,3]

axs[1].plot(x,sampled_function,label="f(x)")
axs[1].plot(x,inverse_fourier,label=r"$F^{-1}[F(f(x))]$")
axs[1].legend()

plt.show()