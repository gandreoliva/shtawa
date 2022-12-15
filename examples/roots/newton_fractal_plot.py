import matplotlib.pyplot as plt
import numpy as np
import sys

figfilename = "data/newton_fractal_plot.png"
colors = np.fromfile(sys.argv[1],dtype=np.int32)
nx = int(sys.argv[2])
ny = int(sys.argv[3])
colors = colors.reshape(nx,ny)

l = float(sys.argv[4])
x = np.linspace(-l,l,nx)
y = np.linspace(-l,l,ny)

plt.xlabel("Re(z)")
plt.ylabel("Im(z)")
plt.pcolormesh(x,y,colors,shading='nearest')
plt.colorbar(label="Solution number")
plt.savefig(figfilename,dpi=200)
print(f" Figure saved in {figfilename}")