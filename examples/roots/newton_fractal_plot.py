import matplotlib.pyplot as plt
import numpy as np
import sys

figfilename = "data/newton_fractal_plot.png"
colors = np.fromfile(sys.argv[1],dtype=np.int32)
colors = colors.reshape(int(sys.argv[2]),int(sys.argv[3]))

plt.imshow(colors)
plt.colorbar()
plt.savefig(figfilename,dpi=200)
print(f" Figure saved in {figfilename}")