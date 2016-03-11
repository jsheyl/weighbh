import numpy as np
import matplotlib.pyplot as plt
import sys

for f in sys.argv[1:]:
    m,x,y,z,vx,vy,vz = np.loadtxt(f,unpack=True,skiprows=1,usecols=[0,1,2,3,4,5,6])
    r2=x*x+y*y+z*z
    r=r2**0.5
    rs=np.sort(r)
    i=np.arange(0,len(rs),1000)
    
    plt.plot(np.log10(rs[i+5]),np.log10((rs[i+9]-rs[i])/rs[i+5]**2))
    
plt.show()
