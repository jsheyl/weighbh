import numpy as np
import matplotlib.pyplot as plt
import sys

for f in sys.argv[1:]:
    m,x,y,z,vx,vy,vz = np.loadtxt(f,unpack=True,skiprows=1,usecols=[0,1,2,3,4,5,6])
    r2=x*x+y*y+z*z
    r=r2**0.5
    rs=np.sort(r)
    step=len(rs)/32
    i=np.arange(0,len(rs),step)
    
    plt.plot(np.log10(rs[i+step/2]),
             -np.log10(rs[i-1+step]-rs[i])-2*np.log10(rs[i+step/2]))

#    plt.plot(np.log10(rs),np.log10(np.linspace(1,len(rs),len(rs))))
    
plt.plot(np.log10(rs),-np.log10(rs)-3*np.log10(rs+1)+1.9)
plt.xlim(-1.5,3)
plt.ylim(-10,3)
plt.show()
