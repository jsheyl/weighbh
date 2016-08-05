import numpy as np

x=np.logspace(-2,2,300)
s=x/13.3
xfunk=np.where(s<1,np.arccosh(1/s)/(1-s*s)**0.5,np.arccos(1/s)/(s*s-1)**0.5)
ii=((2+s*s)*xfunk-3)/(1-s*s)**2
np.savetxt('hern.dat',np.transpose([x,-2.5*np.log10(ii)+1.393240331430417589e+00+13.8,ii,s]))

