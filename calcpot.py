import numpy as np
import matplotlib.pyplot as plt
nrand=1000000
mbh=1e-3
b=np.sqrt(mbh)
rbh=b/(1-b)
print("rbh = %g" % rbh)
rbh2=rbh*rbh
massfrac=np.random.uniform(size=nrand)*mbh
b=np.sqrt(massfrac)
dbh=b/(1-b)
muang=np.random.uniform(low=-1,high=1,size=nrand)
xarr=np.linspace(0,rbh*0.01,101)
pot=0*xarr
for i,x in enumerate(xarr):
    r2=dbh**2+x**2+2*dbh*x*muang
    r2=np.where(r2>rbh2,r2,rbh2)
    pot[i]=np.mean(1/(np.sqrt(r2)+1))

plt.plot(xarr,(pot[0]-pot)/ (0.15*xarr**2/mbh**0.5) )
plt.show()    
