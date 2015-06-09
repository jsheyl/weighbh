import numpy as np

i1=np.random.random(1000)*1.5707963268
i2=np.random.random(1000)*1.5707963268

c3i1=np.cos(i1)**3
c3i2=np.cos(i2)**3
angterm=c3i2*c3i2+c3i1*c3i1-2*c3i1*c3i2*np.cos(0.7853981634)

angterm=np.sort(angterm)

print ('mean = %g ( %g %g %g )' % (np.mean(angterm),angterm[50],angterm[500],angterm[-50]))
