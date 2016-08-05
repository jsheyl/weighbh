import numpy as np

levels=np.loadtxt('ds9_contour_levels')
with open('ds9_contour','r') as cfile:
    x=[]
    y=[]
    i=0
    for line in cfile:
        a=line.split()
        if (len(a)>1):
            x.append(float(a[0]))
            y.append(float(a[1]))
        else:
            xm=np.mean(x)
            ym=np.mean(y)
            x2=np.mean((x-xm)**2)
            y2=np.mean((y-ym)**2)
            xy=np.mean((x-xm)*(y-ym))
            m=np.matrix([[x2, xy],[xy,y2]])
            ee,v=np.linalg.eig(m)
            ee=np.sort(ee)
            print('%10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g %10.8g' % (xm,ym,x2,y2,xy,levels[i],ee[0]**0.5,ee[1]**0.5,(ee[0]*ee[1])**0.25))
            i=i+1
            x=[]
            y=[]
            

