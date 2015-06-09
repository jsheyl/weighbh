import numpy as np
import string

def dollarconv(thestring):
    return(float(string.translate(thestring,None,'$\\')))

convlist = { 5: dollarconv, 6 : dollarconv, 7 : dollarconv }
Name,Morphology,GAL_RA,GAL_DEC,dist,x,y,z = np.loadtxt('grebel.tex',dtype=[('f0',str),('f1',str),('f2',str),('f3',str),('f4',str),('f5',float),('f6',float),('f7',float)], skiprows=13, unpack = True, delimiter='&', converters = convlist, comments='#')

for m in zip(Name,x,y,z):
    print("%s %g %g %g" % m)
