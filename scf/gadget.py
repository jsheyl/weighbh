import numpy as np


def readsnap(filename):
    with open(filename,"r") as f:
        snap={}
        np.fromfile(f,dtype='uint32',count=1)
        snap['npart']=np.fromfile(f,dtype='int32',count=6)
        snap['mass']=np.fromfile(f,dtype='float',count=6)
        snap['time']=np.fromfile(f,dtype='float',count=1)[0]
        snap['redshift']=np.fromfile(f,dtype='float',count=1)[0]
        snap['flag_sfr']=np.fromfile(f,dtype='int32',count=1)[0]
        snap['flag_feedback']=np.fromfile(f,dtype='int32',count=1)[0]
        snap['npartTotal']=np.fromfile(f,dtype='int32',count=6)
        snap['flag_cooling']=np.fromfile(f,dtype='int32',count=1)[0]
        snap['num_files']=np.fromfile(f,dtype='int32',count=1)[0]
        snap['boxsize']=np.fromfile(f,dtype='float',count=1)[0]
        snap['Omega0']=np.fromfile(f,dtype='float',count=1)[0]
        snap['OmegaLambda']=np.fromfile(f,dtype='float',count=1)[0]
        snap['HubbleParam']=np.fromfile(f,dtype='float',count=1)[0]
        np.fromfile(f,dtype='int8',
                    count=256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8)
        np.fromfile(f,dtype='uint32',count=1)

        totpart=np.sum(snap['npart'])

        np.fromfile(f,dtype='uint32',count=1)
        pos=np.fromfile(f,dtype='float32',count=totpart*3)
        snap['Position']=np.reshape(pos,[totpart,3])
        np.fromfile(f,dtype='uint32',count=1)

        np.fromfile(f,dtype='uint32',count=1)
        vel=np.fromfile(f,dtype='float32',count=totpart*3)
        snap['Velocity']=np.reshape(vel,[totpart,3])
        np.fromfile(f,dtype='uint32',count=1)

        np.fromfile(f,dtype='uint32',count=1)
        snap['ID']=np.fromfile(f,dtype='int32',count=totpart)
        np.fromfile(f,dtype='uint32',count=1)
        
        return snap

class Snapshot:
    def __init__(self):
        self.snap = {}
        self.filename = []

    def __init__(self,filename):
        self.snap = readsnap(filename)
        self.filename = filename
    
    def loadsnap(self,filename):
        self.snap = readsnap(filename)
        self.filename = filename
    
    def centerofmass(self):
        pos=self.snap['Position']
        vel=self.snap['Velocity']
        mass=self.snap['mass']
        npart=self.snap['npart']
        cmpos=np.zeros(3)
        cmvel=np.zeros(3)
        totmass=0
        startpos=0
        for i in range(0,5):
          endpos=startpos+npart[i]
          for j in range(0,3):
              cmpos[j]+=mass[i]*np.sum(pos[startpos:endpos,j])
              cmvel[j]+=mass[i]*np.sum(vel[startpos:endpos,j])
          startpos=endpos
        totmass=np.sum(mass*npart)
        cmpos/=totmass
        cmvel/=totmass
        return cmpos, cmvel
        
    def __str__(self):
        return str(self.snap)

        
import sys

def main():
    for a in sys.argv[1:]:
        bob=Snapshot(a)
        cmr, cmv = bob.centerofmass()
        print (a+' %g %g %g %g %g %g %g %g %g' % (tuple(bob.snap['Position'][-1])+
           tuple(cmr)+tuple(bob.snap['Position'][0]))+' %g' % (bob.snap['time']))

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()
        

