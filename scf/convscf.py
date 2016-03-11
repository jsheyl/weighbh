import gadget
import sys

bob=gadget.Snapshot(sys.argv[1])
pos=bob.snap['Position']
vel=bob.snap['Velocity']
mass=bob.snap['mass']
npart=bob.snap['npart']

print('%20d %14.6e' % (sum(npart),bob.snap['time']))
startpos=0
m=[0]
for i in range(0,5):
    endpos=startpos+npart[i]
    m[0]=mass[i]
    for n in range(startpos,endpos):
        print('%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e' %
              (tuple(m)+tuple(pos[n])+tuple(vel[n])))
