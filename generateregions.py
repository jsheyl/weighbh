
output='''
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
panda(12:30:49.404,+12:23:28.11,30.2166,90,1,12",24",1) || # panda=(30.2166 90 180 270 189.968)(12" 24")
panda(12:30:49.404,+12:23:28.11,90,180,1,12",24",1) || # panda=ignore
panda(12:30:49.404,+12:23:28.11,180,270,1,12",24",1) || # panda=ignore
panda(12:30:49.404,+12:23:28.11,270,549.968,1,12",24",1) || # panda=ignore
panda(12:30:49.404,+12:23:28.11,210.217,450,1,12",24",1) # panda=(210.217 90 180 270 9.968)(12" 24")
panda(12:30:49.404,+12:23:28.11,90,180,1,12",24",1) # panda=ignore
panda(12:30:49.404,+12:23:28.11,180,270,1,12",24",1) # panda=ignore
panda(12:30:49.404,+12:23:28.11,270,369.968,1,12",24",1) # panda=ignore
'''

for minval in range(36,80):
    maxval=minval+1
    with open('isoregion_%02d.reg' % minval,'w') as f:
        f.write((output.replace('12"',str(minval)+'"')).replace('24"',str(maxval)+'"'))
