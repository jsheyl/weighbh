
output='''
# Region file format: DS9 version 4.1
panda(12:30:49.404,+12:23:28.11,30.2166,169.968,1,12",24",1) 
panda(12:30:49.404,+12:23:28.11,210.217,369.968,1,12",24",1) 
'''

for minval in range(1,70):
    maxval=minval+1
    with open('isoregion_%02d.reg' % minval,'w') as f:
        f.write((output.replace('12"',str(minval)+'"')).replace('24"',str(maxval)+'"'))
