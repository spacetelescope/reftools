#make the libaray file and image file user input
from __future__ import print_function

dqheader='/user/sosey/wfc3/wf3dq.h'
#image=

f=open(dqheader)
lines=f.readlines()
f.close()

from collections import OrderedDict
dqval=OrderedDict()

#parse the header file
for line in lines:
    if '/*' in line[:2]:
        print(line)

    if '#' in line[:2]:
        cols=line.split()
        start=line.find("/*")
        end=line.find("*/")
        dqval[int(cols[3])] = line[start+2:end]

revkeys=dqval.keys()
revkeys.reverse()

#check out the user image, ask for ver from user
image='/user/lim/ssb/calacs/unittest/pcte_2k/orig/jbk102qpq_blv_tmp.fits'
ext='DQ'
ver=1

from astropy.io import fits
import numpy as np

data=fits.getdata(image,extname=ext,extver=ver)
pixel=np.where(data > 0)
pixlist=np.dstack(pixel)

for loc in pixlist[0][:3]:
    dval=data[loc[0],loc[1]]
    leftover=dval
    print("(%i,%i): %i\n"%(loc[0],loc[1],dval))
    while leftover > 0:
        for dq in revkeys:
            delta=leftover-dq
            if  (delta == 0):
                print("\t%s\n"%(dqval[dq]))
                leftover =0
            elif (delta > 0):
                print("\t%s\n"%(dqval[dq]))
                leftover -= dq





