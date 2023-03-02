import h5py
import numpy as np
import sys

fn = sys.argv[1]
f = h5py.File(fn)
x = np.asarray(f['aux']['fld'], dtype='float64')
y = np.cumsum(x)/np.sum(x)
cutoff = np.argmax(y)
print(cutoff)
# max insert size. 

# cutoff = np.argmax(y > .95)
# 95% CI insert size. If needed. 

## Mean insert size. If needed. 
# fn = sys.argv[1]
# f = h5py.File(fn)
# x = np.asarray(f['aux']['fld'], dtype='float64')
# t = np.arange(0,len(x),1)
# mean = np.sum(x*t)/np.sum(x)
# print(mean)
