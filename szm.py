import numpy as np
from math import *

from ctypes import *
import os
# compilation
# g++ -fPIC -shared -o libszm.so szm.c
lizm = CDLL(os.path.abspath('libszm.so'))

lizm.szm.argtypes = (POINTER(c_ushort),c_int,c_int,
                        POINTER(c_float),POINTER(c_float),
                                     POINTER(c_int))

lizm.szm.restype = None

def szm(img):
    h,w = np.shape(img)
    ns = c_int()
    x = np.zeros(np.size(img), dtype='float32')
    y = np.zeros(np.size(img), dtype='float32')
    lizm.szm(img.ctypes.data_as(POINTER(c_ushort)),
               c_int(w),c_int(h),
               x.ctypes.data_as(POINTER(c_float)),
               y.ctypes.data_as(POINTER(c_float)),
               byref(ns))
    return x[0:ns.value],y[0:ns.value]
    # return np.ctypeslib.as_array(x, shape=(ns.value,)).astype('float64')
    
lizm.szm2.argtypes = (POINTER(c_ushort),c_int,c_int,
                        POINTER(c_float),POINTER(c_float),
                        POINTER(c_float),POINTER(c_float),
                                     POINTER(c_int))

lizm.szm2.restype = None
    
def szm2(img):
    h,w = np.shape(img)
    ns = c_int()
    x = np.zeros(np.size(img), dtype='float32')
    y = np.zeros(np.size(img), dtype='float32')
    l = np.zeros(np.size(img), dtype='float32')
    pa = np.zeros(np.size(img), dtype='float32')
    lizm.szm2(img.ctypes.data_as(POINTER(c_ushort)),
               c_int(w),c_int(h),
               x.ctypes.data_as(POINTER(c_float)),
               y.ctypes.data_as(POINTER(c_float)),
               l.ctypes.data_as(POINTER(c_float)),
               pa.ctypes.data_as(POINTER(c_float)),
               byref(ns))
    return x[0:ns.value],y[0:ns.value],l[0:ns.value],pa[0:ns.value]
    