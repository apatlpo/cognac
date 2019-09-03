
import sys, os
from shutil import copyfile
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import pickle

class wrapper():
    # wrappers around C library

    def __init__(self):
        self._lines = ctypes.cdll.LoadLibrary('../compileSO/Lines.so')
        # nm -gU ../compileSO/Lines.so | grep Line
        # https://pgi-jcns.fz-juelich.de/portal/pages/using-c-from-python.html
        # https://docs.python.org/3/library/ctypes.html

        #print(dir(_lines.LinesInit))
        self._lines.LinesInit.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                          ctypes.c_double, ctypes.c_double)

        self._lines.LinesCalc.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                          ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.c_double)

    def __del__(self):
        self._lines.LinesClose()
        print('MoorDyn closed')

    def LinesInit(self, x, xd, U0, U1):
        #global _lines
        nx=len(x)
        nxd=len(xd)
        x_array_type = ctypes.c_double * nx
        xd_array_type = ctypes.c_double * nxd
        self._lines.LinesInit(x_array_type(), xd_array_type(), ctypes.c_double(U0), ctypes.c_double(U1))

    def LinesCalc(self, x, xd, Flines, t, dt):
        #global _lines
        nx=len(x)
        nxd=len(xd)
        nFlines=len(Flines)
        x_array_type = ctypes.c_double * nx
        xd_array_type = ctypes.c_double * nxd
        Flines_array_type = ctypes.c_double * nFlines
        self._lines.LinesCalc(x_array_type(), xd_array_type(), Flines_array_type(), ctypes.c_double(t), ctypes.c_double(dt))


# utils
def read_output():
    fname='Mooring/Line1.out'
    print('Read: '+fname)
    f = open(fname, 'r')
    #
    header = f.readline().strip().split()
    for line in f:
        d = line.strip().split()
        #print(' '.join(d))
        d = np.array([float(s) for s in d])
        t = d[0]
        x = d[1::3]
        y = d[2::3]
        z = d[3::3]
    f.close()
    return x, y, z, t

def plot_output():
    # look at output
    print('Plot a solution')
    x, y, z, t = read_output()
    for v in [x, y, z]:
        print(v)

    plt.figure()
    plt.plot(x, z)
    plt.grid()
    plt.show()
    plt.ion()


#
if __name__ == '__main__':

    #
    # X = [0., 0., 0., 0., 0., 0.] # plateform position
    # Xd = [0., 0., 0., 0., 0., 0.] # plateform velocity
    X = [0., 0., 0.]  # plateform position
    Xd = [0., 0., 0.]  # plateform velocity
    #
    if len(sys.argv)==1 or sys.argv[1]=='compute':
        # compute solution
        print('Compute a solution')
        #
        if len(sys.argv)>=3:
            u0 = float(sys.argv[2])
        if len(sys.argv)>=4:
            u1 = float(sys.argv[3])
        #
        w = wrapper()
        w.LinesInit(X, Xd, u0, u1)
        del w
        #
        plot_output()

    elif sys.argv[1]=='all'
        pref='run_U_'
        # store lines.txt
        copyfile('Mooring/lines.txt','Mooring/'+pref+'lines.txt')
        #U0 = [0. , 0.,  ,]
        U1 = [.05, .1, .2, .5, 1.]
        u1=.1
        for i, u1 in enumerate(U1):
            w = wrapper()
            w.LinesInit(X, Xd, u1, u1)
            del w
            x, y, z, t = read_output()
            f = open('Mooring/'+pref+'%.03d.p'%i, 'wb')
            pickle.dump([u1,u1,x,y,z],f)
            f.close()



    elif sys.argv[1]=='plot':
        plot_output()

    elif sys.argv[1]=='clean':
        print('Remove output files')
        os.remove('Mooring/Line1.out')
        os.remove('Mooring/Lines.out')

    else:
        print('command line argument should be absent or \'compute\', \'all\' or  \'plot\' ')









