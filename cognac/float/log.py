import numpy as np
import matplotlib.pyplot as plt

class logger():
    ''' Store a log of the float trajectory

    Parameters
    ----------
    var: list of strings
        List containing the name of variables that will be logged

    '''

    def __init__(self, var):
        self.var = var
        for item in var:
            setattr(self, item, np.array([]))

    def store(self, **kwargs):
        ''' Appends variables to the logger database:

        Usage:

        log.store(t=10., v=1.)

        The above line will append values 10. and 1. to variables t and v respectively

        '''
        for item in self.var:
            if item in kwargs:
                setattr(self, item, np.hstack((getattr(self, item), kwargs[item])))

#
def plot_log(log, f, z_target=None, eta=None, title=None):
    t = log.t
    #
    if hasattr(log,'nrg'):
        # extrapolate to a 30d long simulation
        nrg = (log.nrg[-1]-log.nrg[0])/(t[-1]-t[0])
        print( 'Extrapolated energy conssumption: %.1f Wh/day = %.1f Wh/30day'
              %( nrg*86400, nrg*86400*30. ))
    #
    plt.figure(figsize=(12, 10))
    #
    ax = plt.subplot(321)
    ax.plot(t / 60., log.z, label='z')
    if z_target is not None:
        ax.plot(t / 60., z_target(t), color='r', label='target')
        if eta is not None:
            ax.plot(t / 60., z_target(t) + eta(t), color='green', label='target+eta')
    ax.legend(loc=0)
    ax.set_ylabel('z [m]')
    if title is not None:
        ax.set_title(title)
    ax.grid()
    #
    ax = plt.subplot(322)
    if z_target is not None:
        if hasattr(f, 'ctrl'):
            # (x,y) # width # height
            ax.fill_between(t / 60., t * 0. - f.ctrl['dz_nochattering'], t * 0. + f.ctrl['dz_nochattering'],
                            facecolor='orange', alpha=.5)
        ax.plot(t / 60., log.z - z_target(t), label='z-ztarget')
        ax.legend()
        ax.set_ylabel('[m]')
        ax.set_ylim([-2., 2.])
        if title is not None:
            ax.set_title(title)
        ax.grid()
        ax.yaxis.set_label_position('right')
        ax.yaxis.tick_right()
    #
    ax = plt.subplot(323, sharex=ax)
    ax.plot(t / 60., log.w * 1.e2, label='dzdt')
    ax.legend()
    ax.set_ylabel('[cm/s]')
    ax.grid()
    #
    ax = plt.subplot(324)
    ax.plot(t / 60., log.v * 1.e6, '-', label='v')
    # ax.axhline(f.piston.dv*1.e6,ls='--',color='k')
    ax.axhline(f.piston.vol_min * 1.e6, ls='--', color='k')
    ax.axhline(f.piston.vol_max * 1.e6, ls='--', color='k')
    ax.legend()
    ax.set_ylabel('[cm^3]')
    ax.grid()
    ax.yaxis.set_label_position('right')
    ax.yaxis.tick_right()
    #
    ax = plt.subplot(325, sharex=ax)
    ax.plot(t / 60., log.dwdt, label='d2zdt2')
    ax.legend()
    ax.set_xlabel('t [min]')
    ax.set_ylabel('[m/s^2]')
    ax.grid()
    #
    if hasattr(log,'nrg'):
        ax = plt.subplot(326, sharex=ax)
        ax.plot(t / 60., log.nrg, label='nrg')
        ax.legend()
        ax.set_xlabel('t [min]')
        ax.set_ylabel('[Wh]')
        ax.grid()
        ax.yaxis.set_label_position('right')
        ax.yaxis.tick_right()
