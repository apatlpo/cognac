import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from bokeh.palettes import Category10 as palette

s2m = 1./60.

class logger():
    ''' Store a log of the float trajectory

    Parameters
    ----------
    var: list of strings
        List containing the name of variables that will be logged

    '''

    def __init__(self, logs):
        #self.var = var
        #for item in var:
        #    setattr(self, item, np.array([]))
        #
        self._df = pd.DataFrame(columns=['time']+logs)

    def __getitem__(self,key):
        ''' return column of the log
        '''
        if key is not 'time':
            return self._df.set_index('time')[key]
        else:
            return self._df['time']

    def __repr__(self):
        return str(self._df.head())

    def store(self, **kwargs):
        ''' Appends variables to the logger database:
        Usage:
            log.store(time=10., v=1.)
        The above line will append values 10. and 1. to variables t and v respectively
        '''
        #for item in self.var:
        #    if item in kwargs:
        #        setattr(self, item, np.hstack((getattr(self, item), kwargs[item])))
        self._df = self._df.append(kwargs, ignore_index=True)
        #self._df.update(kwargs)

    def cleanup(self):
        ''' Merge duplicated temporal information
        '''
        self._df = self._df.groupby('time').mean().reset_index()

    def plot(self, **kwargs):
        _dkwargs = {'subplots': True, 'grid': True}
        _dkwargs.update(kwargs)
        ax = self._df.set_index('time').plot(**_dkwargs)
        return ax

    def plot_bk(self):

        _d = self._df

        output_notebook()
        TOOLS = 'pan,wheel_zoom,box_zoom,reset,help'
        colors = palette[10]
        
        s1 = figure(tools=TOOLS, plot_width=600, plot_height=300, title=None)
        iters = ((_d[c],c) for c in _d.columns if c is not 'time')
        for c, col in zip(iters, colors):
            s1.line(_d['time']/60, c[0], legend=c[1], color=col, line_width=3)
        show(s1)

#
def plot_logs(log, f, z_target=None, eta=None, title=None):
    #
    plt.figure(figsize=(12, 10))
    #
    state = log['state']
    t = state['time']
    ax = plt.subplot(321)
    ax.plot(t *s2m, state['z'], color='k', label='z')
    if z_target is not None:
        ax.plot(t *s2m, z_target(t), color='b', label='target')
        if eta is not None:
            ax.plot(t *s2m, z_target(t) + eta(t), color='green', label='target+eta')
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
            ax.fill_between(t *s2m, t * 0. - f.ctrl.dz_nochattering,
                            t * 0. + f.ctrl.dz_nochattering,
                            facecolor='orange', alpha=.5)
        ax.plot(t *s2m, state['z'].values - z_target(t), \
                color='k', label='z-ztarget')
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
    ax.plot(t *s2m, state['w'] * 1.e2, color='k', label='dzdt')
    ax.legend()
    ax.set_ylabel('[cm/s]')
    ax.grid()
    #
    ax = plt.subplot(324)
    ax.plot(t *s2m, state['v'] * 1.e6, '-', color='k', label='v')
    # ax.axhline(f.piston.dv*1.e6,ls='--',color='k')
    ax.axhline(f.piston.vol_min * 1.e6, ls='--', color='0.5')
    ax.axhline(f.piston.vol_max * 1.e6, ls='--', color='0.5')
    ax.legend()
    ax.set_ylabel('[cm^3]')
    ax.grid()
    ax.yaxis.set_label_position('right')
    ax.yaxis.tick_right()
    #
    ax = plt.subplot(325, sharex=ax)
    ax.plot(t *s2m, state['dwdt'], color='k', label='d2zdt2')
    ax.legend()
    ax.set_xlabel('t [min]')
    ax.set_ylabel('[m/s^2]')
    ax.grid()
    #
    if 'piston' in log:
        piston = log['piston']
        t = piston['time']
        #
        ax = plt.subplot(326, sharex=ax)
        ax.plot(t *s2m, piston['work'], color='k', label='piston work')
        ax.legend()
        ax.set_xlabel('t [min]')
        ax.set_ylabel('[Wh]')
        ax.grid()
        ax.yaxis.set_label_position('right')
        ax.yaxis.tick_right()
        #
        # extrapolate to a 30d long simulation
        nrg = (piston['work'].iloc[-1]-piston['work'].iloc[0]) \
                /(t.iloc[-1]-t.iloc[0])
        print( 'Extrapolated energy conssumption: %.1e Wh/day = %.1f Wh/30day' \
              %( nrg*86400, nrg*86400*30. ))

def plot_kalman(log, f, V_e=None, gamma_e=None, z_target=None):
    alpha = 0.7
    state, t = log['state'], log['state']['time']*s2m
    k, tk = log['kalman'], log['kalman']['time']*s2m
    #
    fig = plt.figure(figsize=(15,10))
    #
    ax=fig.add_subplot(231)
    #ax.plot(tk, - k['gamma_diag1'])
    ax.fill_between(tk, -k['z'] - np.sqrt(k['gamma_diag0']),
                        -k['z'] + np.sqrt(k['gamma_diag0']),
                    facecolor='orange', alpha=alpha)
    ax.plot(tk,-k['z'], color='r', label ="estimated depth")
    ax.plot(t, state['z'], color='k', label = "real depth")
    if z_target is not None:
        ax.plot(t, z_target(log['state']['time']),
                color='b', label='target')
    ax.set_title("z  [m]")
    #ax.set_xlabel("t (min)")
    ax.grid()
    legend = ax.legend(loc='best', shadow=True, fontsize='medium')
    #
    ax=fig.add_subplot(232)
    ax.fill_between(tk, -k['w'] - np.sqrt(k['gamma_diag1']),
                        -k['w'] + np.sqrt(k['gamma_diag1']),
                    facecolor='orange', alpha=alpha)
    ax.plot(tk,-k['w'], color='r', label ="estimated velocity")
    ax.plot(t,state['w'], color='k', label = "real velocity")
    ax.set_title("velocity [m/s]")
    #ax.set_xlabel("t (min)")
    ax.grid()
    legend = ax.legend(loc='best', shadow=True, fontsize='medium')
    #
    ax=fig.add_subplot(234)
    ax.fill_between(tk, (k['V_e'] - np.sqrt(k['gamma_diag3']))*1e6,
                        (k['V_e'] + np.sqrt(k['gamma_diag3']))*1e6,
                    facecolor='orange', alpha=alpha)
    ax.plot(tk,k['V_e']*1e6, color='r', label ="estimated V_e volume")
    if V_e is not None:
        if isinstance(V_e,float) or V_e.size==1:
            ax.axhline(V_e*1e6, color='k')
        else:
            ax.plot(tk,V_e*1e6, color='k', label = "real V_e")
    ax.set_title("volume V_e [cm^3]")
    ax.set_xlabel("t (min)")
    ax.grid()
    legend = ax.legend(loc='best', shadow=True, fontsize='medium')
    #
    ax=fig.add_subplot(235)
    ax.fill_between(tk, k['gamma_e'] - np.sqrt(k['gamma_diag2']),
                        k['gamma_e'] + np.sqrt(k['gamma_diag2']),
                    facecolor='orange', alpha=alpha)
    ax.plot(tk,k['gamma_e'], color='r',
            label ="estimated equivalent compressibility")
    if gamma_e is not None:
        if isinstance(gamma_e,float) or gamma_e.size==1:
            ax.axhline(gamma_e,color='k')
        else:
            ax.plot(tk,gamma_e, color='k',
                label = "gamma_e")
    else:
        ax.plot(tk,k['gammaV'], color='k',
            label = "float compressibility x float volume")
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.set_title(r"equivalent compressibility $\gamma_e$ [m^2]")
    ax.set_xlabel("t (min)")
    ax.grid()
    legend = ax.legend(loc='best', shadow=True, fontsize='medium')
    #
    ax=fig.add_subplot(233)
    ax.plot(tk,k['dwdt'], color='r', label ="estimated acceleration")
    ax.plot(t,state['dwdt'], color='k', label = "real acceleration")
    ax.set_title("acceleration dw/dt [m.s^-2]")
    ax.set_xlabel("t (min)")
    ax.grid()
    legend = ax.legend(loc='best', shadow=True, fontsize='medium')
    #
    ax=fig.add_subplot(236)
    ax.plot(tk,k['dwdt_diff'], color='tab:blue', label ="acceleration/force difference")
    ax.set_title("acceleration/force difference [m.s^-2]")
    ax.set_xlabel("t (min)")
    ax.grid()
    #legend = ax.legend(loc='best', shadow=True, fontsize='medium')
