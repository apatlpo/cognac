import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import gridspec
import matplotlib.colors as colors
import matplotlib.cm as cmx

from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from bokeh.palettes import Category10 as palette

sec2min = 1./60.

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

    def plot_bk(self, gridded=True, size=None):
        """ Bokeh (dynamic) plot of the log
        """
        if size:
            width, height = size
        elif gridded:
            width, height = (400, 200)
        else:
            width, height = (600, 400)

        output_notebook()
        TOOLS = 'pan,wheel_zoom,box_zoom,reset,help'
        colors = palette[10]

        _d = self._df
                
        if gridded:
            col = 'black'
            def bfig(c, x_range):
                s = figure(tools=TOOLS, 
                           plot_width=width, 
                           plot_height=height,
                           title=c, 
                           x_range=x_range)
                s.line(_d['time']/60, _d[c], legend=c, 
                       color=col, line_width=3)
                return s
            c0 = next(c for c in _d.columns if c is not 'time')
            s = [bfig(c0, None)]
            s = s+[bfig(c, s[0].x_range) for c in _d.columns 
                   if c not in ['time', c0]]
            grid = gridplot(s, ncols=2)
            show(grid)
        else:
            s1 = figure(tools=TOOLS,
                        plot_width=width,
                        plot_height=height,
                        title=None)
            iters = ((_d[c],c) for c in _d.columns if c is not 'time')
            for c, col in zip(iters, colors):
                s1.line(_d['time']/60, c[0], legend=c[1], 
                        color=col, line_width=3)
            show(s1)

#
def plot_logs(logs,
              f,
              z_target=None,
              eta=None,
              title=None,
              figsize=(12, 10)
              ):
    """ Standard plot of a deployment
    
    Parameters:
    -----------
    logs: logger, dict
        logs to plot, a dict of multiple logs may be passed
    f: autonomous_float
    z_target: lambda, optional
        trajectory prescribed
    eta: lambda, optional
        isopycnal displacements
    title: str, optional
    figsize: tuple
    """
    if isinstance(logs, dict) and 'state' not in logs:
        cols = get_cmap_colors(len(logs), cmap='plasma')
        def _plot(ax, log, variable, name, scale=1., offset=0.):
            for key, c in zip(logs, cols):
                t, _log = logs[key][log]['time'], logs[key][log]
                ax.plot(t *sec2min, _log[variable]*scale + offset, 
                        color=c, 
                        label=name+' '+key)
        key = list(logs)[0]
        log0 = logs[key]
    else:
        cols = 'k'
        def _plot(ax, log, variable, name, scale=1., offset=0.):
            t, _log = logs[log]['time'], logs[log]
            ax.plot(t *sec2min, _log[variable]*scale + offset, 
                    color=cols, 
                    label=name)
        log0 = logs
    t = log0['state']['time']
    #
    fig, axes = plt.subplots(nrows=3,ncols=2, sharex=True, figsize=figsize)
    #
    ax = axes[0,0]
    _plot(ax, 'state', 'z', 'z')
    if z_target is not None:
        ax.plot(t *sec2min, z_target(t), color='b', label='target')
        if eta is not None:
            ax.plot(t *sec2min, z_target(t) + eta(t), 
                    color='green', label='target+eta')
    ax.legend(loc=0)
    ax.set_ylabel('z [m]')
    if title is not None:
        ax.set_title(title)
    ax.grid()
    #
    ax = axes[0,1]
    if z_target is not None:
        if hasattr(f, 'ctrl'):
            # (x,y) # width # height
            ax.fill_between(t *sec2min, t * 0. - f.ctrl.dz_nochattering,
                            t * 0. + f.ctrl.dz_nochattering,
                            facecolor='orange', alpha=.5)
        #
        if eta is not None:
            ax.plot(t *sec2min, z_target(t) + eta(t), 
                    color='green', label='target+eta')
        _plot(ax, 'state', 'z', 'z', offset=-z_target(t))
        ax.legend()
        ax.set_ylabel('[m]')
        ax.set_ylim([-2., 2.])
        if title is not None:
            ax.set_title(title)
        ax.grid()
        ax.yaxis.set_label_position('right')
        ax.yaxis.tick_right()        
    #
    ax = axes[1,0]
    _plot(ax, 'state', 'w', 'dzdt', scale=1e2)
    ax.legend()
    ax.set_ylabel('[cm/s]')
    ax.grid()
    #
    ax = axes[1,1]
    _plot(ax, 'state', 'v', 'v', scale=1e6)    
    ax.axhline(f.piston.vol_min * 1.e6, ls='--', color='0.5')
    ax.axhline(f.piston.vol_max * 1.e6, ls='--', color='0.5')
    ax.legend()
    ax.set_ylabel('[cm^3]')
    ax.grid()
    ax.yaxis.set_label_position('right')
    ax.yaxis.tick_right()
    #
    ax = axes[2,0]
    _plot(ax, 'state', 'dwdt', 'd2zdt2')
    ax.legend()
    ax.set_xlabel('t [min]')
    ax.set_ylabel('[m/s^2]')
    ax.grid()
    #
    if 'piston' in log0:
        #
        ax = axes[2,1]
        _plot(ax, 'piston', 'work', 'piston work')
        ax.legend()
        ax.set_xlabel('t [min]')
        ax.set_ylabel('[Wh]')
        ax.grid()
        ax.yaxis.set_label_position('right')
        ax.yaxis.tick_right()
        #
        # extrapolate to a 30d long simulation
        piston = log0['piston']
        t = piston['time']
        nrg = (piston['work'].iloc[-1]-piston['work'].iloc[0]) \
                /(t.iloc[-1]-t.iloc[0])
        print( 'Extrapolated energy conssumption: %.1e Wh/day = %.1f Wh/30day' \
              %( nrg*86400, nrg*86400*30. ))

# ------------------------- kalman specific plots ------------------------------

def _plot_fill(k, t, i, ax, alpha):
    ax.fill_between(t, k['x_%d'%i] - np.sqrt(k['gamma_%d'%i]),
                        k['x_%d'%i] + np.sqrt(k['gamma_%d'%i]),
                    facecolor='orange', alpha=alpha)
    ax.plot(t,k['x_%d'%i], color='r', label ="estimation")

def plot_kalman(log, f, z_target=None, truth=None, figsize=(10,6)):
    alpha = 0.7
    #state, t = log['state'], log['state']['time']*sec2min
    k, tk = log['kalman'], log['kalman']['time']*sec2min
    if truth is not None:
        tt = truth.reset_index()['time']*sec2min
    #
    N = f.kalman.x_hat.size
    cols = 2
    rows = int(np.ceil(N / cols))
    gs = gridspec.GridSpec(rows, cols)
    fig = plt.figure(figsize=figsize)
    for i, name in enumerate(f.kalman.names):
        ax = fig.add_subplot(gs[i])
        _plot_fill(k, tk, i, ax, alpha)
        if truth is not None:
            ax.plot(tt, truth['x_%d'%i], color='k', label = "truth")
        if i==1 and z_target:
            ax.plot(tk, -z_target(log['kalman']['time']),
                    color='b', label='target')
        ax.set_title(name)
        ax.grid()
        legend = ax.legend(loc='best', shadow=True, fontsize='medium')
    fig.tight_layout()  

def plot_kalman_v0(log, f, V_e=None, gamma_e=None, z_target=None):
    alpha = 0.7
    state, t = log['state'], log['state']['time']*sec2min
    k, tk = log['kalman'], log['kalman']['time']*sec2min
    #
    fig = plt.figure(figsize=(15,10))
    #
    ax=fig.add_subplot(231)
    #ax.plot(tk, - k['gamma_1'])
    ax.fill_between(tk, k['z'] - np.sqrt(k['gamma_0']),
                        k['z'] + np.sqrt(k['gamma_0']),
                    facecolor='orange', alpha=alpha)
    ax.plot(tk,k['z'], color='r', label ="estimated depth")
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
    ax.fill_between(tk, k['w'] - np.sqrt(k['gamma_1']),
                        k['w'] + np.sqrt(k['gamma_1']),
                    facecolor='orange', alpha=alpha)
    ax.plot(tk,k['w'], color='r', label ="estimated velocity")
    ax.plot(t,state['w'], color='k', label = "real velocity")
    ax.set_title("velocity [m/s]")
    #ax.set_xlabel("t (min)")
    ax.grid()
    legend = ax.legend(loc='best', shadow=True, fontsize='medium')
    #
    ax=fig.add_subplot(234)
    ax.fill_between(tk, (k['V_e'] - np.sqrt(k['gamma_3']))*1e6,
                        (k['V_e'] + np.sqrt(k['gamma_3']))*1e6,
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
    ax.fill_between(tk, k['gamma_e'] - np.sqrt(k['gamma_2']),
                        k['gamma_e'] + np.sqrt(k['gamma_2']),
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

def get_cmap_colors(Nc, cmap='plasma'):
    """ load colors from a colormap to plot lines
    
    Parameters
    ----------
    Nc: int
        Number of colors to select
    cmap: str, optional
        Colormap to pick color from (default: 'plasma')
    """
    scalarMap = cmx.ScalarMappable(norm=colors.Normalize(vmin=0, vmax=Nc),
                                   cmap=cmap)
    return [scalarMap.to_rgba(i) for i in range(Nc)]
