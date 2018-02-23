
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from netCDF4 import Dataset
import gsw

import matplotlib.pyplot as plt
import cartopy.crs as ccrs


# useful parameters
g=9.81
watth = 1e-6/3.6

#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


#
class autonomous_float():
    
    def __init__(self,**kwargs):
        ''' Float constructor, the float is assumed to be a cylinder
        
        Parameters
        ----------
            a : float radius [m]
            L : float length [m]
            m : float mass [kg]
            gamma : mechanical compressibility [1/dbar]
            alpha : thermal compressibility [1/degC]
            temp0: reference temperature used for thermal compressibility [degC]
            
        '''
        # default parameters
        params = {'a': 0.05, 'L': 0.4, 'gamma': 2.e-6, 'alpha': 7.e-5, 'temp0': 15.}
        params['m']= 1000. * np.pi * params['a']**2 * params['L']
        #
        if 'model' in kwargs:
            if kwargs['model'] is 'ENSTA':
                params = {'a': 0.06, 'L': 0.5, 'gamma': None, 'alpha': 0., 'temp0': 0.}
                params['m'] = 1000. * np.pi * params['a'] ** 2 * params['L']
        #
        params.update(kwargs)
        for key,val in params.items():
            setattr(self,key,val)
        # compute the volume of the float:
        self.V = np.pi*self.a**2*self.L

    def __repr__(self):
        strout='Float parameters: \n'
        strout+='  L     = %.2f m      - float length\n'%(self.L)
        strout+='  a     = %.2f m      - float radius\n'%(self.a)
        strout+='  m     = %.2f kg     - float radius\n'%(self.m)
        strout+='  V     = %.2e cm^3   - float volume\n'%(self.V*1.e6)
        strout+='  gamma = %.2e /dbar  - mechanical compressibility\n'%(self.gamma)
        strout+='  alpha = %.2e /degC  - thermal compressibility\n'%(self.alpha)
        strout+='  temp0 = %.2e  degC  - reference temperature\n'%(self.temp0)
        if hasattr(self,'piston'):
            strout+=str(self.piston)
        return strout
        
    def rho(self,p=None,temp=None,v=None,z=None,waterp=None):
        ''' Returns float density i.e. mass over volume
        '''
        if v is None:
            if hasattr(self,'v'):
                v = self.v
            else:
                v = 0.
        if p is not None and temp is not None:
            return self.m/(self.V*(1.-self.gamma*p+self.alpha*(temp-self.temp0))+v)
        elif z is not None and waterp is not None:
            # assumes thermal equilibrium
            p, tempw = waterp.get_p(self.z), waterp.get_temp(self.z)
            return self.rho(p=p,temp=tempw,v=v)
        else:
            print('You need to provide p/temp or z/waterp')

    def volume(self,**kwargs):
        ''' Returns float volume (V+v)
        '''
        return self.m/self.rho(**kwargs)    
    
    def volume4equilibrium(self,p_eq,temp_eq,rho_eq):
        ''' Find volume that needs to be added in order to be at equilibrium 
            prescribed pressure, temperature, density
        '''
        v = fsolve(lambda x: rho_eq-self.rho(p=p_eq,temp=temp_eq,v=x),0.)[0]
        return v

    def adjust_m(self,p_eq,temp_eq,rho_eq):
        ''' Find mass that needs to be added in order to be at equilibrium 
            prescribed pressure, temperature, density
        '''
        def func(m):
            self.m = m
            return rho_eq-self.rho(p=p_eq,temp=temp_eq)
        m0=self.m
        self.m = fsolve(func,self.m)[0]
        print('%.1f g were added to the float in order to be at equilibrium at %.0f dbar \n'%((self.m-m0)*1.e3,p_eq))

    def init_piston(self,**kwargs):
        self.piston = piston(**kwargs)
        
    def _f0(self,p,rhof,rhow,Lv):
        ''' Compute the vertical force exterted on the float
        '''
        f = -self.m*g
        f += self.m*rhow/rhof*g # Pi0, we ignore DwDt terms for now
        f += -self.m/Lv*np.abs(self.w)*self.w # Pi'
        return f

    def _f(self,z,waterp,Lv,v=None,w=None):
        ''' Compute the vertical force exterted on the float
        '''
        p, tempw = waterp.get_p(z), waterp.get_temp(z)
        rhow = waterp.get_rho(z)
        rhof = self.rho(p=p,temp=tempw,v=v)
        #
        f = -self.m*g
        f += self.m*rhow/rhof*g # Pi0, we ignore DwDt terms for now
        #
        if w is None:
            w=self.w
        f += -self.m/Lv*np.abs(w)*w # Pi'
        return f
    
    def _df(self,z,waterp,Lv):
        ''' Compute gradients of the vertical force exterted on the float
        '''
        df1 = ( self._f(z+5.e-2,waterp,Lv) - self._f(z-5.e-2,waterp,Lv) ) /1.e-1
        df2 = ( self._f(z,waterp,Lv,w=self.w+5.e-3) - self._f(z,waterp,Lv,w=self.w-5.e-3) ) /1.e-2
        df3 = ( self._f(z,waterp,Lv,v=self.v+5.e-5) - self._f(z,waterp,Lv,v=self.v-5.e-5) ) /1.e-4
        return df1, df2, df3
        
    def compute_bounds(self,waterp,zmin,zmax=0.,Lv=None):
        ''' Compute approximate bounds on velocity and acceleration
        '''
        #self.w=0.
        if Lv is None:
            if hasattr(self,'Lv'):
                Lv=self.Lv
            else:
                Lv=self.L
        #
        z=zmax
        if hasattr(self,'piston'):
            v=self.piston.vol_max
        else:
            v=None
        fmax=self._f(z,waterp,Lv,v=v,w=0.) # upward force
        #
        z=zmin
        if hasattr(self,'piston'):
            v=self.piston.vol_min
        else:
            v=None
        fmin=self._f(z,waterp,Lv,v=v,w=0.) # downward force
        #
        afmax = np.amax((np.abs(fmax),np.abs(fmin)))
        wmax = np.sqrt( afmax * Lv/self.m)
        print('Acceleration and velocity bounds (zmin=%.0fm,zmax=%.0fm):' %(zmin,zmax))
        print('fmax/m=%.1e m^2/s, fmin/m= %.1e m^2/s, wmax= %.1f cm/s' %(fmax/self.m,fmin/self.m,wmax*100.) )
        print('For accelerations, equivalent speed reached in 1min:')
        print('  fmax %.1e cm/s, fmin/m= %.1e cm/s' %(fmax/self.m*60.*100.,fmin/self.m*60.*100.) )
        #
        #if hasattr(self,'piston'):
        #    dv = self.piston.dv
        #    p, tempw = waterp.get_p(z), waterp.get_temp(z)
        #    rhow = self.rho(p=p,temp=tempw,v=self.piston.vol_min)
        #    rhof = self.rho(p=p,temp=tempw,v=self.piston.vol_min+dv)
        #    df = self.m*(-1.+rhow/rhof)*g
        #    print('Acceleration after an elementary piston displacement: %.1e m^2/s' %(df[0]/self.m))
        #    print('  corresponding speed and displacement after 1 min: %.1e m/s, %.1e m \n' \
        #          %(df[0]/self.m*60,df[0]/self.m*60**2/2.))
        
        return fmax, fmin, afmax, wmax        
        
    def time_step(self, waterp, T=600., dt_step=1.,
                  z=None, w=None, v=None, t0=0., Lv=None,
                  usepiston=False, z_target=None, 
                  ctrl=None,
                  eta=lambda t: 0.,
                  log=['t','z','w','v','dwdt'], dt_store=60., log_nrg=True, c_friction=2., p_float=1.e5, 
                  **kwargs):
        ''' Time step the float position given initial conditions
        
        Parameters
        ----------
            
        waterp: water profile object
                Contains information about the water profile
        T: float
            Length of the simulation in seconds
        dt_step: float
            Simulation time step
        z: float
            Initial position, 0. at the surface and negative downward
        w: float
            Initial vertical velocity, negative for downward motions
        v: float
            Initial volume adjustement (total volume is V+v)
        t0: float
            Initial time
        Lv: float
            Drag length scale
        usepiston: boolean, default is False
            Turns piston usage
        z_target: function
            Target trajectory as a function of time
        w_target: function
            Target velocity as a function of time
        tau_ctrl: float
            Control time scale
        d3y_ctrl: float
            Sliding control parameter
        dz_nochattering: float
            Width of the window around z_target where control is turned off
        eta: function
            Isopycnal displacement as a function of time
        log: list of strings or False
            List of variables that will logged
        dt_store: float
            Time interval between log storage
        log_nrg: boolean, default is True
            Turns on/off nrg computation and storage
        p_float: float
            Internal float pressure in Pa
        '''
        t=t0
        #
        if z is None:
            if not hasattr(self,'z'):
                z=0.
            else:
                z=self.z
        self.z = z
        #
        if w is None:
            if not hasattr(self,'w'):
                w=0.
            else:
                w=self.w
        self.w = w
        #
        if usepiston:
            if v is not None:
                self.piston.update_vol(v)
            self.v=self.piston.vol
            u = 0.
            if ctrl:
                ctrl_default = {'tau': 60., 'dz_nochattering': 1., 'mode': 'sliding'}
                ctrl_default.update(ctrl)
                self.ctrl = ctrl_default
                print(ctrl)
        elif v is None:
            if not hasattr(self,'v'):
                self.v=0.
        else:
            self.v=v
        v0 = self.v            
        #
        if Lv is None:
            Lv = self.L
        self.Lv = Lv
        #
        if log_nrg:
            if 'nrg' not in log:
                log.append('nrg')
            self.nrg = 0. # Wh
        if log:
            if hasattr(self,'log'):
                delattr(self,'log')
            self.log = logger(log)
        #
        print('Start time stepping for %d min ...'%(T/60.))
        #
        _f=0.
        while t<t0+T:
            #
            # get vertical force on float
            waterp.eta = eta(t) # update isopycnal displacement
            _f = self._f(self.z, waterp, self.Lv)
            # control starts here
            if usepiston and self.ctrl:
                # activate control only if difference between the target and actual vertical 
                # position is more than the dz_nochattering threshold                
                if np.abs(self.z-z_target(t)) > self.ctrl['dz_nochattering']:
                    u = self._control(t, waterp, z_target, _f)
                    #
                    v0 = self.piston.vol
                    self.piston.update(dt_step, u)
                    self.v = self.piston.vol
                # energy integration, 1e4 converts from dbar to Pa
                if log_nrg and (self.v != v0):
                    self.nrg += dt_step * c_friction * np.abs((waterp.get_p(self.z)*1.e4 - p_float)*u) * watth
            # store
            if log:
                if (dt_store is not None) and t_modulo_dt(t, dt_store, dt_step):
                    self.log.store(t=t, z=self.z, w=self.w, v=self.v, dwdt=_f/self.m)
                    if log_nrg:
                        self.log.store(nrg=self.nrg)
            # update variables
            self.z += dt_step*self.w
            self.z = np.amin((self.z,0.))
            self.w += dt_step*_f/self.m
            t+=dt_step
        print('... time stepping done')

    def _control(self, t, waterp, z_target, _f):
        z_t = z_target(t)        
        dz_t = (z_target(t+.05)-z_target(t-.05))/.1
        d2z_t = (z_target(t+.05)-2.*z_target(t)+z_target(t-.05))/.05**2
        #
        #x1=self.z
        x2=self.w
        #x3=self.V+self.v
        #f1=x2
        f2=_f/self.m
        f3=( self.volume(z=self.z+.5, waterp=waterp) - self.volume(z=self.z-.5, waterp=waterp) )/1. *x2 # dVdz*w
        df1, df2, df3 = self._df(self.z, waterp, self.Lv)
        df1, df2, df3 = df1/self.m, df2/self.m, df3/self.m
        #
        if self.ctrl['mode'] is 'sliding':
            d3y = self.ctrl['d3y_ctrl']*self._control_sliding(self.z, self.w, _f/self.m, z_t, dz_t, d2z_t, self.ctrl['tau'])
            u = df1*x2 + df2*f2 + df3*f3 - d3y
            u = -u/df3
        else:
            print('!! mode '+self.ctrl['mode']+' is not implemented, exiting ...')
            sys.exit()
        return u

    def _control_sliding(self,  z, dz, d2z, z_t, dz_t, d2z_t, tau_ctrl):
        ''' Several inputs are required:
        (z_target,w_target,dwdt_target) - describes the trajectory
        tau_ctrl  - a time scale of control 
        '''
        return np.sign( d2z_t - d2z + 2.*(dz_t-dz)/tau_ctrl + (z_t-z)/tau_ctrl**2 )

    
#
def compute_gamma(R,t,material=None,E=None,mu=.35):
    ''' Compute the compressibility to pressure of a cylinder
    
    Parameters
    ----------
    R: float, [m]
        cylinder radius
    t: float, [m]
        cylinder thickness
    E: float, [GPa]
        Young's modulus
    mu: float, []
        Poisson ratio
    
    Returns
    -------
    gamma: float, [1/dbar]
        approximate compressibility of the float
        
    '''
    pmat = {'glass': {'E': 90., 'mu': .25}, 'aluminium': {'E': 70., 'mu': .35},
            'pom': {'E': 3.5, 'mu': .35}, 'polycarbonat':  {'E': 2.2, 'mu': .37}}
    if material is not None:
        if material in pmat:
            E = pmat[material]['E']
            mu = pmat[material]['mu']
        else:
            print('material not in our database')
            sys.exit()
    elif E is None:
        print('You need to provide material or E !')
        sys.exit()
    # convert E to dbar
    E=E*1.e5
    return R*(6.-7.*mu)/2./E/t
    
    
# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------


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

# utils
def plot_float_density(z, f, waterp):
    ''' Plot float density with respect to that of water profile
    
    Parameters
    ----------
    z: ndarray
        depth grid
    f: float object
    waterp: water profile object
    
    '''
    #
    rho_w, p, temp = waterp.get_rho(z), waterp.get_p(z), waterp.get_temp(z)
    #
    plt.figure()
    ax = plt.subplot(111)
    #
    iz = np.argmin(np.abs(z+500))
    rho_f = f.rho(p, temp, v=0.)
    rho_f_vmax=f.rho(p, temp, v=f.piston.vol_max)
    rho_f_vmin=f.rho(p, temp, v=f.piston.vol_min)
    #
    ax.fill_betweenx(z, rho_f_vmax, rho_w, where=rho_f_vmax>=rho_w, facecolor='red', interpolate=True)
    ax.plot(rho_w, z, 'b', label='rho water')
    ax.plot(rho_f_vmax, z, '-+', color='orange', label='rho float vmax', markevery=10)
    ax.plot(rho_f_vmin, z, '--', color='orange', label='rho float vmin')
    ax.legend()
    ax.set_xlabel('[kg/m^3]')
    ax.set_ylim((np.amin(z),0.))
    ax.set_ylabel('z [m]')
    ax.grid()

#
def plot_log(f, z_target=None, eta=None, title=None):
    log = f.log
    t = log.t
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


# build scenarios
def descent(Tmax, zt, f, waterp):
    ''' Contruct trajectory of a descent to some depth
    Parameters
    ----------
    Tmax: float
        Time length of the trajectory in SECONDS
    zt: target depth level
    f: float object
        Used to compute maximum accelerations
    waterp: water profile object
        Used to compute maximum accelerations
        
    '''
    # compute bounds on motions
    fmax, fmin, afmax, wmax = f.compute_bounds(waterp,-500.)
    # build time line
    t = np.arange(0.,Tmax,1.)
    # build trajectory
    z_target = np.zeros_like(t)
    dzdt_target = -t*afmax/2./f.m
    dzdt_target[np.where(-dzdt_target>wmax)]=-wmax
    z_target = np.cumsum(dzdt_target*1.)
    z_target[np.where(z_target<zt)] = zt

    # convert to callable function
    return interp1d(t, z_target, kind='linear', fill_value='extrapolate')


#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


class piston():
    ''' Piston object, facilitate float buoyancy control
    '''
    
    def __init__(self,**kwargs):
        """ Piston object

        Parameters
        ----------
        a: float [m]
            piston radius
        phi: float [rad]
            angle of rotation
        d: float [m]
            current piston displacement
        vol: float [m^3]
            current volume
            vol = d x pi x a^2
        omega: float [rad/s]
            current rotation rate, omega=dphi/dt
        lead: float [m]
            screw lead (i.e. displacement after one screw revolution)
            d = phi/2/pi x lead
        phi_max: float [rad]
            maximum angle of rotation
        phi_min: float [rad]
            minimum angle of rotation
        d_max: float [m]
            screw displacement when piston is fully out
        d_min: float [m]
            screw displacement when piston is fully in
        vol_max: float [m^3]
            volume when piston is fully out
        vol_min: float [m^3]
            volume when piston is fully in
        omega_max: float [rad/s]
            maximum rotation rate
        omega_min: float [rad/s]
            minimum rotation rate

        """
        # default parameters: ENSTA float
        params = {'a': 0.025, 'phi': 0., 'd': 0., 'vol': 0., 'omega': 0., 'lead': 0.0175, \
                  'phi_min': 0., 'd_min': 0., 'd_max': 0.07, 'vol_min': 0., \
                  'omega_max': 124.*2.*np.pi/60., 'omega_min': 12.4*2.*np.pi/60.}
        #
        params.update(kwargs)
        for key,val in params.items():
            setattr(self,key,val)
        # assumes here volumes are given
        #if 'd_min' not in kwargs:
        #    self.d_min = self.vol2d(self.vol_min)
        # (useless as will reset d_min to 0.)
        if 'vol_max' in kwargs:
            if 'd_max' in kwargs:
                self.d_max=kwargs['d_max']
                self.vol_min = self.d2vol_no_volmin(self.d_min)
            self.d_max = self.vol2d(self.vol_max)
            print('Piston max displacement set from max volume')
        elif 'd_max' in kwargs:
            self.vol_max = self.d2vol(self.d_max)
            print('Piston max volume set from max displacement')
        else:
            print('You need to provide d_max or vol_max')
            sys.exit()
        #
        self.phi_max = self.d2phi(self.d_max)
        self.update_dvdt()

    def __repr__(self):
        strout='Piston parameters and state: \n'
        strout+='  a     = %.2f cm        - piston radius\n'%(self.a*1.e2)
        strout+='  phi   = %.2f rad       - present angle of rotation\n'%(self.phi)
        strout+='  d     = %.2f cm        - present piston displacement\n'%(self.d*1.e2)
        strout+='  vol   = %.2f cm^3      - present volume addition\n'%(self.vol*1.e6)
        strout+='  lead  = %.2f cm        - screw lead\n'%(self.lead*1.e2)
        strout+='  phi_max = %.2f deg     - maximum rotation\n'%(self.phi_max*1.e2)
        strout+='  phi_min = %.2f deg     - minimum rotation\n'%(self.phi_min*1.e2)
        strout+='  d_max = %.2f cm        - maximum piston displacement\n'%(self.d_max*1.e2)
        strout+='  d_min = %.2f cm        - minimum piston displacement\n'%(self.d_min*1.e2)
        strout+='  vol_min = %.2f cm^3    - min volume displaced\n'%(self.vol_min*1.e6)
        strout+='  vol_max = %.2f cm^3    - max volume displaced\n'%(self.vol_max*1.e6)
        strout+='  omega_max = %.2f deg/s - maximum rotation rate\n'%(self.omega_max*180./np.pi)
        strout+='  omega_min = %.2f deg/s - minimum rotation rate\n'%(self.omega_min*180./np.pi)
        return strout

#------------------------------------------- update methods -----------------------------------------------

    def update(self, dt, dvdt):
        """ Update piston position given time interval and desired volume change

        Parameters
        ----------
        dt: float
            time interval
        dvdt: float
            desired volume change
        """
        omega=self.dvdt2omega(dvdt)
        self.update_omega(omega)
        self.update_dvdt()
        self.update_phi(dt)

    def update_phi(self,dt):
        self.phi+=self.omega*dt
        self._checkbounds()
        self._bcast_phi()

    def update_vol(self,vol):
        self.phi=self.vol2phi(vol)
        self._checkbounds()
        self._bcast_phi()

    def update_omega(self,omega):
        if np.abs(omega)<self.omega_min:
            self.omega=0.
        else:
            self.omega=np.sign(omega)*np.amin([np.abs(omega),self.omega_max])

    def update_dvdt(self):
        self.dvdt = self.omega2dvdt(self.omega)

    def _bcast_phi(self):
        self.d = self.phi2d(self.phi)
        self.vol = self.phi2vol(self.phi)

#------------------------------------------- conversion methods -----------------------------------------------

    def omega2dvdt(self,omega):
        # /np.pi*np.pi has been simplified
        return omega*self.lead/2.*self.a**2

    def dvdt2omega(self,dvdt):
        # /np.pi*np.pi has been simplified
        return dvdt/(self.lead/2.*self.a**2)

    def phi2d(self,phi):
        return self.d_min+(phi-self.phi_min)/2./np.pi*self.lead

    def phi2vol(self,phi):
        return self.d2vol(self.phi2d(phi))

    def d2phi(self,d):
        return self.phi_min+(d-self.d_min)*2.*np.pi/self.lead

    def d2vol(self,d):
        return self.vol_min+(d-self.d_min)*np.pi*self.a**2

    def d2vol_no_volmin(self,d):
        return self.vol_max+(d-self.d_max)*np.pi*self.a**2
        
    def vol2d(self,vol):
        return self.d_min+(vol-self.vol_min)/(np.pi*self.a**2)

    def vol2phi(self,vol):
        return self.d2phi(self.vol2d(vol))

    def _checkbounds(self):
        self.phi = np.amin([np.amax([self.phi,self.phi_min]),self.phi_max])




#------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------


#
class waterp():
    ''' Data holder for a water column based on the World Ocean Atlas (WOA) climatology, see:
    https://www.nodc.noaa.gov/OC5/woa13/
    
    Parameters
    ----------
    lon: float
        longitude of the selected location
    lat: float
        latitude of the selected location
    
    '''
    
    def __init__(self,lon,lat):
        self._woa=True
        #
        self._tfile = 'woa13_decav_t00_01v2.nc'
        nc = Dataset(self._tfile,'r')
        #
        glon = nc.variables['lon'][:]
        glat = nc.variables['lat'][:]
        ilon = np.argmin(np.abs(lon-glon))
        ilat = np.argmin(np.abs(lat-glat))
        self.lon = glon[ilon]
        self.lat = glat[ilat]
        #
        self.z = -nc.variables['depth'][:]
        self.p = gsw.p_from_z(self.z,self.lat)
        #
        self.temp = nc.variables['t_an'][0,:,ilat,ilon]
        nc.close()
        #
        self._sfile = 'woa13_decav_s00_01v2.nc'
        nc = Dataset(self._sfile,'r')
        self.s = nc.variables['s_an'][0,:,ilat,ilon]
        nc.close()
        # derive absolute salinity and conservative temperature
        self.SA = gsw.SA_from_SP(self.s, self.p, self.lon, self.lat)
        self.CT = gsw.CT_from_t(self.SA, self.temp, self.p)
        # isopycnal displacement
        self.eta=0.

    
    def show_on_map(self):
        if self._woa:
            nc = Dataset(self._tfile,'r')
            glon = nc.variables['lon'][:]
            glat = nc.variables['lat'][:]
            temps = nc.variables['t_an'][0,0,:,:]
            nc.close()
            #
            crs=ccrs.PlateCarree()
            plt.figure(figsize=(10, 5))
            ax = plt.axes(projection=crs)
            hdl = ax.pcolormesh(glon,glat,temps,transform = crs,cmap=plt.get_cmap('CMRmap_r'))
            ax.plot(self.lon,self.lat,'*',markersize=10,markerfacecolor='CadetBlue',markeredgecolor='w',transform=crs)
            ax.coastlines(resolution='110m')
            ax.gridlines()
            plt.colorbar(hdl,ax=ax)
            ax.set_title('sea surface temperature [degC]')
            plt.show()
    
    def __repr__(self):
        if self._woa:
            strout = 'WOA water profile at lon=%.0f, lat=%.0f'%(self.lon,self.lat)
        plt.figure(figsize=(7,5))
        ax = plt.subplot(121)
        ax.plot(self.temp,self.z,'k')
        ax.set_ylabel('z [m]')
        ax.set_title('in situ temperature [degC]')
        plt.grid()
        ax = plt.subplot(122)
        ax.plot(self.s,self.z,'k')
        ax.set_yticklabels([])
        #ax.set_ylabel('z [m]')
        ax.set_title('practical salinity [psu]')
        plt.grid()
        return strout

    def get_temp(self,z):
        ''' get in situ temperature
        '''
        #return interp(self.z, self.temp, z)
        SA = interp(self.z, self.SA, z-self.eta)
        CT = interp(self.z, self.CT, z-self.eta)
        p = self.get_p(z)
        return gsw.conversions.t_from_CT(SA,CT,p)
                        
    def get_s(self,z):
        ''' get practical salinity
        '''
        return interp(self.z, self.s, z-self.eta)

    def get_p(self,z):
        ''' get pressure
        '''
        return interp(self.z, self.p, z)

    def get_theta(self,z):
        ''' get potential temperature
        '''
        SA = interp(self.z, self.SA, z-self.eta)
        CT = interp(self.z, self.CT, z-self.eta)
        return gsw.conversions.pt_from_CT(SA,CT)
    
    def get_rho(self,z, ignore_temp=False):
        #s = self.get_s(z,eta=eta)
        p = self.get_p(z)
        #SA = gsw.SA_from_SP(s, p, self.lon, self.lat)
        SA = interp(self.z, self.SA, z-self.eta)
        #
        #temp = self.get_temp(z,eta=eta)
        #if ignore_temp:
        #    temp[:]=self.temp[0]
        #    print('Uses a uniform temperature in water density computation, temp= %.1f degC' %self.temp[0])
        #CT = gsw.CT_from_t(SA, temp, p)
        CT = interp(self.z, self.CT, z-self.eta)
        if ignore_temp:
            CT[:]=self.CT[0]
            print('Uses a uniform conservative temperature in water density computation, CT= %.1f degC' %self.CT[0])
        return gsw.density.rho(SA, CT, p)



# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------
# utils functions

#
def t_modulo_dt(t, dt, dt_step):
    threshold = 0.25 * dt_step / dt
    if np.abs(t / dt - np.rint(t / dt)) < threshold:
        return True
    else:
        return False

#
def interp(z_in, v_in, z_out):
    return interp1d(z_in, v_in, kind='linear', fill_value='extrapolate')(z_out)

