from math import atan
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve, brentq
import matplotlib.pyplot as plt

import yaml
from pprint import pformat
import gsw

from .log import *
from .regulation import control, kalman_filter

# useful parameters
g=9.81
watth = 1e-6/3.6

#-------------------------- float object with dynamics -------------------------

class autonomous_float():

    def __init__(self,**kwargs):
        ''' Autonomous float object, the float is assumed to be cylinder-shaped

        Parameters
        ----------
            m : float mass [kg]
            r : float radius [m]
            L : float length [m]
            V : float volume [m^3]
            gamma : mechanical compressibility [1/dbar]
            alpha : thermal compressibility [1/degC]
            temp0: reference temperature used for thermal compressibility [degC]
            a : float added mass [no dimension]
            c0, c1 : float drag parameters [no dimension]

            All parameters are optional, default values that of ifremer prototype

        '''
        #
        if 'model' not in kwargs:
            self.model = 'ifremer'
        else:
            self.model = kwargs['model']
        #
        if self.model.lower() == 'ifremer':
            params = {'r': 0.07, 'L': 0.8278,
                      'gamma': 3.78039e-06, 'alpha': 6.98e-5,
                      'temp0': 0., 'a': 1., 'c0': 0., 'c1': 1.}
            params['m'] = 11.630
        elif self.model.lower() == 'ensta':
            params = {'r': 0.06, 'L': 0.5,
                      'gamma': 9.30e-5, 'alpha': 0.,
                      'temp0': 0., 'a': 1., 'c0': 0., 'c1': 1.}
            params['m'] = 9.045
        # avec 16 piles : 11.630, avec 32 piles : 13.315
        # compressibility from Paul Troadec's report (p39)
        #
        params.update(kwargs)
        for key,val in params.items():
            setattr(self,key,val)
        # compute the volume of the float:
        if 'V' not in params:
            self.V = np.pi*self.r**2*self.L

        #auxiliary parameters
        self.rho_cte= self.m / self.V #kg.m^-3
        self.gammaV = self.gamma*self.V #m^2

    def __repr__(self):
        strout='Float parameters: \n'
        strout+='  L     = %.2f m      - float length\n'%(self.L)
        strout+='  r     = %.2f m      - float radius\n'%(self.r)
        strout+='  m     = %.2f kg     - float mass\n'%(self.m)
        strout+='  V     = %.2e cm^3   - float volume\n'%(self.V*1.e6)
        strout+='  rho_cte = m/V = %.2e kg.cm^3   - float baseline density\n'%(self.rho_cte*1.e6)
        strout+='  gamma = %.2e /dbar  - mechanical compressibility\n'%(self.gamma)
        strout+='  gamma x V = %.2e cm^3/dbar  - normalized compressibility\n'%(self.gamma*self.V*1.e6)
        strout+='  alpha = %.2e /degC  - thermal compressibility\n'%(self.alpha)
        strout+='  alpha x V = %.2e cm^3/degC  - normalized thermal compressibility\n'%(self.alpha*self.V*1.e6)
        strout+='  temp0 = %.2e  degC  - reference temperature\n'%(self.temp0)
        strout+='  a = %.2e  (no dimension)  - float added mass\n'%(self.a)
        strout+='  c0 = %.2e  (no dimension)  - float drag parameter 0\n'%(self.c0)
        strout+='  c1 = %.2e  (no dimension)  - float drag parameter 1\n'%(self.c1)
        if hasattr(self,'piston'):
            strout+=str(self.piston)
        return strout

    def rho(self, p=None, temp=None, v=None, z=None, waterp=None):
        ''' Returns float density i.e. mass over volume
        '''
        if v is None:
            if hasattr(self, 'v'):
                v = self.v
            else:
                v = 0.
        if p is not None and temp is not None:
            _V = self.V*(1.-self.gamma*p+self.alpha*(temp-self.temp0)) + v
            return self.m/_V
        elif z is not None and waterp is not None:
            # assumes thermal equilibrium
            p, tempw = waterp.get_p(z), waterp.get_temp(z)
            return self.rho(p=p, temp=tempw, v=v)
        else:
            print('You need to provide p/temp or z/waterp')

    def volume(self, **kwargs):
        ''' Returns float volume (V+v)
        '''
        return self.m/self.rho(**kwargs)

    def volume4equilibrium(self, p_eq, temp_eq, rho_eq):
        ''' Find volume that needs to be added in order to be at equilibrium
            prescribed pressure, temperature, density
        '''
        _f = lambda x: rho_eq - self.rho(p=p_eq, temp=temp_eq, v=x)
        v = fsolve(_f, 0.)[0]
        return v

    def z4equilibrium(self, waterp, z0=-10):
        ''' Find depth where float is at equilibrium
        '''
        _f = lambda z: waterp.get_rho(z) - self.rho(z=z, waterp=waterp)
        z = fsolve(_f, z0)[0]
        return z

    def adjust_m(self, p_eq, temp_eq, rho_eq, offset=0.):
        ''' Find mass that needs to be added in order to be at equilibrium
            prescribed pressure, temperature, density

            Parameters
            ----------
            p_eq: Float
                pressure for equilibrium
            temp_eq: Float
                temperature for equilibrium
            rho_eq: float
                density for equilibrium
            offset: float
                offset in grams to add to the float mass

        '''
        def _f(m):
            self.m = m
            return rho_eq - self.rho(p=p_eq, temp=temp_eq)
        m0 = self.m
        self.m = fsolve(_f, self.m)[0]
        self.m += offset*1e-3
        self.rho_cte= self.m / self.V #kg.m^-3
        print('%.1f g '%((self.m-m0)*1.e3+offset) + \
              ' were added to the float in order to be at equilibrium' + \
              ' at %.0f dbar \n'%(p_eq))

    def init_piston(self,**kwargs):
        ''' wrapper around piston initializer
        '''
        self.piston = piston(self.model, **kwargs)
        self.v = self.piston.vol

    def piston_update_vol(self, vol=None):
        ''' wrapper around piston volume update
        '''
        if vol is None:
            vol = self.piston.vol
        else:
            self.piston.update_vol(vol)
        #
        self.v = vol

    def set_piston4equilibrium(self, p_eq, temp_eq, rho_eq):
        ''' Adjust piston to be at equilibrium at a given pressure,
            temperature and density
        '''
        self.v = self.piston.vol
        #
        def _f(vol):
            self.piston_update_vol(vol)
            return rho_eq - self.rho(p=p_eq, temp=temp_eq)
        vol = brentq(_f, self.piston.vol_min, self.piston.vol_max)
        self.piston_update_vol(vol)
        print('Piston reset for equilibrium : vol=%.1e cm^3  ' % (vol*1e6))
        return vol

    def compute_force(self, z, w, waterp, Lv, v=None, sum=True):
        ''' Compute the vertical force exterted on the float
        We assume thermal equilibrium
        Drag is quadratic

        Parameters
        ----------

        z: float
            Depth of the Float
        waterp: water profile object
            Water profile used to compute the upward buoyancy force and water
            vertical flow
        Lv: float
            Length of the float used for drag
        v: float
            Piston volume
        w: float
            float vertical velocity (positive upward)

        '''
        # gravity
        f_b = -self.m*g
        # upward buoyancy
        p, tempw = waterp.get_p(z), waterp.get_temp(z)
        rhow = waterp.get_rho(z)
        # for latter storage
        self.w_tempw = tempw
        self.w_rhow = rhow
        #
        rhof = self.rho(p=p,temp=tempw,v=v)
        f_b += self.m*rhow/rhof*g
        # drag
        cd = self.c1/(2*Lv) * np.abs(w - waterp.detadt)
        if self.c0 != 0:
            N2 = waterp.get_N2(z)
            if N2>0:
                N = np.sqrt(N2)
            else:
                N = 0.
            cd += N * self.c0
        f_d = -self.m * cd * (w - waterp.detadt)
        # we ignore DwDt terms for now
        if sum:
            return f_b+f_d
        else:
            return f_b+f_d, f_b, f_d

    def compute_dforce(self, z, w, waterp, Lv):
        ''' Compute gradients of the vertical force exterted on the float
        '''
        df1 = ( self.compute_force(z+5.e-2, w, waterp, Lv, v=self.v)
                    - self.compute_force(z-5.e-2,0.,waterp,Lv, v=self.v) ) /1.e-1
        df2 = ( self.compute_force(z, w+5.e-3, waterp, Lv, v=self.v)
                    - self.compute_force(z, w-5.e-3, waterp, Lv, v=self.v) ) /1.e-2
        df3 = ( self.compute_force(z, w, waterp, Lv, v=self.v+5.e-5)
                    - self.compute_force(z, w, waterp, Lv, v=self.v-5.e-5) ) /1.e-4
        return df1, df2, df3

    def compute_bounds(self,waterp,zmin,zmax=0.):
        ''' Compute approximate bounds on velocity and acceleration
        '''
        z=zmax
        if hasattr(self,'piston'):
            v=self.piston.vol_max
        else:
            v=None
        fmax=self.compute_force(z, 0., waterp, self.L, v=v) # upward force
        #
        z=zmin
        if hasattr(self,'piston'):
            v=self.piston.vol_min
        else:
            v=None
        fmin=self.compute_force(z, 0., waterp, self.L, v=v) # downward force
        #
        afmax = np.amax((np.abs(fmax),np.abs(fmin)))
        wmax = np.sqrt( afmax * self.m*2*self.L / self.c1)
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
        
    def get_drag_velocity(self, dm, Lv=None):
        """ From a mass offset compute the drag velocity
        
        Parameters:
        -----------
        dm: float
            mass offset [kg]
        """
        if not Lv:
            if hasattr(self, 'Lv'):
                Lv = self.Lv
            else:
                Lv = self.L
        return np.sqrt( g*dm/self.m * 2*Lv / self.c1)

    def time_step(self, waterp, T=600., dt_step=1.,
                  z=None, w=None, v=None, t0=0., Lv=None,
                  ctrl=None, z_target=None,
                  kalman=None,
                  eta=lambda t: 0.,
                  log=True,
                  dt_log=10.,
                  log_nrg=True, p_float=1.e5,
                  verbose=0,
                  **kwargs):
        ''' Time step the float position given initial conditions

        Parameters
        ----------

        waterp: water profile object
                Contains information about the water profile
        T: float
            Length of the simulation in seconds [s]
        dt_step: float
            Simulation time step [s]
        z: float
            Initial position, 0. at the surface and negative downward [m]
        w: float
            Initial vertical velocity, negative for downward motions [m.s^-1]
        v: float
            Initial volume adjustement (total volume is V+v) [m^3]
        t0: float
            Initial time [t]
        Lv: float
            Drag length scale [m]
        z_target: function
            Target trajectory as a function of time [m]
        w_target: function
            Target velocity as a function of time [m.^s-1]
        ctrl: dict
            Contains control parameters
        eta: function
            Isopycnal displacement as a function of time
        log: boolean
            List of variables that will logged, default to True
        dt_log: float
            Time interval between log storage
        log_nrg: boolean, default is True
            Turns on/off nrg computation and storage
        p_float: float [Pa]
            Internal float pressure in Pa
        '''
        t=t0
        #
        if z is not None:
            self.z = z
        elif not hasattr(self,'z'):
            self.z = 0.
        #
        if w is not None:
            self.w = w
        elif not hasattr(self,'w'):
            self.w = 0.
        self.dwdt = 0.
        #
        if Lv is not None:
            self.Lv = Lv
        else:
            self.Lv = self.L
        # log init
        _log = {'state':['z','w','v','dwdt','w_temp','w_rho']}
        _log['dynamics'] = ['acceleration','buoyancy','drag']
        # kalman initialisation
        if kalman:
            self.init_kalman(kalman, self.w, self.z, verbose)
            _log['kalman'] = ['z','w', 'V_e', 'gamma_e', \
                              'dwdt', 'dwdt_diff'] + \
                              ['gamma_diag%d'%i for i in range(4)]
            print(self.kalman)
        #
        if ctrl:
            self.init_control(ctrl, v, dt_step)
            print(self.ctrl)
            _log['piston'] = ['u', 'work']
            _log['control'] = []
            u = 0.
        elif v is None:
            if not hasattr(self,'v'):
                self.v = 0.
        else:
            self.v = v
        v0 = self.v
        #
        if log:
            self.log = {logname:logger(vars) for logname, vars in _log.items()}
            _log_now = False
        else:
            self.log = None
        if log and ctrl:
            # used to update contributions to u
            self.ctrl.log = self.log['control']
            #self.ctrl._log_now = self._log_now
        #
        print('Start time stepping for %d min ...'%(T/60.))
        #
        _force=0.
        while t<t0+T:
            #
            if log and (dt_log is not None) and t_modulo_dt(t, dt_log, dt_step):
                _log_now = True
            else:
                _log_now = False
            # get vertical force on float
            waterp.update_eta(eta, t) # update isopycnal displacement
            _force, _force_b, _force_d = \
                    self.compute_force(self.z, self.w, waterp,
                                       self.Lv, v=self.v,
                                       sum=False)
            #
            # update kalman state estimation
            if kalman and t_modulo_dt(t, self.kalman.dt, dt_step):
                self.kalman.update_kalman(self.v, self.z)
            #
            # control starts here
            if ctrl:
                # activate control only if difference between the target and actual vertical
                # position is more than the dz_nochattering threshold
                _ctrl_now = ( (np.abs(self.z-z_target(t)) > self.ctrl.dz_nochattering) \
                          and t_modulo_dt(t, self.ctrl.dt, dt_step) )
                if _ctrl_now:
                    u = self.get_control(z_target, t, _log_now)
                if self.ctrl.continuous:
                    _dt = dt_step
                elif _ctrl_now:
                    _dt = self.ctrl.dt
                else:
                    _dt = 0.
                if self.ctrl.mode == 'kalman_feedback1' and _dt>0:
                    u = (u-self.piston.vol)/_dt
                _v0 = self.piston.vol
                self.piston.update(_dt, u)
                self.v = self.piston.vol
                # energy integration, 1e4 converts from dbar to Pa
                if self.v != _v0:
                    self.piston_work += _dt \
                                * np.abs((waterp.get_p(self.z)*1.e4 - p_float)*u) \
                                * watth /self.piston.efficiency
            #
            # log data
            if _log_now:
                self.log['state'].store(time=t, z=self.z, w=self.w, v=self.v,
                                        dwdt=_force/(1+self.a)/self.m,
                                        w_temp=self.w_tempw,
                                        w_rho=self.w_rhow)
                self.log['dynamics'].store(time=t,
                                           acceleration=_force/(1+self.a)/self.m,
                                           buoyancy=_force_b/(1+self.a)/self.m,
                                           drag=_force_d/(1+self.a)/self.m,
                                           )
                if ctrl:
                    _info = {'time': t, 'u': u, 'work': self.piston_work}
                    self.log['piston'].store(**_info)
                if kalman:
                    _k = self.kalman
                    _dwdt = _k.A_coeff*(_k.x_hat[3] \
                                        -_k.x_hat[2]*_k.x_hat[1] \
                                        + self.v) \
                            +_k.B_coeff*_k.x_hat[0]*np.abs(_k.x_hat[0])
                    #
                    self.log['kalman'].store(time=t,
                               z=_k.x_hat[1], w=_k.x_hat[0],
                               V_e=_k.x_hat[3], gamma_e=_k.x_hat[2],
                               **{'gamma_diag%d'%i: _k.gamma[i,i] for i in range(4)},
                               dwdt=_dwdt,
                               dwdt_diff=_force/(1+self.a)/self.m - _dwdt)
            #
            # update variables
            self.z += dt_step*self.w
            self.z = np.amin((self.z,0.))
            self.w += dt_step*_force/(1+self.a)/self.m
            self.dwdt = _force/(1+self.a)/self.m
            t+=dt_step

        # cleanup logs
        #for key, item in self.log.items():
        #    item.cleanup()

        print('... time stepping done')


    def init_control(self, ctrl, v, dt_step):
        if v is not None:
            self.piston.update_vol(v)
            self.piston.reset_phi_float()
        self.v = self.piston.vol
        self.piston_work = 0.
        if ctrl:
            ctrl_default = {'dt': dt_step, 'dz_nochattering': 0.}
            if ctrl['mode'] == 'sliding':
                ctrl_default.update({'tau': 60., 'mode': 'sliding',
                                     'waterp': waterp, 'Lv': self.L})
            elif ctrl['mode'] == 'pid':
                ctrl_default.update({'error':0.,'integral': 0.})
            elif 'feedback' in ctrl['mode']:
                ctrl_default.update({'tau': 3.25,  # Set the root of feed-back regulation # s assesed by simulation
                                     'nu': 0.03*2./np.pi, # Set the limit speed : 3cm/s assesed by simulation
                                     'delta': 0.11, #length scale that defines the zone of influence around the target depth, assesed by simulation
                                     'gamma': self.gamma, #mechanical compressibility [1/dbar]
                                     'm': self.m, 'a': self.a,
                                     'Lv': self.Lv, 'c1': self.c1,
                                     'gammaV': self.gammaV,
                                     'rho_cte': self.rho_cte,
                                     })
            if 'kalman' in ctrl['mode']:
                _k = self.kalman
                ctrl_default.update({
                                     'm': _k.m, 'a': _k.a,
                                     'Lv': _k.Lv, 'c1': _k.c1,
                                     'gammaV': _k.gamma_e0,
                                     'rho_cte': _k.rho_cte,
                                     })
            #
            ctrl_default.update(ctrl)
            self.ctrl = control(**ctrl_default)
            #print('Control parameters:') # should be moved into control class
            #for key, val in ctrl_default.items():
            #    if key not in ['waterp','f']:
            #        print('  '+key+' = '+str(val))

    def get_control(self, z_target, t, log):
        _c = self.ctrl
        if _c.mode=='sliding':
            u = _c.get_u_sliding(z_target, t, self.z, self.w, self)
        elif _c.mode=='pid':
            u = _c.get_u_pid(z_target, t, self.z)
        elif _c.mode=='feedback2':
            u = _c.get_u_feedback2(z_target, t, self.z, self.w,
                                   self.dwdt, self.gammaV, log)
        elif _c.mode=='kalman_feedback1':
            _k = self.kalman
            u = _c.get_u_feedback1(z_target, t, -_k.x_hat[1], -_k.x_hat[0],
                                   _k.x_hat[3], _k.x_hat[2], log)
        elif _c.mode=='kalman_feedback2':
            _k = self.kalman
            _dwdt = _k.A_coeff*(_k.x_hat[3] \
                                -_k.x_hat[2]*_k.x_hat[1] \
                                +self.v) \
                    +_k.B_coeff*_k.x_hat[0]*np.abs(_k.x_hat[0])
            u = _c.get_u_feedback2(z_target, t, -_k.x_hat[1], -_k.x_hat[0],
                                  _dwdt, _k.x_hat[2], log)
        else:
            print('%s is not a valid control method'%_c.mode)
            return
        return u

    def init_kalman(self, kalman, w, z, verbose):
        _params = {'m': self.m, 'a':self.a, 'rho_cte': self.rho_cte,
                   'c1':self.c1, 'Lv': self.L,
                   'gamma_e0': self.gammaV, 'V_e0': 0.,
                   'verbose': verbose}
        if type(kalman) is dict:
            _params.update(kalman)
        _x0 = [-w, -z, _params['gamma_e0'], _params['V_e0']]
        self.kalman = kalman_filter(_x0, **_params)

    def plot_logs(self, **kwargs):
        ''' wrapper around plot_logs
        '''
        plot_logs(self.log, self, **kwargs)


#---------------------------- piston -------------------------------------------
class piston():
    ''' Piston object, facilitate float buoyancy control
    '''

    def __init__(self, model='ifremer', **kwargs):
        """ Piston object

        Parameters
        ----------
        r: float [m]
            piston radius
        phi: float [rad]
            angle of rotation
        d: float [m]
            current piston displacement
        vol: float [m^3]
            current volume
            vol = d x pi x r^2
        omega: float [rad/s]
            current rotation rate, omega=dphi/dt
        lead: float [m]
            screw lead (i.e. displacement after one screw revolution)
            d = phi/2/pi x lead
        tick_per_turn: [no dimension]
            number of notches on the thumbwheel of the piston
        d_increment: [m]
            smallest variation of translation motion for the piston
        #increment_error: [no dimension]
        #    coefficient measuring the accuracy on the smallest variation of
        #    translation motion for the piston (coefficient >= 1)
        vol_increment: [m^3]
            smallest variation of volume possible for the piston
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
        #omega_min: float [rad/s] (not relevant anymore with d_increment and co)
        #    minimum rotation rate
        efficiency: float [<1]
            Piston efficiency, i.e. mechanical work produced over electrical
            work supplied

        """

        if model.lower() == 'ifremer':
            # default parameters: IFREMER float
            params = {'r': 0.0195/2, 'phi': 0., 'd': 0., 'vol': 0.,
                      'omega': 0., 'lead': 1,
                      'phi_min': 0., 'd_min': 0., 'd_max': 0.09,'vol_min': 0.,
                      'translation_max': 0.12,
                      'translation_min': 0.03,
                      'efficiency': .1,
                      'd_increment': 0.12*4./5600.}
        elif model.lower() == 'ensta':
            # default parameters: ENSTA float
            params = {'r': 0.025, 'phi': 0., 'd': 0., 'vol': 0.,
                      'omega': 0., 'lead': 0.00175, 'tick_per_turn': 48,
                      'phi_min': 0., 'd_min': 0., 'd_max': 0.07,
                      'vol_max': 1.718e-4, 'vol_min': 0.,
                      'omega_max': 60./48*2.*np.pi,
                      'efficiency':.1}
            self.d_increment = params['lead']/params['tick_per_turn']

            #translation_max = d_increment*(225 pulses par seconde)  (vitesse de translation max en m/s)
            #translation_min fixe arbitrairement pour l'instant

            #d_increment le 4 vient du facteur 4 de la roue codeuse de thomas
            #d_increment = 0.12/5600 ou 0.090/4200 selon la prise en compte ou non du gros piston

            #dmax = 0.102 ancienne valeur pour IFREMER

            #verifier si importance lead et angles lors de la regulation, ecraser parametres redondants
            #vol_max = 0.090*np.pi*(0.0195/2)**2+0.030*np.pi*(0.080/2)**2 = 1.777e-4

    	#48 encoches
    	#vitesse max de 60 encoches par seconde
        #
        params.update(kwargs)
        for key,val in params.items():
            setattr(self,key,val)
        # assumes here volumes are given
        #if 'd_min' not in kwargs:
        #    self.d_min = self.vol2d(self.vol_min)
        # (useless as will reset d_min to 0.)
        if 'vol_max' in params:
            if 'd_max' in params:
                self.d_max=params['d_max']
                self.vol_min = self.d2vol_no_volmin(self.d_min)
                print('Piston vol_min reset from d_min, d_max and vol_max')
            else:
                self.d_max = self.vol2d(self.vol_max)
                print('Piston max displacement set from max volume')
        elif 'd_max' in params:
            self.vol_max = self.d2vol(self.d_max)
            print('Piston max volume set from max displacement')
        else:
            print('You need to provide d_max or vol_max')
            sys.exit()
        if 'vol' not in kwargs:
            self.vol = self.vol_max
        if 'translation_max' in params:
            self.omega_max = params['translation_max']*2.*np.pi/self.lead
        if 'translation_min' in params:
            self.omega_min = params['translation_min']*2.*np.pi/self.lead
        #
        self.phi_max = self.d2phi(self.d_max)
        self.update_dvdt()
        #
        self.u_max = self.omega2dvdt(self.omega_max)
        #
        self.phi_float = self.phi
        self.phi_increment = self.d2phi(self.d_min+self.d_increment) \
                              -self.d2phi(self.d_min)
        # not sure if vol_increment is still used
        self.vol_increment = self.d_increment*np.pi*self.r**2 #*self.increment_error
        self.tick_to_volume = self.vol_increment

    def __repr__(self):
        strout='Piston parameters and state: \n'
        strout+='  r     = %.2f cm        - piston radius\n'%(self.r*1.e2)
        #strout+='  phi   = %.2f rad       - present angle of rotation\n'%(self.phi)
        strout+='  d     = %.2f cm        - present piston displacement\n'%(self.d*1.e2)
        strout+='  vol   = %.2f cm^3      - present volume addition\n'%(self.vol*1.e6)
        #strout+='  lead  = %.2f cm        - screw lead\n'%(self.lead*1.e2)
        #strout+='  tick_per_turn  = %.2f no dimension        - number of notches on the thumbwheel of the piston\n'%(self.tick_per_turn)
        strout+='  d_increment  = %.2e mm        - smallest variation of translation motion for the piston\n'%(self.d_increment*1e3)
        strout+='  vol_increment  = %.2e cm^3        - smallest variation of volume possible for the piston\n'%(self.vol_increment*1e6)
        #strout+='  phi_max = %.2f deg     - maximum rotation\n'%(self.phi_max*1.e2)
        #strout+='  phi_min = %.2f deg     - minimum rotation\n'%(self.phi_min*1.e2)
        strout+='  d_max = %.2f mm        - maximum piston displacement\n'%(self.d_max*1.e3)
        strout+='  d_min = %.2f mm        - minimum piston displacement\n'%(self.d_min*1.e3)
        strout+='  vol_min = %.2f cm^3    - min volume displaced\n'%(self.vol_min*1.e6)
        strout+='  vol_max = %.2f cm^3    - max volume displaced\n'%(self.vol_max*1.e6)
        #strout+='  omega_max = %.2f deg/s - maximum rotation rate\n'%(self.omega_max*180./np.pi)
        #strout+='  omega_min = %.2f deg/s - minimum rotation rate\n'%(self.omega_min*180./np.pi)
        strout+='  u_max = %.2f cm^3/s - maximum volume rate of change\n'%(self.u_max*1e6)
        strout+='  efficiency = %.2f - mechanical work produced / electrical work supplied\n'%(self.efficiency)
        return strout

#-------------------------- update methods -------------------------------------

    def update(self, dt, dvdt):
        """ Update piston position given time interval and a volume rate of change

        Parameters
        ----------
        dt: float
            time interval
        dvdt: float
            desired volume rate of change
        """
        omega=self.dvdt2omega(dvdt)
        self.update_omega(omega)
        self.update_dvdt()
        self.update_phi(dt)

    def update_phi(self,dt):
        self.phi_float += self.omega*dt
        dphi = self.phi_float - self.phi
        self.phi += np.trunc(dphi/self.phi_increment)*self.phi_increment
        self._checkbounds()
        self._bcast_phi()

    def update_vol(self,vol):
        self.phi = self.vol2phi(vol)
        self._checkbounds()
        self._bcast_phi()

    def update_d(self,d):
        self.phi = self.d2phi(d)
        self._checkbounds()
        self._bcast_phi()

    def update_omega(self, omega):
        #if np.abs(omega)<self.omega_min:
        #    self.omega=0.
        #else:
        self.omega=np.sign(omega)*min([np.abs(omega),self.omega_max])

    def update_dvdt(self):
        self.dvdt = self.omega2dvdt(self.omega)

    def _bcast_phi(self):
        self.d = self.phi2d(self.phi)
        self.vol = self.phi2vol(self.phi)

    def _checkbounds(self):
        self.phi = min([max([self.phi,self.phi_min]),self.phi_max])
        self.phi_float = min([max([self.phi_float,self.phi_min]),self.phi_max])

    def reset_phi_float(self):
        self.phi_float = self.phi

#--------------------- conversion methods --------------------------------------

    def omega2dvdt(self,omega):
        #to compute omega for the ENSTA float:
        #the motor of the piston needs 48 notches to complete a full rotation
        #it can reach until 30 rotations a seconde
        #so omega = 2*pi*30/48 = 3.9 rad/s
        return omega*self.lead/2.*self.r**2

    def dvdt2omega(self,dvdt):
        return dvdt/(self.lead/2.*self.r**2)

    def phi2d(self,phi):
        return self.d_min+(phi-self.phi_min)/2./np.pi*self.lead

    def phi2vol(self,phi):
        return self.d2vol(self.phi2d(phi))

    def d2phi(self,d):
        return self.phi_min+(d-self.d_min)*2.*np.pi/self.lead

    def d2vol(self,d):
        return self.vol_min+(d-self.d_min)*np.pi*self.r**2

    def d2vol_no_volmin(self,d):
        return self.vol_max+(d-self.d_max)*np.pi*self.r**2

    def vol2d(self,vol):
        return self.d_min+(vol-self.vol_min)/(np.pi*self.r**2)

    def vol2phi(self,vol):
        return self.d2phi(self.vol2d(vol))

# ----------------------------- balast -----------------------------------------

_bvariables = ['mass_in_water', 'mass_in_air',
               'water_temperature',
               'piston_displacement']

class balast(object):

    # in db
    pressure = 1.
    # brest:
    lon, lat = 4.5, 38.5

    def __init__(self, file=None):
        ''' Class handling the balasting procedure
        '''
        self._d = {}
        if file is not None:
            self.load(file)

    def __repr__(self):
        return pformat(self._d, indent=2, width=1)

    def add_balast(self, name, comments='None', **kwargs):
        ''' Add a mass measurement

        Parameters
        ----------
        name: str
            Name of the measurement with for example date, location
        mass_in_water: float
            in grammes
        mass_in_air: float
            in grammes
        water_temperature: float
            in degC
        water_salinity: float
            in psu
        piston_displacement: float
            in cm
        comments: str
            additional comments
        '''
        if not all(v in kwargs for v in _bvariables):
            print('All following variables should be passed as kwargs:')
            print(_bvariables)
            print('Abort')
            return
        self._d[name] = {v: kwargs[v] for v in _bvariables}
        _d = self._d[name]
        # derive absolute salinity and conservative temperature
        if 'water_conductivity' in kwargs:
            _SP = gsw.SP_from_C(kwargs['water_conductivity'],
                                 kwargs['water_temperature'],
                                 self.pressure
                                 )
            _d['water_salinity'] = _SP
            print('Practical salinity derived from conductivity: SP=%.2f'%_SP)
        else:
            _SP = kwargs['water_salinity']
        SA = gsw.SA_from_SP(_SP,self.pressure, self.lon, self.lat)
        CT = gsw.CT_from_t(SA, kwargs['water_temperature'], self.pressure)
        rho_water = gsw.density.rho(SA, CT, self.pressure) # in situ density
        # store derived variables
        _d['SA'] = SA
        _d['CT'] = CT
        _d['rho_water'] = rho_water
        _d['V'] = _d['mass_in_air']/1e3/rho_water
        _d['comments'] = comments

    def compute_mass_adjustment(self, f, w=None, verbose=False, **kwargs):
        """
        \delta_m = -V \delta \rho_w - \rho_w \delta V
        """
        if w is not None:
            z = -self.pressure
            water_temperature = w.get_temp(z)
            water_salinity = w.get_s(z)
            lon, lat = w.lon, w.lat
        else:
            water_temperature = kwargs['temperature']
            if 'salinity' in kwargs:
                water_salinity = kwargs['salinity']
            elif 'conductivity' in kwargs:
                water_salinity = gsw.SP_from_C(kwargs['conductivity'],
                                     kwargs['temperature'],
                                     self.pressure
                                     )
            lon, lat = kwargs['lon'], kwargs['lat']
            if isinstance(lon,str):
                lon = ll_conv(lon)
            if isinstance(lat,str):
                lat = ll_conv(lat)
        SA = gsw.SA_from_SP(water_salinity, self.pressure, lon, lat)
        CT = gsw.CT_from_t(SA, water_temperature, self.pressure)
        rho_water = gsw.density.rho(SA, CT, self.pressure) # new in situ density
        #
        for name, b in self._d.items():
            delta_rho_water = rho_water - b['rho_water']
            #
            f.m = b['mass_in_air']*1e-3  # need to convert to kg
            f.piston.update_d(b['piston_displacement']*1e-2) # need to convert to m
            f.piston_update_vol()
            # reset volume and temperature reference
            f.V = b['V'] - f.v # reset volume to match volume infered from balasting
            f.temp0 = b['water_temperature']
            #
            Vb = f.volume(p=self.pressure, temp=b['water_temperature'])
            V = f.volume(p=self.pressure, temp=water_temperature)
            delta_V = V - Vb
            #
            delta_m0 = -b['mass_in_water']*1e-3
            delta_m1 = V * delta_rho_water + rho_water * delta_V
            delta_m = delta_m0 + delta_m1
            print('-- According to balast %s, you need add %.1f g'
                    %(name, delta_m*1e3))
            print('(if the weight is external to the float, remember this must be the value in water)')
            #
            if verbose:
                print('Independent mass correction (balast -mass_in_water): %.1f [g]'%(delta_m0*1e3))
                print('  Water change mass correction: %.1f [g]'%(delta_m1*1e3))
                delta_m_t = - gsw.density.alpha(SA, CT, self.pressure) \
                               * (CT-b['CT']) * b['rho_water'] *  V
                print('    T contribution: %.1f [g]'%(delta_m_t*1e3))
                delta_m_s = gsw.density.beta(SA, CT, self.pressure) \
                              * (SA-b['SA']) * b['rho_water'] * V
                print('    S contribution: %.1f [g]'%(delta_m_s*1e3))
                print('  New/old in situ density: %.2f, %.2f  [kg/m^3] '%(rho_water, b['rho_water']))
                print('    difference: %.1f [kg/m^3]'%(delta_rho_water))
                print('  New/old float volume: %.2e, %.2e  [m^3] '%(V, b['V']))
                print('    difference: %.1f [cm^3]'%(delta_V*1e6))

    def store(self, file):
        _d = dict(self._d)
        _v = ['SA','CT','rho_water', 'V', 'water_salinity']
        for name, b in _d.items():
            for v in _v:
                b[v]=float(b[v])
        with open(file, 'w') as ofile:
            documents = yaml.dump(self._d, ofile)
        print('Balasting data stored in %s'%file)

    def load(self, file):
        with open(file) as ifile:
            d = yaml.full_load(ifile)
        # should check that content of d is valid
        self._d.update(d)
        print('File %s loaded'%file)

# ----------------------------- utils ------------------------------------------

def ll_conv(ll):
    """ Returns lon or lat as floating degree and vice versa
    string format should be '43N09.007'
    This piece of code should be elsewhere
    """
    if isinstance(ll,str):
        if 'N' in ll:
            _ll = ll.split('N')
            sign = 1.
        elif 'S' in ll:
            _ll = ll.split('S')
            sign = -1.
        elif 'E' in ll:
            _ll = ll.split('E')
            sign = 1.
        elif 'W' in ll:
            _ll = ll.split('W')
            sign = -1.
        return sign*(float(_ll[0])+float(_ll[1])/60.)
    elif isinstance(ll,float):
        return '%dX%.6f'%(int(ll),(ll-int(ll))*60.)

def plot_float_density(z, f, waterp, mid=False, ax=None):
    ''' Plot float density with respect to that of water profile

    Parameters
    ----------
    z: ndarray
        depth grid
    f: float object
    waterp: water profile object
    mid: boolean, True for curve at piston mid displacement
    '''
    #
    rho_w, p, temp = waterp.get_rho(z), waterp.get_p(z), waterp.get_temp(z)
    #
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    #
    #rho_f = f.rho(p=p, temp=temp, v=0.)
    rho_f_vmax = f.rho(p=p, temp=temp, v=f.piston.vol_max)
    rho_f_vmin = f.rho(p=p, temp=temp, v=f.piston.vol_min)
    #
    ax.fill_betweenx(z, rho_f_vmax, rho_w, where=rho_f_vmax>=rho_w, facecolor='red', interpolate=True)
    ax.plot(rho_w, z, 'b', label='rho water')
    ax.plot(rho_f_vmax, z, '-+', color='orange', label='rho float vmax', markevery=10)
    ax.plot(rho_f_vmin, z, '--', color='orange', label='rho float vmin')
    if mid:
        #rho_f_vmid=f.rho(p=p, temp=temp, v=(f.piston.vol_max+f.piston.vol_min)*.5)
        rho_f_vmid=f.rho(p=p, temp=temp, v=mid)
        ax.plot(rho_f_vmid, z, '--', color='grey', label='rho float vmid')
    ax.legend()
    ax.set_xlabel('[kg/m^3]')
    ax.set_ylim((np.amin(z),0.))
    ax.set_ylabel('z [m]')
    ax.grid()
    iz = np.argmin(np.abs(z))
    ax.set_title('extra mass @surface, piston out: %.1f g' \
                    %( ( (f.V+f.piston.vol_max) * rho_w[iz] - f.m)*1e3 ) )
    #
    y_annotation = ax.get_ylim()[1]-.1*(ax.get_ylim()[1]-ax.get_ylim()[0])
    ax.annotate('',
                xy=(ax.get_xticks()[1], y_annotation),
                xytext=(ax.get_xticks()[2], y_annotation),
                arrowprops=dict(arrowstyle="<->"))
    ax.text(ax.get_xticks()[1]*.5+ax.get_xticks()[2]*.5,
            y_annotation,
            '%.1f g'%(f.V*(ax.get_xticks()[2]-ax.get_xticks()[1])*1e3),
            {'ha': 'center', 'va': 'bottom'})
    #        ax.get_ylim()[1]-10,
    return ax

#
def plot_float_volume(z, f, waterp, ax=None):
    ''' Plot float volume with respect to depth

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
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    #
    #iz = np.argmin(np.abs(z+500))
    v = f.volume(p=p, temp=temp, v=0.)
    vmax = f.volume(p=p, temp=temp, v=f.piston.vol_max)
    vmin = f.volume(p=p, temp=temp, v=f.piston.vol_min)
    #
    ax.plot(vmax*1e6, z, '-+', color='orange', label='vmax float', markevery=10)
    ax.plot(vmin*1e6, z, '--', color='orange', label='vmin float')
    ax.legend()
    ax.set_xlabel('[cm^3]')
    ax.set_ylim((np.amin(z),0.))
    ax.set_ylabel('z [m]')
    ax.grid()
    return ax

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

def t_modulo_dt(t, dt, dt_step):
    threshold = 0.25 * dt_step / dt
    if np.abs(t / dt - np.rint(t / dt)) < threshold:
        return True
    else:
        return False

# build scenarios
def descent(Tmax, zt, f=None, waterp=None, wmax=None, zstart=0):
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
    # build time line
    t = np.arange(0.,Tmax,1.)
    # compute bounds on motions
    if f is not None and waterp is not None:
        fmax, fmin, afmax, wmax = f.compute_bounds(waterp,zt)
        dzdt_target = -t*afmax/2./f.m
        dzdt_target[np.where(-dzdt_target>wmax)] = -np.abs(wmax)
    elif wmax is not None:
        dzdt_target = t*0. - np.abs(wmax)
    # build trajectory
    z_target = zstart + np.cumsum(dzdt_target*1.)
    z_target[np.where(z_target<zt)] = zt
    # convert to callable function
    return interp1d(t, z_target, kind='linear', fill_value='extrapolate')
