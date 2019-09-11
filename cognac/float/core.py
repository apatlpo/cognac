from math import atan
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve, brentq
import matplotlib.pyplot as plt

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
            params = {'r': 0.07, 'L': 0.8278, 'gamma': 3.78039e-06, 'alpha': 0.,
                      'temp0': 0., 'a': 1., 'c0': 0., 'c1': 1.}
            params['m'] = 11.630
        elif self.model.lower() == 'ensta':
            params = {'r': 0.06, 'L': 0.5, 'gamma': 9.30e-5, 'alpha': 0.,
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
        strout+='  rho_cte = m/V = %.2e kg.cm^3   - float baseline density\n'%(self.m/self.V*1.e6)
        strout+='  gamma = %.2e /dbar  - mechanical compressibility\n'%(self.gamma)
        strout+='  alpha = %.2e /degC  - thermal compressibility\n'%(self.alpha)
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
        v = fsolve(_f,0.)[0]
        return v

    def z4equilibrium(self, waterp):
        ''' Find depth where float is at equilibrium
        '''
        _f = lambda z: waterp.get_rho(z) - self.rho(z=z, waterp=waterp)
        z = fsolve(_f, 0.)[0]
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
        self.m = fsolve(_f,self.m)[0]
        self.m += offset*1e-3
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

    def compute_force(self, z, w, waterp, Lv, v=None):
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
        f = -self.m*g
        # upward buoyancy
        p, tempw = waterp.get_p(z), waterp.get_temp(z)
        rhow = waterp.get_rho(z)
        rhof = self.rho(p=p,temp=tempw,v=v)
        f += self.m*rhow/rhof*g
        # drag
        if self.c0 != 0:
            print(' !!! linear drag coefficient not implemented yet')
            return None
        f += -self.m*self.c1/(2*Lv) * np.abs(w - waterp.detadt) * (w - waterp.detadt)
        # we ignore DwDt terms for now
        return f

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

    def init_kalman(self, kalman, w, z, gammaE, Ve, verbose):
        _params = {'m': self.m, 'a':self.a, 'rho_cte': self.rho_cte,
                   'c1':self.c1, 'Lv': self.L, 'gammaV': self.gammaV,
                   'verbose': verbose}
        if type(kalman) is dict:
            _params.update(kalman)
        _x0 = [-w, -z, gammaE, Ve]
        self.kalman = kalman_filter(_x0, **_params)

    def time_step(self, waterp, T=600., dt_step=1.,
                  z=None, w=None, v=None, Ve=None, t0=0., Lv=None,
                  piston=False, z_target=None, gammaE=None,
                  ctrl=None,
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
        Ve: float
            Volume offset (total volume is Ve+v) [m^3]
        t0: float
            Initial time [t]
        Lv: float
            Drag length scale [m]
        piston: boolean, default is False
            Turns piston usage on and off [no dimension]
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
        if z is None:
            if not hasattr(self,'z'):
                z=0.
            else:
                z=self.z
        self.z = z

        #
        if Ve is None:
            if not hasattr(self,'Ve'):
                Ve=0.
            else:
                Ve=self.Ve
        self.Ve = Ve

        #
        if w is None:
            if not hasattr(self,'w'):
                w=0.
            else:
                w=self.w
        self.w = w
        self.dwdt = 0.
        #
        _log = {'state':['z','w','v','dwdt']}
        #
        if gammaE is None:
            gammaE = self.gammaV
        #kalman initialisation
        if kalman:
            self.init_kalman(kalman, w, z, gammaE, Ve, verbose)
            _log['kalman'] = ['Ve', 'gammaV',
                              'z_kalman','w_kalman', 'Ve_kalman',
                              'dwdt_kalman', 'gammaE'] + \
                              ['gamma_diag%d'%i for i in range(4)]
        #
        if piston:
            if v is not None:
                self.piston.update_vol(v)
            self.v=self.piston.vol
            self.piston_work = 0.
            u = 0.
            if ctrl:
                ctrl_default={'dt_ctrl': dt_step, 'dz_nochattering': 0.}
                if ctrl['mode'] == 'sliding':
                    ctrl_default = {'tau': 60., 'mode': 'sliding',
                                    'waterp': waterp, 'Lv': self.L, }
                elif ctrl['mode'] == 'pid':
                    ctrl_default = {'error':0.,'integral': 0.}
                elif ctrl['mode'] == 'feedback':
                    ctrl_default = {'tau': 3.25,  # Set the root of feed-back regulation # s assesed by simulation
                                    'nu': 0.03*2./np.pi, # Set the limit speed : 3cm/s assesed by simulation
                                    'delta': 0.11, #length scale that defines the zone of influence around the target depth, assesed by simulation
                                    'gamma': self.gamma, #mechanical compressibility [1/dbar]
                                    'L': self.L, 'c1': self.c1,'m': self.m,
                                    'gammaV': self.gammaV,
                                    'rho_cte': self.rho_cte,
                                    'a': self.a,
                                    'waterp': waterp}
                #
                elif ctrl['mode'] == 'kalman_feedback':
                    ctrl_default= {'tau': 3.25,  # Set the root of feed-back regulation # s assesed by simulation
                                   'nu': 0.03*2./np.pi, # Set the limit speed, 3cm/s assesed by simulation
                                   'delta': 0.11, #length scale that defines the zone of influence around the target depth, assesed by simulation
                                   'kalman': self.kalman}
                #
                #print(ctrl_default)
                ctrl_default.update(ctrl)
                ctrl = ctrl_default
                self.ctrl = ctrl
                for key, val in ctrl_default.items():
                    if key not in ['waterp','f']:
                        print(' ctrl: '+key+' = '+str(val))
                #
                _log['piston'] = ['u','work']
        elif v is None:
            if not hasattr(self,'v'):
                self.v = 0.
        else:
            self.v = v
        v0 = self.v
        #
        if Lv is None:
            Lv = self.L
        self.Lv = Lv
        #
        if log:
            self.log = {logname:logger(vars) for logname, vars in _log.items()}
        #
        print('Start time stepping for %d min ...'%(T/60.))
        #
        _force=0.
        while t<t0+T:
            #
            # get vertical force on float
            waterp.update_eta(eta, t) # update isopycnal displacement
            _force = self.compute_force(self.z, self.w, waterp,
                                        self.Lv, v=self.v)
            #
            # update kalman state estimation
            if kalman and t_modulo_dt(t, self.kalman.dt, dt_step):
                self.kalman.update_kalman(u, self.v, self.z)
            #
            # control starts here
            if piston and ctrl and t_modulo_dt(t, ctrl['dt_ctrl'], dt_step):
                # activate control only if difference between the target and actual vertical
                # position is more than the dz_nochattering threshold
                if np.abs(self.z-z_target(t)) > ctrl['dz_nochattering']:
                    if verbose>0:
                        print('[-w, -z, -dwdt, gammaV, Ve]',
                              [-self.w, -self.z, -self.dwdt, self.gammaV,
                               self.Ve])
                    u = control(self.z, z_target, ctrl, t=t, w=self.w,
                                dwdt=self.dwdt, v=self.v) #, f=ctrl['f'])
                    #
                    v0 = self.piston.vol
                    self.piston.update(dt_step, u)
                    self.v = self.piston.vol
                # energy integration, 1e4 converts from dbar to Pa
                if self.v != v0:
                    self.piston_work += dt_step \
                                * np.abs((waterp.get_p(self.z)*1.e4 - p_float)*u) \
                                * watth /self.piston.efficiency

            # store log
            if log and (dt_log is not None) and t_modulo_dt(t, dt_log, dt_step):
                self.log['state'].store(time=t, z=self.z, w=self.w, v=self.v,
                                        dwdt=_force/self.m)
                if piston:
                    _info = {'time': t, 'u': u, 'work': self.piston_work}
                    self.log['piston'].store(**_info)
                if kalman:
                    _k = self.kalman
                    # Ve, gammae
                    _Ve = _force/(g*self.rho_cte) - self.gammaV * self.z - self.v
                    #_gammae =
                    _dwdt = -_k.A_coeff* \
                        (_k.x_hat[2] + _k.x_hat[3] -_k.gammaV*_k.x_hat[1]) \
                        -_k.B_coeff*abs(_k.x_hat[0])*_k.x_hat[0]
                    []
                    self.log['kalman'].store(time=t,
                               z_kalman=_k.x_hat[1], w_kalman=_k.x_hat[0],
                               gammaE_kalman=_k.x_hat[2], Ve_kalman=_k.x_hat[3],
                               **{'gamma_diag%d'%i: _k.gamma[i,i] for i in range(4)},
                               dwdt_kalman = _dwdt,
                               Ve=_Ve, gammaV=self.gammaV)

            # update variables
            self.z += dt_step*self.w
            self.z = np.amin((self.z,0.))
            self.w += dt_step*_force/(1+self.a)/self.m
            self.dwdt = _force/(1+self.a)/self.m
            t+=dt_step
        print('... time stepping done')

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
        increment_error: [no dimension]
            coefficient measuring the accuracy on the smallest variation of
            translation motion for the piston (coefficient >= 1)
        vol_error: [m^3]
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
        omega_min: float [rad/s]
            minimum rotation rate
        efficiency: float [<1]
            Piston efficiency, i.e. mechanical work produced over electrical work supplied

        """

        if model.lower() == 'ifremer':
            # default parameters: IFREMER float
            params = {'r': 0.0195/2, 'phi': 0., 'd': 0., 'vol': 0.,
                      'omega': 0., 'lead': 1,
                      'phi_min': 0., 'd_min': 0., 'd_max': 0.09,'vol_min': 0.,
                      'translation_max': 0.12/5600.*225.,
                      'translation_min': 0.12/5600.*10.,
                      'efficiency':.1,
                      'd_increment' : 0.12/5600.,'increment_error' : 10}
        elif model.lower() == 'ensta':
            # default parameters: ENSTA float
            params = {'r': 0.025, 'phi': 0., 'd': 0., 'vol': 0.,
                      'omega': 0., 'lead': 0.00175, 'tick_per_turn': 48,
                      'phi_min': 0., 'd_min': 0., 'd_max': 0.07,
                      'vol_max': 1.718e-4,'vol_min': 0.,
                      'omega_max': 60./48*2.*np.pi, 'omega_min': 0.,
                      'efficiency':.1, 'increment_error' : 1}
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
        self.vol_error = self.d_increment*((self.r)**2)*np.pi*self.increment_error

    def __repr__(self):
        strout='Piston parameters and state: \n'
        strout+='  r     = %.2f mm        - piston radius\n'%(self.r*1.e3)
        strout+='  phi   = %.2f rad       - present angle of rotation\n'%(self.phi)
        strout+='  d     = %.2f mm        - present piston displacement\n'%(self.d*1.e3)
        strout+='  vol   = %.2f cm^3      - present volume addition\n'%(self.vol*1.e6)
        strout+='  lead  = %.2f cm        - screw lead\n'%(self.lead*1.e2)
        #strout+='  tick_per_turn  = %.2f no dimension        - number of notches on the thumbwheel of the piston\n'%(self.tick_per_turn)
        strout+='  d_increment  = %.2f m        - smallest variation of translation motion for the piston\n'%(self.d_increment)
        strout+='  vol_error  = %.2e m^3        - smallest variation of volume possible for the piston\n'%(self.vol_error)
        strout+='  phi_max = %.2f deg     - maximum rotation\n'%(self.phi_max*1.e2)
        strout+='  phi_min = %.2f deg     - minimum rotation\n'%(self.phi_min*1.e2)
        strout+='  d_max = %.2f mm        - maximum piston displacement\n'%(self.d_max*1.e3)
        strout+='  d_min = %.2f mm        - minimum piston displacement\n'%(self.d_min*1.e3)
        strout+='  vol_min = %.2f cm^3    - min volume displaced\n'%(self.vol_min*1.e6)
        strout+='  vol_max = %.2f cm^3    - max volume displaced\n'%(self.vol_max*1.e6)
        strout+='  omega_max = %.2f deg/s - maximum rotation rate\n'%(self.omega_max*180./np.pi)
        strout+='  omega_min = %.2f deg/s - minimum rotation rate\n'%(self.omega_min*180./np.pi)
        strout+='  efficiency = %.2f - mechanical work produced / electrical work supplied\n'%(self.efficiency)
        return strout

#-------------------------- update methods -------------------------------------

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
        self.phi = self.vol2phi(vol)
        self._checkbounds()
        self._bcast_phi()

    def update_d(self,d):
        self.phi = self.d2phi(d)
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

#--------------------- conversion methods --------------------------------------

    def omega2dvdt(self,omega):
        # /np.pi*np.pi has been simplified
        #to compute omega for the ENSTA float:
        #the motor of the piston needs 48 notches to complete a full rotation
        #it can reach until 30 rotations a seconde
        #so omega = 2*pi*30/48 = 3.9 rad/s
        return omega*self.lead/2.*self.r**2

    def dvdt2omega(self,dvdt):
        # /np.pi*np.pi has been simplified
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

    def _checkbounds(self):
        self.phi = np.amin([np.amax([self.phi,self.phi_min]),self.phi_max])


# ----------------------------- utils ------------------------------------------

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
    ax.annotate('',
                xy=(ax.get_xticks()[1], ax.get_ylim()[1]-10),
                xytext=(ax.get_xticks()[2], ax.get_ylim()[1]-10),
                arrowprops=dict(arrowstyle="<->"))
    ax.text(ax.get_xticks()[1]*.5+ax.get_xticks()[2]*.5,
            ax.get_ylim()[1]-10,
            '%.1f g'%(f.V*(ax.get_xticks()[2]-ax.get_xticks()[1])*1e3),
            {'ha': 'center', 'va': 'bottom'})
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
def descent(Tmax, zt, f, waterp, zstart = 0):
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
    #z_target = np.zeros_like(t)
    dzdt_target = -t*afmax/2./f.m
    dzdt_target[np.where(-dzdt_target>wmax)]=-wmax
    z_target = zstart + np.cumsum(dzdt_target*1.)
    z_target[np.where(z_target<zt)] = zt

    # convert to callable function
    return interp1d(t, z_target, kind='linear', fill_value='extrapolate')
