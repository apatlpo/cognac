import numpy as np

# useful parameters
g=9.81

from .seabot import load_config, load_config_from_log, get_kf_parameters

#------------------------------- controls --------------------------------------

class control(object):

    def __init__(self, **kwargs):
        for key, item in kwargs.items():
            setattr(self, key, item)
        if not hasattr(self, 'continuous'):
            self.continuous = False
        if 'feedback1' in self.mode:
            self.lbd1 = 1./self.tau
            self.A_coeff = g*self.rho_cte/((self.a+1)*self.m)
            self.B_coeff = self.c1/(2*self.Lv*(1+self.a))
        elif 'feedback2' in self.mode:
            self.lbd1 = 2/self.tau
            self.lbd2 = 1/self.tau**2
            self.A_coeff = g*self.rho_cte/((self.a+1)*self.m)
            self.B_coeff = self.c1/(2*self.Lv*(1+self.a))
        if "pid" in self.mode:
            self.error= 0.
            self.integral = 0.
            self.derivative = 0.

    def __repr__(self):
        _core_params = ['dt', 'dz_nochattering',
                        'tau', 'nu', 'delta',
                        'Kp','Ki','Kd',
                        'continuous',
                        ]
        strout='Control parameters: \n'
        strout+='  mode = %s \n'%(self.mode)
        for p in _core_params:
            if hasattr(self, p):
                strout+='  %s = %.2e \n'%(p, getattr(self, p))
        return strout

    def get_u_sliding(self, z_target, t, z, w, f):
        z_t = z_target(t)
        dz_t = (z_target(t+.05)-z_target(t-.05))/.1
        d2z_t = (z_target(t+.05)-2.*z_target(t)+z_target(t-.05))/.05**2
        #
        #x1=self.z
        x2 = w
        #x3=self.V+self.v
        #
        #f1=x2
        f2 = f.compute_force(z, self.waterp, self.Lv)/f.m
        f3 = ( f.volume(z=z+.5, waterp=self.waterp)
             - f.volume(z=z-.5, waterp=self.waterp) )/1. *x2 # dVdz*w
        df1, df2, df3 = f.compute_dforce(f.z, self.waterp, self.Lv)
        df1, df2, df3 = df1/f.m, df2/f.m, df3/f.m
        # (z, dz, d2z, z_t, dz_t, d2z_t, tau_ctrl)
        # (z, w, f2, z_t, dz_t, d2z_t,self.tau)
        d3y = (self.d3y_ctrl
                * np.sign( d2z_t - f2 + 2.*(dz_t-w)/self.tau
                            + (z_t-z)/self.tau**2 )
              )
        #
        u = df1*x2 + df2*f2 + df3*f3 - d3y
        u = -u/df3
        return u

    def get_u_pid(self, z_target, t, z, log):
        error = z_target(t) - z
        self.integral += error*self.dt
        self.derivative = (error - self.error)/self.dt
        self.error = error
        u = self.Kp*self.error \
            + self.Ki*self.integral \
            + self.Kd*self.derivative
        if log:
            self.log.store(time=t,
                           error=self.error,
                           derivative=self.derivative,
                           integral=self.integral,
                           u=u,
                           )
        return u

    def get_u_feedback1(self, z_target, t, z, w, V,
                        gamma, log
                        ):
        u = _control_feedback1(self.lbd1, self.nu, self.delta,
                              z, w, z_target(t), V, gamma,
                              self.A_coeff, self.B_coeff)
        if log:
            self.log.store(time=t, u=sum(u),\
                           **{'u%d'%i: u[i] for i in range(len(u))})
        return sum(u)

    def get_u_feedback2(self, z_target, t, z, w, dwdt,
                        gamma1, log,
                        gamma2=0., c1=1.
                        ):
        u = _control_feedback2(self.lbd1, self.lbd2,
                               self.nu, self.delta,
                               z, w, dwdt,
                               z_target(t),
                               gamma1,
                               gamma2,
                               c1,
                               self.A_coeff,
                               self.B_coeff
                               )
        if log:
            self.log.store(time=t, u=sum(u),\
                           **{'u%d'%i: u[i] for i in range(len(u))})
        return sum(u)

def _control_feedback1(lbd1,
                       nu, delta,
                       z, dz,
                       z_t,
                       V,
                       gamma,
                       A,
                       B,
                       ):
    ''' Control feedback of the float position
    Parameters
    ----------

    lbd1: float
        float control parameter 1 [s^-1]
    nu: float
        Travel velocity when the float is far from the target position [m.s^-1]
    delta: float
        length scale that defines the zone of influence around the target depth [m]
    z: float
        Position of the float, 0. at the surface and negative downward [m]
    dz: float
        Vertical velocity of the float, negative for downward motions [m.s^-1]
    z_target: float
        Target depth [m]
    V: float
        equivalent volume [m^3]
    gamma: float
        equivalent compressibility [m^2]
    A, B: float
        A = g rho / [(1+a)m], B = c_1/[2(1+a)L]
    '''
    #
    x0 = -dz
    x1 = -z
    x1bar = -z_t
    x2 = gamma
    x3 = V
    #
    e = x1bar - x1
    D = 1 + e**2/delta**2
    #
    y = x0 - nu*np.arctan(e/delta)
    #
    return ((1/A)*lbd1*y,
            (1/A)*nu/delta*x0/D,
            -(1/A)*B*np.abs(x0)*x0,
            -x3,
            x2*x1
           )

def _control_feedback2(lbd1, lbd2,
                       nu, delta,
                       z, dz, d2z,
                       z_t,
                       gamma1,
                       gamma2,
                       c1,
                       A,
                       B,
                       ):
    ''' Control feedback of the float position
    Parameters
    ----------

    lbd1: float
        float control parameter 1 [s^-1]
    lbd2: float
        float control parameter 2 [s^-2]
    nu: float
        Travel velocity when the float is far from the target position [m.s^-1]
    delta: float
        length scale that defines the zone of influence around the target depth [m]
    z: float
        Position of the float, 0. at the surface and negative downward [m]
    dz: float
        Vertical velocity of the float, negative for downward motions [m.s^-1]
    d2z: float
        Vertical acceleration of the float, negative for downward accelerations [m.s^-2]
    z_target: float
        Target depth [m]
    gamma: float
        equivalent compressibility [m^2]
    A, B: float
        A = g rho / [(1+a)m], B = c_1/[2(1+a)L]

    '''
    #
    x0 = -dz
    dx0 = -d2z
    x1 = -z
    x1bar = -z_t
    #
    e = x1bar - x1
    D = 1 + e**2/delta**2
    #
    y = x0 - nu*np.arctan(e/delta)
    dy = dx0 + nu*x0/(delta*D)
    #
    return ((1/A)*lbd1*dy,
            (1/A)*lbd2*y,
            (1/A)*nu/delta*(dx0*D + 2*e*x0**2/delta**2)/(D**2),
            (1/A)*2*B*c1*abs(x0)*dx0,
            gamma1*x0+gamma1*x0**2,
           )


#------------------------------- kalman filter ---------------------------------

class kalman_core(object):
    ''' Kalman filter for float state estimation
    '''

    def __init__(self, **params_in):
        # check if mandatory parameters are here
        _required_params = ['m','a','rho_cte',
                            'x_init', 'gamma_init',
                            ]
        assert all([v in params_in for v in _required_params]), \
                'Kalman filter parameters should include: ' \
                +', '.join(_required_params)
        #
        params = {'dt': 1., 'depth_error': 1e-3, 'verbose': 0}
        params.update(params_in)
        # initial state covariance
        params['gamma'] = np.diag(params['gamma_init'])
        # dynamical noise covariance
        if 'sqrt' in params and params['sqrt']:
            params['sqrt'] = np.sqrt(params['dt'])
        else:
            params['sqrt'] = 1.
        if 'gamma_alpha' in params:
            params['gamma_alpha'] = np.diag(params['gamma_alpha'])
        elif 'gamma_alpha_scaled' in params:
            params['gamma_alpha'] = (params['dt']**2
                                    * np.diag(params['gamma_alpha_scaled'])
                                    )
        # observation noise covariance
        if 'gamma_beta' not in params:
            params['gamma_beta'] = np.diag([params['depth_error']**2])
        # set parameters as attributes
        for key,val in params.items():
            setattr(self,key,val)
        # constant coefficients
        self.A_coeff = g*self.rho_cte/((self.a+1)*self.m)
        # init x and operators
        self.init_x()
        self.init_operators()

    def __repr__(self):
        strout='Kalman filter: \n'
        strout+='  dt     = %.2f s     - filter time step\n'%(self.dt)
        strout+='  x_hat   = [' \
                +', '.join(['{:.2e}'.format(_x) for _x in self.x_hat]) \
                +'] - kalman state \n'
        strout+='  sqrt(diag(gamma)): '
        strout+=np.array2string(np.sqrt(np.diag(self.gamma)),
                                formatter={'float_kind':lambda x: "%.2e" % x})
        strout+='\n  sqrt(gamma_alpha) / dt: '
        strout+=np.array2string(np.sqrt(np.diag(self.gamma_alpha))/self.dt,
                                formatter={'float_kind':lambda x: "%.2e" % x})
        strout+='\n  sqrt(gamma_beta): '
        strout+=np.array2string(np.sqrt(np.diag(self.gamma_beta)),
                                formatter={'float_kind':lambda x: "%.2e" % x})
        return strout

    def init_x(self):
        """ Initialize state vector
        """
        self.x_hat = np.array(self.x_init)
        self.Nx = self.x_hat.size

    def init_operators(self):
        """ Initialize operators: linearized dynamics and observation
        Only depth is observed here
        """
        self.A = np.eye(self.Nx)
        self.C = np.zeros((1,self.Nx))
        self.C[0,1] = 1

    def step(self, v, z):
        # update dynamical operator
        self.update_A()
        # generate synthetic observations
        y = self.gen_obs(z)
        if self.verbose>0:
            print("x kalman", self.x_hat)
        self.x_hat, self.gamma = self.core(self.x_hat, self.gamma,
                                           v, y,
                                           self.A, self.gamma_alpha,
                                           self.C, self.gamma_beta
                                           )

    def core(self, x0, gamma0,
             v, y,
             A, gamma_alpha,
             C, gamma_beta):
        xup, Gup = self.correct(x0, gamma0, y, C, gamma_beta)
        return self.predict(xup, Gup, v, A, gamma_alpha)

    def correct(self, x0, gamma0, y, C, gamma_beta):
        S = C @ gamma0 @ C.T + self.gamma_beta
        K = gamma0 @ C.T @ np.linalg.inv(S)
        ytilde = np.array(y) - C @ x0
        Gup = (np.eye(len(x0))- K @ C) @ gamma0
        xup = x0 + K@ytilde
        #
        self.S = S
        self.K = K
        self.ytilde = ytilde
        #
        if self.verbose>1:
            print('K', K)
            print('ytilde', ytilde)
        return xup, Gup

    def predict(self, xup, Gup, v, A, gamma_alpha):
        x1 = xup + self.f(xup, v)*self.dt
        gamma1 = (A @ Gup @ A.T)
        gamma1 += self.sqrt * gamma_alpha
        return x1, gamma1

    def gen_obs(self, z):
        # build observations
        y_depth = -z \
                 + np.random.normal(loc=0.0, scale=self.depth_error)
        return [y_depth]

    def init_log(self, variables=[]):
        self.log_variables =   ['x_%d'%i for i in range(self.Nx)] \
                             + ['gamma_%d'%i for i in range(self.Nx)] \
                             + ['dwdt', 'dwdt_diff'] \
                             + variables

    def get_log(self, t, v, dwdt):
        _dwdt = -self.f(self.x_hat, v)[0]
        out ={'time': t,
              **{'x_%d'%i: self.x_hat[i] for i in range(self.Nx)},
              **{'gamma_%d'%i: self.gamma[i,i] for i in range(self.Nx)},
              'dwdt': _dwdt,
              'dwdt_diff': dwdt - _dwdt}
        return out

class kalman_v0(kalman_core):
    ''' Kalman filter for float state estimation

    State variables are:
        - downward velocity
        - depth
        - equivalent compressibility
        - equivalent volume
    '''

    def __init__(self, **params_in):
        super().__init__(**params_in)
        _required_params = ['c1', 'Lv', 'a']
        assert all([v in params_in for v in _required_params]), \
                'Kalman filter parameters should include: ' \
                +', '.join(_required_params)
        # state vector
        self.names = ['speed', 'depth', 'gamma_e', 'V_e']
        # useful constants:
        self.B_coeff = self.c1/(2*self.Lv*(1+self.a))
        # initialize log variables
        self.init_log(['z','w', 'V_e', 'gamma_e'])

    def update_A(self):
        x, dt, A, B = self.x_hat, self.dt, self.A_coeff, self.B_coeff
        self.A[0,0] = 1 - dt*B*2.*abs(x[0])
        self.A[0,1] = dt*A*x[2]
        self.A[0,2] = dt*A*x[1]
        self.A[0,3] = -dt*A
        self.A[1,0] = dt

    def f(self, x, v):
        dx = np.array(x)
        dx[0] = -self.A_coeff*(x[3] - x[2]*x[1] + v) \
                -self.B_coeff*x[0]*np.abs(x[0])
        dx[1] = x[0]
        dx[2] = 0.0
        dx[3] = 0.0
        return dx

    def get_log(self, t, v, dwdt):
        _dwdt = -self.f(self.x_hat, v)[0]
        out ={'time': t,
              **{'x_%d'%i: self.x_hat[i] for i in range(4)},
              **{'gamma_%d'%i: self.gamma[i,i] for i in range(4)},
              'z': -self.x_hat[1],
              'w': -self.x_hat[0],
              'V_e':self.x_hat[3],
              'gamma_e': self.x_hat[2],
              'dwdt': _dwdt,
              'dwdt_diff': dwdt - _dwdt}
        return out

class kalman_v1(kalman_core):
    ''' Kalman filter for float state estimation

    State variables are:
        - downward velocity
        - depth
        - equivalent compressibility
        - equivalent compressibility square
        - equivalent volume
        - drag coefficient
    '''

    def __init__(self, **params_in):
        super().__init__(**params_in)
        _required_params = ['c1', 'Lv', 'a']
        assert all([v in params_in for v in _required_params]), \
                'Kalman filter parameters should include: ' \
                +', '.join(_required_params)
        # state vector
        self.names = ['velocity', 'depth', 'offset', 'chi', 'chi2', 'cz']
        # useful constants:
        #Cf = np.pi*(ph['diam_collerette']/2.)**2
        #B = 0.5*ph['rho']*Cf/ph['m'] = 0.5 m pi rc^2 / (pi r^2 L) /m = 0.5 (rc/r)^2 /L
        self.B_coeff = self.c1/(2*self.Lv*(1+self.a))
        #
        # initialize log variables
        self.init_log()

    def to_seabot(self, f):
        """ Return parameters as scaled for seabot
        Could do a from_seabot ...
        """
        tick_to_volume = f.piston.vol_increment
        #
        _sqrt = 1. if self.sqrt else np.sqrt(self.dt)
        _d = np.sqrt(np.diag(self.gamma_alpha/_sqrt))
        print(' gamma_alpha_velocity: {:.2e}'.format(_d[0]))
        print(' gamma_alpha_depth: {:.2e}'.format(_d[1]))
        print(' gamma_alpha_offset: {:.2e}'.format(_d[2]/tick_to_volume))
        print(' gamma_alpha_chi: {:.2e}'.format(_d[3]/tick_to_volume))
        print(' gamma_alpha_chi2: {:.2e}'.format(_d[4]/tick_to_volume))
        print(' gamma_alpha_cz: {}'.format(_d[5]))
        #
        _d = np.sqrt(self.gamma_init)
        print(' gamma_init_velocity: {:.2e}'.format(_d[0]))
        print(' gamma_init_depth: {:.2e}'.format(_d[1]))
        print(' gamma_init_offset: {:.2e}'.format(_d[2]/tick_to_volume))
        print(' gamma_init_chi: {:.2e}'.format(_d[3]/tick_to_volume))
        print(' gamma_init_chi2: {:.2e}'.format(_d[4]/tick_to_volume))
        print(' gamma_init_cz: {:.2e}'.format(_d[5]))
        #
        _d = np.sqrt(self.gamma_beta[0])
        print(' gamma_beta_depth: {:.2e}'.format(_d[0]))
        #
        print(' init_chi: {:.2e}'.format(self.x_init[3]/tick_to_volume))
        print(' init_chi2: {:.2e}'.format(self.x_init[4]/tick_to_volume))

    def compare_with_seabot(self, f, cfg=None, log=None, update=False):
        """ Read seabot config parameters and update kalman filter
        """
        tick_to_volume = f.piston.vol_increment

        assert any([cfg, log]), 'A path to a config or a log file is needed'
        if cfg:
            cfg = load_config(file)
        elif log:
            cfg = load_config_from_log(log)

        # compare frequencies:
        dt_seabot = 1./cfg['kalman']['frequency']
        print('dt: {:3.3e}s vs {:3.3e}s (seabot)'.format(self.dt, dt_seabot) )

        # compare A, B parameters
        ph = cfg['physics']
        A = ph['g']*ph['rho']/ph['m']
        Cf = np.pi*(ph['diam_collerette']/2.)**2
        B = 0.5*ph['rho']*Cf/ph['m']
        print('A_coeff: {:3.3e} vs {:3.3e} (seabot)'.format(self.A_coeff, A) )
        print('B_coeff: {:3.3e} vs {:3.3e} (seabot)'.format(self.B_coeff, B) )

        # compare covariances
        gamma_alpha, gamma_init, gamma_beta = get_kf_parameters(
                                        cfg,
                                        tick_to_volume,
                                        dt_seabot
                                        )
        print('gamma_alpha         : '
                +' ,'.join(['{:2.2e}'.format(g) for g
                                    in np.diag(self.gamma_alpha)]))
        print('gamma_alpha (seabot): '
                +' ,'.join(['{:2.2e}'.format(g) for g in gamma_alpha]))
        print('gamma_init          : '
                +' ,'.join(['{:2.2e}'.format(g) for g in self.gamma_init]))
        print('gamma_init (seabot) : '
                +' ,'.join(['{:2.2e}'.format(g) for g in gamma_init]))
        print('gamma_beta          : {:2.2e}'.format(float(self.gamma_beta)))
        print('gamma_beta (seabot) : {:2.2e}'.format(float(gamma_beta)))

        # update covariances
        if update:
            print('Update kalman filter with seabot values')
            self.gamma_alpha = np.diag(gamma_alpha)
            self.gamma_init = np.diag(gamma_init)
            self.gamma_beta = np.diag([gamma_beta])

        # print in terms of seabot inputs
        print('--- in terms of seabot input parameters (with corresponding dt):')
        _g_alpha, _g_init, _g_beta = _seabot_input_convert(tick_to_volume,
                                                           self.gamma_alpha,
                                                           self.gamma_init,
                                                           self.gamma_beta,
                                                           self.dt)
        _g_alpha_sb, _g_init_sb, _g_beta_sb = _seabot_input_convert(tick_to_volume,
                                                                    gamma_alpha,
                                                                    gamma_init,
                                                                    gamma_beta,
                                                                    dt_seabot)
        print('gamma_alpha         : '
                +' ,'.join(['{:2.2e}'.format(g) for g in _g_alpha]))
        print('gamma_alpha (seabot): '
                +' ,'.join(['{:2.2e}'.format(g) for g in _g_alpha_sb]))
        print('gamma_init          : '
                +' ,'.join(['{:2.2e}'.format(g) for g in _g_init]))
        print('gamma_init (seabot) : '
                +' ,'.join(['{:2.2e}'.format(g) for g in _g_init_sb]))
        print('gamma_beta          : {:2.2e}'.format(float(_g_beta)))
        print('gamma_beta (seabot) : {:2.2e}'.format(float(_g_beta_sb)))

    def update_A(self):
        x, dt, A, B = self.x_hat, self.dt, self.A_coeff, self.B_coeff
        self.A[0,0] = 1 - dt*2.*B*abs(x[0])*x[5]
        self.A[0,1] = dt*A*(x[2]+2.*x[4]*x[1])
        self.A[0,2] = -dt*A
        self.A[0,3] = dt*A*x[1]
        self.A[0,4] = dt*A*x[1]**2
        self.A[0,5] = -B*abs(x[0])*x[0]
        self.A[1,0] = dt

    def f(self, x, v):
        dx = np.zeros_like(x)
        dx[0] = -self.A_coeff*(v + x[2] - (x[3]*x[1]+x[4]*x[1]**2)) \
                -self.B_coeff*x[5]*x[0]*np.abs(x[0])
        dx[1] = x[0]
        dx[2] = 0.0
        dx[3] = 0.0
        dx[4] = 0.0
        dx[5] = 0.0
        return dx

def _gamma_to_1d(gamma):
    if isinstance(gamma, list):
        return np.array(gamma)
    elif isinstance(gamma, np.ndarray):
        if len(gamma.shape)==2:
            return np.diag(gamma)
        else:
            return gamma

def _seabot_input_convert(tick_to_volume,
                          gamma_alpha,
                          gamma_init,
                          gamma_beta,
                          dt,
                          x_init=None,
                         ):
    """ Convert gamma_alpha, gamma_init, gamma_beta in terms of seabot input
    parameters
    Valid for kalman filter version v1
    """
    sqrt = np.sqrt(dt)
    _gamma_alpha = np.sqrt(_gamma_to_1d(gamma_alpha)/sqrt)
    _gamma_init = np.sqrt(_gamma_to_1d(gamma_init))
    for i in [2,3,4]:
        _gamma_alpha[i] = _gamma_alpha[i]/tick_to_volume
        _gamma_init[i] = _gamma_init[i]/tick_to_volume
    _gamma_beta = np.sqrt(gamma_beta)
    return _gamma_alpha, _gamma_init, _gamma_beta

class kalman_profile(kalman_core):
    ''' Kalman filter for float state estimation

    State variables are:
        - downward velocity
        - depth
        - drag coefficient
        - equivalent compressibility profile
    '''

    def __init__(self, **params_in):
        super().__init__(**params_in)
        _required_params = ['z', 'Lv', 'a']
        assert all([v in params_in for v in _required_params]), \
                'Kalmanf filter parameters should include: ' \
                +', '.join(_required_params)
        # state vector
        self.names = ['speed', 'depth', 'gamma_e', 'V_e']
        self.x_hat = np.array([-self.dzdt, -self.z, self.gamma_e0, self.V_e0])
        #
        # linearized dynamical operator
        self.A = np.eye(4)
        self.B_coeff = self.c1/(2*self.Lv*(1+self.a))
        self.A += self.dt * \
                 np.array([[-self.B_coeff*abs(self.x_hat[0]),
                            self.A_coeff*self.x_hat[2],
                            self.A_coeff*self.x_hat[1], -self.A_coeff],
                           [1., 0., 0, 0],
                           [0, 0, 0., 0],
                           [0, 0, 0, 0.]])
        # observation operator
        self.C = np.array([[0, 1, 0, 0.]])
        # log variables
        self.log_variables = ['z','w', 'V_e', 'gamma_e', \
                              'dwdt', 'dwdt_diff'] + \
                              ['gamma_diag%d'%i for i in range(4)]

    def update_A(self):
        self.A[0,0] = 1 - self.dt*self.B_coeff*np.abs(self.x_hat[0])
        self.A[0,1] = 1 + self.dt*self.A_coeff*self.x_hat[2]
        self.A[0,2] = 1 + self.dt*self.A_coeff*self.x_hat[1]

    def f(self, x, v):
        dx = np.array(x)
        dx[0] = -self.A_coeff*(x[3] - x[2]*x[1] + v) \
                -self.B_coeff*x[0]*np.abs(x[0])
        dx[1] = x[0]
        dx[2] = 0.0
        dx[3] = 0.0
        return dx

    def get_log(self, t, v, dwdt):
        _dwdt = -self.f(self.x_hat, v)[0]
        out ={'time': t,
              **{'x_%d'%i: self.x_hat[i] for i in range(4)},
              **{'gamma_%d'%i: self.gamma[i,i] for i in range(4)},
              'z': -self.x_hat[1],
              'w': -self.x_hat[0],
              'V_e':self.x_hat[3],
              'gamma_e': self.x_hat[2],
              'dwdt': _dwdt,
              'dwdt_diff': dwdt - _dwdt}
        return out

#----------------------------- feedback parameters -----------------------------
# Functions to guide the estimation of feedback parameters

def omega2dvdt(omega=12.4*2.*np.pi/60., lead=0.0175, r_piston=0.025):
    '''
    Function computing the piston flow u
    parameters:
        omega: float [rad/s]
            current rotation rate, omega=dphi/dt
            for ENSTA float, omega_max = 124.*2.*np.pi/60.,
            omega_min = 12.4*2.*np.pi/60.
        lead: float [m]
            screw lead (i.e. displacement after one screw revolution)
            d = phi/2/pi x lead
        r_piston: float [m]
            piston radius
    '''
    return omega*lead/2.*r_piston**2

def zf(t, params):

    '''
    Function computing the float position depending on time and float parameters
    for initial conditions zf = 0 and vf = 0 at the beginning
    '''
    rho_w = 997 #kg.m^3
    g = 9.81 #m.s^-2
    return (params['u']*g*rho_w*t**3) /6 /params['m'] /(1+params['a'])

def vf(t, params):
    '''
    Function computing the float speed depending on time and float parameters
    for initial conditions zf = 0 and vf = 0 at the beginning
    '''
    rho_w = 997 #kg.m^3
    g = 9.81 #m.s^-2
    return (params['u']*g*rho_w*t**2) / (2*params['m']*(1+params['a']))

def tv(v, params):
    '''
    Function computing the time necessary for the float to reach the speed v
    '''
    rho_w = 997 #kg.m^3
    g = 9.81 #m.s^-2
    return np.sqrt(2*v*params['m']*(1+params['a'])/(g*rho_w*params['u']))

def zv(v, params):
    '''
    Function computing the distance necessary for the float to reach the speed v
    '''
    return zf(tv(v,params),params)
