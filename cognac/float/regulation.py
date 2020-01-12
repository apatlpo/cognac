import numpy as np

# useful parameters
g=9.81

#------------------------------- controls --------------------------------------

class control(object):

    def __init__(self, **kwargs):
        for key, item in kwargs.items():
            setattr(self, key, item)
        if not hasattr(self, 'continuous'):
            self.continuous = True
        if 'feedback' in self.mode:
            self.ldb1 = 2/self.tau
            self.ldb2 = 1/self.tau**2
            self._A = g*self.rho_cte/((self.a+1)*self.m)
            self._B = self.c1/(2*self.Lv*(1+self.a))

    def __repr__(self):
        _core_params = ['dt', 'dz_nochattering', 'tau', 'nu', 'delta',
                        'Kp','Ki','Kd', 'continuous']
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

    def get_u_pid(self, z_target, t, z):
        error = z_target(t) - z
        self.integral += error*self.dt
        self.derivative = (error - self.error)/self.dt
        self.error = error
        u = self.Kp*self.error \
            + self.Ki*self.integral \
            + self.Kd*self.derivative
        return u

    def get_u_feedback1(self, z_target, t, z, w, V, gamma, log):
        u = _control_feedback1(self.ldb1, self.nu, self.delta,
                              z, w, z_target(t), V, gamma,
                              self._A, self._B)
        if log:
            self.log.store(time=t, u=sum(u),\
                           **{'u%d'%i: u[i] for i in range(len(u))})
        return sum(u)

    def get_u_feedback2(self, z_target, t, z, w, dwdt, gamma, log):
        u = _control_feedback2(self.ldb1, self.ldb2, self.nu, self.delta,
                              z, w, dwdt, z_target(t), gamma,
                              self._A, self._B)
        if log:
            self.log.store(time=t, u=sum(u),\
                           **{'u%d'%i: u[i] for i in range(len(u))})
        return sum(u)

def _control_feedback1(lbd1, nu, delta, z, dz, z_t, V, gamma,
                      A, B):
    ''' Control feedback of the float position
    Parameters
    ----------

    ldb1: float
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
    #
    e = x1bar - x1
    D = 1 + e**2/delta**2
    #
    y = x0 - nu*np.arctan(e/delta)
    #
    return (1/A)*lbd1*y, \
            (1/A)*nu/delta*x0/D, \
            -(1/A)*2*B*np.abs(x0)*x0, \
            -V, gamma*x0

def _control_feedback2(lbd1, lbd2, nu, delta, z, dz, d2z, z_t, gamma,
                      A, B):
    ''' Control feedback of the float position
    Parameters
    ----------

    ldb1: float
        float control parameter 1 [s^-1]
    ldb2: float
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
    if x0 > 0:
        _sign = 1.
    else:
        _sign = -1.
    return (1/A)*lbd1*dy, (1/A)*lbd2*y, \
            (1/A)*nu/delta*(dx0*D + 2*e*x0**2/delta**2)/(D**2), \
            (1/A)*_sign*2*B*x0*dx0, gamma*x0

#    return (1/A)*(lbd1*dy + lbd2*y\
#               + nu/delta*(dx0*D + 2*e*x0**2/delta**2)/(D**2)\
#               + _sign*2*B*x0*dx0) + gamma*x0


#------------------------------- kalman filter ---------------------------------

class kalman_filter(object):
    ''' Kalman filter for float state estimation
    State vector is:
     [downward velocity, depth, equivalent compressibility, equivalent volume]
    '''

    def __init__(self, x0, **params_in):
        # default parameters
        params = {'dt': 1., 'depth_error': 1e-3, 'verbose': 0}
        # check if mandatory parameters are here, bad form ...
        for key in ['m', 'a','rho_cte','c1','Lv']:
            assert key in params_in
        #
        params.update(params_in)
        # initial state covariance
        if 'gamma' in params:
            params['gamma'] = np.diag(params['gamma'])
        # dynamical noise covariance
        if 'gamma_alpha_scaled' in params:
            params['gamma_alpha'] = params['dt']**2 * np.diag(params['gamma_alpha_scaled'])
        elif 'gamma_alpha' in params:
            params['gamma_alpha'] = np.diag(params['gamma_alpha'])
        # observation noise covariance
        if 'gamma_beta' not in params:
            params['gamma_beta'] = np.diag([params['depth_error']**2])
        # set parameters as attributes
        for key,val in params.items():
            setattr(self,key,val)
        # state vector
        self.x_hat = np.array(x0)
        # coefficients
        self.A_coeff = g*self.rho_cte/((self.a+1)*self.m)
        self.B_coeff = self.c1/(2*self.Lv*(1+self.a))
        # linearized dynamical operator
        self.A = np.eye(4)
        self.A += self.dt * \
                 np.array([[-self.B_coeff*abs(self.x_hat[0]),
                            self.A_coeff*self.x_hat[2],
                            self.A_coeff*self.x_hat[1], -self.A_coeff],
                           [1., 0., 0, 0],
                           [0, 0, 0., 0],
                           [0, 0, 0, 0.]])
        # observation operator
        self.C = np.array([[0, 1, 0, 0.]])

    def __repr__(self):
        strout='Kalman filter: \n'
        strout+='  dt     = %.2f s     - filter time step\n'%(self.dt)
        strout+='  x_hat   = [%.2e,%.2e,%.2e,%.2e] - kalman state \n'%tuple(self.x_hat)
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

    def gen_obs(self, z):
        # build observations
        y_depth = -z \
                 + np.random.normal(loc=0.0, scale=self.depth_error)
        return [y_depth]

    def update_kalman(self, v, z):
        # update state
        self.A[0,0] = 1 - self.dt*self.B_coeff*np.abs(self.x_hat[0])
        self.A[0,1] = 1 + self.dt*self.A_coeff*self.x_hat[2]
        self.A[0,2] = 1 + self.dt*self.A_coeff*self.x_hat[1]
        y = self.gen_obs(z)
        if self.verbose>0:
            print("x kalman", self.x_hat)
        (self.x_hat, self.gamma) = self.kalman(self.x_hat, self.gamma, v, y,
                                               self.A)
        if self.verbose>1:
            print("x0 iteration", self.x_hat)
            print('x_hat', self.x_hat)
            print('z', z)
            print('gamma', self.gamma)

    def kalman(self, x0, gamma0, v, y, A):
        xup, Gup = self.kalman_correc(x0, gamma0, y)
        x1, gamma1 = self.kalman_predict(xup, Gup, v, A)
        return x1, gamma1

    def kalman_predict(self, xup, Gup, v, A):
        gamma1 = (A @ Gup @ A.T)
        gamma1 += self.gamma_alpha
        x1 = xup + self.f(xup, v)*self.dt
        return x1, gamma1

    def kalman_correc(self, x0, gamma0, y):
        C = self.C
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

    def f(self, x, v):
        dx = np.array(x)
        dx[0] = -self.A_coeff*(x[3] - x[2]*x[1] + v) \
                -self.B_coeff*x[0]*np.abs(x[0])
        dx[1] = x[0]
        dx[2] = 0.0
        dx[3] = 0.0
        return dx

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
