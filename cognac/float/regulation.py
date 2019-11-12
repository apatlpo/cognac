import numpy as np

# useful parameters
g=9.81

#------------------------------- controls --------------------------------------

def control(z, z_target, ctrl, t=None, w=None, f=None, dwdt=None, v=None):
    ''' Implements the control of the float position
    '''
    z_t = z_target(t)
    dz_t = (z_target(t+.05)-z_target(t-.05))/.1
    d2z_t = (z_target(t+.05)-2.*z_target(t)+z_target(t-.05))/.05**2
    #
    if ctrl['mode'] == 'sliding':
        #x1=self.z
        x2 = w
        #x3=self.V+self.v
        #f1=x2
        f2 = f.compute_force(z, ctrl['waterp'], ctrl['Lv'])/f.m
        f3 = ( f.volume(z=z+.5, waterp=ctrl['waterp'])
             - f.volume(z=z-.5, waterp=ctrl['waterp']) )/1. *x2 # dVdz*w
        df1, df2, df3 = f.compute_dforce(f.z, ctrl['waterp'], ctrl['Lv'])
        df1, df2, df3 = df1/f.m, df2/f.m, df3/f.m
        #
        d3y = ctrl['d3y_ctrl']*control_sliding(z, w, f2, z_t, dz_t, d2z_t, ctrl['tau'])
        u = df1*x2 + df2*f2 + df3*f3 - d3y
        u = -u/df3

    elif ctrl['mode'] == 'pid':
        error = z_t - z
        ctrl['integral'] += error*ctrl['dt']
        ctrl['derivative'] = (error - ctrl['error'])/ctrl['dt']
        ctrl['error'] = error
        u = ctrl['Kp']*ctrl['error'] \
            + ctrl['Ki']*ctrl['integral'] \
            + ctrl['Kd']*ctrl['derivative']

    elif ctrl['mode'] == 'feedback':
        ctrl['ldb1'] = 2/ctrl['tau'] # /s
        ctrl['ldb2'] = 1/ctrl['tau']**2 # /s^2
        #f2 = f.compute_force(z, ctrl['waterp'], ctrl['L'])/f.m
        u = control_feedback(z, w, dwdt, z_t, ctrl['nu'], ctrl['gammaV'], ctrl['L'], ctrl['c1'],
                             ctrl['m'], ctrl['rho_cte'], ctrl['a'], ctrl['waterp'],
                             ctrl['ldb1'], ctrl['ldb2'], ctrl['delta'])

    elif ctrl['mode'] == 'kalman_feedback':
        u = control_kalman_feedback(z_t, v, ctrl)

    else:
        print('!! mode '+ctrl['mode']+' is not implemented, exiting ...')
        sys.exit()
    return u

def control_sliding(z, dz, d2z, z_t, dz_t, d2z_t, tau_ctrl):
    ''' Several inputs are required:
    (z_target,w_target,dwdt_target) - describes the trajectory
    tau_ctrl  - a time scale of control
    '''
    return np.sign( d2z_t - d2z + 2.*(dz_t-dz)/tau_ctrl + (z_t-z)/tau_ctrl**2 )

def control_feedback(z, dz, d2z, z_t, nu, gammaV, L, c1, m, rho, a, waterp,
                     lbd1, lbd2, delta):

    ''' Control feedback of the float position
    Parameters
    ----------

    z: float
        Position of the float, 0. at the surface and negative downward [m]
    dz: float
        Vertical velocity of the float, negative for downward motions [m.s^-1]
    d2z: float
        Vertical acceleration of the float, negative for downward accelerations [m.s^-2]
    z_target: float
        Target depth [m]
    nu: float
        Travel velocity when the float is far from the target position [m.s^-1]
    gammaV: float
        Float mechanical compressibility x float volume [m^3/dbar]
    L: float
        Float length [m]
    c1: float
        Float drag parameter
    m: float
        Float mass [kg]
    rho: float
        Float constant density [kg.m^-3]
    a: float
        Float added mass [no dimension]
    waterp: water profile object
            Contains information about the water profile
    ldb1: float
        float control parameter 1 [s^-1]
    ldb2: float
        float control parameter 2 [s^-2]
    delta: float
        length scale that defines the zone of influence around the target depth [m]
    '''
    A = g*rho/((a+1)*m)
    B = c1/(2*L*(1+a))
    x1 = -dz
    dx1 = -d2z
    x2 = -z
    x2bar = -z_t
    e = x2bar - x2
    D = 1 + (e**2)/(delta**2)
    y = x1 - nu*np.arctan(e/delta)
    dy = dx1 + nu*x1/(delta*D)

    if dz < 0:
        return (1/A)*(lbd1*dy + lbd2*y\
               + nu/delta*(dx1*D + 2*e*x1**2/delta**2)/(D**2)\
               + 2*B*x1*dx1) + gammaV*x1
    else: #dz >= 0 not differentiable at value 0 : critical value
        return (1/A)*(lbd1*dy + lbd2*y\
               + nu/delta*(dx1*D + 2*e*x1**2/delta**2)/(D**2)\
               - 2*B*x1*dx1) + gammaV*x1

def control_kalman_feedback(depth_target, v, ctrl):
    ''' Control feedback of the float position with kalman filter
    Parameters
    ----------
    z: float
        Position of the float, 0. at the surface and negative downward [m]
    '''
    kalman = ctrl['kalman']
    x_control = kalman.x_hat
    return kalman.control(x_control, v, depth_target, ctrl)


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
        for key in ['m', 'a','rho_cte','c1','Lv','gammaV']:
            assert key in params_in
        #
        params.update(params_in)
        # derived parameters
        _dt = params['dt']
        #if 'vel_rms' not in params:
        #    params['vel_rms'] = params['depth_rms']/_dt
        #
        # initial state covariance
        if 'gamma' in params:
            params['gamma'] = np.diag(params['gamma'])
        #else:
        #    # old and should probably not be used
        #    params['gamma'] = np.diag([params['vel_rms']**2,
        #                               params['depth_rms']**2,
        #                               params['gamma_alpha_gammaE']**2,
        #                               params['piston_volume_error']**2])
        # dynamical noise covariance
        if 'gamma_alpha_scaled' in params:
            params['gamma_alpha'] = _dt**2 * np.diag(params['gamma_alpha_scaled'])
        elif 'gamma_alpha' in params:
            params['gamma_alpha'] = np.diag(params['gamma_alpha'])
        #else:
        #    # old and should probably not be used
        #    params['gamma_alpha'] = np.diag([params['vel_rms']**2,
        #                                params['depth_rms']**2,
        #                                params['gamma_alpha_gammaE']**2,
        #                                params['piston_volume_error']**2])
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
        strout='Kalman filter object: \n'
        strout+='  dt     = %.2f s     - filter time step\n'%(self.dt)
        strout+='  x_hat   = [%.2e,%.2e,%.2e,%.2e] - kalman state \n'%tuple(self.x_hat)
        strout+=np.array2string(self.gamma, formatter={'float_kind':lambda x: "%.2e" % x})
        return strout

    def gen_obs(self, z, scale = 1.0):
        # build observations
        y_depth = -z + np.random.normal(loc=0.0,
                    scale=np.sqrt(self.gamma_beta[0,0]))
        return [y_depth]

    def update_kalman(self, u, v, z):
        # update state
        self.A[0,0] = 1 - self.dt*self.B_coeff*np.abs(self.x_hat[0])
        self.A[0,1] = 1 + self.dt*self.A_coeff*self.x_hat[2]
        self.A[0,2] = 1 + self.dt*self.A_coeff*self.x_hat[1]
        y = self.gen_obs(z)
        if self.verbose>0:
            print("x kalman", self.x_hat)
        (self.x_hat, self.gamma) = self.kalman(self.x_hat, self.gamma, u, v, y,
                                              self.A)
        if self.verbose>1:
            print("x0 iteration", self.x_hat)
            print('x_hat', self.x_hat)
            print('u', u)
            print('z', z)
            print('gamma', self.gamma)

    def kalman(self,x0,gamma0,u, v, y, A):
        xup, Gup = self.kalman_correc(x0, gamma0, y)
        x1, gamma1 = self.kalman_predict(xup, Gup, u, v, A)
        return x1, gamma1

    def kalman_predict(self, xup, Gup, u, v, A):
        gamma1 = (A @ Gup @ A.T)
        gamma1 += self.gamma_alpha
        x1 = xup + self.f(xup, u, v)*self.dt
        return x1, gamma1

    def kalman_correc(self,x0,gamma0,y):
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

    def f(self,x, u, v):
        dx = np.array(x)
        dx[0] = -self.A_coeff*(x[3] - x[2]*x[1] + v) \
                -self.B_coeff*x[0]*np.abs(x[0])
        dx[1] = x[0]
        dx[2] = 0.0
        dx[3] = 0.0
        return dx

    def control(self, x, v, depth_target, ctrl):
        l1 = 2/ctrl['tau'] # /s
        l2 = 1/ctrl['tau']**2 # /s^2
        nu = ctrl['nu'] # Set the limit speed : 3cm/s # m.s^-1 assesed by simulation
        delta = ctrl['delta'] #length scale that defines the zone of influence around the target depth, assesed by simulation

        e = -depth_target - x[1]
        y = x[0] - nu*np.arctan(e/delta)

        dx1 = -self.A_coeff*(x[3] - x[2]*x[1] + v) \
              -self.B_coeff*x[0]*np.abs(x[0])

        D = 1. + (e/delta)**2
        dy = dx1 + nu*x[0]/(delta*D)

        if x[0] > 0:
            return (1/self.A_coeff)*(l1*dy + l2*y \
                    + nu/delta*(dx1*D + 2*e*x[0]**2/delta**2)/(D**2) \
                    + 2*self.B_coeff*x[0]*dx1) + x[2]*x[0]
        else:
            return (1/self.A_coeff)*(l1*dy + l2*y \
                    + nu/delta*(dx1*D + 2*e*x[0]**2/delta**2)/(D**2) \
                    - 2*self.B_coeff*x[0]*dx1) + x[2]*x[0]


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
