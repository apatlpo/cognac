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
        # add tests: if w is None, if f is None ...
        #x1=self.z
        x2 = w
        #x3=self.V+self.v
        #f1=x2
        f2 = f._f(z, ctrl['waterp'], ctrl['Lv'])/f.m
        f3 = ( f.volume(z=z+.5, waterp=ctrl['waterp'])
             - f.volume(z=z-.5, waterp=ctrl['waterp']) )/1. *x2 # dVdz*w
        df1, df2, df3 = f._df(f.z, ctrl['waterp'], ctrl['Lv'])
        df1, df2, df3 = df1/f.m, df2/f.m, df3/f.m
        #
        d3y = ctrl['d3y_ctrl']*control_sliding(z, w, f2, z_t, dz_t, d2z_t, ctrl['tau'])
        u = df1*x2 + df2*f2 + df3*f3 - d3y
        u = -u/df3

    elif ctrl['mode'] == 'pid':
        error = z_t - z
        ctrl['integral'] += error*ctrl['dt_ctrl']
        ctrl['derivative'] = (error - ctrl['error'])/ctrl['dt_ctrl']
        ctrl['error'] = error
        u = ctrl['Kp']*ctrl['error'] + ctrl['Ki']*ctrl['integral'] + ctrl['Kd']*ctrl['derivative']

    elif ctrl['mode'] == 'feedback':
        ctrl['ldb1'] = 2/ctrl['tau'] # /s
        ctrl['ldb2'] = 1/ctrl['tau']**2 # /s^2
        #f2 = f._f(z, ctrl['waterp'], ctrl['L'])/f.m
        u = control_feedback(z, w, dwdt, z_t, ctrl['nu'], ctrl['gammaV'], ctrl['L'], ctrl['c1'],
                             ctrl['m'], ctrl['rho'], ctrl['a'], ctrl['waterp'],
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

class Kalman(object):
    ''' Kalman filter for float state estimation
    '''

    def __init__(self, x0, **params):

        for key,val in params.items():
            setattr(self,key,val)

        self.x_hat = np.array(x0)
        #self.u = 0

        self.A_coeff = g*self.rho/((self.a+1)*self.m)
        self.B_coeff = self.c1/(2*self.L*(1+self.a))

        self.A = self.dt * \
                 np.array([[-self.B_coeff*abs(self.x_hat[0]), self.A_coeff*self.x_hat[2], self.A_coeff*self.x_hat[1], -self.A_coeff],
                           [1., 0., 0, 0],
                           [0, 0, 0., 0],
                           [0, 0, 0, 0.]])
        self.A += np.eye(4)
        self.C = np.array([[0, 1, 0, 0.]])

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
            print("x0 initial", self.x_hat)
        (self.x_hat, self.gamma) = self.kalman(self.x_hat, self.gamma, u, v, y,
                                              self.A)
        if self.verbose>0:
            print("x0 iteration", self.x_hat)
            print('x_hat', self.x_hat)
            print('u', u)
            print('z', z)
            print('gamma', self.gamma)


    def kalman(self,x0,gamma0,u, v, y, A):
        xup,Gup = self.kalman_correc(x0, gamma0, y)
        x1,gamma1=self.kalman_predict(xup, Gup, u, v, A)
        return x1, gamma1

    def kalman_predict(self, xup, Gup, u, v, A):
        gamma1 = (A @ Gup @ A.T)
        gamma1 += self.gamma_alpha
        x1 = xup + self.f(xup, u, v)*self.dt
        return x1, gamma1

    def kalman_correc(self,x0,gamma0,y):
        C = self.C
        if self.verbose>0:
            print(C.shape, self.gamma_beta.shape, gamma0.shape)
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
        if self.verbose>0:
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
        y = x[0] - nu*atan(e/delta)

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
#Functions necesary to estimate parameters for feedback regulation

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
