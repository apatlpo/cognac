import os
from glob import glob
import yaml

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


def find_bags(bag_dir):
    """find log files and sort them properly
    Assumes names of the type: '20200618_m6.log'
    where m6 means the 6th mission
    """
    return sorted(glob(bag_dir + "pd_*"))


def find_logs(bag_dir):
    """find log files and sort them properly
    Assumes names of the type: '20200618_m6.log'
    where m6 means the 6th mission
    """
    files = glob(bag_dir + "*.log")
    k_sorted = [
        int(f.split("/")[-1].split("_")[-1].replace(".log", "")[1:]) for f in files
    ]
    return [f for i, f in sorted(zip(k_sorted, files))]


def load_bag(bag_dir):
    """Load data from seabot ros bags
    Data is in a directory with csv files
    Bags must have been preprocessed with cognac/float/bag.py

    Parameters:
    -----------
    bag_dir: str
        Path to directory where bag processed
    """
    files = glob(os.path.join(bag_dir, "*.csv"))
    bdata = {
        f.split("/")[-1].replace(".csv", ""): pd.read_csv(f, parse_dates=[0])
        for f in files
    }
    # set date as index:
    for d in bdata:
        bdata[d] = (
            bdata[d]
            .rename(columns={"time": "seconds", "Unnamed: 0": "time"})
            .set_index("time")
        )
    return bdata


def load_config(file):
    """Load config from yaml config file"""
    # should add arguments to overwrite config file content
    with open(file, "r") as stream:
        cfg = yaml.load(stream)
    return cfg


def load_config_from_log(file):
    """Load config from ros stdout files created with commands such as:
    roslaunch seabot mission.launch >20200618_m1.log 2>&1
    """
    params = {}
    with open(file) as f:
        for line in f:
            if line[:4] == " * /":
                line = line.strip(" \t\n\r")
                value = line.split(":")[1]
                if "[" in value:
                    value = list(
                        map(
                            lambda v: float(v),
                            value.replace("[", "").replace("]", "").split(","),
                        )
                    )
                else:
                    try:
                        value = float(value)
                    except:
                        if "True" in value:
                            value = True
                        elif "False" in value:
                            value = False
                #
                keys = line.split(":")[0].split("/")[1:]
                _p = params
                for i, k in enumerate(keys):
                    if i == len(keys) - 1:
                        _p[k] = value
                    elif k not in _p:
                        _p[k] = {}
                        _p = _p[k]
                    else:
                        _p = _p[k]
    return params


def compare_config(d1, d2, path=""):
    """Compare two configurations"""
    for k in d1:
        if k not in d2:
            print(path, ":")
            print(k + " as key not in d2", "\n")
        else:
            if type(d1[k]) is dict:
                if path == "":
                    path = k
                else:
                    path = path + "->" + k
                compare_config(d1[k], d2[k], path)
            else:
                if d1[k] != d2[k]:
                    print(path, ":")
                    print(" - ", k, " : ", d1[k])
                    print(" + ", k, " : ", d2[k])


def resample_join(dt, *V):
    """Resample and join on a common time line pandas Series or dataframes

    Parameters:
    -----------
    dt: str
        "rule" in pandas terminology, e.g. '1s'
    *V: pandas.Series, pandas.DataFrame
        Time series to combine
    """
    return pd.concat(
        [d.resample(dt).mean() for d in V],
        join="inner",
        axis=1,
    )


def compute_u(
    f, df, root_regulation, limit_velocity, approach_velocity, decomposed=False
):
    """Compute relaxation feedback command

    Parameters:
    -----------
    f: cognac.float.autonomous_float
    df: pandas.DataFrame
        must contain the following variables:
        velocity, depth, position, offset, chi, chi2, cz
    root_regulation: float
        feedback root parameter (Hz)
    limit_velocity: float
        m/s
    approach_velocity: float
        m
    decomposed: boolean, optional
        returns a decomposed output if True
        Default is False
    """

    set_point = df["set_point"]

    A, B = f.ctrl._A * (1 + f.a), f.ctrl._B
    # print(A) # 741.6 check
    # print(B) # 3.6 check
    s = root_regulation

    x1 = df["velocity"]
    x2 = df["depth"]
    x3 = -df["position"] * f.piston.vol_increment
    x4 = df["offset"]
    x5 = df["chi"]
    x6 = df["chi2"]
    x7 = df["cz"]

    beta = limit_velocity * 2.0 / np.pi
    alpha = approach_velocity
    gamma = beta / alpha

    e = (set_point - x2) / alpha
    y = x1 - gamma * np.arctan(e)
    D = 1 + e**2
    dx1 = -A * (x3 + x4 - (x5 * x2 + x6 * x2**2)) - B * x7 * abs(x1) * x1
    # double dx1 = -A*(x3+x4-(x5*x2+x6*pow(x2,2)))-B*x7*abs(x1)*x1;
    dy = dx1 + gamma * x1 / D

    u = (
        (-2.0 * s * dy / A).rename("u_dy"),
        (s**2 * y / A).rename("u_y"),
        (gamma * (dx1 * D + 2.0 * e * x1**2 / alpha**2) / (D**2) / A).rename(
            "u_gamma"
        ),
        (-2.0 * B * x7 * abs(x1) * dx1 / A).rename("u_B"),
        (x1 * (x5 + 2.0 * x6 * x2)).rename("u_chi"),
    )
    # return (-2.*s*dy+pow(s,2)*y+ gamma*(dx1*D+2.*e*pow(x1,2)/pow(alpha,2))/(pow(D,2))-2.*B*x7*abs(x1)*dx1)/A+x1*(x5+2.*x6*x2);

    if decomposed:
        return {
            "y": (x1.rename("x1"), (-gamma * np.arctan(e)).rename("gamma")),
            "dy": (dx1.rename("dx1"), (-gamma * x1 / D).rename("gamma")),
            "dx1": (
                (-A * x3).rename("v"),
                (-A * x4).rename("offset"),
                (A * (x5 * x2 + x6 * x2**2)).rename("chi"),
                (-B * x7 * abs(x1) * x1).rename("drag"),
            ),
            "u": u,
        }

    else:
        return u


def get_depth_profile(df, threshold=0.02, dz=0.1):
    """select ascent and bin by depth"""
    z = -df["depth"]
    # dt = df.index.to_series().diff().mean().total_seconds()
    dt = df.index.to_series().diff().apply(lambda x: x.total_seconds())
    # compute speed of descent
    dzdt = z.diff() / dt
    #
    df = df[dzdt > threshold]  # ascending values only
    df["z"] = z
    # dzdt.plot(subplots=True)
    #
    z_rounded = np.round(z / dz).rename("bins")
    df = df.groupby(by=z_rounded).mean()
    #
    df = df.drop(columns="depth").reset_index().drop(columns="bins").set_index("z")
    return df


class kalman(object):
    """Kalman filter for float state estimation

        State vector is:
        - x[0] = downward velocity
        - x[1] = depth
        - x[2] = equivalent volume - offset
        - x[3] = equivalent compressibility - chi
        - x[4] = equivalent compressibility2 - chi2
        - x[5] = drag coefficient - cz

    Parameters:
    -----------
    x_init: list, array
    gamma_init: list, array
    gamma_alpha: list, array
    gamma_beta: list, array
    dt: float
        Kalman filter time step
    A_coeff: float
        =g*rho/m
    B_coeff: float
        =0.5*rho*Cf/m
        Cf = M_PI*pow(diam_collerette/2.0, 2);
    sqrt: boolean, optional
        Rescaling of gamma_alpha with sqrt(dt).
        Default is True
    tick_to_volume: float, optional
        Default is 1
    verbose: int, optional
        Turn verbosity on
    """

    def __init__(
        self,
        x_init,
        gamma_init,
        gamma_alpha,
        gamma_beta,
        dt,
        A_coeff,
        B_coeff,
        sqrt=True,
        tick_to_volume=1.0,
        verbose=0,
    ):
        #
        self.verbose = verbose
        self.dt = dt
        self.A_coeff = A_coeff
        self.B_coeff = B_coeff
        self.sqrt = sqrt
        # initial state covariance
        self.gamma_init = np.diag(gamma_init)
        self.gamma = self.gamma_init
        # dynamical noise covariance
        self.gamma_alpha = np.diag(gamma_alpha)
        # observation noise covariance
        if isinstance(gamma_beta, float):
            self.gamma_beta = np.diag([gamma_beta])
        else:
            self.gamma_beta = np.diag(gamma_beta)
        # rescale gamma_alpha to match seabot:
        if sqrt:
            self.gamma_alpha = self.gamma_alpha * np.sqrt(dt)
        self.tick_to_volume = tick_to_volume
        # state vector
        self.names = ["velocity", "depth", "offset", "chi", "chi2", "cz"]
        self.init_x(x_init)
        # linearized dynamical operator
        self.A = self.compute_A(self.x_hat)
        # observation operator
        self.C = np.array([[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]])

    def init_x(self, x_init):
        """Initialize state vector"""
        self.x_hat = np.array(x_init)

    def simulate(self, df):
        """run the kalman filter over observations"""
        out = pd.DataFrame(columns=self.names)
        for i, d in df.iterrows():
            self.step(-d["v"], y=d["depth"])
            out.loc[i] = self.x_hat
        return out

    def gen_obs(self, z):
        """generate synthetic observations based on prescribed covariances"""
        y_depth = -z + np.random.normal(loc=0.0, scale=self.depth_error)
        return [y_depth]

    def compute_A(self, x):
        A = np.eye(6)
        _A = np.zeros((6, 6))
        cA = self.A_coeff
        cB = self.B_coeff
        _A[0, 0] = -2.0 * cB * abs(x[0]) * x[5]
        _A[0, 1] = cA * (x[3] + 2.0 * x[4] * x[1])
        _A[0, 2] = -cA
        _A[0, 3] = x[1] * cA
        _A[0, 4] = x[1] ** 2 * cA
        _A[0, 5] = -cB * abs(x[0]) * x[0]
        _A[1, 0] = 1.0
        A += _A * self.dt
        return A

    def step(self, v, z=None, y=None):
        # update arrays
        self.A = self.compute_A(self.x_hat)
        if y is None:
            y = self.gen_obs(z)
        if self.verbose > 0:
            print("x kalman", self.x_hat)
        (self.x_hat, self.gamma) = self.core(self.x_hat, self.gamma, v, y, self.A)
        if self.verbose > 1:
            print("x0 iteration", self.x_hat)
            print("x_hat", self.x_hat)
            print("z", z)
            print("gamma", self.gamma)

    def core(self, x0, gamma0, v, y, A):
        xup, Gup = self.correct(x0, gamma0, y)
        x1, gamma1 = self.predict(xup, Gup, v, A)
        return x1, gamma1

    def predict(self, xup, Gup, v, A):
        gamma1 = A @ Gup @ A.T
        gamma1 += self.gamma_alpha
        x1 = xup + self.f(xup, v) * self.dt
        return x1, gamma1

    def correct(self, x0, gamma0, y):
        C = self.C
        S = C @ gamma0 @ C.T + self.gamma_beta
        K = gamma0 @ C.T @ np.linalg.inv(S)
        ytilde = np.array(y) - C @ x0
        Gup = (np.eye(len(x0)) - K @ C) @ gamma0
        xup = x0 + K @ ytilde
        #
        self.S = S
        self.K = K
        self.ytilde = ytilde
        #
        if self.verbose > 1:
            print("K", K)
            print("ytilde", ytilde)
        return xup, Gup

    def f(self, x, v):
        dx = np.array(x)
        dx[0] = -self.A_coeff * (
            v + x[2] - (x[3] * x[1] + x[4] * x[1] ** 2)
        ) - self.B_coeff * x[5] * x[0] * np.abs(x[0])
        dx[1] = x[0]
        dx[2] = 0.0
        dx[3] = 0.0
        dx[4] = 0.0
        dx[5] = 0.0
        return dx


def get_kf_parameters(cfg, tick2volume, dt, **kwargs):
    """Compute covariances from seabot configuration
    Kalman filter version 1
    """
    _cfg = cfg["kalman"]
    square = lambda v, scale: [v**2 * scale for v in v]
    sqrt = np.sqrt(dt)
    #
    x = ["velocity", "depth", "offset", "chi", "chi2", "cz"]
    _params = ["gamma_alpha_" + v for v in x]
    gamma_alpha = [float(_cfg[p]) if p not in kwargs else kwargs[p] for p in _params]
    #
    _params = ["gamma_init_" + v for v in x]
    gamma_init = [float(_cfg[p]) if p not in kwargs else kwargs[p] for p in _params]
    #
    for i in [2, 3, 4]:
        gamma_alpha[i] = gamma_alpha[i] * tick2volume
        gamma_init[i] = gamma_init[i] * tick2volume
    #
    gamma_alpha = [v**2 * sqrt for v in gamma_alpha]
    gamma_init = [v**2 for v in gamma_init]
    gamma_beta = float(_cfg["gamma_beta_depth"]) ** 2
    return gamma_alpha, gamma_init, gamma_beta
