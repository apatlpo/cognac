import os, sys
from glob import glob
import yaml
from pprint import pformat

import numpy as np
import xarray as xr

from scipy.interpolate import interp1d
from scipy.optimize import fsolve, brentq

import matplotlib.pyplot as plt

import gsw

from .log import *
from .regulation import *

# useful parameters
g = 9.81  # m/s^2
watth = 1e-6 / 3.6
dbar2Pa = 1.0e4
cm3 = 1e-6  # cm^3 to m^3 conversion
rho0 = 1030

# -------------------------- float object with dynamics -------------------------


class autonomous_float:
    def __init__(self, model="ifremer", **kwargs):
        """Autonomous float object, the float is assumed to be cylinder-shaped

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

        """
        #
        params = dict()
        if model is None:
            model = "custom"
        self.model = model
        #
        if model.lower() == "ifremer":
            params = {
                "r": 0.07,
                "L": 0.8278,
                "gamma": 3.78039e-06,
                "alpha": 6.98e-5,
                "temp0": 20.0,
                "a": 1.0,
                "drag_model": "quadratic",
                "c1": 1.0,
            }
            params["m"] = 11.630
            # avec 16 piles : 11.630, avec 32 piles : 13.315
            # compressibility from Paul Troadec's report (p39)
        elif model.lower() == "seabot_v0":
            params = {
                "r": 0.06,
                "L": 0.5,
                "gamma": 9.30e-5,
                "alpha": 0.0,
                "temp0": 0.0,
                "a": 1.0,
                "drag_model": "quadratic",
                "c1": 1.0,
            }
            params["m"] = 9.045
        elif model.lower() == "seabot":
            m = 12 # kg
            params = {
                "r": 0.06,
                "L": 1.,
                "gamma": 3.5e-06,
                "alpha": 1.e-4,
                "temp0": 0.0,
                "a": 1.0,
                "drag_model": "polynomial",
                "drag_order": 5,
                "drag_c0": .4588/m,
                "drag_c1": .0388/m,
                "drag_c2": 29.224/m,
                "drag_c3": -1.0928/m,
                "drag_c4": -119.9/m,
            }
            params["m"] = m
        elif model.lower() == "minion":
            params = {
                "r": 0.09 / 2.0,
                "L": 0.4,
                "gamma": 1e-1 * 1e5 / (2.85 * 1e9),
                "alpha": 12e-6,
                "temp0": 0.0,
                "a": 1.0,
                "drag_model": "quadratic",
                "drag_c1": 1.0,
            }
        #
        params.update(kwargs)
        # compute the volume of the float:
        if "V" not in params:
            if "m" in params:
                # derive volume from mass rather that from geometrical dimensions
                V = params["m"]/rho0
            else:
                V = np.pi * params["r"]**2 * params["L"]
            params["V"] = V
        if "m" not in params:
            m = params["V"] * rho0
            print(
                "infer mass from volume with rho=1030 kg/m3, m = {:.3f} kg".format(
                    m
                )
            )
            params["m"] = m

        # auxiliary parameters
        params["rho_cte"] = params["m"] / params["V"]  # kg.m^-3
        params["gammaV"] = params["gamma"] * params["V"]  # m^2
        # fill in drag parameters if need be:
        if params["drag_model"] == "polynomial":
            for i in range(params["drag_order"]):
                if f"drag_c{i}" not in params:
                    params[f"drag_c{i}"] = 0.

        self.params = params
        self._set_params_as_attrs() # will have to be turned off eventually

    def __repr__(self):
        strout = "Float parameters: \n"
        strout += "  L     = %.2f m      - float length\n" % (self["L"])
        strout += "  r     = %.2f m      - float radius\n" % (self["r"])
        strout += "  m     = %.2f kg     - float mass\n" % (self["m"])
        strout += "  V     = %.2e cm^3   - float volume\n" % (self["V"] / cm3)
        strout += "  rho_cte = m/V = %.2e kg.cm^3   - float baseline density\n" % (
            self["rho_cte"] / cm3
        )
        strout += "  gamma = %.2e /dbar  - mechanical compressibility\n" % (self["gamma"])
        strout += "  gamma x V = %.2e cm^3/dbar  - normalized compressibility\n" % (
            self["gamma"] * self["V"] / cm3
        )
        strout += "  alpha = %.2e /degC  - thermal compressibility\n" % (self["alpha"])
        strout += (
            "  alpha x V = %.2e cm^3/degC  - normalized thermal compressibility\n"
            % (self["alpha"] * self["V"] / cm3)
        )
        strout += "  temp0 = %.2e  degC  - reference temperature\n" % (self["temp0"])
        strout += "  a = %.2e  (no dimension)  - float added mass\n" % (self["a"])
        strout += "  drag_model = %s  0\n" % (self["drag_model"])
        if self["drag_model"]=="quadratic":
            strout += "  drag_c1 = %.2e  (no dimension)  - float drag parameter\n" % (self["drag_c1"])
        elif self["drag_model"]=="polynomial":
            for i in range(self["drag_order"]):
                _c = self[f"drag_c{i}"]
                strout += f"  drag_c{i} = {_c:.2e}  (??)  - float drag parameter\n"
        elif self["drag_model"]=="internal_wave":
            strout += "  drag_c0 = %.2e  (??)  - float drag parameter\n" % (self["drag_c0"])
        if hasattr(self, "piston"):
            strout += str(self.piston)
        return strout
    
    def __getitem__(self, key):
        """ to be implemented to access parameters """
        return self.params[key]
    
    def _set_params_as_attrs(self):
        for key, val in self.params.items():
            setattr(self, key, val)

    def rho(
        self,
        p=None,
        temp=None,
        v=None,
        z=None,
        waterp=None,
        m=None,
        air=None,
    ):
        """Returns float density i.e. mass over volume
        Includes the piston volume

        Parameters
        ----------
        p: float, optional
            pressure in dbar
        temp: float, optional
            temperature in degC
        v: float, optional
            piston volume in m^3
        z: float, optional
            depth in meters, negative downward
        waterp: cognac.water.waterp, optional
            water profile
        m: float
            float mass in kg
        air: tuple
            (v_air, t_air) volume of air (m^3) and temperature (degC) at the surface
        """
        V = self["V"]
        #
        if v is None:
            if hasattr(self, "v"):
                v = self.v
            else:
                v = 0.0
        V += v
        #
        if m is None:
            m = self["m"]
        #
        if waterp is not None:
            water = waterp.at(z=z)
        if p is None:
            p = water["pressure"]
        if temp is None:
            # assumes thermal equilibrium
            temp = water["temperature"]
        #
        if air is not None:
            p_surf, v_surf, temp_surf = 10, air[0], air[1]
            ## taking temperature into account
            # p v_air / temp = p(surface) v_air(surface) / temp(surface)
            #   = n R = m/M R
            # p(surface) = 10 dbar
            # temp_abs_zero = 273.15 # degC
            # R = 8.31446262 # m2 kg s-2 K-1 mol-1, gaz constant
            # M_air = 28.97e-3 # kg/mol, molar mass
            # m_air = dbar2Pa*p_surf * v_surf/(temp_abs_zero+temp_surf)  *M_air/R
            # v_air = v_surf * p_surf/(p_surf+p) \
            #    * (temp+temp_abs_zero)/(temp_surf+temp_abs_zero)
            ## ignoring temperature
            rho_air = 1.225  # kg/m^3   =  M_air/R
            m_air = v_surf * rho_air
            v_air = v_surf * p_surf / (p_surf + p)
            m += m_air  # this must ridiculously small
            self._m_air = m_air
            V += v_air
        #
        V += self["V"] * (-self["gamma"] * p + self["alpha"] * (temp - self["temp0"]))
        self._volume = V  # for log purposes
        return m / V

    def volume(self, **kwargs):
        """Returns float volume (V+v)"""
        return self["m"] / self.rho(**kwargs)

    def volume4equilibrium(self, p_eq, temp_eq, rho_eq):
        """Find volume that needs to be added in order to be at equilibrium
        prescribed pressure, temperature, density
        """
        _f = lambda x: rho_eq - self.rho(p=p_eq, temp=temp_eq, v=x)
        v = fsolve(_f, 0.0)[0]
        return v

    def z4equilibrium(self, waterp, z0=-10):
        """Find depth where float is at equilibrium"""
        _f = lambda z: waterp.at(z=z)["rho"] - self.rho(z=z, waterp=waterp)
        z = fsolve(_f, z0)[0]
        return z

    def adjust_m(self, p_eq, temp_eq, rho_eq, piston=True, offset=0.0):
        """Find mass that needs to be added in order to be at equilibrium
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

        """

        v = 0.
        if piston and hasattr(self, "piston"):
            v = self.piston.vol
        def _f(m):
            return rho_eq - self.rho(p=p_eq, temp=temp_eq, m=m, v=v)

        m0 = self["m"]
        m = fsolve(_f, self["m"])[0]
        m += offset * 1e-3
        self.params["m"] = m
        self.params["rho_cte"] = m / self["V"]  # kg.m^-3
        self._set_params_as_attrs()
        print(
            "%.1f g " % ((m - m0) * 1.0e3)
            + " were added to the float in order to be at equilibrium"
            + " at %.0f dbar \n" % (p_eq)
        )

    def get_gamma(self, v_air=None, p=None):
        """compute the float mechanical compressibility, potentially including the effect of trapped air

                Parameters
                ----------
        v       v_air: float, np.array, optional
                    Volume of air at the surface in m^3
                p: float, optional
                    Pressure in dbar where the compressibility has to be computed
                    Required if v_air is passed

                Returns
                -------
                gamma: float, np.array
                    Float compressibility in 1/dbar
        """
        gamma = self["gamma"]
        if v_air is not None:
            assert p is not None, "Need to pass pressure if v_air is provided"
            p_atmo = 10  # dbar
            gamma += v_air / self["V"] * p_atmo / (p + p_atmo) ** 2
        return gamma

    def init_piston(self, **kwargs):
        """wrapper around piston initializer"""
        self.piston = piston(self.model, **kwargs)
        self.v = self.piston.vol

    def piston_update_vol(self, vol=None):
        """wrapper around piston volume update"""
        if vol is None:
            vol = self.piston.vol
        else:
            self.piston.update_vol(vol)
        #
        self.v = vol

    def set_piston4equilibrium(self, p_eq, temp_eq, rho_eq):
        """Adjust piston to be at equilibrium at a given pressure,
        temperature and density
        """
        self.v = self.piston.vol
        #
        def _f(vol):
            self.piston_update_vol(vol)
            return rho_eq - self.rho(p=p_eq, temp=temp_eq)

        vol = brentq(_f, self.piston.vol_min, self.piston.vol_max)
        self.piston_update_vol(vol)
        print("Piston reset for equilibrium : vol=%.1e cm^3  " % (vol / cm3))
        return vol

    def compute_force(self, z, w, water, v=None, sum=True, **kwargs):
        """Compute the vertical force exterted on the float
        We assume thermal equilibrium
        Drag is quadratic

        Parameters
        ----------

        z: float
            Depth of the Float
        water: water profile object, dict
            Water profile used to compute the upward buoyancy force and water
            vertical flow
        Lv: float
            Length of the float used for drag
        v: float
            Piston volume
        w: float
            float vertical velocity (positive upward)

        """
        # gravity
        f_b = -self["m"] * g
        # upward buoyancy
        if not isinstance(water, dict):
            water = water.at(z=z)
        p, temp_w, rho_w = water["pressure"], water["temperature"], water["rho"]
        #
        rho_f = self.rho(p=p, temp=temp_w, v=v)
        f_b += self["m"] * rho_w / rho_f * g
        # drag
        f_d = self.get_drag_force(w, water["detadt"], **kwargs)
        # water accelerative terms, ignoring nonlinear accelerative terms including w_r dw/dz
        f_w = self["m"] * (1 + self.a) * water["d2etadt2"]
        #
        f_total = f_b + f_d + f_w
        if sum:
            return f_total
        else:
            return f_total, f_b, f_d, f_w

    def get_cd(self, w, w_water, Lv=None, N2=None):
        """return drag coefficient such that the drag force is -cd*(w-w_water)"""
        model = self["drag_model"]
        if model=="quadratic":
            assert Lv is not None, "Lv must be provided with quadratic drag"
            cd = 0.5 * 0.82 * self["drag_c1"] / Lv * np.abs(w - w_water)
            # long cylinder: cd = 0.82
            # https://en.wikipedia.org/wiki/Drag_coefficient
            # we need c1 = 1 and Lv = L
        elif model=="polynomial":
            u = w - w_water
            cd = 0
            for i in range(self["drag_order"]):
                cd += self[f"drag_c{i}"]* u**i
        elif model=="internal_wave":
            assert N2 is not None, "pass N2 if c0 is not 0 (iwave drag)"
            if N2 > 0:
                N = np.sqrt(N2)
            else:
                N = 0.0
            cd += N * self["drag_c0"]
        return cd

    def get_drag_force(self, w, w_water, *args, **kwargs):
        cd = self.get_cd(w, w_water, *args, **kwargs)
        f_d = -self["m"] * cd * (w - w_water)
        # m = rho x V = rho x A x L
        # cd = c1/Lv x U
        # F = rho x A x L x c1/Lv x U^2
        #   = c1 x m /Lv x U^2
        #   ~ c1 x rho x A x U^2
        # if the drag force is given in Newtons, for instance: f(u) = C x u^3
        # then one need to enforce: c1 x m /Lv x  = 
        return f_d

    def compute_dforce(self, z, w, waterp, **kwargs):
        """Compute gradients of the vertical force exterted on the float"""
        df1 = (
            self.compute_force(z + 5.0e-2, w, waterp, v=self.v, **kwargs)
            - self.compute_force(z - 5.0e-2, 0.0, waterp, v=self.v, **kwargs)
        ) / 1.0e-1
        df2 = (
            self.compute_force(z, w + 5.0e-3, waterp, v=self.v, **kwargs)
            - self.compute_force(z, w - 5.0e-3, waterp, v=self.v, **kwargs)
        ) / 1.0e-2
        df3 = (
            self.compute_force(z, w, waterp, v=self.v + 5.0e-5, **kwargs)
            - self.compute_force(z, w, waterp, v=self.v - 5.0e-5, **kwargs)
        ) / 1.0e-4
        return df1, df2, df3

    def compute_bounds(self, waterp, zmin, zmax=0.0, **kwargs):
        """Compute approximate bounds on velocity and acceleration"""
        z = zmax
        if hasattr(self, "piston"):
            v = self.piston.vol_max
        else:
            v = None
        fmax = self.compute_force(z, 0.0, waterp, v=v, **kwargs)  # upward force
        #
        z = zmin
        if hasattr(self, "piston"):
            v = self.piston.vol_min
        else:
            v = None
        fmin = self.compute_force(z, 0.0, waterp, v=v, **kwargs)  # downward force
        #
        afmax = np.amax((np.abs(fmax), np.abs(fmin)))
        wmax = np.sqrt(afmax * self["m"] * 2 * self.L / self.c1)
        print(
            "Acceleration and velocity bounds (zmin=%.0fm,zmax=%.0fm):" % (zmin, zmax)
        )
        print(
            "fmax/m=%.1e m^2/s, fmin/m= %.1e m^2/s, wmax= %.1f cm/s"
            % (fmax / self["m"], fmin / self["m"], wmax * 100.0)
        )
        print("For accelerations, equivalent speed reached in 1min:")
        print(
            "  fmax %.1e cm/s, fmin/m= %.1e cm/s"
            % (fmax / self["m"] * 60.0 * 100.0, fmin / self["m"] * 60.0 * 100.0)
        )
        #
        # if hasattr(self,'piston'):
        #    dv = self.piston.dv
        #    p, tempw = waterp.get_p(z), waterp.get_temp(z)
        #    rhow = self.rho(p=p,temp=tempw,v=self.piston.vol_min)
        #    rhof = self.rho(p=p,temp=tempw,v=self.piston.vol_min+dv)
        #    df = self["m"]*(-1.+rhow/rhof)*g
        #    print('Acceleration after an elementary piston displacement: %.1e m^2/s' %(df[0]/self["m"]))
        #    print('  corresponding speed and displacement after 1 min: %.1e m/s, %.1e m \n' \
        #          %(df[0]/self["m"]*60,df[0]/self["m"]*60**2/2.))
        return fmax, fmin, afmax, wmax

    def get_terminal_velocity(self, dm):
        """ from a mass offset compute the terminal velocity that balances drag

        Parameters:
        -----------
        dm: float
            mass offset [kg]
        """        
        def _f(w):
            cd = self.get_cd(w, 0.)
            return g * dm + self["m"] * cd * w 
        if isinstance(dm, float) or isinstance(dm, int):
            w0 = 0.
            return fsolve(_f, w0)[0]
        else:
            w0 = np.zeros_like(dm)
            wt = fsolve(_f, w0)
            return wt, (1+self["a"])/self.get_cd(wt, 0.)

    def get_transfer_functions(
        self,
        waterp,
        w,
        omega,
        omega_num=100,
        v_air=0.0,
        cd=None,
        cd_kwargs={},
        Ap=0.0,
    ):
        """get transfer function with respect to moving isopycnal and piston
        
        Parameters
        ----------
        waterp: ufloat.water.waterp
            Water profile
        w: float
            typical vertical velocity
        omega: tuple, array-like
            min/max frequencies in Hz
        omega_num: int, optional
            size of frequency vector
        v_air: float, optional
            volume of air in m3
        cd: float, optional
            drag coefficient
        Lv: float, optional
            length scale used for the computation of the drag coefficient, see get_cd
        Ap: np.array?
            piston dispalcements ... to be implemented really
        """

        # massage water data
        water = waterp.at(extra=True)
        p = water["pressure"]
        z = water["z"]
        rho = water["rho"]
        drho_dz = water["drho_dz"]
        dt_dz = water["dtemperature_dz"]
        N2 = water["N2"]
        dtheta_dz = water["dtheta_dz"]
        Gamma = water["adiabatic_lapse_rate"]
        if not isinstance(water, dict):
            # miss the opportunity for broadcasting ...
            p = p.data
            z = z.data
            rho = rho.data
            drho_dz = drho_dz.data
            dt_dz = dt_dz.data
            N2 = N2.data
            dtheta_dz = dtheta_dz.data
            Gamma = Gamma.data

        gamma = self.get_gamma(v_air=v_air, p=p)
        M2 = g * (-drho_dz / rho - self.alpha * dt_dz - gamma)
        L2 = g * self.alpha * dtheta_dz
        L2b = g * self.alpha * (dt_dz + Gamma)

        if cd is None:
            cd = self.get_cd(w, 0.0, **cd_kwargs)

        if isinstance(omega, tuple):
            omega_hz = np.logspace(*omega, num=omega_num)  # Hz = 1/s
        else:
            omega_hz = omega
        omega = omega_hz * 2 * np.pi  # rad/s

        ds = xr.Dataset(
            None,
            coords=dict(z=("z", z), omega=("omega", omega_hz)),
        )
        ds["_omega"] = ("omega", omega)
        ds["v_air"] = v_air
        ds["M2"] = ("z", M2)
        ds["N2"] = ("z", N2)
        ds["L2"] = ("z", L2)
        ds["L2b"] = ("z", L2b)

        _denum = ds.M2 - (1 + self.a) * ds._omega**2 - 1j * ds._omega * cd
        ds["H_w^r"] = (ds.M2 - ds.N2 + ds.L2) / _denum
        ds["H_p^r"] = -g * Ap / self["V"] / _denum

        ds["H_w^f"] = 1 - ds["H_w^r"]
        ds["H_p^f"] = -ds["H_p^r"]

        ds["H_f^r"] = ds["H_w^r"] / ds["H_w^f"]
        ds["H_p^rbis"] = ds["H_p^r"] / (1 - ds["H_w^r"])

        return ds

    def time_step(
        self,
        waterp,
        T=600.0,
        dt_step=1.0,
        z=None,
        w=None,
        v=None,
        t0=0.0,
        ctrl=None,
        z_target=None,
        kalman=None,
        eta=lambda t: 0.0,
        log=True,
        dt_log=10.0,
        p_float=1.0e5,
        restart=False,
        verbose=0,
        **kwargs,
    ):
        """Time step the float position given initial conditions

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
        z_target: function
            Target trajectory as a function of time [m]
        w_target: function
            Target velocity as a function of time [m.^s-1]
        ctrl: dict
            Contains control parameters
        kalman: boolean, dict
            Activates kalman state estimate.
        eta: function
            Isopycnal displacement as a function of time
        log: boolean
            List of variables that will logged, default to True
        dt_log: float
            Time interval between log storage
        p_float: float [Pa]
            Internal float pressure in Pa (energy conssumption)
        restart:
        verbose:
        **kwargs: passed to compute_force
        """
        if restart:
            t0 = self._t
        t = t0
        #
        if restart:
            assert hasattr(self, "z"), "Restart but attribute z not existing !?"
            pass
        elif z is not None:
            self.z = z
        elif not hasattr(self, "z"):
            print("Default initial position used: z=0")
            self.z = 0.0
        #
        if restart:
            assert hasattr(self, "w"), "Restart but attribute w not existing !?"
            pass
        elif w is not None:
            self.w = w
        elif not hasattr(self, "w"):
            print("Default initial vertical velocity used: w=0")
            self.w = 0.0
        self.dwdt = 0.0
        #
        #if Lv is not None:
        #    self.Lv = Lv
        #else:
        #    self.Lv = self.L
        # log init
        logs = ["state", "dynamics", "water"]
        # _log = {'state':['z','w','v','dwdt', 'volume',
        #                 'water_temperature','water_rho'
        #                 ]
        #        }
        # _log['dynamics'] = ['acceleration','buoyancy','drag']
        # kalman initialisation
        if kalman:
            if not hasattr(self, "kalman") or isinstance(kalman, dict):
                self.init_kalman(kalman, self.w, self.z, verbose)
            # _log['kalman'] = self.kalman.log_variables
            logs.append("kalman")
            print(self.kalman)
        #
        if ctrl:
            self.init_control(ctrl, v, dt_step)
            if verbose > 0:
                print(self.ctrl)
            logs.append("piston")
            logs.append("control")
            # _log['piston'] = ['u', 'work']
            # _log['control'] = []
            u = 0.0
        elif v is None:
            if not hasattr(self, "v"):
                self.v = 0.0
        else:
            self.v = v
        v0 = self.v
        #
        if log:
            if not restart:
                # self.log = {logname:logger(vars) for logname, vars in _log.items()}
                self.log = {key: logger() for key in logs}
            _log_now = False
        else:
            self.log = None
        if log and ctrl:
            # used to update contributions to u
            self.ctrl.log = self.log["control"]
            # self.ctrl._log_now = self._log_now
        #
        print("Start time stepping for %d min ..." % (T / 60.0))
        #
        _force = 0.0
        while t < t0 + T:
            #
            if log and (dt_log is not None) and t_modulo_dt(t, dt_log, dt_step):
                _log_now = True
            else:
                _log_now = False
            # get vertical force on float
            waterp.update_eta(eta, t)  # update isopycnal displacement
            water = waterp.at(z=self.z)
            _force, _force_b, _force_d, _force_w = self.compute_force(
                self.z, self.w, water, v=self.v, sum=False, 
                **kwargs,
            )
            #
            # update kalman state estimation
            if kalman and t_modulo_dt(t, self.kalman.dt, dt_step):
                self.kalman.step(self.v, self.z)
            #
            # control starts here
            if ctrl:
                # activate control only if difference between the target and actual vertical
                # position is more than the dz_nochattering threshold
                ctrl_now = (
                    np.abs(self.z - z_target(t)) > self.ctrl.dz_nochattering
                ) and t_modulo_dt(t, self.ctrl.dt, dt_step)
                if ctrl_now:
                    u = self.get_control(z_target, t, water, _log_now)
                    _dt = self.ctrl.dt
                else:
                    _dt = 0.0
                if self.ctrl.continuous:
                    # apply control at each model timestep
                    _dt = dt_step
                _v0 = self.piston.vol
                if "feedback1" in self.ctrl.mode:
                    if ctrl_now:
                        self.piston.update_vol(u)
                else:
                    # u units = m^3/s
                    self.piston.update(_dt, u)
                self.v = self.piston.vol
                # energy integration, 1e4 converts from dbar to Pa
                if self.v != _v0:
                    _dv = self.v - _v0
                    self.piston_work += (
                        abs((water["pressure"] * 1.0e4 - p_float) * _dv)
                        * watth
                        / self.piston.efficiency
                    )
            #
            # log data
            if _log_now:
                self.log["state"].log(
                    time=t,
                    z=self.z,
                    w=self.w,
                    v=self.v,
                    volume=self._volume,
                    dwdt=_force / (1 + self.a) / self["m"],
                )
                self.log["water"].log(
                    time=t,
                    temperature=water["temperature"],
                    salinity=water["salinity"],
                    rho=water["rho"],
                    eta=water["eta"],
                    detadt=water["detadt"],
                    d2etadt2=water["d2etadt2"],
                )
                self.log["dynamics"].log(
                    time=t,
                    acceleration=_force / (1 + self.a) / self["m"],
                    buoyancy=_force_b / (1 + self.a) / self["m"],
                    drag=_force_d / (1 + self.a) / self["m"],
                    water=_force_w / (1 + self.a) / self["m"],
                )
                if ctrl:
                    _info = {"time": t, "u": u, "work": self.piston_work}
                    self.log["piston"].log(**_info)
                if kalman:
                    #
                    _dwdt = _force / (1 + self.a) / self["m"]
                    _klog = self.kalman.get_log(t, self.v, _dwdt)
                    self.log["kalman"].log(**_klog)
            #
            # update variables
            self.dwdt = _force / (1 + self.a) / self["m"]
            self.w += dt_step * self.dwdt
            self.z += dt_step * self.w
            self.z = np.amin((self.z, 0.0))
            t += dt_step

        self._t = t  # for restarts

        # finalize logs
        for _, _log in self.log.items():
            _log.finalize()

        print("... time stepping done")

    def init_control(self, ctrl, v, dt_step):
        if v is not None:
            self.piston.update_vol(v)
            self.piston.reset_phi_float()
        self.v = self.piston.vol
        self.piston_work = 0.0
        if ctrl:
            ctrl_default = {"dt": dt_step, "dz_nochattering": 0.0}
            if ctrl["mode"] == "sliding":
                ctrl_default.update(
                    {"tau": 60.0, "mode": "sliding", "waterp": waterp, "Lv": self.L}
                )
            elif ctrl["mode"] == "pid":
                ctrl_default.update(Kp=0.0, Ki=0.0, Kd=0.0)
            elif ctrl["mode"] == "pid_position":
                ctrl_default.update(dvdt=0.0, Kp=0.0, Ki=0.0, Kd=0.0)
            elif "feedback" in ctrl["mode"]:
                ctrl_default.update(
                    {
                        "tau": 3.25,  # Set the root of feed-back regulation # s assesed by simulation
                        "nu": 0.03
                        * 2.0
                        / np.pi,  # Set the limit speed : 3cm/s assesed by simulation
                        "delta": 0.11,  # length scale that defines the zone of influence around the target depth, assesed by simulation
                        "gamma": self.gamma,  # mechanical compressibility [1/dbar]
                        "m": self["m"],
                        "a": self.a,
                        "Lv": self.Lv,
                        "c1": self.c1,
                        "gammaV": self.gammaV,
                        "rho_cte": self.rho_cte,
                        "perturbation": 0.0,
                    }
                )
            if "kalman" in ctrl["mode"]:
                _k = self.kalman
                ctrl_default.update(
                    {
                        "m": _k.m,
                        "a": _k.a,
                        "Lv": _k.Lv,
                        "c1": _k.c1,
                        "rho_cte": _k.rho_cte,
                    }
                )
            #
            ctrl_default.update(ctrl)
            self.ctrl = control(**ctrl_default)

    def get_control(self, z_target, t, waterp, log, **kwargs):
        _c = self.ctrl
        if _c.mode == "sliding":
            u = _c.get_u_sliding(z_target, t, self.z, self.w, self)
        elif _c.mode == "pid":
            u = _c.get_u_pid(z_target, t, self.z, log)
        elif _c.mode == "pid_position":
            # this is not proper PID ...
            d_pid = _c.get_u_pid(z_target, t, self.z, log)
            # u = np.sign(d_pid - self.piston.d) * _c.dvdt
            u = np.sign(d_pid) * _c.dvdt
        elif _c.mode == "feedback1":
            _Ve = self.compute_force(self.z, 0.0, waterp, v=0.0, sum=True, **kwargs) / (
                g * self.rho_cte
            )
            _Ve += _c.perturbation  # sensitivity tests
            u = _c.get_u_feedback1(z_target, t, self.z, self.w, _Ve, 0.0, log)
        elif _c.mode == "feedback2":
            dwdt = self.dwdt
            dwdt += _c.perturbation  # sensitivity tests
            u = _c.get_u_feedback2(z_target, t, self.z, self.w, dwdt, self.gammaV, log)
        elif _c.mode == "kalman_feedback1":
            _k = self.kalman
            assert (
                _k.version == "v0"
            ), "Not implemented for kalman version different than v0"
            u = _c.get_u_feedback1(
                z_target, t, -_k.x_hat[1], -_k.x_hat[0], _k.x_hat[3], _k.x_hat[2], log
            )
        elif _c.mode == "kalman_feedback2":
            _k = self.kalman
            _dwdt = -_k.f(_k.x_hat, self.v)[0]
            if _k.version == "v0":
                _gamma1 = _k.x_hat[2]
                _gamma2 = 0.0
                _c1 = 1.0
            elif _k.version == "v1":
                _gamma1 = _k.x_hat[3]
                _gamma2 = _k.x_hat[4]
                _c1 = _k.x_hat[5]
            u = _c.get_u_feedback2(
                z_target,
                t,
                -_k.x_hat[1],
                -_k.x_hat[0],
                _dwdt,
                _gamma1,
                log,
                gamma2=_gamma2,
                c1=_c1,
            )
        else:
            print("%s is not a valid control method" % _c.mode)
            return
        return u

    def init_kalman(self, kalman, w, z, verbose):
        _params = {
            "version": "v0",
            "m": self["m"],
            "a": self.a,
            "rho_cte": self.rho_cte,
            "c1": self.c1,
            "Lv": self.L / 4,
            "x_init": [-w, -z, self.gammaV, 0.0],
            "verbose": verbose,
        }
        if isinstance(kalman, dict):
            _params.update(kalman)
        self.kalman = eval("kalman_" + _params["version"])(**_params)

    def plot_logs(self, **kwargs):
        """wrapper around plot_logs"""
        plot_logs(self.log, self, **kwargs)

    def store_logs(self, data_dir, label, **kwargs):
        """ store logs in netcdf files """
        for log_label, log in self.log.items():
            file = os.path.join(data_dir, label+"_"+log_label+".nc")
            log.store(file, **kwargs)
            print(f"log {log_label} stored in {file}")

    def load_logs(self, data_dir, label):
        """ reload logs from netcdf files """
        log_files = sorted(glob(os.path.join(data_dir, label+"_*.nc")))
        log = {}
        for file in log_files:
            log_label = file.split("_")[-1].replace(".nc", "")
            log[log_label] = logger(file=file)
        self.log = log


# ---------------------------- piston -------------------------------------------
class piston:
    """Piston object, facilitate float buoyancy control"""

    def __init__(self, model="ifremer", **kwargs):
        """Piston object

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
        params = {
            "phi": 0.0,
            "d": 0.0,
            "vol": 0.0,
            "omega": 0.0,
            "phi_min": 0.0,
            "d_min": 0.0,
            "vol_min": 0.0,
        }
        if model is None:
            model = "custom"
        if model.lower() == "ifremer":
            # default parameters: IFREMER float
            params.update(
                r=0.0195 / 2,
                lead= 1,
                d_max= 0.09,
                translation_max= 0.01,  # m/s !!! wild guess
                translation_min= 0.001,  # m/s !!! wild guess
                efficiency= 0.1,
                d_increment= 0.12 * 4.0 / 5600.0,
                d_offset= 0.03,  # big piston offset
            )
        elif model.lower() == "seabot_v0":
            # default parameters: ENSTA float
            params.update(
                r=0.025,
                lead= 0.00175,
                tick_per_turn= 48,
                d_max= 0.07,
                vol_max= 1.718e-4,
                omega_max= 60.0 / 48 * 2.0 * np.pi,
                efficiency= 0.1,
                d_offset= 0.0,
            )
            self.d_increment = params["lead"] / params["tick_per_turn"]

            # translation_max = d_increment*(225 pulses par seconde)  (vitesse de translation max en m/s)
            # translation_min fixe arbitrairement pour l'instant

            # d_increment le 4 vient du facteur 4 de la roue codeuse de thomas
            # d_increment = 0.12/5600 ou 0.090/4200 selon la prise en compte ou non du gros piston

            # dmax = 0.102 ancienne valeur pour IFREMER

            # verifier si importance lead et angles lors de la regulation, ecraser parametres redondants
            # vol_max = 0.090*np.pi*(0.0195/2)**2+0.030*np.pi*(0.080/2)**2 = 1.777e-4

            # 48 encoches
            # vitesse max de 60 encoches par seconde
        elif model.lower() == "seabot":
            # default parameters: ENSTA float
            params.update(
                r= 4.5e-2*0.5,
                lead= 1e-3, # m/turn
                tick_per_turn= 8192,
                d_max= 11e-2, # m
                omega_max= 38.0 / 60 * 2.0 * np.pi,
                efficiency= 0.1,
            )
            self.d_increment = params["lead"] / params["tick_per_turn"]

        #
        params.update(kwargs)
        for key, val in params.items():
            setattr(self, key, val)
        self.A = np.pi * self.r**2 # piston area
        # assumes here volumes are given
        # if 'd_min' not in kwargs:
        #    self.d_min = self.vol2d(self.vol_min)
        # (useless as will reset d_min to 0.)
        if "vol_max" in params:
            if "d_max" in params:
                self.d_max = params["d_max"]
                self.vol_min = self.d2vol_no_volmin(self.d_min)
                print("Piston vol_min reset from d_min, d_max and vol_max")
            else:
                self.d_max = self.vol2d(self.vol_max)
                print("Piston max displacement set from max volume")
        elif "d_max" in params:
            self.vol_max = self.d2vol(self.d_max)
            print("Piston max volume set from max displacement")
        else:
            print("You need to provide d_max or vol_max")
            sys.exit()
        if "vol" not in kwargs:
            self.vol = self.vol_max
        if "translation_max" in params:
            self.omega_max = params["translation_max"] * 2.0 * np.pi / self.lead
        if "translation_min" in params:
            self.omega_min = params["translation_min"] * 2.0 * np.pi / self.lead
        #
        self.phi_max = self.d2phi(self.d_max)
        self.update_dvdt()
        #
        self.u_max = self.omega2dvdt(self.omega_max)
        #
        self.phi_float = self.phi
        self.phi_increment = self.d2phi(self.d_min + self.d_increment) - self.d2phi(
            self.d_min
        )
        # not sure if vol_increment is still used
        self.vol_increment = (
            self.d_increment * self.A
        )  # *self.increment_error
        self.tick_to_volume = self.vol_increment

    def __repr__(self):
        strout = "Piston parameters and state: \n"
        strout += "  r     = %.2f cm        - piston radius\n" % (self.r * 1.0e2)
        # strout+='  phi   = %.2f rad       - present angle of rotation\n'%(self.phi)
        strout += "  d     = %.2f cm        - present piston displacement\n" % (
            self.d * 1.0e2
        )
        strout += "  vol   = %.2f cm^3      - present volume addition\n" % (
            self.vol * 1.0e6
        )
        # strout+='  lead  = %.2f cm        - screw lead\n'%(self.lead*1.e2)
        # strout+='  tick_per_turn  = %.2f no dimension        - number of notches on the thumbwheel of the piston\n'%(self.tick_per_turn)
        strout += (
            "  d_increment  = %.2e mm        - smallest variation of translation motion for the piston\n"
            % (self.d_increment * 1e3)
        )
        strout += (
            "  vol_increment  = %.2e cm^3    - smallest variation of volume possible for the piston\n"
            % (self.vol_increment / cm3)
        )
        # strout+='  phi_max = %.2f deg     - maximum rotation\n'%(self.phi_max*1.e2)
        # strout+='  phi_min = %.2f deg     - minimum rotation\n'%(self.phi_min*1.e2)
        strout += "  d_max = %.2f mm        - maximum piston displacement\n" % (
            self.d_max * 1.0e3
        )
        strout += "  d_min = %.2f mm        - minimum piston displacement\n" % (
            self.d_min * 1.0e3
        )
        strout += "  vol_min = %.2f cm^3    - min volume displaced\n" % (
            self.vol_min * 1.0e6
        )
        strout += "  vol_max = %.2f cm^3    - max volume displaced\n" % (
            self.vol_max * 1.0e6
        )
        # strout+='  omega_max = %.2f deg/s - maximum rotation rate\n'%(self.omega_max*180./np.pi)
        # strout+='  omega_min = %.2f deg/s - minimum rotation rate\n'%(self.omega_min*180./np.pi)
        strout += "  u_max = %.2f cm^3/s - maximum volume rate of change\n" % (
            self.u_max / cm3
        )
        strout += (
            "  efficiency = %.2f - mechanical work produced / electrical work supplied\n"
            % (self.efficiency)
        )
        return strout

    # -------------------------- update methods -------------------------------------

    def update(self, dt, dvdt):
        """Update piston position given time interval and a volume rate of change

        Parameters
        ----------
        dt: float
            time interval in seconds
        dvdt: float
            desired volume rate of change in m^3/s
        """
        omega = self.dvdt2omega(dvdt)
        self.update_omega(omega)
        self.update_dvdt()
        self.update_phi(dt)

    def update_phi(self, dt):
        self.phi_float += self.omega * dt
        dphi = self.phi_float - self.phi
        self.phi += np.trunc(dphi / self.phi_increment) * self.phi_increment
        self._checkbounds()
        self._bcast_phi()

    def update_vol(self, vol):
        self.phi = self.vol2phi(vol)
        self._checkbounds()
        self._bcast_phi()

    def update_d(self, d):
        self.phi = self.d2phi(d)
        self._checkbounds()
        self._bcast_phi()

    def update_omega(self, omega):
        # if np.abs(omega)<self.omega_min:
        #    self.omega=0.
        # else:
        self.omega = np.sign(omega) * min([np.abs(omega), self.omega_max])

    def update_dvdt(self):
        self.dvdt = self.omega2dvdt(self.omega)

    def _bcast_phi(self):
        self.d = self.phi2d(self.phi)
        self.vol = self.phi2vol(self.phi)

    def _checkbounds(self):
        self.phi = min([max([self.phi, self.phi_min]), self.phi_max])
        self.phi_float = min([max([self.phi_float, self.phi_min]), self.phi_max])

    def reset_phi_float(self):
        self.phi_float = self.phi

    # --------------------- conversion methods --------------------------------------

    def omega2dvdt(self, omega):
        # to compute omega for the ENSTA float:
        # the motor of the piston needs 48 notches to complete a full rotation
        # it can reach until 30 rotations a seconde
        # so omega = 2*pi*30/48 = 3.9 rad/s
        return omega * self.lead / 2.0 * self.r**2

    def dvdt2omega(self, dvdt):
        return dvdt / (self.lead / 2.0 * self.r**2)

    def phi2d(self, phi):
        return self.d_min + (phi - self.phi_min) / 2.0 / np.pi * self.lead

    def phi2vol(self, phi):
        return self.d2vol(self.phi2d(phi))

    def d2phi(self, d):
        return self.phi_min + (d - self.d_min) * 2.0 * np.pi / self.lead

    def d2vol(self, d):
        return self.vol_min + (d - self.d_min) * self.A

    def d2vol_no_volmin(self, d):
        return self.vol_max + (d - self.d_max) * self.A

    def vol2d(self, vol):
        return self.d_min + (vol - self.vol_min) / self.A

    def vol2phi(self, vol):
        return self.d2phi(self.vol2d(vol))


# ----------------------------- ballast -----------------------------------------

_bvariables = [
    "mass_in_water",
    "mass_in_air",
    "water_temperature",
    "piston_displacement",
]


class ballast(object):

    # in db
    pressure = 1.0
    # brest:
    lon, lat = 4.5, 38.5

    def __init__(self, file=None):
        """Class handling the ballasting procedure"""
        self._d = {}
        if file is not None:
            self.load(file)

    def __repr__(self):
        return pformat(self._d, indent=2, width=1)

    def add_ballast(
        self,
        name,
        comments="None",
        **kwargs,
    ):
        """Add a mass measurement

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
        water_salinity: float, optional
            in psu, one of water salinity or conductivity must be passed
        water_conductivity: float, optional
            in S, one of water salinity or conductivity must be passed
        piston_displacement: float, str
            positive out and in m or 'max' for all out
        comments: str
            additional comments
        """
        assert all(v in kwargs for v in _bvariables), print(
            "All following variables should be passed as kwargs: "
            + ", ".join(_bvariables)
        )
        if "water_conductivity" in kwargs:
            _SP = gsw.SP_from_C(
                kwargs["water_conductivity"], kwargs["water_temperature"], self.pressure
            )
            _d["water_salinity"] = _SP
            print("Practical salinity derived from conductivity: SP=%.2f" % _SP)
        elif "water_salinity" in kwargs:
            _SP = kwargs["water_salinity"]
        else:
            print("salinity or conductivity must be provided")
        self._d[name] = {v: kwargs[v] for v in _bvariables}
        self._d[name]["water_salinity"] = _SP
        # derive absolute salinity and conservative temperature
        SA = gsw.SA_from_SP(_SP, self.pressure, self.lon, self.lat)
        CT = gsw.CT_from_t(SA, kwargs["water_temperature"], self.pressure)
        rho_water = gsw.density.rho(SA, CT, self.pressure)  # in situ density
        # store derived variables
        _d = self._d[name]
        _d["SA"] = SA
        _d["CT"] = CT
        _d["rho_water"] = rho_water
        _d["V"] = (_d["mass_in_air"] - _d["mass_in_water"]) / 1e3 / rho_water
        _d["comments"] = comments

    def compute_mass_adjustment(self, f, w=None, verbose=False, **kwargs):
        """
        Compute the mass adjustement to perform prior to underwater float deployment
            \delta_m = -V(new) \rho_w(new) - V(ballast) \rho_w(ballast)
        Or its approximation:
            \delta_m = -V \delta \rho_w - \rho_w \delta V

        !! this method automatically adjust the mass and volume of the underwater float object !!

        Parameters
        ----------
        f: cognac.ufloat.core.ufloat
            float object
        w: cognac.ufloat.water.waterp
            water profile
        verbose: boolean, optional
            print extra information
        """
        if w is not None:
            z = -self.pressure
            water_temperature = w.get_temp(z)
            water_salinity = w.get_s(z)
            lon, lat = w.lon, w.lat
        else:
            water_temperature = kwargs["temperature"]
            if "salinity" in kwargs:
                water_salinity = kwargs["salinity"]
            elif "conductivity" in kwargs:
                water_salinity = gsw.SP_from_C(
                    kwargs["conductivity"], kwargs["temperature"], self.pressure
                )
            lon, lat = kwargs["lon"], kwargs["lat"]
            if isinstance(lon, str):
                lon = ll_conv(lon)
            if isinstance(lat, str):
                lat = ll_conv(lat)
        SA = gsw.SA_from_SP(water_salinity, self.pressure, lon, lat)
        CT = gsw.CT_from_t(SA, water_temperature, self.pressure)
        rho_water = gsw.density.rho(SA, CT, self.pressure)  # new in situ density
        print("New in situ density: {:.2f} kg/m3".format(rho_water))
        #
        for name, b in self._d.items():
            delta_rho_water = rho_water - b["rho_water"]
            # reset float mass
            f.m = b["mass_in_air"] * 1e-3  # need to convert from g to kg
            if b["piston_displacement"] == "max":
                _piston_displacement = f.piston.d_max
            else:
                _piston_displacement = b["piston_displacement"]  # m
            f.piston.update_d(_piston_displacement)
            f.piston_update_vol()
            # reset volume and temperature reference
            f.V = b["V"] - f.v  # reset volume to match volume infered from ballasting
            f.temp0 = b["water_temperature"]
            #
            Vb = f.volume(p=self.pressure, temp=b["water_temperature"])
            V = f.volume(p=self.pressure, temp=water_temperature)
            delta_V = V - Vb
            #
            delta_m0 = -b["mass_in_water"] * 1e-3
            # delta_m1 = V * delta_rho_water + rho_water * delta_V # approximation of more accurate form:
            delta_m1 = V * rho_water - Vb * b["rho_water"]
            delta_m = delta_m0 + delta_m1
            #
            print(
                "-- According to ballast %s, you need add %.1f g"
                % (name, delta_m * 1e3)
            )
            print(
                "! If the weight is external to the float, remember this must be the value in water"
            )
            #
            if verbose:
                print(
                    "Independent mass correction (ballast -mass_in_water): %.1f [g]"
                    % (delta_m0 * 1e3)
                )
                print("  Water change mass correction: %.1f [g]" % (delta_m1 * 1e3))
                delta_m_t = (
                    -gsw.density.alpha(SA, CT, self.pressure)
                    * (CT - b["CT"])
                    * b["rho_water"]
                    * V
                )
                print("    T contribution: %.1f [g]" % (delta_m_t * 1e3))
                delta_m_s = (
                    gsw.density.beta(SA, CT, self.pressure)
                    * (SA - b["SA"])
                    * b["rho_water"]
                    * V
                )
                print("    S contribution: %.1f [g]" % (delta_m_s * 1e3))
                print(
                    "  New/old in situ density: %.2f, %.2f  [kg/m^3] "
                    % (rho_water, b["rho_water"])
                )
                print("    difference: %.1f [kg/m^3]" % (delta_rho_water))
                print("  New/old float volume: %.5e, %.5e  [m^3] " % (V, b["V"]))
                print("    difference: %.1f [cm^3]" % (delta_V / cm3))
            #
            return delta_m

    def store(self, file):
        _d = dict(self._d)
        _v = ["SA", "CT", "rho_water", "V", "water_salinity"]
        for name, b in _d.items():
            for v in _v:
                b[v] = float(b[v])
        with open(file, "w") as ofile:
            documents = yaml.dump(self._d, ofile)
        print("ballasting data stored in %s" % file)

    def load(self, file):
        with open(file) as ifile:
            d = yaml.full_load(ifile)
        # should check that content of d is valid
        self._d.update(d)
        print("File %s loaded" % file)


# ----------------------------- utils ------------------------------------------


def ll_conv(ll):
    """Returns lon or lat as floating degree and vice versa
    string format should be '43N09.007'
    This piece of code should be elsewhere
    """
    if isinstance(ll, str):
        if "N" in ll:
            _ll = ll.split("N")
            sign = 1.0
        elif "S" in ll:
            _ll = ll.split("S")
            sign = -1.0
        elif "E" in ll:
            _ll = ll.split("E")
            sign = 1.0
        elif "W" in ll:
            _ll = ll.split("W")
            sign = -1.0
        return sign * (float(_ll[0]) + float(_ll[1]) / 60.0)
    elif isinstance(ll, float):
        return "%dX%.6f" % (int(ll), (ll - int(ll)) * 60.0)


def plot_float_density(
        z, f, waterp, 
        mid=False, 
        ax=None, 
        v_air=None, 
        show_no_air=False,
        xlim=None,
    ):
    """Plot float density with respect to that of water profile

    Parameters
    ----------
    z: ndarray
        depth grid
    f: cognac.ufloat.core.float
        underwater float
    waterp: cognac.ufloat.water.waterp
        water profile
    mid: boolean,
        True for curve at piston mid displacement
    ax: matplotlib.pyplot.axes
    v_air: float
        Volume of air at the surface in m^3
    """
    #
    w = waterp.at(z=z)
    rho_w, p, temp = w["rho"], w["pressure"], w["temperature"]
    if v_air is not None:
        air = (v_air, waterp.at(z=-1.0)["temperature"])
    else:
        air = None
    #
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    if xlim is not None:
        ax.set_xlim(*xlim)
    lw = 4
    #
    # rho_f = f.rho(p=p, temp=temp, v=0.)
    piston = hasattr(f, "piston")
    if piston:
        rho_f_vmax = f.rho(p=p, temp=temp, v=f.piston.vol_max, air=air)
        rho_f_vmin = f.rho(p=p, temp=temp, v=f.piston.vol_min, air=air)
        ax.fill_betweenx(
            z,
            rho_f_vmax,
            rho_w,
            where=rho_f_vmax >= rho_w,
            facecolor="red",
            interpolate=True,
        )
        if air is not None and show_no_air:
            rho_f_vmax_no_air = f.rho(p=p, temp=temp, v=f.piston.vol_max)
            rho_f_vmin_no_air = f.rho(p=p, temp=temp, v=f.piston.vol_min)
    else:
        rho_f = f.rho(p=p, temp=temp, air=air)
        if air is not None and show_no_air:
            rho_f_no_air = f.rho(p=p, temp=temp)
    #
    ax.plot(rho_w, z, "b", lw=lw, label="water")
    if piston:
        ax.plot(
            rho_f_vmax, z, "-", color="orange", lw=lw, label="float v_max", markevery=10
        )
        ax.plot(rho_f_vmin, z, "--", color="orange", lw=lw, label="float v_min")
        if air is not None and show_no_air:
            ax.plot(
                rho_f_vmax_no_air,
                z,
                "-",
                color="orange",
                lw=lw / 2,
                label="float v_max - no air",
                markevery=10,
            )
            ax.plot(
                rho_f_vmin_no_air,
                z,
                "--",
                color="orange",
                lw=lw / 2,
                label="float v_min - no air",
            )
    else:
        ax.plot(rho_f, z, "-", color="orange", lw=lw, label="float")
        if air is not None and show_no_air:
            ax.plot(
                rho_f_no_air, z, "-", color="orange", lw=lw / 2, label="float - no air"
            )

    if piston and mid:
        # rho_f_vmid=f.rho(p=p, temp=temp, v=(f.piston.vol_max+f.piston.vol_min)*.5)
        rho_f_vmid = f.rho(p=p, temp=temp, v=mid, air=air)
        ax.plot(rho_f_vmid, z, "-", color="grey", lw=lw, label="float v_mid")
        if air is not None and show_no_air:
            rho_f_vmid_no_air = f.rho(p=p, temp=temp, v=mid)
            ax.plot(rho_f_vmid_no_air, z, "--", color="grey", lw=lw/2, label="float v_mid - no air")

    ax.legend()
    ax.set_xlabel("[kg/m^3]")
    ax.set_ylim((np.amin(z), 0.0))
    ax.set_ylabel("z [m]")
    ax.grid()
    iz = np.argmin(np.abs(z))
    if piston:
        delta_mass_surf = f.m * (rho_w[iz] / rho_f_vmax[iz] - 1.0) * 1e3
        title = "extra mass @surface (piston full out): %.1f g" % (delta_mass_surf)
        if air:
            title = title + f"\n volume air @surface = {v_air/cm3} cm3"
        ax.set_title(title)
    #
    y_annotation = ax.get_ylim()[1] - 0.1 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.annotate(
        "",
        xy=(ax.get_xticks()[1], y_annotation),
        xytext=(ax.get_xticks()[2], y_annotation),
        arrowprops=dict(arrowstyle="<->"),
    )
    ax.text(
        ax.get_xticks()[1] * 0.5 + ax.get_xticks()[2] * 0.5,
        y_annotation,
        "%.1f g" % (f.V * (ax.get_xticks()[2] - ax.get_xticks()[1]) * 1e3),
        {"ha": "center", "va": "bottom"},
    )
    #        ax.get_ylim()[1]-10,
    return ax


#
def plot_float_volume(z, f, waterp, ax=None):
    """Plot float volume with respect to depth

    Parameters
    ----------
    z: ndarray
        depth grid
    f: float object
    waterp: water profile object

    """
    #
    rho_w, p, temp = waterp.get_rho(z), waterp.get_p(z), waterp.get_temp(z)
    #
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    lw = 4
    #
    # iz = np.argmin(np.abs(z+500))
    v = f.volume(p=p, temp=temp, v=0.0)
    vmax = f.volume(p=p, temp=temp, v=f.piston.vol_max)
    vmin = f.volume(p=p, temp=temp, v=f.piston.vol_min)
    #
    ax.plot(
        vmax / cm3, z, "-", color="orange", lw=lw, label="float v_max", markevery=10
    )
    ax.plot(vmin / cm3, z, "--", color="orange", lw=lw, label="float v_min")
    ax.legend()
    ax.set_xlabel("[cm^3]")
    ax.set_ylim((np.amin(z), 0.0))
    ax.set_ylabel("z [m]")
    ax.grid()
    return ax


def plot_equilibrium_volume(f, w, zlim=(-100, 0)):
    """Plot the piston volume required to be at equilibrium as a function
    of depth
    """
    z = np.arange(zlim[0], zlim[1], 5.0)
    rho_w, p, temp = w.get_rho(z), w.get_p(z), w.get_temp(z)
    v = np.array([f.volume4equilibrium(p, t, rho) for p, t, rho in zip(p, temp, rho_w)])

    fig, ax = plt.subplots()

    ax = plt.subplot(111)
    ax.plot(v / cm3, z, linewidth=4)
    ax.set_xlabel("equilibrium v [cm^3]")
    ax.set_ylabel("z [m]")
    ax.grid()

    ax1 = ax.twiny()
    p = f.piston
    tick2volume = p.vol_increment
    v2tick = lambda v: (p.vol_max - v / cm3) / tick2volume + p.d_offset / p.d_increment
    ax1.set_xlim((v2tick(l) for l in ax.get_xlim()))
    _ = ax1.set_xlabel("equilibrium tick")
    return fig, ax, ax1


def compute_gamma(R, t, material=None, E=None, mu=0.35):
    """Compute the compressibility to pressure of a cylinder

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

    """
    pmat = {
        "glass": {"E": 90.0, "mu": 0.25},
        "aluminium": {"E": 70.0, "mu": 0.35},
        "pom": {"E": 3.5, "mu": 0.35},
        "polycarbonat": {"E": 2.2, "mu": 0.37},
    }
    if material is not None:
        if material in pmat:
            E = pmat[material]["E"]
            mu = pmat[material]["mu"]
        else:
            print("material not in our database")
            sys.exit()
    elif E is None:
        print("You need to provide material or E !")
        sys.exit()
    # convert E to dbar
    E = E * 1.0e5
    return R * (6.0 - 7.0 * mu) / 2.0 / E / t


def t_modulo_dt(t, dt, dt_step):
    threshold = 0.25 * dt_step / dt
    if np.abs(t / dt - np.rint(t / dt)) < threshold:
        return True
    else:
        return False


# build scenarios
def descent(Tmax, zt, f=None, waterp=None, wmax=None, zstart=0):
    """Contruct trajectory of a descent to some depth
    Parameters
    ----------
    Tmax: float
        Time length of the trajectory in SECONDS
    zt: target depth level
    f: float object
        Used to compute maximum accelerations
    waterp: water profile object
        Used to compute maximum accelerations

    """
    # build time line
    t = np.arange(0.0, Tmax, 1.0)
    # compute bounds on motions
    if f is not None and waterp is not None:
        fmax, fmin, afmax, wmax = f.compute_bounds(waterp, zt)
        dzdt_target = -t * afmax / 2.0 / f.m
        dzdt_target[np.where(-dzdt_target > wmax)] = -np.abs(wmax)
    elif wmax is not None:
        dzdt_target = t * 0.0 - np.abs(wmax)
    # build trajectory
    z_target = zstart + np.cumsum(dzdt_target * 1.0)
    z_target[np.where(z_target < zt)] = zt
    # convert to callable function
    return interp1d(t, z_target, kind="linear", fill_value="extrapolate")
