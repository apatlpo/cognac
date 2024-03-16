import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import hvplot.pandas  # noqa
import hvplot

import gsw


# ----------------------------------- core --------------------------------------

#
class waterp:
    """Data holder for a water column based on either:
     - in situ temperature, salinity, pressure profiles
     - the World Ocean Atlas (WOA) climatology, see:
        https://www.nodc.noaa.gov/OC5/woa18/

    Should base this on pandas or xarray !!

    Parameters
    ----------
    pressure: np.ndarray, optional
        pressure in dbar
    temperature: np.ndarray, optional
        in situ temperature in degC
    salinity: np.ndarray, optional
        salinity in PSU
    lon: float, optional
        longitude of the selected location
    lat: float, optional
        latitude of the selected location

    """

    def __init__(
        self,
        pressure=None,
        temperature=None,
        salinity=None,
        lon=None,
        lat=None,
        name=None,
        interp_kwargs=None,
        **kwargs
    ):

        self._pts, self._woa = False, False
        self._interp_kwargs = {
            "kind": "quadratic",
            "fill_value": None,
            "bounds_error": False,
        }
        if interp_kwargs:
            self._interp_kwargs.update(**interp_kwargs)

        args = [pressure, temperature, salinity, lon, lat]
        if all([a is not None for a in args]):
            self._load_from_pts(pressure, temperature, salinity, lon, lat, name)
        elif all([lon, lat]):
            self._load_from_woa(lon, lat, name)
        else:
            print("Inputs missing")

        # init isopycnal displacement and velocity
        self.update_eta(0.0, detadt=0.0, d2etadt2=0.0)

    def _load_from_pts(self, pressure, temperature, salinity, lon, lat, name):
        self._pts = True
        #
        if all(
            [isinstance(i, float) for i in [pressure, temperature, salinity, lon, lat]]
        ):
            # pressure = np.array([-1, 1e4])
            # temperature = np.array([temperature, temperature])
            # salinity = np.array([salinity, salinity])
            # lon = np.array([lon,lon])
            # lat = np.array([lat,lat])
            pressure = [-1, 1e4]
            temperature = [temperature, temperature]
            salinity = [salinity, salinity]
            lon = [lon, lon]
            lat = [lat, lat]
        if isinstance("pressure", pd.Series):
            pressure = pressure.values
        if isinstance("temperature", pd.Series):
            temperature = temperature.values
        if isinstance("salinity", pd.Series):
            salinity = salinity.values
        if not isinstance(lon, (int, float)):
            lon = np.mean(lon)
        if not isinstance(lat, (int, float)):
            lat = np.mean(lat)
        #
        self.lon, self.lat = lon, lat
        #
        self.update_background_profile(temperature, salinity, p=pressure)
        #
        if name is None:
            self.name = "Provided water profile at lon=%.0f, lat=%.0f" % (
                self.lon,
                self.lat,
            )

    def _load_from_woa(self, lon, lat, name):
        self._woa = True
        #
        self._tfile = "woa18_A5B7_t00_01.nc"
        nc = Dataset(self._tfile, "r")
        #
        glon = nc.variables["lon"][:]
        glat = nc.variables["lat"][:]
        ilon = np.argmin(np.abs(lon - glon))
        ilat = np.argmin(np.abs(lat - glat))
        self.lon = glon[ilon]
        self.lat = glat[ilat]
        #
        z = -nc.variables["depth"][:].data
        #
        temp = nc.variables["t_an"][0, :, ilat, ilon]
        nc.close()
        #
        self._sfile = "woa18_A5B7_s00_01.nc"
        nc = Dataset(self._sfile, "r")
        s = nc.variables["s_an"][0, :, ilat, ilon]
        nc.close()
        #
        self.update_background_profile(temp, s, z=z)
        #
        if name is None:
            self.name = "WOA water profile at lon=%.0f, lat=%.0f" % (self.lon, self.lat)

    def update_background_profile(self, temp, s, p=None, z=None):
        """update eos from (s, temp, p, lon, lat)"""
        if p is None:
            p = gsw.p_from_z(z, self.lat)
        if z is None:
            z = gsw.z_from_p(p, self.lat)
        # derive absolute salinity and conservative temperature
        SA = gsw.SA_from_SP(s, p, self.lon, self.lat)
        CT = gsw.CT_from_t(SA, temp, p)
        bg = dict(
            temperature=temp,
            salinity=s,
            SA=SA,
            CT=CT,
            pressure=p,
            z=z,
        )
        bg = xr.Dataset({k: ("z", v) for k, v in bg.items()}).set_coords(
            ["SA", "CT", "z", "pressure"]
        )
        self._append_derived_properties(bg, True)
        self.bg = bg

    def now(self, z=None, dz=1.0, extra=False, **kwargs):
        """get current water properties, i.e. potentially displaced vertically compared to
        the background profile (isopycnal displacement eta)

        Parameters
        ----------
        z: float, np.array, xr.DataArray, optional
            Vertical coordinate where water properties are thought
            Use the background profile coordinate by default
        dz: float, optional
            If z is a float, used for vertical derivative computation
        **kwargs: passed to interp
        """
        _kwargs = dict(method="linear", kwargs={"fill_value": "extrapolate"})
        _kwargs.update(**kwargs)
        single = False
        if z is None:
            z = self.bg["z"].values
        elif isinstance(z, (int, float)):
            z = np.array([z - dz, z, z + dz])
            single = True
        bg = self.bg.reset_coords()
        pr = bg[["SA", "CT"]].interp(z=(z - self.eta), **_kwargs)
        pr = pr.assign_coords(z_bg=pr.z)  # store isopycnal original depth
        pr["z"] = z  # need to reset z
        pr["pressure"] = bg["pressure"].interp(z=z, **_kwargs)  # need to reset pressure
        pr["eta"] = self.eta
        pr["detadt"] = self.detadt
        pr["d2etadt2"] = self.d2etadt2
        self._append_derived_properties(pr, extra)
        if single:
            pr = pr.isel(z=1).reset_coords()
            pr = {k: float(v) for k, v in pr.items()}
            pr["z"] = z[1]
        return pr

    def _append_derived_properties(self, ds, extra):
        """derive useful water properties from SA, CT, pressure"""
        SA, CT, p = ds["SA"], ds["CT"], ds["pressure"]
        if "temperature" not in ds:
            ds["temperature"] = gsw.t_from_CT(SA, CT, p)
        if "salinity" not in ds:
            ds["salinity"] = gsw.SP_from_SA(SA, p, self.lon, self.lat)
        # derive standard variables / could add units ...
        ds["rho"] = gsw.density.rho(SA, CT, p)  # kg/m^3
        if extra:
            ds["theta"] = gsw.conversions.pt_from_CT(SA, CT)  # degC
            ds["sigma0"] = gsw.density.sigma0(SA, CT)  # kg/m^3
            ds["kappa"] = gsw.density.kappa(SA, CT, p) * 1e4  # 1/dbar
            ds["alpha"] = gsw.density.alpha(SA, CT, p)  # 1/degC
            ds["beta"] = gsw.density.beta(SA, CT, p)  # kg/g
            ds["adiabatic_lapse_rate"] = (
                gsw.adiabatic_lapse_rate_from_CT(SA, CT, p) * 1e4
            )  # degC/dbar
            # derive vertical derivate variables
            ds["dtemperature_dz"] = ds["temperature"].differentiate(
                "z"
            )  # self.dz(bg["temperature"], z),
            ds["dtheta_dz"] = ds["theta"].differentiate("z")
            ds["drho_dz"] = ds["rho"].differentiate("z")
            ds["dsigma0_dz"] = ds["sigma0"].differentiate("z")
            g = 9.81  # m/s^2
            ds["N2"] = -g * (
                ds["drho_dz"] / ds["rho"] + ds["kappa"]
            )  # it's key that kappa is in decibar here

    def show_on_map(self, zdiff=None):
        if self._woa:
            ds = xr.open_dataset(self._tfile, decode_times=False).squeeze()
            #
            if zdiff is not None:
                toplt = ds["t_an"].sel(depth=zdiff[0]) - ds["t_an"].sel(depth=zdiff[1])
                title = "t(%dm) - t(%dm) [degC]" % (zdiff[0], zdiff[1])
                print("t(0m) - t(500m) global maximum = %.1f degC" % toplt.max())
            else:
                toplt = ds["t_an"].sel(depth=0)
                title = "sea surface temperature [degC]"
            #
            crs = ccrs.PlateCarree()
            plt.figure(figsize=(10, 5))
            ax = plt.axes(projection=crs)
            hdl = ax.pcolormesh(
                ds.lon, ds.lat, toplt, transform=crs, cmap=plt.get_cmap("CMRmap_r")
            )
            ax.plot(
                self.lon,
                self.lat,
                "*",
                markersize=10,
                markerfacecolor="CadetBlue",
                markeredgecolor="w",
                transform=crs,
            )
            ax.coastlines(resolution="110m")
            ax.gridlines()
            plt.colorbar(hdl, ax=ax)
            ax.set_title(title)
            plt.show()
        else:
            print("No map to show")

    def _plot_matplotlib(self, **kwargs):
        df = self.bg
        fig, axes = plt.subplots(1, 3, figsize=(9, 5), sharey=True)
        ax = axes[0]
        ax.plot(df["temperature"], df["z"], "k")
        ax.set_ylabel("z [m]")
        ax.set_title("in situ temperature [degC]")
        ax.grid()
        #
        ax = axes[1]
        ax.plot(df["salinity"], df["z"], "k")
        ax.set_title("practical salinity [psu]")
        ax.grid()
        #
        ax = axes[2]
        ax.plot(df["rho"], df["z"], "k")
        ax.set_title("density [kg/m3]")
        ax.grid()
        return fig, axes

    def _plot_hv(self, **kwargs):
        # df = self._get_df()
        df = self.now(extra=True).to_dataframe()[["temperature", "salinity", "rho"]]
        p = df.hvplot(
            grid=True,
            subplots=True,
            shared_axes=False,
            width=200,
            invert=True,
        )
        return p

    def plot(self, type="hv", **kwargs):
        # print(self)
        if type == "hv":
            return self._plot_hv(**kwargs)
        elif type == "matplotlib":
            return self._plot_matplotlib(**kwargs)

    def __repr__(self):
        return self.name

    def __getitem__(self, key):
        if self.no is None:
            self.no = self.now()
        v = self.now[key]
        if v["z"].size == 1:
            v = float(v)
        return v

    def update_eta(self, eta, t=None, detadt=0.0, d2etadt2=0.0):
        """Update isopycnal diplacement and water velocity given a function
        for isopycnal displacement and time

        Parameters
        ----------
        eta: func
            Isopycnal as a function time
        t: float
            Time in seconds
        """
        if t is not None:
            self.eta = eta(t)
            _dt = 0.1  # must be small compared to isopycnal fluctuations temporal scale
            detadt = (eta(t + _dt) - eta(t - _dt)) / _dt / 2
            d2etadt2 = (eta(t + _dt) - 2 * eta(t) + eta(t - _dt)) / _dt**2
        else:
            assert isinstance(eta, (int, float)), "eta must be a number"
            self.eta = eta
        self.detadt = detadt
        self.d2etadt2 = d2etadt2
        # self.no = self.now()
        self.no = None
