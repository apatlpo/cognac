
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from netCDF4 import Dataset
import gsw

import matplotlib.pyplot as plt
import cartopy.crs as ccrs


# 
g=9.81

#
def t_modulo_dt(t,dt,threshold):
    if np.abs(t/dt-np.rint(t/dt)) < threshold:
        return True
    else:
        return False

def interp(z_in,v_in,z_out):
    return interp1d(z_in,v_in,kind='linear',fill_value='extrapolate')(z_out)

#
class autonomous_float():
    
    def __init__(self,**kwargs):
        # default parameters
        params = {'a': 0.05, 'L': 0.4, 'gamma': 2.e-6, 'alpha': 7.e-5, 'temp0': 15.}
        params['m']= 1000. * np.pi * params['a']**2 * params['L']
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
        strout+='  temp0 = %.2e /degC  - reference temperature\n'%(self.temp0)
        if hasattr(self,'piston'):
            strout+=str(self.piston)
        return strout
    
    def rho(self,p,temp,v=None):
        ''' Returns float density i.e. mass over volume
        '''
        if v is None:
            if hasattr(self,'v'):
                v = self.v
            else:
                v = 0.
        return self.m/(self.V*(1.-self.gamma*p+self.alpha*(temp-self.temp0))+v)

    def adjust_balast(self,p_eq,temp_eq,rho_eq):
        ''' Find volume that needs to be added in order to be at equilibrium 
            prescribed pressure, temperature, density
        '''
        v = fsolve(lambda x: rho_eq-self.rho(p_eq,temp_eq,v=x),0.)
        return v

    def adjust_m(self,p_eq,temp_eq,rho_eq):
        ''' Find mass that needs to be added in order to be at equilibrium 
            prescribed pressure, temperature, density
        '''
        def func(m):
            self.m = m
            return rho_eq-self.rho(p_eq,temp_eq)
        m0=self.m
        self.m = fsolve(func,self.m)
        print('%.1f g were added to the float in order to be at equilibrium at %.0f dbar \n'%((self.m-m0)*1.e3,p_eq))
    
    def time_step(self,waterp,T=600.,dt_step=1.,dt_plt=None,dt_store=60., \
                  z=0.,w=0.,v=None,piston=False,t0=0.,Lv=None, \
                  z_target=0., w_target=0., dwdt_target=0., tau_ctrl=60., dt_ctrl=None, **kwargs):
        ''' Time step the float position given initial conditions
        '''
        t=t0
        self.z = z
        self.w = w
        if piston:
            self.v=self.piston.v
        elif v is None:
            if not hasattr(self,'v'):
                self.v=0.
        else:
            self.v=v
        #
        if Lv is None:
            Lv=self.L
        self.Lv=Lv
        #
        if dt_store is not None:
            threshold_store = 0.5* dt_step/dt_store
        #
        if dt_ctrl is not None:
            threshold_ctrl = 0.5* dt_step/dt_ctrl
        else:
            dt_ctrl = dt_step
            threshold_ctrl = 0.5* dt_step/dt_ctrl
        #
        if hasattr(self,'X'):
            delattr(self,'X')
        #
        print('Start time stepping for %d min ...'%(T/60.))
        #
        _f=0.
        while t<t0+T:
            #
            p, tempw = waterp.get_p(self.z), waterp.get_temp(self.z)
            rhow = waterp.get_rho(self.z)
            rhof = self.rho(p,tempw)
            # store
            if (dt_store is not None) and t_modulo_dt(t,dt_store,threshold_store):
                ndata = [t,self.z,self.w,_f,rhof,rhow]
                if not hasattr(self,'X'):
                    self.X=np.array(ndata)
                else:
                    self.X=np.vstack((self.X,ndata))
                #print('Data stored at t=%.0f'%t)
                #print('z=%.2f'%self.z)
            # plot
            # get vertical force on float
            _f = self._f(p,rhof,rhow,Lv)            
            # decide wether piston will be adjusted
            if piston and t_modulo_dt(t,dt_ctrl,threshold_ctrl):
                if self._control(self.z,self.w,_f/self.m,z_target,w_target,dwdt_target,tau_ctrl)>0:
                    self.piston.push()
                else:
                    self.piston.pull()
            # update
            self.z += dt_step*self.w
            self.z = np.amin((self.z,0.))
            self.w += dt_step*_f/self.m
            if piston:
                self.v = self.piston.v
            t+=dt_step
        print('... time stepping done')
        
    def _f(self,p,rhof,rhow,Lv):
        ''' Compute the vertical force exterted on the float
        '''
        f = -self.m*g
        f += self.m*rhow/rhof*g # Pi0, we ignore DwDt terms for now
        f += -self.m/Lv*np.abs(self.w)*self.w # Pi'
        return f
    
    def init_piston(self,vmin,vmax,dv=None,v=None):
        if dv is None:
            dv = (vmax-vmin)/10.
        self.piston = piston(vmin,vmax,dv,v=v)
    
    def compute_bounds(self,waterp,zmin,zmax=0.,Lv=None):
        # compute approximate bounds on velocity and acceleration
        self.w=0.
        if Lv is None:
            if hasattr(self,'Lv'):
                Lv=self.Lv
            else:
                Lv=self.L
        #
        z=zmax
        p, tempw = waterp.get_p(z), waterp.get_temp(z)
        rhow = waterp.get_rho(z)
        if hasattr(self,'piston'):
            rhof = self.rho(p,tempw,v=self.piston.vmax)
        else:
            rhof = self.rho(p,tempw)
        fmax=self._f(p,rhof,rhow,Lv) # upward force
        #
        z=zmin
        p, tempw = waterp.get_p(z), waterp.get_temp(z)
        rhow = waterp.get_rho(z)
        if hasattr(self,'piston'):
            rhof = self.rho(p,tempw,v=self.piston.vmin)
        else:
            rhof = self.rho(p,tempw)
        fmin=self._f(p,rhof,rhow,Lv) # downward force
        #
        afmax = np.amax((np.abs(fmax),np.abs(fmin)))
        wmax = np.sqrt( afmax * Lv/self.m)
        print('Acceleration and velocity bounds (zmin=%.0fm,zmax=%.0fm):' %(zmin,zmax))
        print('fmax/m=%.1e m^2/s, fmin/m= %.1e m^2/s, wmax= %.1f cm/s' %(fmax/self.m,fmin/self.m,wmax*100.) )
        print('For accelerations, equivalent speed reached in 1min:')
        print('  fmax %.1e cm/s, fmin/m= %.1e cm/s' %(fmax/self.m*60.*100.,fmin/self.m*60.*100.) )
        return fmax, fmin, afmax, wmax
        
    def _control(self,z,w,dwdt,z_target,w_target,dwdt_target,tau_ctrl):
        ''' Several inputs are required:
        (z_target,w_target,dwdt_target) - describes the trajectory
        tau_ctrl  - a time scale of control 
        '''
        return np.sign( (dwdt_target - dwdt)*tau_ctrl**2 + 2*(w_target-w)*tau_ctrl + z_target-z )
        #return np.sign( z_target-z )

    
class piston():
    
    def __init__(self,vmin,vmax,dv,v=None):
        self.vmin=vmin
        self.vmax=vmax
        self.dv=dv
        if v is None:
            self.v = vmax
        else:
            self.v = v

    def __repr__(self):
        strout='Piston parameters and state: \n'
        strout+='  vmin = %.2f cm^3      - min volume\n'%(self.vmin*1.e6)
        strout+='  vmax = %.2f cm^3      - max volume\n'%(self.vmax*1.e6)
        strout+='  dv   = %.2f cm^3      - volume increments\n'%(self.dv*1.e6)
        strout+='  v    = %.2f cm^3      - present volume addition\n'%(self.v*1.e6)
        return strout
    
    def push(self,dv=None):
        if dv is None:
            self.v+=self.dv
        
    def pull(self,dv=None):
        if dv is None:
            self.v+=-self.dv

#
class waterp():
    
    def __init__(self,lon=None,lat=None):
        if lon is not None and lat is not None:
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
            #

    
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

    def get_temp(self,z,eta=0.):
        ''' get in situ temperature
        '''
        return interp(self.z,self.temp,z)

    def get_s(self,z,eta=0.):
        ''' get practical salinity
        '''
        return interp(self.z,self.s,z)

    def get_p(self,z,eta=0.):
        ''' get pressure
        '''
        return interp(self.z,self.p,z)

    def get_theta(self,z,eta=0.):
        ''' get potential temperature
        '''
        pass
    
    def get_rho(self,z,eta=0.,ignore_temp=False):
        s = self.get_s(z,eta=eta)
        p = self.get_p(z,eta=eta)
        SA = gsw.SA_from_SP(s, p, self.lon, self.lat)
        #
        temp = self.get_temp(z,eta=eta)
        if ignore_temp:
            temp[:]=self.temp[0]
            print('Uses a uniform temperature in water density computation, temp= %.1f degC' %self.temp[0])
        CT = gsw.CT_from_t(SA, temp, p)
        #
        return gsw.density.rho(SA, CT, p)



