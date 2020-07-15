import numpy as np

_default = {'deployment': {'depth': 500,
                           'T':30, 
                           'drho':5,
                           },
            'hull': {'length':.8, 
                     'radius':.08,
                     'thickness':0.005,
                     'density':2700,
                     },
            'piston': {'length':.2, 
                       'radius':.019/2., 
                       'density':2700,
                       'efficiency':.1,
                       'speed':.2/3600.,
                       },
            'electronics': {'volume':.2*.05*.05, 
                            'mass':200., 
                            'c':1.,
                            },
            'batteries': {'mass':3.,
                          'density':96e-3/(np.pi*.0332**2*.0615),
                          'e':32.4*3600/96e-3,
                          }
            }

rho0 = 1030 # kg/m3
g = 9.81 # m/s^2

#class part(object):
#    
#    def __init__(self, name, attrs):

class float(object):
    
    def __init__(self,
                 deployment={},
                 hull={},
                 piston={},
                 electronics={},
                 batteries={},
                 ):
        for k, v in _default.items():
            setattr(self, k, v)
        self.parts = list(_default.keys())
        self.update_deployment(**deployment)
        self.update_hull(**hull)
        self.update_piston(**piston)
        self.update_electronics(**electronics)
        self.update_batteries(**batteries)
    
    def __getitem__(self,item):
        assert item in self.parts
        return getattr(self,item)
    
    def __repr__(self):
        out = []
        for c in self.parts:
            out.append('-- {}'.format(c))
            for k, v in self[c].items():
                out.append(' {}: {:.2e}'.format(k,v))
        return '\n'.join(out)

    def update(self, **kwargs):
        for p in self.parts:
            eval('self.update_'+p+'(**kwargs['+p+'])')
    
    def update_deployment(self, 
                          **kwargs
                          ):
        """ deployment characteristics
        
        Parameters:
        -----------
        depth: float
            Deployment depth [m]
        T: float
            Deployment length [days]
        drho: float
            Water density differences encountered [kg/m^3]
        """
        for k in _default['deployment']:
            if k in kwargs and kwargs[k] is not None:
                self.deployment[k] = kwargs[k]

    def update_hull(self, 
                    **kwargs
                    ):
        """ deployment
        
        Parameters:
        -----------
        length: float
            Hull length [m]
        radius: float
            Hull radius [m]
        thickness: float
            Hull thickness [m]
        density: float
            Hull material density [kg/m^3]
        """
        for k in _default['hull']:
            if k in kwargs and kwargs[k] is not None:
                self.hull[k] = kwargs[k]
        #
        h = self.hull
        h['volume'] = np.pi*h['radius']**2*h['length']
        h['mass'] = 2.*np.pi*h['radius']*h['length']*h['thickness'] \
                    *h['density']

    def update_piston(self, 
                    **kwargs
                    ):
        """ piston
        
        Parameters:
        -----------
        length: float
            piston length [m]
        radius: float
            piston radius [m]
        density: float
            density [kg/m^3]
        efficiency: float
            piston energetic efficiency [1]
        speed: float
            piston average speed [m/s]
        """
        for k in _default['piston']:
            if k in kwargs and kwargs[k] is not None:
                self.piston[k] = kwargs[k]
        #p['length'] = np.min(self.hull('length'), p['length'])
        #p['radius'] = np.min(self.hull('radius'), p['radius'])
        #
        p = self.piston
        p['volume'] = np.pi*p['radius']**2*p['length']
        p['mass'] = p['volume']*p['density']

        
    def update_electronics(self,
                          **kwargs
                          ):
        """ electronics
        
        Parameters:
        -----------
        volume: float
            electronics volume [m^3]
        mass: float
            electronics mass [kg]
        c: float
            electric consumption [watts]
        """
        for k in _default['electronics']:
            if k in kwargs and kwargs[k] is not None:
                self.electronics[k] = kwargs[k]

    def update_batteries(self,
                         **kwargs
                         ):
        """ electronic
        
        Parameters:
        -----------
        mass: float
            battery mass [kg]
        density: float
            material density [kg/m^3]
        e: float
            enertical density [J/kg]
        
        alcaline: 1.2V*1Ah=1.2Wh pour 136g
        lithium: 3.6V*9Ah=32.4Wh pour 96g
        """
        for k in _default['batteries']:
            if k in kwargs and kwargs[k] is not None:
                self.batteries[k] = kwargs[k]

    def get_hull_mass(self):
        h = self.hull
        return 2*np.pi*h['radius']*h['e']*h['length']*h['density']
    
    def get_gamma(self):
        """ Compute minimum volume ratio that must be displace by the piston 
        """
        return self.deployment['drho']/rho0
        
    def get_piston_radius(self, gamma=None):
        if gamma is None:
            gamma = self.get_gamma()
        h, p = self.hull, self.piston
        return np.sqrt(h['length']/p['length'] * gamma) * h['radius']

    def get_piston_length(self):
        h, p = self.hull, self.piston
        gamma = self.get_gamma()
        return h['length'] * gamma * h['radius']**2/p['radius']**2
    
    def get_piston_consumption(self):
        d, h, p = self.deployment, self.hull, self.piston
        return (rho0 * g * d['depth']
                * np.pi * p['radius']**2 
                * p['speed']/p['efficiency'])
                
    def get_battery_mass(self, cp=None):
        if cp is None:
            cp = self.get_piston_consumption()
        d, h, p, e, b = (self.deployment, self.hull, self.piston,
                         self.electronics, self.batteries
                         )
        return d['T']*86400.*(e['c']+cp)/b['e']

    def get_volume(self, gamma=None, mb=None):
        if gamma is None:
            gamma = self.get_gamma()
        if mb is None:
            mb = self.get_battery_mass()
        d, h, p, e, b = (self.deployment, self.hull, self.piston,
                         self.electronics, self.batteries
                         )
        volume = ( (mb+e['mass'])
                   /(rho0 
                     - 2*h['thickness']/h['radius']*h['density'] 
                     - p['density']*gamma
                    )
                 )
        if isinstance(volume, float):
            if volume<0:
                return np.NaN
        elif isinstance(volume, np.ndarray):
            volume[np.where(volume<0)] = np.NaN
            return volume

if __name__=='__main__':
    
    d = float()
    print(d)
    
    gamma = d.get_gamma()
    print('gamma={:.2e}'.format(gamma))
    
    # need modifications in overleaf !!!
    ap = d.get_piston_radius()
    print('piston radius required: {:.2e} m'
            .format(ap))

    #Lp = d.get_piston_length(gamma)
    #print('piston length required: {:.2e} m'
    #        .format(Lp))
    
    cp = d.get_piston_consumption()
    print('piston consumption: {:.2e} W'
            .format(cp))
            
    mb = d.get_battery_mass()
    print('battery mass: {:.2e} kg'
            .format(mb))
            
    v = d.get_volume()
    print('volume: {:.2e} cm^3'
            .format(v*1e6))
            
    d = float(hull={'radius':np.arange(.05,1.,.01)})
    d.update_piston(radius=d.get_piston_radius())
    print(d.get_volume())

