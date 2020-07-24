import numpy as np
#import pandas as pd

import param
from bokeh.plotting import figure
from bokeh.layouts import gridplot

import panel as pn
#pn.extension('plotly')

from graphviz import Graph

rho0 = 1030 # kg/m3
g = 9.81 # m/s^2

# https://en.wikipedia.org/wiki/D_battery
D_volume = np.pi*(33.2/2.)**2 * 61.5*1e-9
D_mass_lithium = 96e-3
D_mass_alcaline = 136e-3

_items = ['deployment','hull', 'piston', 'battery', 'electronics']
_colors = {'deployment': 'orange', 
           'hull': 'cadetblue',
           'piston': 'gray',
           'electronics': 'salmon',
           'battery': 'orange',
           'total': 'black',
           }

class ufloat(param.Parameterized):
    
    d_depth = param.Number(500, bounds=(10,4000), step=10, 
                           label='depth [m]', doc="Deployment depth [m]")
    d_T = param.Number(30, bounds=(1,300), step=1, 
                       label='T [days]', doc="Deployment time length [days]")
    d_delta_rho = param.Number(5, bounds=(1,30), step=.5, 
                               label='delta_rho [kg/m3]', doc="Deployment density delta [kg/m^3]")
    
    #h_length = param.Number(0.8, bounds=(.1,2.), step=.02, 
    #                        label='length [m]', doc="Hull length [m]")
    h_radius = param.Range(default=(1., 50.), bounds=(1.,100.), 
                           label='radius [cm]', doc="Hull radius [cm]")
    h_thickness = param.Number(.5, bounds=(.1,2.), step=.1, 
                               label='thickness [cm]', doc="Hull thickness [cm]")
    h_density = param.Number(2700, bounds=(1000,3000), step=50., 
                             label='density [kg/m3]', doc="Hull density [kg/m^3]")

    p_length = param.Number(.2, bounds=(.1,1.), step=.01, 
                            label='length',doc="Piston length [m]")
    #radius = param.Number(.0095, bounds=(.001,1.), step=.001, doc="Piston radius [m]") # dependent variable
    p_density = param.Number(2700, bounds=(1000,3000), step=50., 
                             label='density [kg/m3]', doc="Piston density [kg/m^3]")
    p_efficiency = param.Number(.1, bounds=(.01,1.), step=.01, 
                                label='efficiency [1]', doc="Piston energetic efficiency [1]")
    p_speed = param.Number(.1, bounds=(.01,10), step=.01, label='speed [m/h]',
                           doc="Piston speed [m/h]")

    e_volume = param.Number(250, bounds=(10,1000), step=10, 
                            label='volume [cm3]', doc="Electronics volume [cm3]")
    e_mass = param.Number(.400, bounds=(.1,1.), step=.05,
                          label='mass [kg]', doc="Electronics mass [kg]")
    e_c = param.Number(.1, bounds=(.01,1.), step=.01, 
                       label='conssumption [W]', doc="Electronics conssumption [W]")

    #mass = param.Number(3., bounds=(.1,20.), step=.1, doc="Battery mass [kg]") # dependent variable
    b_cell = param.Selector(objects=["D",], label='size', doc="Battery size")
    b_lithium = param.Boolean(True, label='lithium', doc="Battery type")

    _params = ['d_depth', 'd_T', 'd_delta_rho',
               'h_radius', 'h_thickness', 'h_density', # 'h_length', 
               'p_length', 'p_density', 'p_efficiency', 'p_speed',
               'e_volume', 'e_mass', 'e_c', 
               'b_cell', 'b_lithium',
              ]
    _scales = {'d_depth': 1., 
               'd_T': 86400., 
               'd_delta_rho': 1., #'h_length': 1.,
               'h_radius': 1e-2, 
               'h_thickness': 1e-2,
               'h_density': 1,
               'p_length': 1.,
               'p_density': 1., 
               'p_efficiency': 1.,
               'p_speed': 1/3600.,
               'e_volume': 1e-6, 
               'e_mass': 1., 
               'e_c': 1., 
               'b_cell': None, 
               'b_lithium': None,
               }
    #_parts = ['deployment', 'hull', 'piston', 'electronics', 'battery']
    _mapping = {p[0]: p for p in _items}
    
    def get_part(self, item):
        assert item in self._params
        return self._mapping[item.split('_')[0]]
    
    def get_key_and_part(self, item):
        assert item in self._params
        item_split = item.split('_')
        return '_'.join(item_split[1:]), self._mapping[item_split[0]]

    def __init__(self, **params):
        super(ufloat, self).__init__(**params)

        # will mirror data in dicts:
        self.D = {i: {} for i in _items}
        self.D['total'] = {}

        # init figures
        self.renderers = {v: {} for v in ['volume', 'length', 'mass']}
        line = lambda f, c, l: f.line([],[], 
                                      color=c,
                                      line_width=3,
                                      legend_label=l
                                      )
        _tools = 'pan,wheel_zoom,box_zoom,reset'
        f_volume = figure(y_axis_label='volume [liters]',
                          tools=_tools
                          )
        _r = self.renderers['volume']
        for i in ['hull', 'piston', 'battery', 'electronics', 'total']:
            _r[i] = line(f_volume, _colors[i], i)
        #
        f_length = figure(y_axis_label='length [cm]',
                          tools=_tools,
                          x_range=f_volume.x_range
                          )
        _r = self.renderers['length']
        for i in ['hull', 'piston']:
            _r[i] = line(f_length, _colors[i], i+' radius')
        _r['hull_length'] = f_length.line([],[],
                                          color=_colors[i], 
                                          line_dash='4 4',
                                          legend_label='hull length')
        #
        f_mass = figure(y_axis_label='mass [kg/m3]',
                        tools=_tools,
                        x_range=f_volume.x_range
                        )
        _r = self.renderers['mass']
        for i in ['battery', 'hull', 'total']:
            _r[i] = line(f_mass, _colors[i], i)
        #
        f_length.xaxis.axis_label = 'hull radius [cm]'
        f_mass.xaxis.axis_label = 'hull radius [cm]'
        #x_range=s1.x_range,
        self.figures = gridplot([[f_volume, f_mass], [f_length,None]], 
                                plot_width=400, 
                                plot_height=300,
                                #toolbar_location=None,
                                )
        #
        self.update()
    
    #@param.depends('b_cell', 'b_lithium')
    def _update_battery(self):
        """ J/kg
        alcaline: 1.2V*1Ah=1.2Wh pour 136g
        lithium: 3.6V*9Ah=32.4Wh pour 96g
        """
        if self.b_cell is "D":
            if self.b_lithium:
                density = D_mass_lithium/D_volume
                edensity = 32.4*3600/D_mass_lithium
            else:
                density = D_mass_alcaline/D_volume
                edensity = 1.2*3600/D_mass_alcaline
        self.D['battery']['density'] = density
        self.D['battery']['edensity'] = edensity

    @param.depends(*_params)
    def variables_view(self):
        return self.figures

    def _update_dict_params(self):
        # push and rescale data in central dict
        for p in self._params:
            _k, _p = self.get_key_and_part(p)
            _d = self.D[_p]
            if p=='h_radius':
                _d['radius'] = np.linspace(self.h_radius[0], 
                                           self.h_radius[1], 100) /1e2
            elif self._scales[p]:
                _d[_k] = getattr(self, p) * self._scales[p]

    def _update_line(self, renderer_key, item, variable, x, scale):
        _d = self.D[item][variable]
        if isinstance(_d, float):
            _d = np.ones_like(x)*_d
        (self.renderers[renderer_key][item]
         .data_source
         .data
         .update({'x': x*1e2, 
                  'y': _d*scale,
                  }
                )
         )

    @param.depends(*_params, watch=True)
    def update(self):
        
        # update data containers
        self._update_dict_params()
        self._update_battery()
        
        d, h, p, b, e = (self.D[i] for i in _items)
        t = self.D['total']

        # renormalizations
        d_T = self.d_T*86400.
        h_radius = np.linspace(self.h_radius[0], self.h_radius[1], 100) /1e2
        h_thickness = self.h_thickness/1e2
        p_speed = self.p_speed/3600.
        e_volume = self.e_volume*1e-6

        # solve for volume first
        gamma = d['delta_rho']/rho0
        sigma = (gamma/p['length'] 
                 *rho0*g*d['depth'] 
                 *p['speed']/p['efficiency']
                 )
        h['volume'] = ( (d['T']*e['c']/b['edensity'] + e['mass'])
                        /( rho0 
                          - 2*h['thickness']/h['radius']*h['density']
                          - p['density'] * gamma
                          - d['T']*sigma/b['edensity']
                          )
                       )
        h['volume'][np.where(h['volume']<0)] = np.NaN
        
        # propagate volume 
        h['length'] = h['volume']/(np.pi*h['radius']**2)
        
        # piston radius
        #gamma = self.d_delta_rho/rho0
        #p_radius = np.sqrt(self.h_length/self.p_length * gamma) * h_radius
        p['radius'] = np.sqrt(h['length']/p['length'] * gamma) * h['radius']
        
        # piston conssumption
        p['c'] =  (rho0 * g * d['depth']
                   * np.pi * p['radius']**2 
                   * p['speed']/p['efficiency']
                   )
    
        # battery mass
        #self.update_battery()
        b['mass'] = d['T']*(e['c']+p['c'])/b['edensity']
        b['volume'] = b['mass']/b['density']
        
        
        # other parameters
        h['mass'] = h['volume'] * h['density']
        p['volume'] = np.pi*p['radius']**2*p['length']
        p['mass'] = p['density']*p['volume']
        t['mass'] = h['mass'] + p['mass'] + b['mass'] + e['mass']
        t['volume'] = p['volume'] + b['volume'] + e['volume']

        # update plots
        for _p in ['hull', 'piston', 'battery', 'electronics', 'total']:
            self._update_line('volume', _p, 'volume', h['radius'], 1e3)
        #
        for _p in ['hull', 'piston', ]:
            self._update_line('length', _p, 'radius', h['radius'], 1e2)
        (self.renderers['length']['hull_length']
         .data_source.data.update({'x': h['radius']*1e2, 
                                   'y': self.D['hull']['length']*1e2}
                                  )
         )
        #
        for _p in ['hull', 'battery', 'total']:
            self._update_line('mass', _p, 'mass', h['radius'], 1)

    def panel(self):
        w = self._widgets_panel()
        return pn.Column(pn.Row(w[0], w[1], w[2]),
                         pn.Row(w[3],w[4]),
                         self.variables_view,
                        )
    
    def _widgets_panel(self):
        return (pn.Column('### {}'.format('Deployment'),
                          pn.panel(self.param.d_depth),
                          pn.panel(self.param.d_T),
                          pn.panel(self.param.d_delta_rho)
                         ),
                pn.Column('### {}'.format('Hull'),
                          pn.panel(self.param.h_radius), #pn.panel(self.param.h_length),
                          pn.panel(self.param.h_thickness),
                          pn.panel(self.param.h_density),
                         ),
               pn.Column('### {}'.format('Piston'),
                          pn.panel(self.param.p_length),
                          pn.panel(self.param.p_density),
                          pn.panel(self.param.p_efficiency),
                          pn.panel(self.param.p_speed),
                        ),
               pn.Column('### {}'.format('Electronics'),
                          pn.panel(self.param.e_volume),
                          pn.panel(self.param.e_mass),
                          pn.panel(self.param.e_c),
                        ),
               pn.Column('### {}'.format('Battery'),
                          pn.panel(self.param.b_cell),
                          pn.panel(self.param.b_lithium),
                        ),
               )
    

_part_colors = {'deployment': 'orange', 
                'hull': 'cadetblue',
                'piston': 'lightgreen',
                'electronics': 'salmon',
                'battery': 'lightgrey',
               }

def build_graph(f, name='float design'):
        
    # central graph
    g = Graph(name)
    #g.attr(compound='true') # to make subgraph: https://github.com/xflr6/graphviz/blob/master/examples/notebook.ipynb
    g.attr(rankdir='RL', size='15,15')        

    params = [p for p in f._params if p not in ['b_lithium', 'b_cell']]
    for p in params:
        part = f.get_part(p)
        color = _part_colors[part]
        g.attr('node', shape='ellipse', style='filled', color=color)
        g.node(part+' '+'_'.join(p.split('_')[1:]))
    
    # gamma
    #var = 'gamma'
    #d.append(var)
    #g.attr('node', shape='diamond', style='filled', color=self.deployment.color)
    #g.node(var)
    #g.edge(var, 'deployment delta_rho')

    # piston radius
    var = 'piston radius'
    g.attr('node', shape='diamond', style='filled', color=_part_colors['piston'])
    g.node(var)
    for v in ['hull length', 'hull radius', 'piston length', 'deployment delta_rho']:
        g.edge(var, v)

    # piston conssumption
    var = 'piston conssumption'
    g.attr('node', shape='diamond', style='filled', color=_part_colors['piston'])
    g.node(var)
    for v in ['deployment depth', 'piston radius', 'piston speed', 'piston efficiency']:
        g.edge(var, v)

    # battery mass
    var = 'battery mass'
    g.attr('node', shape='diamond', style='filled', color=_part_colors['battery'])
    g.node(var)
    for v in ['piston conssumption', 'battery edensity', 'deployment T', 'electronics c']:
        g.edge(var, v)

    # volume
    var = 'hull volume'
    g.attr('node', shape='diamond', style='filled', color=_part_colors['hull'])
    g.node(var)
    for v in ['deployment delta_rho',  'hull density', 'hull radius', 'hull thickness',
              'piston density',
              'battery mass', 
              'electronics mass',
             ]:
        g.edge(var, v)
    
    return g