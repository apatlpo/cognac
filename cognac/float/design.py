import numpy as np
#import pandas as pd

import param
from bokeh.plotting import figure

import panel as pn
#pn.extension('plotly')

from graphviz import Graph

rho0 = 1030 # kg/m3
g = 9.81 # m/s^2

class ufloat(param.Parameterized):
        
    d_depth = param.Number(500, bounds=(10,4000), step=10, 
                           label='depth', doc="Deployment depth [m]")
    d_T = param.Number(30, bounds=(1,300), step=1, 
                       label='T', doc="Deployment time length [days]")
    d_delta_rho = param.Number(5, bounds=(1,30), step=.5, 
                               label='delta_rho', doc="Deployment density delta [kg/m^3]")
    
    h_length = param.Number(0.8, bounds=(.1,2.), step=.02, 
                            label='length', doc="Hull length [m]")
    h_radius = param.Range(default=(1., 50.), bounds=(1.,100.), 
                           label='radius', doc="Hull radius [cm]")
    h_thickness = param.Number(.5, bounds=(.1,2.), step=.1, 
                               label='thickness', doc="Hull thickness [cm]")
    h_density = param.Number(2700, bounds=(1000,3000), step=50., 
                             label='density', doc="Hull density [kg/m^3]")

    p_length = param.Number(.2, bounds=(.1,1.), step=.01, 
                            label='length',doc="Piston length [m]")
    #radius = param.Number(.0095, bounds=(.001,1.), step=.001, doc="Piston radius [m]") # dependent variable
    p_density = param.Number(2700, bounds=(1000,3000), step=50., 
                             label='density', doc="Piston density [kg/m^3]")
    p_efficiency = param.Number(.1, bounds=(.01,1.), step=.01, 
                                label='efficiency', doc="Piston energetic efficiency [1]")
    p_speed = param.Number(.1, bounds=(.01,10), step=.01, label='speed',
                           doc="Piston speed [m/h]")

    e_volume = param.Number(250, bounds=(10,1000), step=10, 
                            label='volume', doc="Electronics volume [cm^3]")
    e_mass = param.Number(.400, bounds=(.1,1.), step=.05,
                          label='mass', doc="Electronics mass [kg]")
    e_c = param.Number(.1, bounds=(.01,1.), step=.01, 
                       label='conssumption', doc="Electronics conssumption [W]")

    #mass = param.Number(3., bounds=(.1,20.), step=.1, doc="Battery mass [kg]") # dependent variable
    b_cell = param.Selector(objects=["D",], label='size', doc="Battery size")
    b_lithium = param.Boolean(True, label='lithium', doc="Battery type")

    _params = ['d_depth', 'd_T', 'd_delta_rho',
               'h_length', 'h_radius', 'h_thickness', 'h_density',
               'p_length', 'p_density', 'p_efficiency', 'p_speed',
               'e_volume', 'e_mass', 'e_c', 
               'b_cell', 'b_lithium',
              ]
    _parts = ['deployment', 'hull', 'piston', 'electronics', 'battery']
    _mapping = {p[0]: p for p in _parts}
    
    def get_part(self, item):
        assert item in self._params
        return self._mapping[item.split('_')[0]]

    def __init__(self, **params):
        super(ufloat, self).__init__(**params)
                
        # init figures
        self.volume_figure=figure(plot_width=800, plot_height=300)
        #self.figure = figure(x_range=(-1, 10), y_range=(-1, 10))
        self.volume_renderer = self.volume_figure.line([], [],
                                                       line_width=3, 
                                                       legend_label='hull volume [cm^3]',
                                                      )
        #
        self.piston_figure=figure(plot_width=800, plot_height=300,
                                  x_range=self.volume_figure.x_range)
        self.piston_renderer = self.piston_figure.line([], [],
                                                       line_width=3,
                                                       legend_label='piston radius [cm]'
                                                      )
        #
        self.battery_figure=figure(plot_width=800, plot_height=300,
                                   x_range=self.volume_figure.x_range)
        self.battery_renderer = self.battery_figure.line([], [],
                                                         line_width=3, 
                                                         legend_label='battery mass [kg]',
                                                        )
        self.battery_figure.xaxis.axis_label = 'hull radius'
        #x_range=s1.x_range, 
        #from bokeh.layouts import gridplot
        #p = gridplot([[s1, s2, s3]], toolbar_location=None)
        #
        self.update()
    
    #@param.depends('b_cell', 'b_lithium')
    def b_edensity(self):
        if self.b_cell is "D":
            if self.b_lithium:
                return 32.4*3600/96e-3
            else:
                return 1.2*3600/136e-3

    @param.depends(*_params)
    def volume_view(self):
        return self.volume_figure

    @param.depends(*_params)
    def piston_view(self):
        return self.piston_figure
    
    @param.depends(*_params)
    def battery_view(self):
        return self.battery_figure
        
    @param.depends(*_params, watch=True)
    def update(self):
        #xs, ys = self._get_coords()
        
        # renormalizations
        d_T = self.d_T*86400.
        h_radius = np.linspace(self.h_radius[0], self.h_radius[1], 100) /1e2
        h_thickness = self.h_thickness/1e2
        p_speed = self.p_speed/3600.
        
        # piston radius
        gamma = self.d_delta_rho/rho0
        p_radius = np.sqrt(self.h_length/self.p_length * gamma) * h_radius
        
        # piston conssumption
        p_c = (rho0 * g * self.d_depth
                * np.pi * p_radius**2 
                * p_speed/self.p_efficiency)
    
        # battery mass
        b_mass = d_T*(self.e_c+p_c)/self.b_edensity()
        
        # hull volume
        h_volume = ( (b_mass+self.e_mass)
                     /(rho0 
                       - 2*h_thickness/h_radius*self.h_density
                       - self.p_density*gamma
                     )
                   )
        h_volume[np.where(h_volume<0)] = np.NaN

        self.volume_renderer.data_source.data.update({'x': h_radius, 'y': h_volume*1e6})
        self.piston_renderer.data_source.data.update({'x': h_radius, 'y': p_radius*1e2})
        self.battery_renderer.data_source.data.update({'x': h_radius, 'y': b_mass})
            
    def panel(self):
        w = self._widgets_panel()
        return pn.Column(pn.Row(w[0], w[1], w[2]),
                         pn.Row(w[3],w[4]),
                         self.volume_view,
                         self.piston_view,
                         self.battery_view,
                        )
    
    def _widgets_panel(self):
        return (pn.Column('### {}'.format('Deployment'),
                          pn.panel(self.param.d_depth),
                          pn.panel(self.param.d_T),
                          pn.panel(self.param.d_delta_rho)
                         ),
                pn.Column('### {}'.format('Hull'),
                          pn.panel(self.param.h_length),
                          pn.panel(self.param.h_radius),
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