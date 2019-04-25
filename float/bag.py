#!/usr/bin/env python
# coding: utf-8

# load float log

import os
import pandas as pd
#from datetime import *
from load_data import *


bagfile = '/export/home/ahoudevi/log_float/2019-04-16-12-47-31_0.bag'
startTime = load_bag(bagfile)

startDate = datetime.datetime.fromtimestamp(startTime.to_time())

out_dir = './'+bagfile.split('/')[-1].split('.')[0]
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    
    
def store(topic, time, L):
    df = pd.DataFrame(index=[startDate+datetime.timedelta(t) for t in time],
                      data={l: globals()[l] for l in L})
    df.to_pickle(out_dir+'/'+topic+'.p')
    

#
store('piston_position', time_piston_position, ['piston_position'])

#
store('piston_state', time_piston_state,
      ['piston_state_position', 'piston_state_switch_out', 'piston_state_switch_in',
       'piston_state_state', 'piston_state_motor_on', 'piston_state_enable_on',
       'piston_state_position_set_point', 'piston_state_motor_speed'])



