#!/usr/bin/env python
# coding: utf-8

# load float log
bagfile = '/export/home/ahoudevi/log_float/2019-04-16-12-47-31_0.bag'

#import load_data
#from rosbag_pandas import bag_to_dataframe
#bg = bag_to_dataframe('/export/home/ahoudevi/log_float/2019-04-16-12-47-31_0.bag')

from load_data import *
load_bag(bagfile)


