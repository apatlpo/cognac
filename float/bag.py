#!/usr/bin/env python
# coding: utf-8

# load float log

import os, sys
import pandas as pd
from glob import glob
#from datetime import *
#from load_data import *
import collections


Datasets = [
     'rosoutData', 'rosoutAggData', 'pistonStateData', 'pistonSetPointData',
	 'imuData', 'magData', 'eulerData', 'pistonVelocityData', 'pistonDistanceData', 
	 'pistonSpeedData', 'batteryData', 'sensorExtData', 'sensorIntData', 'engineData',
	 'engineCmdData', 'fixData', 'temperatureData', 'batteryFusionData', 'sensorIntFusionData',
	 'depthFusionData', 'poseFusionData', 'kalmanData', 'regulationData', 'regulationHeadingData',
	 'regulationHeadingSetPointData', 'missionData', 'safetyData', 'safetyDebugData', 
	 'safetyCpu', 'iridiumStatusData', 'iridiumSessionData', 'regulationWaypointData', 
	 'imuDebugData'
	 ]

def get_datadict():
	datadict = collections.OrderedDict()
	for d in Datasets:
		datadict[d] = eval(d[0].capitalize()+d[1:]+'()')
	return datadict

#--------------------- start copy paste from seabot/tools/rosbag/load_data.py ----------

import rospy
import rosbag
import yaml
import numpy as np
import datetime
#from numba import jit
import time

####################### Driver #######################

class SeabotData(object):
    def __init__(self, topic_name="", bag=None):
        if bag==None:
            self.nb_elements = 0
        else:
            self.nb_elements = bag.get_message_count(topic_name)
        self.time = np.empty([self.nb_elements])
        self.k=0
        self.topic_name=topic_name
    def add_time(self, t, startTime):
        _t = round((t-startTime).to_sec(), 3)
        if (_t>86400*60 or _t<1e-3) and self.k>1: 
        	# deployments should be shorter than 2 months (filter out bad values)
        	print('bad time values found')
        	_t = 2*self.time[self.k]-self.time[self.k-1]
        self.time[self.k] = _t
        self.k=self.k+1
    def to_pandas(self, startDate):
    	""" convert to dataframe
    	"""
    	attrs = [a for a in dir(self) if isinstance(getattr(self,a), np.ndarray)]
    	data = {a: [_a for _a,_t in zip(getattr(self,a), self.time)] 
    			for a in attrs}
    	index = []
    	for t in self.time:
    		# treats corrupted times
    		try:
    			index.append(startDate+pd.to_timedelta(t, unit='s'))
    		except:
    			index.append(startDate+pd.to_timedelta(-10, unit='s'))
    	df = pd.DataFrame(index=index, data=data)
    	return df

class RosoutData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/rosout", bag)
        self.level = np.empty([self.nb_elements], dtype=np.unicode_)
        self.name = np.empty([self.nb_elements], dtype=np.unicode_)
        self.msg = np.empty([self.nb_elements], dtype=np.unicode_)
        self.file = np.empty([self.nb_elements], dtype=np.unicode_)
        self.function = np.empty([self.nb_elements], dtype=np.unicode_)
        self.line = np.empty([self.nb_elements], dtype=np.unicode_)

class RosoutAggData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/rosout_agg", bag)
        self.level = np.empty([self.nb_elements], dtype=np.unicode_)
        self.name = np.empty([self.nb_elements], dtype=np.unicode_)
        self.msg = np.empty([self.nb_elements], dtype=np.unicode_)
        self.file = np.empty([self.nb_elements], dtype=np.unicode_)
        self.function = np.empty([self.nb_elements], dtype=np.unicode_)
        self.line = np.empty([self.nb_elements], dtype=np.unicode_)

class PistonStateData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/piston/state", bag)
        self.position = np.empty([self.nb_elements], dtype=np.float32)
        self.switch_out = np.empty([self.nb_elements], dtype=np.uint8)
        self.switch_in = np.empty([self.nb_elements], dtype=np.uint8)
        self.state = np.empty([self.nb_elements], dtype=np.uint16)
        self.motor_on = np.empty([self.nb_elements], dtype=np.uint8)
        self.enable_on = np.empty([self.nb_elements], dtype=np.uint8)
        self.position_set_point = np.empty([self.nb_elements], dtype=np.float32)
        self.motor_speed = np.empty([self.nb_elements], dtype=np.float32)

class PistonSetPointData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/piston/position", bag)
        self.position = np.empty([self.nb_elements], dtype=np.uint16)

class ImuData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/imu", bag)
        self.acc_x = np.empty([self.nb_elements])
        self.acc_y = np.empty([self.nb_elements])
        self.acc_z = np.empty([self.nb_elements])
        self.gyro_x = np.empty([self.nb_elements])
        self.gyro_y = np.empty([self.nb_elements])
        self.gyro_z = np.empty([self.nb_elements])

class ImuDebugData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/imu_debug", bag)
        self.accelValid = np.empty([self.nb_elements], dtype=np.uint8)
        self.fusionPoseValid = np.empty([self.nb_elements], dtype=np.uint8)
        self.fusionQPoseValid = np.empty([self.nb_elements], dtype=np.uint8)
        self.gyroValid = np.empty([self.nb_elements], dtype=np.uint8)
        self.compassValid = np.empty([self.nb_elements], dtype=np.uint8)
        self.readValid = np.empty([self.nb_elements], dtype=np.uint8)

class MagData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/mag", bag)
        self.x = np.empty([self.nb_elements])
        self.y = np.empty([self.nb_elements])
        self.z = np.empty([self.nb_elements])

class EulerData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/euler", bag)
        self.x = np.empty([self.nb_elements])
        self.y = np.empty([self.nb_elements])
        self.z = np.empty([self.nb_elements])

class PistonVelocityData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/piston/velocity", bag)
        self.velocity = np.empty([self.nb_elements])

class PistonDistanceData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/piston/distance_travelled", bag)
        self.distance = np.empty([self.nb_elements], dtype=np.uint64)

class PistonSpeedData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/piston/speed", bag)
        self.speed_in = np.empty([self.nb_elements],dtype=np.uint16)
        self.speed_out = np.empty([self.nb_elements],dtype=np.uint16)

class BatteryData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/power/battery", bag)
        self.b1 = np.empty([self.nb_elements], dtype=np.float32)
        self.b2 = np.empty([self.nb_elements], dtype=np.float32)
        self.b3 = np.empty([self.nb_elements], dtype=np.float32)
        self.b4 = np.empty([self.nb_elements], dtype=np.float32)

class SensorExtData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/sensor_external", bag)
        self.pressure = np.empty([self.nb_elements], dtype=np.float32)
        self.temperature = np.empty([self.nb_elements], dtype=np.float32)

class SensorIntData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/sensor_internal", bag)
        self.pressure = np.empty([self.nb_elements], dtype=np.float32)
        self.temperature = np.empty([self.nb_elements], dtype=np.float32)
        self.humidity = np.empty([self.nb_elements], dtype=np.float32)

class EngineData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/thruster/engine", bag)
        self.left = np.empty([self.nb_elements], dtype=np.uint8)
        self.right = np.empty([self.nb_elements], dtype=np.uint8)

class EngineCmdData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/thruster/cmd_engine", bag)
        self.linear = np.empty([self.nb_elements], dtype=np.float32)
        self.angular = np.empty([self.nb_elements], dtype=np.float32)

class FixData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/fix", bag)
        self.status = np.empty([self.nb_elements])
        self.latitude = np.empty([self.nb_elements])
        self.longitude = np.empty([self.nb_elements])
        self.altitude = np.empty([self.nb_elements])
        self.track = np.empty([self.nb_elements])
        self.speed = np.empty([self.nb_elements])
        self.gdop = np.empty([self.nb_elements])
        self.pdop = np.empty([self.nb_elements])
        self.hdop = np.empty([self.nb_elements])
        self.vdop = np.empty([self.nb_elements])
        self.tdop = np.empty([self.nb_elements])
        self.err = np.empty([self.nb_elements])
        self.err_horz = np.empty([self.nb_elements])
        self.err_vert = np.empty([self.nb_elements])
        self.err_track = np.empty([self.nb_elements])
        self.err_speed = np.empty([self.nb_elements])
        self.err_time = np.empty([self.nb_elements])

class TemperatureData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/driver/sensor_temperature", bag)
        self.temperature = np.empty([self.nb_elements], dtype=np.float32)

####################### Fusion #######################

class BatteryFusionData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/fusion/battery", bag)
        self.b1 = np.empty([self.nb_elements], dtype=np.float32)
        self.b2 = np.empty([self.nb_elements], dtype=np.float32)
        self.b3 = np.empty([self.nb_elements], dtype=np.float32)
        self.b4 = np.empty([self.nb_elements], dtype=np.float32)

class SensorIntFusionData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/fusion/sensor_internal", bag)
        self.pressure = np.empty([self.nb_elements], dtype=np.float32)
        self.temperature = np.empty([self.nb_elements], dtype=np.float32)
        self.humidity = np.empty([self.nb_elements], dtype=np.float32)

class DepthFusionData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/fusion/depth", bag)
        self.depth = np.empty([self.nb_elements], dtype=np.float32)
        self.velocity = np.empty([self.nb_elements], dtype=np.float32)
        self.zero_depth_pressure = np.empty([self.nb_elements], dtype=np.float32)

class PoseFusionData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/fusion/pose", bag)
        self.north = np.empty([self.nb_elements], dtype=np.float64)
        self.east = np.empty([self.nb_elements], dtype=np.float64)
        self.heading = np.empty([self.nb_elements], dtype=np.float64)

class KalmanData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/fusion/kalman", bag)
        self.depth = np.empty([self.nb_elements], dtype=np.float64)
        self.volume = np.empty([self.nb_elements], dtype=np.float64)
        self.velocity = np.empty([self.nb_elements], dtype=np.float64)
        self.offset = np.empty([self.nb_elements], dtype=np.float64)
        self.chi = np.empty([self.nb_elements], dtype=np.float64)
        self.chi2 = np.empty([self.nb_elements], dtype=np.float64)
        self.cz = np.empty([self.nb_elements], dtype=np.float64)
        self.offset_total = np.empty([self.nb_elements], dtype=np.float64)
        self.cov_depth = np.empty([self.nb_elements], dtype=np.float64)
        self.cov_velocity = np.empty([self.nb_elements], dtype=np.float64)
        self.cov_offset = np.empty([self.nb_elements], dtype=np.float64)
        self.cov_chi = np.empty([self.nb_elements], dtype=np.float64)
        self.cov_chi2 = np.empty([self.nb_elements], dtype=np.float64)
        self.cov_cz = np.empty([self.nb_elements], dtype=np.float64)
        self.valid = np.empty([self.nb_elements], dtype=np.uint8)

####################### Regulation #######################

class RegulationData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/regulation/debug", bag)
        self.u = np.empty([self.nb_elements], dtype=np.float64)
        self.y = np.empty([self.nb_elements], dtype=np.float64)
        self.dy = np.empty([self.nb_elements], dtype=np.float64)
        self.set_point = np.empty([self.nb_elements], dtype=np.float32)
        self.mode = np.empty([self.nb_elements], dtype=np.uint8)

class RegulationHeadingData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/regulation/debug_heading", bag)
        self.set_point = np.empty([self.nb_elements], dtype=np.float64)
        self.error = np.empty([self.nb_elements], dtype=np.float64)
        self.p_var = np.empty([self.nb_elements], dtype=np.float64)
        self.d_var = np.empty([self.nb_elements], dtype=np.float64)
        self.command = np.empty([self.nb_elements], dtype=np.float64)
        self.command_limit = np.empty([self.nb_elements], dtype=np.float64)

class RegulationHeadingSetPointData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/regulation/heading_set_point", bag)
        self.set_point = np.empty([self.nb_elements], dtype=np.float64)

class RegulationWaypointData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/regulation/waypoint_debug", bag)
        self.yaw_set_point = np.empty([self.nb_elements], dtype=np.float64)
        self.yaw_error = np.empty([self.nb_elements], dtype=np.float64)
        self.distance_error = np.empty([self.nb_elements], dtype=np.float64)
        self.enable_regulation = np.empty([self.nb_elements], dtype=np.uint8)
        self.hysteresis_inside = np.empty([self.nb_elements], dtype=np.uint8)
        self.angular = np.empty([self.nb_elements], dtype=np.float64)
        self.angular_limit = np.empty([self.nb_elements], dtype=np.float64)
        self.valid_time = np.empty([self.nb_elements], dtype=np.uint8)

####################### Mission #######################

class MissionData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/mission/set_point", bag)
        self.north = np.empty([self.nb_elements], dtype=np.float64)
        self.east = np.empty([self.nb_elements], dtype=np.float64)
        self.depth = np.empty([self.nb_elements], dtype=np.float32)
        self.limit_velocity = np.empty([self.nb_elements], dtype=np.float32)
        self.approach_velocity = np.empty([self.nb_elements], dtype=np.float32)
        self.mission_enable = np.empty([self.nb_elements], dtype=np.uint8)
        self.enable_thrusters = np.empty([self.nb_elements], dtype=np.uint8)
        self.waypoint_number = np.empty([self.nb_elements], dtype=np.uint16)
        self.wall_time = np.empty([self.nb_elements], dtype=np.uint64)
        self.time_to_next_waypoint = np.empty([self.nb_elements], dtype=np.uint64)
        self.seafloor_landing = np.empty([self.nb_elements], dtype=np.uint8)

####################### Safety #######################

class SafetyData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/safety/safety", bag)
        self.published_frequency = np.empty([self.nb_elements], dtype=np.uint8)
        self.depth_limit = np.empty([self.nb_elements], dtype=np.uint8)
        self.batteries_limit = np.empty([self.nb_elements], dtype=np.uint8)
        self.depressurization = np.empty([self.nb_elements], dtype=np.uint8)
        self.seafloor = np.empty([self.nb_elements], dtype=np.uint8)
        self.piston = np.empty([self.nb_elements], dtype=np.uint8)

class SafetyDebugData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/safety/debug", bag)
        self.flash = np.empty([self.nb_elements], dtype=np.uint8)
        self.ratio_p_t = np.empty([self.nb_elements], dtype=np.float32)
        self.ratio_delta = np.empty([self.nb_elements], dtype=np.float32)
        self.volume = np.empty([self.nb_elements], dtype=np.float32)
        self.volume_delta = np.empty([self.nb_elements], dtype=np.float32)
        self.zero_depth = np.empty([self.nb_elements], dtype=np.uint8)

class SafetyCpu(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/safety/cpu", bag)
        self.cpu = np.empty([self.nb_elements], dtype=np.float32)
        self.ram = np.empty([self.nb_elements], dtype=np.float32)

####################### Iridium #######################

class IridiumStatusData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/iridium/status", bag)
        self.service = np.empty([self.nb_elements], dtype=np.uint8)
        self.signal_strength = np.empty([self.nb_elements], dtype=np.uint8)
        self.antenna = np.empty([self.nb_elements], dtype=np.uint8)

class IridiumSessionData(SeabotData):
    def __init__(self, bag=None):
        SeabotData.__init__(self, "/iridium/session", bag)
        self.mo = np.empty([self.nb_elements], dtype=np.int8)
        self.momsn = np.empty([self.nb_elements], dtype=np.int8)
        self.mt = np.empty([self.nb_elements], dtype=np.int8)
        self.mtmsn = np.empty([self.nb_elements], dtype=np.int8)
        self.waiting = np.empty([self.nb_elements], dtype=np.int8)

startTime = 0.0
end_time = 0.0

startTime = 0.0
end_time = 0.0

########################################################
####################### Function #######################

#@jit#(nopython=True)
def load_bag(filename, rosoutData, rosoutAggData, pistonStateData, pistonSetPointData, 
			 imuData, magData, eulerData, pistonVelocityData, pistonDistanceData, 
			 pistonSpeedData, batteryData, sensorExtData, sensorIntData, engineData, 
			 engineCmdData, fixData, temperatureData, batteryFusionData, sensorIntFusionData, 
			 depthFusionData, poseFusionData, kalmanData, regulationData, regulationHeadingData, 
			 regulationHeadingSetPointData, missionData, safetyData, safetyDebugData, safetyCpu, 
			 iridiumStatusData, iridiumSessionData, regulationWaypointData, imuDebugData,
			 verbose=True):
    start_time_process = time.time()
    bag = rosbag.Bag(filename, 'r')

    if verbose:
    	print(bag)

    startTime = rospy.Time.from_sec(bag.get_start_time())# + rospy.Duration(600)
    end_time = rospy.Time.from_sec(bag.get_end_time())# + rospy.Duration(100)

    rosoutData.__init__(bag)
    rosoutAggData.__init__(bag)
    pistonStateData.__init__(bag)
    pistonSetPointData.__init__(bag)
    imuData.__init__(bag)
    magData.__init__(bag)
    eulerData.__init__(bag)
    pistonVelocityData.__init__(bag)
    pistonDistanceData.__init__(bag)
    pistonSpeedData.__init__(bag)
    batteryData.__init__(bag)
    sensorExtData.__init__(bag)
    sensorIntData.__init__(bag)
    engineData.__init__(bag)
    engineCmdData.__init__(bag)
    fixData.__init__(bag)
    temperatureData.__init__(bag)
    batteryFusionData.__init__(bag)
    sensorIntFusionData.__init__(bag)
    depthFusionData.__init__(bag)
    poseFusionData.__init__(bag)
    kalmanData.__init__(bag)
    regulationData.__init__(bag)
    regulationHeadingData.__init__(bag)
    regulationHeadingSetPointData.__init__(bag)
    missionData.__init__(bag)
    safetyData.__init__(bag)
    safetyDebugData.__init__(bag)
    safetyCpu.__init__(bag)
    iridiumStatusData.__init__(bag)
    iridiumSessionData.__init__(bag)
    regulationWaypointData.__init__(bag)
    imuDebugData.__init__(bag)

    for topic, msg, t in bag.read_messages(start_time=startTime, end_time=end_time):
        if(topic==pistonSetPointData.topic_name):
            pistonSetPointData.position[pistonSetPointData.k] = msg.position
            pistonSetPointData.add_time(t,startTime)

        elif(topic==pistonStateData.topic_name):
            pistonStateData.position[pistonStateData.k] = msg.position
            pistonStateData.switch_out[pistonStateData.k] = msg.switch_out
            pistonStateData.switch_in[pistonStateData.k] = msg.switch_in
            pistonStateData.state[pistonStateData.k] = msg.state
            pistonStateData.motor_on[pistonStateData.k] = msg.motor_on
            pistonStateData.enable_on[pistonStateData.k] = msg.enable_on
            pistonStateData.position_set_point[pistonStateData.k] = msg.position_set_point
            pistonStateData.motor_speed[pistonStateData.k] = msg.motor_speed
            pistonStateData.add_time(t,startTime)

        elif(topic==rosoutData.topic_name):
            rosoutData.level[rosoutData.k] = msg.level
            rosoutData.name[rosoutData.k] = msg.name
            rosoutData.msg[rosoutData.k] = msg.msg
            rosoutData.file[rosoutData.k] = msg.file
            rosoutData.function[rosoutData.k] = msg.function
            rosoutData.line[rosoutData.k] = msg.line
            rosoutData.add_time(t,startTime)
        elif(topic==rosoutAggData.topic_name):
            rosoutAggData.level[rosoutAggData.k] = msg.level
            rosoutAggData.name[rosoutAggData.k] = msg.name
            rosoutAggData.msg[rosoutAggData.k] = msg.msg
            rosoutAggData.file[rosoutAggData.k] = msg.file
            rosoutAggData.function[rosoutAggData.k] = msg.function
            rosoutAggData.line[rosoutAggData.k] = msg.line
            rosoutAggData.add_time(t,startTime)

        elif(topic==batteryData.topic_name):
            batteryData.b1[batteryData.k] = msg.battery1
            batteryData.b2[batteryData.k] = msg.battery2
            batteryData.b3[batteryData.k] = msg.battery3
            batteryData.b4[batteryData.k] = msg.battery4
            batteryData.add_time(t,startTime)

        elif(topic==batteryFusionData.topic_name):
            batteryFusionData.b1[batteryFusionData.k] = msg.battery1
            batteryFusionData.b2[batteryFusionData.k] = msg.battery2
            batteryFusionData.b3[batteryFusionData.k] = msg.battery3
            batteryFusionData.b4[batteryFusionData.k] = msg.battery4
            batteryFusionData.add_time(t,startTime)

        elif(topic==sensorExtData.topic_name):
            sensorExtData.pressure[sensorExtData.k] = msg.pressure
            sensorExtData.temperature[sensorExtData.k] = msg.temperature
            sensorExtData.add_time(t,startTime)

        elif(topic==sensorIntData.topic_name):
            sensorIntData.pressure[sensorIntData.k] = msg.pressure
            sensorIntData.temperature[sensorIntData.k] = msg.temperature
            sensorIntData.humidity[sensorIntData.k] = msg.humidity
            sensorIntData.add_time(t,startTime)

        elif(topic==sensorIntFusionData.topic_name):
            sensorIntFusionData.pressure[sensorIntFusionData.k] = msg.pressure
            sensorIntFusionData.temperature[sensorIntFusionData.k] = msg.temperature
            sensorIntFusionData.humidity[sensorIntFusionData.k] = msg.humidity
            sensorIntFusionData.add_time(t,startTime)

        elif(topic==engineData.topic_name):
            engineData.left[engineData.k] = msg.left
            engineData.right[engineData.k] = msg.right
            engineData.add_time(t,startTime)

        elif(topic==engineCmdData.topic_name):
            engineCmdData.linear[engineCmdData.k] = msg.linear
            engineCmdData.angular[engineCmdData.k] = msg.angular
            engineCmdData.add_time(t,startTime)

        elif(topic==depthFusionData.topic_name):
            depthFusionData.depth[depthFusionData.k] = msg.depth
            depthFusionData.velocity[depthFusionData.k] = msg.velocity
            depthFusionData.zero_depth_pressure[depthFusionData.k] = msg.zero_depth_pressure
            depthFusionData.add_time(t,startTime)

        elif(topic==regulationData.topic_name):
            regulationData.u[regulationData.k] = msg.u
            regulationData.y[regulationData.k] = msg.y
            regulationData.dy[regulationData.k] = msg.dy
            regulationData.set_point[regulationData.k] = msg.piston_set_point
            regulationData.mode[regulationData.k] = msg.mode
            regulationData.add_time(t,startTime)

        elif(topic==poseFusionData.topic_name):
            poseFusionData.east[poseFusionData.k] = msg.east
            poseFusionData.north[poseFusionData.k] = msg.north
            if hasattr(msg, 'heading'):
                poseFusionData.heading[poseFusionData.k] = msg.heading
            poseFusionData.add_time(t,startTime)

        elif(topic==regulationHeadingData.topic_name):
            regulationHeadingData.error[regulationHeadingData.k] = msg.error
            regulationHeadingData.p_var[regulationHeadingData.k] = msg.p_var
            regulationHeadingData.d_var[regulationHeadingData.k] = msg.d_var
            regulationHeadingData.command[regulationHeadingData.k] = msg.command
            regulationHeadingData.command_limit[regulationHeadingData.k] = msg.command_limit
            regulationHeadingData.set_point[regulationHeadingData.k] = msg.set_point
            regulationHeadingData.add_time(t,startTime)

        elif(topic==regulationHeadingSetPointData.topic_name):
            regulationHeadingSetPointData.set_point[regulationHeadingSetPointData.k] = msg.data
            regulationHeadingSetPointData.add_time(t,startTime)

        elif(topic==missionData.topic_name):
            missionData.north[missionData.k] = msg.north
            missionData.east[missionData.k] = msg.east
            missionData.depth[missionData.k] = msg.depth
            if hasattr(msg, 'limit_velocity'):
                missionData.limit_velocity[missionData.k] = msg.limit_velocity
            else:
                missionData.limit_velocity[missionData.k] = 0
            if hasattr(msg, 'approach_velocity'):
                missionData.approach_velocity[missionData.k] = msg.approach_velocity
            else:
                missionData.approach_velocity[missionData.k] = 0
            missionData.mission_enable[missionData.k] = msg.mission_enable
            missionData.enable_thrusters[missionData.k] = msg.enable_thrusters
            missionData.waypoint_number[missionData.k] = msg.waypoint_number
            missionData.wall_time[missionData.k] = msg.wall_time
            missionData.time_to_next_waypoint[missionData.k] = msg.time_to_next_waypoint
            missionData.seafloor_landing[missionData.k] = msg.seafloor_landing
            missionData.add_time(t,startTime)

        elif(topic==fixData.topic_name):
            fixData.status[fixData.k] = msg.status
            fixData.latitude[fixData.k] = msg.latitude
            fixData.longitude[fixData.k] = msg.longitude
            fixData.altitude[fixData.k] = msg.altitude
            fixData.track[fixData.k] = msg.track
            fixData.speed[fixData.k] = msg.speed
            fixData.gdop[fixData.k] = msg.gdop
            fixData.pdop[fixData.k] = msg.pdop
            fixData.hdop[fixData.k] = msg.hdop
            fixData.vdop[fixData.k] = msg.vdop
            fixData.tdop[fixData.k] = msg.tdop
            fixData.err[fixData.k] = msg.err
            fixData.err_horz[fixData.k] = msg.err_horz
            fixData.err_vert[fixData.k] = msg.err_vert
            fixData.err_track[fixData.k] = msg.err_track
            fixData.err_speed[fixData.k] = msg.err_speed
            fixData.err_time[fixData.k] = msg.err_time
            fixData.add_time(t,startTime)

        elif(topic==pistonVelocityData.topic_name):
            pistonVelocityData.velocity[pistonVelocityData.k] = msg.velocity
            pistonVelocityData.add_time(t,startTime)

        elif(topic==pistonDistanceData.topic_name):
            pistonDistanceData.distance[pistonDistanceData.k] = msg.distance
            pistonDistanceData.add_time(t,startTime)

        elif(topic==pistonStateData.topic_name):
            pistonStateData.speed_in[pistonStateData.k] = msg.speed_in
            pistonStateData.speed_out[pistonStateData.k] = msg.speed_out
            pistonStateData.add_time(t,startTime)

        elif(topic==magData.topic_name):
            if(type(msg).__name__ == "_sensor_msgs__MagneticField"):
                magData.x[magData.k] = msg.magnetic_field.x
                magData.y[magData.k] = msg.magnetic_field.y
                magData.z[magData.k] = msg.magnetic_field.z          
            else:
                magData.x[magData.k] = msg.x
                magData.y[magData.k] = msg.y
                magData.z[magData.k] = msg.z
            magData.add_time(t,startTime)

        elif(topic==imuData.topic_name):
            imuData.acc_x[imuData.k] = msg.linear_acceleration.x
            imuData.acc_y[imuData.k] = msg.linear_acceleration.y
            imuData.acc_z[imuData.k] = msg.linear_acceleration.z
            imuData.gyro_x[imuData.k] = msg.angular_velocity.x
            imuData.gyro_y[imuData.k] = msg.angular_velocity.y
            imuData.gyro_z[imuData.k] = msg.angular_velocity.z
            imuData.add_time(t, startTime)

        elif(topic==imuDebugData.topic_name):
            if(msg.accelValid):
                imuDebugData.accelValid[imuDebugData.k] = 1.0
            else:
                imuDebugData.accelValid[imuDebugData.k] = 0.0
            if(msg.fusionPoseValid):
                imuDebugData.fusionPoseValid[imuDebugData.k] = 1.0
            else:
                imuDebugData.fusionPoseValid[imuDebugData.k] = 0.0
            if(msg.fusionQPoseValid):
                imuDebugData.fusionQPoseValid[imuDebugData.k] = 1.0
            else:
                imuDebugData.fusionQPoseValid[imuDebugData.k] = 0.0
            if(msg.gyroValid):
                imuDebugData.gyroValid[imuDebugData.k] = 1.0
            else:
                imuDebugData.gyroValid[imuDebugData.k] = 0.0
            if(msg.compassValid):
                imuDebugData.compassValid[imuDebugData.k] = 1.0
            else:
                imuDebugData.compassValid[imuDebugData.k] = 0.0
            if(msg.readValid):
                imuDebugData.readValid[imuDebugData.k] = 1.0
            else:
                imuDebugData.readValid[imuDebugData.k] = 0.0
            imuDebugData.add_time(t, startTime)

        elif(topic==eulerData.topic_name):
            eulerData.x[eulerData.k] = msg.x
            eulerData.y[eulerData.k] = msg.y
            eulerData.z[eulerData.k] = msg.z
            eulerData.add_time(t, startTime)

        elif(topic==temperatureData.topic_name):
            temperatureData.temperature[temperatureData.k] = msg.temperature
            temperatureData.add_time(t, startTime)

        elif(topic==kalmanData.topic_name):
            kalmanData.depth[kalmanData.k] = msg.depth
            kalmanData.velocity[kalmanData.k] = msg.velocity
            kalmanData.offset[kalmanData.k] = msg.offset
            kalmanData.chi[kalmanData.k] = msg.chi
            kalmanData.chi2[kalmanData.k] = msg.chi2
            kalmanData.cz[kalmanData.k] = msg.cz
            kalmanData.offset_total[kalmanData.k] = msg.offset_total
            kalmanData.cov_depth[kalmanData.k] = msg.variance[0]
            kalmanData.cov_velocity[kalmanData.k] = msg.variance[1]
            kalmanData.cov_offset[kalmanData.k] = msg.variance[2]
            kalmanData.cov_chi[kalmanData.k] = msg.variance[3]
            kalmanData.cov_chi2[kalmanData.k] = msg.variance[4]
            kalmanData.cov_cz[kalmanData.k] = msg.variance[5]
            if hasattr(msg, 'valid'):
                if(msg.valid):
                    kalmanData.valid[kalmanData.k] = 1
                else:
                    kalmanData.valid[kalmanData.k] = 0
            kalmanData.add_time(t, startTime)


        elif(topic==safetyData.topic_name):
            if(msg.published_frequency):
                safetyData.published_frequency[safetyData.k] = 1
            else:
                safetyData.published_frequency[safetyData.k] = 0
            if(msg.depth_limit):
                safetyData.depth_limit[safetyData.k] = 1
            else:
                safetyData.depth_limit[safetyData.k] = 0
            if(msg.batteries_limit):
                safetyData.batteries_limit[safetyData.k] = 1
            else:
                safetyData.batteries_limit[safetyData.k] = 0
            if(msg.depressurization):
                safetyData.depressurization[safetyData.k] = 1
            else:
                safetyData.depressurization[safetyData.k] = 0
            if(msg.seafloor):
                safetyData.seafloor[safetyData.k] = 1
            else:
                safetyData.seafloor[safetyData.k] = 0
            if hasattr(msg, 'piston'):
                if(msg.piston):
                    safetyData.piston[safetyData.k] = 1
                else:
                    safetyData.piston[safetyData.k] = 0
            safetyData.add_time(t, startTime)

        elif(topic==safetyCpu.topic_name):
            safetyCpu.cpu[safetyCpu.k] = msg.cpu
            safetyCpu.ram[safetyCpu.k] = msg.ram
            safetyCpu.add_time(t, startTime)

        elif(topic==safetyDebugData.topic_name):
            if(msg.flash):
                safetyDebugData.flash[safetyDebugData.k] = 1
            else:
                safetyDebugData.flash[safetyDebugData.k] = 0
            safetyDebugData.ratio_p_t[safetyDebugData.k] = msg.ratio_p_t
            safetyDebugData.ratio_delta[safetyDebugData.k] = msg.ratio_delta
            safetyDebugData.volume[safetyDebugData.k] = msg.volume
            safetyDebugData.volume_delta[safetyDebugData.k] = msg.volume_delta
            if hasattr(msg, 'zero_depth'):
                if(msg.zero_depth):
                    safetyDebugData.zero_depth[safetyDebugData.k] = 1
                else:
                    safetyDebugData.zero_depth[safetyDebugData.k] = 0
            safetyDebugData.add_time(t, startTime)

        elif(topic==iridiumStatusData.topic_name):
            iridiumStatusData.service[iridiumStatusData.k] = msg.service
            iridiumStatusData.signal_strength[iridiumStatusData.k] = msg.signal_strength
            iridiumStatusData.antenna[iridiumStatusData.k] = msg.antenna
            iridiumStatusData.add_time(t, startTime)

        elif(topic==iridiumSessionData.topic_name):
            iridiumSessionData.mo[iridiumSessionData.k] = msg.mo
            iridiumSessionData.momsn[iridiumSessionData.k] = msg.momsn
            iridiumSessionData.mt[iridiumSessionData.k] = msg.mt
            iridiumSessionData.mtmsn[iridiumSessionData.k] = msg.mtmsn
            iridiumSessionData.waiting[iridiumSessionData.k] = msg.waiting
            iridiumSessionData.add_time(t, startTime)

        elif(topic==regulationWaypointData.topic_name):
            regulationWaypointData.yaw_set_point[regulationWaypointData.k] = msg.yaw_set_point
            regulationWaypointData.yaw_error[regulationWaypointData.k] = msg.yaw_error
            regulationWaypointData.distance_error[regulationWaypointData.k] = msg.distance_error
            if(msg.enable_regulation):
                regulationWaypointData.enable_regulation[regulationWaypointData.k] = 1
            else:
                regulationWaypointData.enable_regulation[regulationWaypointData.k] = 0
            if(msg.hysteresis_inside):
                regulationWaypointData.hysteresis_inside[regulationWaypointData.k] = 1
            else:
                regulationWaypointData.hysteresis_inside[regulationWaypointData.k] = 0
            if(msg.valid_time):
                regulationWaypointData.valid_time[regulationWaypointData.k] = 1
            else:
                regulationWaypointData.valid_time[regulationWaypointData.k] = 0
            regulationWaypointData.angular[regulationWaypointData.k] = msg.angular
            regulationWaypointData.angular_limit[regulationWaypointData.k] = msg.angular_limit
            regulationWaypointData.add_time(t, startTime)

    bag.close()

    # if(len(time_regulation_depth_set_point)>0):
    #   time_regulation_depth_set_point.append((end_time-startTime).to_sec())
    #   regulation_depth_set_point.append(regulation_depth_set_point[-1])

    # Data Analysis
    if(len(magData.x)>0) and verbose:
        print("compass_min = ", min(magData.x), min(magData.y), min(magData.z))
        print("compass_max = ", max(magData.x), max(magData.y), max(magData.z))
        print("acc_min = ", min(imuData.acc_x), min(imuData.acc_y), min(imuData.acc_z))
        print("acc_max = ", max(imuData.acc_x), max(imuData.acc_y), max(imuData.acc_z))

        print("gyro_mean = ", max(imuData.gyro_x), max(imuData.gyro_y), max(imuData.gyro_z))

    done_time_process = time.time()
    if verbose:
    	print("  LOAD TIME: %.1f"%(done_time_process-start_time_process))

    return startTime

#--------------------- end copy paste from seabot/tools/rosbag/load_data.py ----------


if __name__ == '__main__':

	if len(sys.argv)==2:
		bagfiles = [sys.argv[1]]
	else:
		bagfiles = glob('*.bag')

	for bagfile in bagfiles:

		print('---')
		print('start processing %s'%bagfile)

		data = get_datadict()
		startTime = load_bag(bagfile, *data.values(), verbose=False)
		startDate = pd.to_datetime(startTime.to_time(), unit='s')
		print('  start time: %s'%startDate) # UTC

		out_dir = './pd_'+bagfile.split('/')[-1].split('.')[0]
		print('  store pandas files in %s'%out_dir)
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)
		
		data['pistonStateData'].to_pandas(startDate)
		data_pd = {d: data[d].to_pandas(startDate) for d in data}

		for d in data_pd:
			data_pd[d].to_csv(out_dir+'/'+d.replace('Data','')+'.csv')
