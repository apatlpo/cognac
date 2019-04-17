import rospy
import rosbag
import yaml
import numpy as np
import datetime

####################### Driver #######################

# /driver/piston/position (set_point)
time_piston_position = []
piston_position = []

# /driver/piston/state
time_piston_state = []
piston_state_position = []
piston_state_position = []
piston_state_switch_out = []
piston_state_switch_in = []
piston_state_state = []
piston_state_motor_on = []
piston_state_enable_on = []
piston_state_position_set_point = []
piston_state_motor_speed = []

# /driver/piston/velocity
time_piston_velocity = []
piston_velocity = []

# /driver/piston/distance_travelled
time_piston_distance_travelled = []
piston_distance_travelled = []

# /driver/piston/speed
time_piston_speed = []
piston_speed_in = []
piston_speed_out = []

# /driver/power/battery
time_battery = []
battery1 = []
battery2 = []
battery3 = []
battery4 = []

# /driver/sensor_external
time_sensor_external = []
sensor_external_pressure = []
sensor_external_temperature = []

# /driver/sensor_internal
time_sensor_internal = []
sensor_internal_pressure = []
sensor_internal_temperature = []
sensor_internal_humidity = []

# /driver/thruster/engine
time_engine = []
engine_left = []
engine_right = []

# /driver/thruster/cmd_engine
time_cmd_engine = []
cmd_engine_linear = []
cmd_engine_angular = []

# /driver/fix
time_fix = []
fix_status = []
fix_latitude = []
fix_longitude = []
fix_altitude = []
fix_track = []
fix_speed = []
fix_gdop = []
fix_pdop = []
fix_hdop = []
fix_vdop = []
fix_tdop = []
fix_err = []
fix_err_horz = []
fix_err_vert = []
fix_err_track = []
fix_err_speed = []
fix_err_time = []

# /driver/mag
time_mag = []
mag_x = []
mag_y = []
mag_z = []

# /driver/imu
time_imu = []
acc_x = []
acc_y = []
acc_z = []
gyro_x = []
gyro_y = []
gyro_z = []

# /driver/euler
time_euler = []
euler_x = []
euler_y = []
euler_z = []

# /driver/sensor_temperature
time_sensor_temperature = []
sensor_temperature = []

####################### Fusion #######################

# /fusion/battery
time_fusion_battery = []
fusion_battery1 = []
fusion_battery2 = []
fusion_battery3 = []
fusion_battery4 = []

# /fusion/sensor_internal
time_fusion_sensor_internal = []
sensor_fusion_internal_pressure = []
sensor_fusion_internal_temperature = []

# /fusion/depth
time_fusion_depth = []
fusion_depth = []
fusion_velocity = []

# /fusion/pose
time_fusion_pose = []
fusion_pose_north = []
fusion_pose_east = []

# /fusion/kalman
time_kalman = []
kalman_depth = []
kalman_volume = []
kalman_velocity = []
kalman_offset = []
kalman_chi = []
kalman_cov_depth = []
kalman_cov_velocity = []
kalman_cov_volume = []
kalman_cov_offset = []
kalman_cov_chi = []

####################### Regulation #######################

# /regulation/debug
time_regulation_debug = []
regulation_u = []
regulation_y = []
regulation_dy = []
regulation_piston_set_point = []
regulation_mode = []

# /regulation/debug_heading
time_regulation_heading = []
regulation_heading_error = []
regulation_heading_p_var = []
regulation_heading_d_var = []
regulation_heading_command = []
regulation_heading_command_limit = []
regulation_heading_set_point = []

# /regulation/heading_set_point
time_regulation_set_point = []
regulation_set_point = []

####################### Mission #######################

time_mission = []
mission_north = []
mission_east = []
mission_depth = []
mission_velocity_depth = []
mission_mission_enable = []
mission_depth_only = []
mission_waypoint_number = []
mission_wall_time = []
mission_time_to_next_waypoint = []

####################### Safety #######################

# /safety/safety
time_safety = []
safety_published_frequency = []
safety_depth_limit = []
safety_batteries_limit = []
safety_depressurization = []
safety_seafloor = []

# /safety/debug
time_safety_debug = []
safety_debug_flash = []
safety_debug_ratio_p_t = []
safety_debug_ratio_delta = []
safety_debug_volume = []
safety_debug_volume_delta = []

####################### Iridium #######################

# /iridium/status
time_iridium_status = []
iridium_status_service = []
iridium_status_signal_strength = []
iridium_status_antenna = []

# /iridium/session
time_iridium_session = []
iridium_session_mo = []
iridium_session_momsn = []
iridium_session_mt = []
iridium_session_mtmsn = []
iridium_session_waiting = []

########################################################
####################### Function #######################

def load_bag(filename):

	bag = rosbag.Bag(filename, 'r')

	print(bag)

	startTime = rospy.Time.from_sec(bag.get_start_time())# + rospy.Duration(600)
	end_time = rospy.Time.from_sec(bag.get_end_time())# + rospy.Duration(100)

	for topic, msg, t in bag.read_messages(start_time=startTime, end_time=end_time):
		if(topic=="/driver/piston/position"):
			if(len(time_piston_position)>0):
				time_piston_position.append((t-startTime).to_sec())
				piston_position.append(piston_position[-1])	
			time_piston_position.append((t-startTime).to_sec())
			piston_position.append(msg.position)

		elif(topic=="/driver/piston/state"):
			if(len(time_piston_state)>0):
				time_piston_state.append((t-startTime).to_sec())
				piston_state_position.append(msg.position)
				piston_state_switch_out.append(piston_state_switch_out[-1])
				piston_state_switch_in.append(piston_state_switch_in[-1])
				piston_state_state.append(piston_state_state[-1])
				piston_state_motor_on.append(piston_state_motor_on[-1])
				piston_state_enable_on.append(piston_state_enable_on[-1])
				piston_state_position_set_point.append(piston_state_position_set_point[-1])
				piston_state_motor_speed.append(piston_state_motor_speed[-1])

			time_piston_state.append((t-startTime).to_sec())
			piston_state_position.append(msg.position)
			piston_state_switch_out.append(msg.switch_out)
			piston_state_switch_in.append(msg.switch_in)
			piston_state_state.append(msg.state)
			piston_state_motor_on.append(msg.motor_on)
			piston_state_enable_on.append(msg.enable_on)
			piston_state_position_set_point.append(msg.position_set_point)
			piston_state_motor_speed.append(msg.motor_speed)

		elif(topic=="/driver/power/battery"):
			time_battery.append((t-startTime).to_sec())
			battery1.append(msg.battery1)
			battery2.append(msg.battery2)
			battery3.append(msg.battery3)
			battery4.append(msg.battery4)

		elif(topic=="/fusion/battery"):
			time_fusion_battery.append((t-startTime).to_sec())
			fusion_battery1.append(msg.battery1)
			fusion_battery2.append(msg.battery2)
			fusion_battery3.append(msg.battery3)
			fusion_battery4.append(msg.battery4)

		elif(topic=="/driver/sensor_external"):
			time_sensor_external.append((t-startTime).to_sec())
			sensor_external_pressure.append(msg.pressure)
			sensor_external_temperature.append(msg.temperature)

		elif(topic=="/driver/sensor_internal"):
			time_sensor_internal.append((t-startTime).to_sec())
			sensor_internal_pressure.append(msg.pressure)
			sensor_internal_temperature.append(msg.temperature)
			sensor_internal_humidity.append(msg.humidity)

		elif(topic=="/fusion/sensor_internal"):
			time_fusion_sensor_internal.append((t-startTime).to_sec())
			sensor_fusion_internal_pressure.append(msg.pressure)
			sensor_fusion_internal_temperature.append(msg.temperature)

		elif(topic=="/driver/thruster/engine"):
			time_engine.append((t-startTime).to_sec())
			engine_left.append(msg.left)
			engine_right.append(msg.right)

		elif(topic=="/driver/thruster/cmd_engine"):
			time_cmd_engine.append((t-startTime).to_sec())
			cmd_engine_linear.append(msg.linear)
			cmd_engine_angular.append(msg.angular)

		elif(topic=="/fusion/depth"):
			time_fusion_depth.append((t-startTime).to_sec())
			fusion_depth.append(msg.depth)
			fusion_velocity.append(msg.velocity)

		elif(topic=="/regulation/debug"):
			time_regulation_debug.append((t-startTime).to_sec())
			regulation_u.append(msg.u)
			regulation_y.append(msg.y)
			regulation_dy.append(msg.dy)
			regulation_piston_set_point.append(msg.piston_set_point)
			regulation_mode.append(msg.mode)

		elif(topic=="/fusion/pose"):
			time_fusion_pose.append((t-startTime).to_sec())
			fusion_pose_east.append(msg.east)
			fusion_pose_north.append(msg.north)

		elif(topic=="/regulation/debug_heading"):
			time_regulation_heading.append((t-startTime).to_sec())
			regulation_heading_error.append(msg.error)
			regulation_heading_p_var.append(msg.p_var)
			regulation_heading_d_var.append(msg.d_var)
			regulation_heading_command.append(msg.command)
			regulation_heading_command_limit.append(msg.command_limit)
			regulation_heading_set_point.append(msg.set_point)

		elif(topic=="/regulation/heading_set_point"):
			time_regulation_set_point.append((t-startTime).to_sec())
			regulation_set_point.append(msg.data)

		elif(topic=="/mission/set_point"):
			time_mission.append((t-startTime).to_sec())
			mission_north.append(msg.north)
			mission_east.append(msg.east)
			mission_depth.append(msg.depth)
			mission_velocity_depth.append(msg.velocity_depth)
			mission_mission_enable.append(msg.mission_enable)
			mission_depth_only.append(msg.depth_only)
			mission_waypoint_number.append(msg.waypoint_number)
			mission_wall_time.append(msg.wall_time)
			mission_time_to_next_waypoint.append(msg.time_to_next_waypoint)			

		elif(topic=="/driver/fix"):
			time_fix.append((t-startTime).to_sec())
			fix_status.append(msg.status)
			fix_latitude.append(msg.latitude)
			fix_longitude.append(msg.longitude)
			fix_altitude.append(msg.altitude)
			fix_track.append(msg.track)
			fix_speed.append(msg.speed)
			fix_gdop.append(msg.gdop)
			fix_pdop.append(msg.pdop)
			fix_hdop.append(msg.hdop)
			fix_vdop.append(msg.vdop)
			fix_tdop.append(msg.tdop)
			fix_err.append(msg.err)
			fix_err_horz.append(msg.err_horz)
			fix_err_vert.append(msg.err_vert)
			fix_err_track.append(msg.err_track)
			fix_err_speed.append(msg.err_speed)
			fix_err_time.append(msg.err_time)

		elif(topic=="/driver/piston/velocity"):
			time_piston_velocity.append((t-startTime).to_sec())
			piston_velocity.append(msg.velocity)

		elif(topic=="/driver/piston/distance_travelled"):
			time_piston_distance_travelled.append((t-startTime).to_sec())
			piston_distance_travelled.append(msg.distance)

		elif(topic=="/driver/piston/speed"):
			time_piston_speed.append((t-startTime).to_sec())
			piston_speed_in.append(msg.speed_in)
			piston_speed_out.append(msg.speed_out)

		elif(topic=="/driver/mag"):
			time_mag.append((t-startTime).to_sec())
			if(type(msg).__name__ == "_sensor_msgs__MagneticField"):
				mag_x.append(msg.magnetic_field.x)
				mag_y.append(msg.magnetic_field.y)
				mag_z.append(msg.magnetic_field.z)				
			else:
				mag_x.append(msg.x)
				mag_y.append(msg.y)
				mag_z.append(msg.z)

		elif(topic=="/driver/imu"):
			time_imu.append((t-startTime).to_sec())
			acc_x.append(msg.linear_acceleration.x)
			acc_y.append(msg.linear_acceleration.y)
			acc_z.append(msg.linear_acceleration.z)
			gyro_x.append(msg.angular_velocity.x)
			gyro_y.append(msg.angular_velocity.y)
			gyro_z.append(msg.angular_velocity.z)

		elif(topic=="/driver/euler"):
			time_euler.append((t-startTime).to_sec())
			euler_x.append(msg.x)
			euler_y.append(msg.y)
			euler_z.append(msg.z)

		elif(topic=="/driver/sensor_temperature"):
			time_sensor_temperature.append((t-startTime).to_sec())
			sensor_temperature.append(msg.temperature)

		elif(topic=="/fusion/kalman"):
			time_kalman.append((t-startTime).to_sec())
			kalman_depth.append(msg.depth)
			kalman_velocity.append(msg.velocity)
			kalman_offset.append(msg.offset)
			kalman_chi.append(msg.chi)
			kalman_cov_depth.append(msg.covariance[0])
			kalman_cov_velocity.append(msg.covariance[1])
			kalman_cov_offset.append(msg.covariance[2])
			kalman_cov_chi.append(msg.covariance[3])

		elif(topic=="/safety/safety"):
			time_safety.append((t-startTime).to_sec())
			if(msg.published_frequency):
				safety_published_frequency.append(1)
			else:
				safety_published_frequency.append(0)
			
			if(msg.depth_limit):
				safety_depth_limit.append(1)
			else:
				safety_depth_limit.append(0)
			if(msg.batteries_limit):
				safety_batteries_limit.append(1)
			else:
				safety_batteries_limit.append(0)
			if(msg.depressurization):
				safety_depressurization.append(1)
			else:
				safety_depressurization.append(0)
			if(msg.seafloor):
				safety_seafloor.append(1)
			else:
				safety_seafloor.append(0)

		elif(topic=="/safety/debug"):
			time_safety_debug.append((t-startTime).to_sec())
			if(msg.flash):
				safety_debug_flash.append(1)
			else:
				safety_debug_flash.append(0)
			safety_debug_ratio_p_t.append(msg.ratio_p_t)
			safety_debug_ratio_delta.append(msg.ratio_delta)
			safety_debug_volume.append(msg.volume)
			safety_debug_volume_delta.append(msg.volume_delta)

		elif(topic=="/iridium/status"):
			time_iridium_status.append((t-startTime).to_sec())
			iridium_status_service.append(msg.service)
			iridium_status_signal_strength.append(msg.signal_strength)
			iridium_status_antenna.append(msg.antenna)

		elif(topic=="/iridium/session"):
			time_iridium_session.append((t-startTime).to_sec())
			iridium_session_mo.append(msg.mo)
			iridium_session_momsn.append(msg.momsn)
			iridium_session_mt.append(msg.mt)
			iridium_session_mtmsn.append(msg.mtmsn)
			iridium_session_waiting.append(msg.waiting)

	bag.close()

	# if(len(time_regulation_depth_set_point)>0):
	# 	time_regulation_depth_set_point.append((end_time-startTime).to_sec())
	# 	regulation_depth_set_point.append(regulation_depth_set_point[-1])


	# Data Analysis
	print("compass_min = ", min(mag_x), min(mag_y), min(mag_z))
	print("compass_max = ", max(mag_x), max(mag_y), max(mag_z))
	print("acc_min = ", min(acc_x), min(acc_y), min(acc_z))
	print("acc_max = ", max(acc_x), max(acc_y), max(acc_z))

	print("gyro_mean = ", max(gyro_x), max(gyro_y), max(gyro_z))

	# if(len(time_fix)>0):
	# 	import gpxpy
	# 	import gpxpy.gpx

	# 	gpx = gpxpy.gpx.GPX()
	# 	last_fix_time = 0.

	# 	for i in range(len(fix_latitude)):
	# 		if(abs(last_fix_time-time_fix[i])>30.):
	# 			last_fix_time = time_fix[i]
	# 			# Create first track in our GPX:
	# 			gpx_track = gpxpy.gpx.GPXTrack()
	# 			gpx.tracks.append(gpx_track)

	# 			# Create first segment in our GPX track:
	# 			gpx_segment = gpxpy.gpx.GPXTrackSegment()
	# 			gpx_track.segments.append(gpx_segment)

	# 		# Create points:
	# 		if(fix_status[i]==3):
	# 			gpx_segment.points.append(gpxpy.gpx.GPXTrackPoint(latitude=fix_latitude[i], 
	# 				longitude=fix_longitude[i],
	# 				elevation=fix_altitude[i],
	# 				time=datetime.datetime.fromtimestamp(time_fix[i]+startTime),
	# 				horizontal_dilution=fix_hdop[i]
	# 				))

	# 	file = open(filename+".gpx","w") 
	# 	file.write(gpx.to_xml()) 
	# 	file.close() 

	