
# workspaceFlotteur compilation

```
cd ~/Cognac/flotteur/workspaceFlotteur
catkin_make
# does not work because dependencies are not installed yet
#pip install catkin_pkg
sudo apt-get install libgps-dev
#
git clone https://github.com/tbesemer/RTIMULib2.git
cd /export/home/localuser/Cognac/flotteur/RTIMULib2/RTIMULib/
mkdir _build
cd _build
cmake .. -DCMAKE_INSTALL_PREFIX=../_install
mkdir ../_install
make
#
# modify src/seabot_driver/i2c_imu/CMakeLists.txt
# http://wiki.ros.org/catkin/CMakeLists.txt#Include_Paths_and_Library_Paths
include_directories(
  ${catkin_INCLUDE_DIRS}, ~/Cognac/flotteur/RTIMULib2/RTIMULib/_install/include
)
link_directories(~/Cognac/flotteur/RTIMULib2/RTIMULib/_install/lib)
#
# note: conda env interferes with ROS catkin installs, you need to turn it off
# to clean up catkin_make:
# rm -rf build/ devel/
sudo apt-get install --install-suggests libi2c-dev
sudo apt-get install libproj-dev
catkin_make
#
import yaml
import smtplib
#
sudo apt install python-pyqtgraph
```


# ROS tutorial notes

[Tutorial website](http://wiki.ros.org/ROS/Tutorials)

## creating workspace

```
source /opt/ros/kinetic/setup.bash
mkdir -p ~/catkin_ws/src
cd ~/catkin_ws/
catkin_make
```

### load workspace

```
source devel/setup.bash
# check
echo $ROS_PACKAGE_PATH
```

*worskpace vs packages* 

- Nodes: A node is an executable that uses ROS to communicate with other nodes.
- Messages: ROS data type used when subscribing or publishing to a topic.
- Topics: Nodes can publish messages to a topic as well as subscribe to a topic to receive messages.
- Master: Name service for ROS (i.e. helps nodes find each other)
- rosout: ROS equivalent of stdout/stderr
- roscore: Master + rosout + parameter server (parameter server will be introduced later) 


``` 
roscore
rosnode list 
rosnode info /rosout
rosrun turtlesim turtlesim_node
```

### topics:

```
rosrun rqt_graph rqt_graph &
rostopic -h
rostopic echo /turtle1/cmd_vel
rostopic list -h
rostopic list -v
```

A topic type is defined by the message type published on it. The type of the message sent on a topic can be determined using rostopic type. 
```
rostopic type /turtle1/cmd_vel
```

```
rostopic pub [topic] [msg_type] [args]
rostopic pub -1 /turtle1/cmd_vel geometry_msgs/Twist -- '[2.0, 0.0, 0.0]' '[0.0, 0.0, 1.8]'
rostopic pub /turtle1/cmd_vel geometry_msgs/Twist -r 1 -- '[2.0, 0.0, 0.0]' '[0.0, 0.0, -1.8]'
```

```
rostopic hz reports the rate at which data is published. 
rostopic hz /turtle1/pose
```

### services:

Services are another way that nodes can communicate with each other. 
Services allow nodes to send a request and receive a response. 

```
rosservice list
rosservice call [service] [args]
rosservice call /clear
rosservice type /spawn | rossrv show
rosservice call /spawn 2 2 0.2 ""
```

rosparam allows you to store and manipulate data on the ROS Parameter Server. The Parameter Server can store integers, floats, boolean, dictionaries, and lists.

```
rosparam list
rosparam set [param_name]
rosparam get [param_name]
rosparam set /background_r 150
rosservice call /clear
rosparam get /background_g 
```

roslaunch starts nodes as defined in a launch file. 
```
roslaunch [package] [filename.launch]
```

```
cd /export/home/localuser/Cognac/flotteur/ros/catkin_ws
source devel/setup.bash
roscd beginner_tutorials
mkdir launch
cd launch
# edit launch file, http://wiki.ros.org/ROS/Tutorials/UsingRqtconsoleRoslaunch
roslaunch beginner_tutorials turtlemimic.launch
rostopic pub /turtlesim1/turtle1/cmd_vel geometry_msgs/Twist -r 1 -- '[2.0, 0.0, 0.0]' '[0.0, 0.0, -1.8]'
rqt    (Plugins > Introspection > Node Graph:)
```

### msg/srv

msg: msg files are simple text files that describe the fields of a ROS message. They are used to generate source code for messages in different languages.

srv: an srv file describes a service. It is composed of two parts: a request and a response. 

msg files are stored in the msg directory of a package, and srv files are stored in the srv directory. 

```
rosmsg show beginner_tutorials/Num
```
