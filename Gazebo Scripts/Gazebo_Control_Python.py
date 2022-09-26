#!/usr/bin/env python

import rospy 
from tf.transformations import quaternion_from_euler,euler_from_quaternion

# from std_msgs.msg import Header 
from std_msgs.msg import String,Float64
from sensor_msgs.msg import Imu,JointState
from gazebo_msgs.msg import ModelStates 


import sensor_msgs as mm 

import numpy as np 

class Server:
    def __init__(self):
        self.Model_Data = None
        self.IMU_Data = None
        self.Angular_Velocity = None 
    def IMU_callback(self, msg):
        self.IMU_Data = msg
        Angles = euler_from_quaternion([
                                self.IMU_Data.orientation.x,
                                self.IMU_Data.orientation.y,
                                self.IMU_Data.orientation.z,
                                self.IMU_Data.orientation.w])

        self.roll = Angles[0]
        self.pitch = Angles[1]
        self.yaw = Angles[2]
        self.roll_dot = self.IMU_Data.angular_velocity.x
        self.pitch_dot = self.IMU_Data.angular_velocity.y
        self.yaw_dot = self.IMU_Data.angular_velocity.z

        self.control()

    def modelState_callback(self, msg):
        self.Model_Data = msg
        # print("Pos")
        self.X = self.Model_Data.pose[1].position.x
        self.Y = self.Model_Data.pose[1].position.y
        self.Z = self.Model_Data.pose[1].position.z
        self.X_dot = self.Model_Data.twist[1].linear.x
        self.Y_dot = self.Model_Data.twist[1].linear.y
        self.Z_dot = self.Model_Data.twist[1].linear.z 
        self.V = np.sqrt(self.X_dot**2+self.Y_dot**2)
        self.L = np.sqrt(self.X**2+self.Y**2)

        self.control()

    def Angular_velocity_callback(self,msg):
        self.Angular_Velocity = msg
        self.control()
    def jointState_callback(self,msg):
        self.joint_data  = msg
        # print(msg)

    def control(self):
        if self.Model_Data is not None and self.IMU_Data is not None:
            self.states = np.transpose([self.L,self.yaw,self.pitch,self.V,self.yaw_dot,self.pitch_dot])
            print("States")
            print(np.transpose([self.L,np.rad2deg(self.yaw),np.rad2deg(self.pitch),self.V,self.yaw_dot,self.pitch_dot]))
            K = [[-0.5477,-0.5477,-22.4256,-1.5618,-0.7051,-4.4735],
                [-0.5477,0.5477,-22.4256,-1.5618,0.7051,-4.4735]]
            K = np.asarray(K)
            Torques = np.dot(K,self.states)
            print("Torques")
            print(Torques)
            left_pub.publish(Torques[0])
            right_pub.publish(Torques[1])

    
if __name__=="__main__":

    left_pub = rospy.Publisher('/teeterbot/left_torque_cmd', Float64,queue_size=10)
    right_pub = rospy.Publisher('/teeterbot/right_torque_cmd', Float64,queue_size=10)

    rospy.init_node("Teerbot",anonymous=True)
    server = Server()
    rospy.Subscriber("/gazebo/model_states",ModelStates,server.modelState_callback)
    # rospy.Subscriber("/gazebo/joint_states",JointState,server.jointState_callback)
    rospy.Subscriber("/teeterbot/imu",Imu,server.IMU_callback)
    

    rospy.spin()