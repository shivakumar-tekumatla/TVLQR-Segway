clc;close all;
% g) RR Bot control using Gazebo 
fprintf("------------Feedback Control Law Implementation in Gazebo----------------\n");

rosshutdown;

rosinit;
left_torque = rospublisher('/teeterbot/left_torque_cmd');%,'std_msgs/Float64');
right_torque = rospublisher('/teeterbot/right_torque_cmd');%,'std_msgs/Float64');
JointStates = rossubscriber('/teeterbot/joint_states');
imu_data    = rossubscriber("/teeterbot/imu");
modelStates = rossubscriber("/gazebo/model_states");

Tl = rosmessage(left_torque);
Tr = rosmessage(right_torque);
joints = rosmessage(JointStates);
model = rosmessage(modelStates);
client = rospublisher("/gazebo/set_model_state");

% client = rossvcclient("/gazebo/set_model_state");
req = rosmessage(client);
req.ModelName = 'teeterbot'; 
%% Initial Conditions 
init_theta = 0;
init_psi = 0;
init_l =0;
eul = [init_psi init_theta  0];
quat= eul2quat(eul,"ZYX");
req.Pose.Orientation.X = quat(2);
req.Pose.Orientation.Y = quat(3);
req.Pose.Orientation.Z = quat(4);
req.Pose.Orientation.W = quat(1);
req.Pose.Position.X = init_l;
req.Pose.Position.Y = 0;
req.Pose.Position.Z = 0.2;
send(client,req);

%% Simulation 
t=0;
tic;
theta=[init_theta];
psi = [init_psi];
l = [init_l];
l_dot=[0];
theta_dot=[0];
psi_dot=[0];
time =[0];
input1=[0];
input2=[0];

send(client,req);
disp("Pausing");
delta_t_array=[0];
% K =[-0.7071   -0.7071  -14.5424   -1.8846   -0.8687   -3.2223;
%    -0.7071    0.7071  -14.5424   -1.8846    0.8687   -3.2223];
% 
% % pause(1);
% K =[5   0  1    1  2   5;
%     0   2  1   1   2   5];

% eigen = [-1+2i -1-2i -1+1i -1-1i -2+1i -2-1i];
% K = place(A_eq,B_eq,eigen);
% r = rateControl(0.5);
while t<10
    t=toc;
    delta_t = t-time(end); %Computing time step
%     delta_t_array = [delta_t_array,delta_t];
    time = [time,t];
    
    IMU = receive(imu_data);
      modelData = receive(modelStates);

      Position = modelData.Pose(2).Position;
      X = Position.X;
      Y = Position.Y;
      temp_l = sqrt(X^2+Y^2);
      
      Orientation = modelData.Pose(2).Orientation;
      Angles = quat2eul([IMU.Orientation.W  IMU.Orientation.X IMU.Orientation.Y IMU.Orientation.Z],"ZYX");
%       disp([X,Y]);
      temp_psi = Angles(1);
      temp_theta = Angles(2);
      if temp_psi<0
          temp_psi=mod(temp_psi,pi);
      end
      if temp_theta<0
          temp_theta=mod(temp_theta,pi);
      end
      %Velocities
      temp_theta_dot =IMU.AngularVelocity.Y;%modelData.Twist(2).Angular.Y;%(temp_theta-theta(end))/delta_t;
      temp_psi_dot =IMU.AngularVelocity.Z;%modelData.Twist(2).Angular.Z;%(temp_psi-psi(end))/delta_t;
      temp_l_dot = sqrt((modelData.Twist(2).Linear.X)^2+(modelData.Twist(2).Linear.Y)^2);%temp_l-l(end))/delta_t;
      
      U = -K*[temp_l;temp_psi;temp_theta;temp_l_dot;temp_psi_dot;temp_theta_dot];
      disp(U);
      
%     disp(modelData);
%     disp(jointData.Velocity);
    Tl.Data = U(1);
    Tr.Data = U(2);
%     input1=[input1,U(1)];
%     input2=[input2,U(2)];
    send(left_torque,Tl);
    send(right_torque,Tr);
    %       disp(rad2deg(Angles));
%       l = [l,temp_l];
%       theta = [theta,temp_theta];
%       psi = [psi,temp_psi];
%       l_dot = [l_dot,temp_l_dot];
%       theta_dot=[theta_dot,temp_theta_dot];
%       psi_dot = [psi_dot,temp_psi_dot];
      disp("Positions");
      disp([temp_l,rad2deg(temp_psi),rad2deg(temp_theta)]);
      disp("Velocities");
      disp([temp_l_dot,temp_psi_dot,temp_theta_dot]);
%     waitfor(r);
end

Tl.Data = 0;
Tr.Data = 0;
send(left_torque,Tl);
send(left_torque,Tr);% disconnectfromroscore

rosshutdown; 
%% Plotting
y = [l' psi' theta' l_dot' psi_dot' theta_dot'];
figure(1);
plot(time,y,'linewidth',2);
title('Time vs l, Psi, Theta, V,Psi-dot,Theta-dot')
lgd = legend('l','Psi','Theta','V','Psi-dot','Theta-dot');
lgd.FontSize = 14;

figure(2);
plot(time,[input1' input2'],'linewidth',2);
title('Generalized Forces - Time vs Left Torque and Right Torque')
lgd = legend('Left Torque','Right Torque');
lgd.FontSize = 14;

%}
