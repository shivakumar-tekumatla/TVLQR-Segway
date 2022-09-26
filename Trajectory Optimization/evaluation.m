clc;clear;close all ;
%% Reaing the optimal trajectories 
M= csvread("OptTrajectory.csv");
figure(1)
time=M(:,1);
l=M(:,2);
psi=M(:,3);
theta=M(:,4);v=M(:,5);psi_dot=M(:,6);theta_dot=M(:,7);Tl_act=M(:,8);Tr_act=M(:,9);
%% Plotting actuals 
plot(time,M(:,2:7),'linewidth',2);
title('Actual States')
lgd = legend('l','Psi','Theta','V','Psi-dot','Theta-dot');
lgd.FontSize = 14;

figure(2)

plot(time,M(:,8:9),'linewidth',2);
title('Actual Torques')
lgd = legend('Left Torque',"Right Torque");
lgd.FontSize = 14;

%% Fitting polynomial for torques 
% Tl_poly = polyfit(t,Tl,30);
% Tr_poly = polyfit(t,Tr,30);

%% Simulation 
[t,y] = ode45(@(t,y) Segway(t,y,Tl_act,Tr_act,time), [0,max(time)],[0.5,pi/3,pi/6,0,0,0]);
figure(3)
plot(t,[interp1(time,Tl_act,t),interp1(time,Tr_act,t)],'linewidth',2);
title('Simulated Torques')
lgd = legend('Left Torque',"Right Torque");
lgd.FontSize = 14;

figure(4);
plot(t,y,'linewidth',2);
title('Simulated States')
lgd = legend('l','Psi','Theta','V','Psi-dot','Theta-dot');
lgd.FontSize = 14;

%% ODE Solver  
 function dz  = Segway(t,z,Tl_act,Tr_act,time)

dz = zeros(6,1);
z = num2cell(z);

[l , psi , theta ,v, psi_dot,theta_dot] = deal(z{:}) ; 

if abs(theta) > 2*pi
    theta = mod(theta, 2*pi);
end 

if abs(psi) > 2*pi
    psi = mod(psi, 2*pi);
    
end 


Tl = interp1(time,Tl_act,t);

Tr = interp1(time,Tr_act,t);
    

disp([t,l,psi,theta,v,psi_dot,theta_dot,Tl,Tr]);
dz(1) = v; 
dz(2) = psi_dot; 
dz(3) = theta_dot;
dz(4) = double((365577050296978040197816623*((2000*Tl)/123 + (2000*Tr)/123 + cos(psi)*((293751*cos(psi)*sin(theta)*psi_dot^2)/500000 + (293751*cos(theta)*sin(psi)*psi_dot*theta_dot)/250000 + (519*v*sin(psi)*psi_dot)/200 + (293751*cos(psi)*sin(theta)*theta_dot^2)/500000) - sin(psi)*((519*psi_dot*v*cos(psi))/200 - (293751*psi_dot^2*sin(psi)*sin(theta))/500000 - (293751*theta_dot^2*sin(psi)*sin(theta))/500000 + (293751*psi_dot*theta_dot*cos(psi)*cos(theta))/250000)))/(281474976710656000000000000000*((189734489104131602862666827337*cos(psi)^2)/56294995342131200000000000000000 + (189734489104131602862666827337*sin(psi)^2)/56294995342131200000000000000000 - (1305476114865129*cos(psi)^4*cos(theta)^2)/1000000000000000000 - (1305476114865129*cos(theta)^2*sin(psi)^4)/1000000000000000000 - (1305476114865129*cos(psi)^2*cos(theta)^2*sin(psi)^2)/500000000000000000 + 24163993013218192887687/140737488355328000000000000)) + (15129*((293751*cos(psi)^2*cos(theta))/500000 + (293751*cos(theta)*sin(psi)^2)/500000)*(Tl + Tr - (288169731*sin(theta))/50000000 - (3010681571494894472289*psi_dot^2*cos(theta)*sin(theta))/8796093022208000000000))/(4000000*((189734489104131602862666827337*cos(psi)^2)/56294995342131200000000000000000 + (189734489104131602862666827337*sin(psi)^2)/56294995342131200000000000000000 - (1305476114865129*cos(psi)^4*cos(theta)^2)/1000000000000000000 - (1305476114865129*cos(theta)^2*sin(psi)^4)/1000000000000000000 - (1305476114865129*cos(psi)^2*cos(theta)^2*sin(psi)^2)/500000000000000000 + 24163993013218192887687/140737488355328000000000000)));
dz(5) = double(-(15129*((110*Tl)/41 - (110*Tr)/41 + v*((293751*psi_dot*cos(psi)^2*sin(theta))/500000 + (293751*psi_dot*sin(psi)^2*sin(theta))/500000) + (3010681571494894472289*psi_dot*theta_dot*cos(theta)*sin(theta))/4398046511104000000000))/(4000000*((3041890139668318641*cos(theta)^2)/144115188075855872000000 + (77939795071369559382053049*sin(theta)^2)/56294995342131200000000000000 + 3412042821/80000000000000)));
dz(6) = double(- (((7851951*cos(psi)^2)/800000000 + (7851951*sin(psi)^2)/800000000 + 1/2000)*(Tl + Tr - (288169731*sin(theta))/50000000 - (3010681571494894472289*psi_dot^2*cos(theta)*sin(theta))/8796093022208000000000))/((189734489104131602862666827337*cos(psi)^2)/56294995342131200000000000000000 + (189734489104131602862666827337*sin(psi)^2)/56294995342131200000000000000000 - (1305476114865129*cos(psi)^4*cos(theta)^2)/1000000000000000000 - (1305476114865129*cos(theta)^2*sin(psi)^4)/1000000000000000000 - (1305476114865129*cos(psi)^2*cos(theta)^2*sin(psi)^2)/500000000000000000 + 24163993013218192887687/140737488355328000000000000) - (15129*((293751*cos(psi)^2*cos(theta))/500000 + (293751*cos(theta)*sin(psi)^2)/500000)*((2000*Tl)/123 + (2000*Tr)/123 + cos(psi)*((293751*cos(psi)*sin(theta)*psi_dot^2)/500000 + (293751*cos(theta)*sin(psi)*psi_dot*theta_dot)/250000 + (519*v*sin(psi)*psi_dot)/200 + (293751*cos(psi)*sin(theta)*theta_dot^2)/500000) - sin(psi)*((519*psi_dot*v*cos(psi))/200 - (293751*psi_dot^2*sin(psi)*sin(theta))/500000 - (293751*theta_dot^2*sin(psi)*sin(theta))/500000 + (293751*psi_dot*theta_dot*cos(psi)*cos(theta))/250000)))/(4000000*((189734489104131602862666827337*cos(psi)^2)/56294995342131200000000000000000 + (189734489104131602862666827337*sin(psi)^2)/56294995342131200000000000000000 - (1305476114865129*cos(psi)^4*cos(theta)^2)/1000000000000000000 - (1305476114865129*cos(theta)^2*sin(psi)^4)/1000000000000000000 - (1305476114865129*cos(psi)^2*cos(theta)^2*sin(psi)^2)/500000000000000000 + 24163993013218192887687/140737488355328000000000000)));
end