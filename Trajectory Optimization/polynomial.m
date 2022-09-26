clc;clear;close all ;
M= csvread("OptTrajectory.csv");

t=M(:,1);
l=M(:,2);
psi=M(:,3);
theta=M(:,4);v=M(:,5);psi_dot=M(:,6);theta_dot=M(:,7);Tl=M(:,8);Tr=M(:,9);
% disp(M);

% l_poly = polyfit(t,l,20);
% psi_poly = polyfit(t,psi,20);
% theta_poly = polyfit(t,theta,20);
% 
% v_poly = polyfit(t,v,20);
% psi_dot_poly = polyfit(t,psi_dot,20);
% theta_dot_poly = polyfit(t,theta_dot,20);
Tl_poly = polyfit(t,Tl,30);
Tr_poly = polyfit(t,Tr,30);

plot(t,polyval(Tr_poly,t),'bo',t,Tr,'r-')


