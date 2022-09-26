clc;clear;close all;
%% Symbols 
syms m Jw r mb h Jx Jz Jw a r mw Jy l l_dot l_ddot theta theta_dot  psi psi_dot  g T Tl Tr v

%% Equations of Motions 
% Before applyting the non-holonomic constraints 
Mq1 = [      m                               0;
             0                               m;
       -mb*h*sin(psi)*sin(theta)       mb*h*cos(psi)*sin(theta);
       mb*h*cos(psi)*cos(theta)        mb*h*sin(psi)*cos(theta);
             0                               0;
             0                               0];
         
 Mq2 = [    -mb*h*sin(psi)*sin(theta)            mb*h*cos(psi)*cos(theta)      0        0;
             mb*h*cos(psi)*sin(theta)            mb*h*sin(psi)*cos(theta)      0        0;
       (Jx+m*h^2)*sin(theta)^2+Jz*cos(theta)^2+2*mw*a^2     0                  0        0;
                       0                         Jy+mb*h^2+2*Jw               Jw        Jw;
                       0                               Jw                     Jw        0;
                       0                               Jw                      0        Jw];
 Mq = [Mq1 Mq2];
                   
% Non Holonomic Constraints 
H =[-sin(psi)  cos(psi)      0           0           0           0;
       1          0    -a*cos(psi)   -r*cos(psi)  -r*cos(psi)    0;
       0          1    -a*sin(psi)   -r*sin(psi)  -r*sin(psi)    0;
       1          0     a*cos(psi)   -r*cos(psi)     0        -r*cos(psi);
       0          1     a*sin(psi)   -r*sin(psi)     0        -r*sin(psi)];
   
 
Bq = [mb*h*cos(psi)*sin(theta)*psi_dot^2+2*mb*h*sin(psi)*cos(theta)*psi_dot*theta_dot+mb*h*sin(theta)*cos(psi)*theta_dot^2;
      mb*h*sin(theta)*sin(psi)*theta_dot^2-2*mb*h*cos(theta)*cos(psi)*theta_dot*psi_dot+mb*h*sin(psi)*sin(theta)*psi_dot^2;
      -2*(Jx-Jz+mb*h^2)*cos(theta)*sin(theta)*psi_dot*theta_dot;
      (Jx-Jz+mb*h^2)*psi_dot^2*cos(theta)*sin(theta)+mb*g*h*sin(theta);
      0;
      0];
  
 Fq = [0 0 0 0 1 0;
       0 0 0 0 0 1]';
 u  =[Tl;Tr];
 
 Phi = [cos(psi)   0   0;
        sin(psi)   0   0;
          0        1   0;
          0        0   1;
          1/r    -a/r -1;
          1/r    a/r   -1];
 p_dot = [v;
          psi_dot;
          theta_dot];
 Phi_dot = [-sin(psi)*psi_dot    0 0;
             cos(psi)*psi_dot    0 0;
                  0              0 0;
                  0              0 0;
                  0              0 0;
                  0              0 0];
 Phi_trans = transpose(Phi);
 % After applyting the non-holonomic constraints 
 Mp = Phi_trans*Mq*Phi;
 Bp =  Phi_trans*(Bq-Mq*Phi_dot*p_dot);
 Fp =  Phi_trans*Fq;
 
RHS = Bp + Fp*u;

p_ddot = inv(Mp)*RHS;

%% Physical and Mechanical Parameters
 a = 0.165  ; % 2*a is the width of the segway 
 h = 0.254  ; % Center of mass 
 g = 9.81   ; % Acceleration due to gravity 
 mb = 2.313 ; % Mass of the body
 mw = 0.141 ; % Mass of the wheel 
 r  = 0.0615; %radius of the wheel 
 lb = 0.0508; % Length along X direction
 wb = 0.1524; % Along Y direction
 hb = 0.9462; % Along Z direction 
 m  = mb+2*mw;
 Jbx = mb*h^2/2; %MoI of body wrt X-Axis 
 
 Jby = mb*h^2/2; %MoI of body wrt Y-Axis 
 
 Jbz = mb*h^2;   %MoI of body wrt Z-Axis 
 
 Jw  = 0.00025;%mw*r^2;  %MoI of wheel in axial direction
 
 Jwt = Jw/3;    %MoI of wheel in perpendicular direction
 
 Jx = m*(hb^2+wb^2)/12; %Jbx + 2*Jwt;
 Jy = m*(hb^2+lb^2)/12; %Jby;
 Jz = m*(lb^2+wb^2)/12; %Jbz + 2*Jwt;

%% State Space Model 
v_dot      = subs(p_ddot(1));
psi_ddot   = subs(p_ddot(2));
theta_ddot = subs(p_ddot(3));

p_dot = sym('p_dot' , [6,1]);

 p_dot(1) = v;
 p_dot(2) = psi_dot;
 p_dot(3) = theta_dot;
 p_dot(4) = v_dot;
 p_dot(5) = psi_ddot;
 p_dot(6) = theta_ddot;
 
 %% Jacobian Linearization for system 1 
 A = jacobian(p_dot,[l,psi,theta,v,psi_dot,theta_dot]);
 B = jacobian(p_dot,[Tl,Tr]);
 
 A_eq = double(subs(A,[l,psi,theta,v,psi_dot,theta_dot,Tl,Tr],[0,0,0,0,0,0,0,0]));
 B_eq = double(subs(B,[l,psi,theta,v,psi_dot,theta_dot,Tl,Tr],[0,0,0,0,0,0,0,0]));
 
 
%  A1_eq = double(subs(A1,[l,theta,l_dot,theta_dot,T],[0,0,0,0,0]));
%  B1_eq = double(subs(B1,[l,theta,l_dot,theta_dot,T],[0,0,0,0,0]));

%% Optimized trajectory 
traj= csvread("OptTrajectory.csv");
% These are the actual/desired trajectories obtained using trajectory
% optimization 
time=traj(:,1);
l_act=traj(:,2);
psi_act=traj(:,3);
theta_act=traj(:,4);
v_act=traj(:,5);
psi_dot_act=traj(:,6);
theta_dot_act=traj(:,7);
Tl_act=traj(:,8);
Tr_act=traj(:,9);

%% LQR Design 

 Q = eye(size(p_dot,1));
 R = eye(size(u,1));
 Qf =  Q.*2;
 
%% Calling ODE function 
[t,y] = ode45(@(t,y) Segway(t,y,A,B,Q,R,Qf,time,l_act,psi_act,theta_act,v_act,psi_dot_act,theta_dot_act,Tl_act,Tr_act), [0,max(time)],[0.6,pi/4,pi/6,0,0,0]);

figure(1);
plot(t,y,'linewidth',2);
title('Time vs l, Psi, Theta, V,Psi-dot,Theta-dot')
lgd = legend('l','Psi','Theta','V','Psi-dot','Theta-dot');
lgd.FontSize = 14;

figure(2);
subplot(3,2,1)      
plot(t,y(:,1),"b-",time,l_act,"r--",'linewidth',2);   
lgd = legend('l','l-desired');
lgd.FontSize = 14;
title('Time vs l')
xlabel("Time in Seconds")
ylabel("Distance in Meters")

subplot(3,2,2)       
plot(t,y(:,2),"b-",time,psi_act,"r--",'linewidth',2);   
lgd = legend('Psi','Psi-desired');
lgd.FontSize = 14;     
title('Time vs Psi ')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(3,2,3)       
plot(t,y(:,3),"b-",time,theta_act,"r--",'linewidth',2);   
lgd = legend('Theta','Theta-desired');
lgd.FontSize = 14;     
title('Time vs Theta')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(3,2,4)       
plot(t,y(:,4),"b-",time,v_act,"r--",'linewidth',2);   
lgd = legend('V','V-desired');
lgd.FontSize = 14;       
title('Time vs l-Dot')
xlabel("Time in Seconds")
ylabel("Velocity in m/s")

subplot(3,2,5)       
plot(t,y(:,5),"b-",time,psi_dot_act,"r--",'linewidth',2);   
lgd = legend('Psi-dot','Psi-dot-desired');
lgd.FontSize = 14;      
title('Time vs Psi-dot ')
xlabel("Time in Seconds")
ylabel("Anglular Velocity in Radians per second")

subplot(3,2,6)       
plot(t,y(:,6),"b-",time,theta_dot_act,"r--",'linewidth',2);   
lgd = legend('Theta-dot','Theta-dot-desired');
lgd.FontSize = 14;    
title('Time vs Theta-dot')
xlabel("Time in Seconds")
ylabel("Angular Velocity in Radians per second")

%% Linearize 
function [Al,Bl] = linearize(A,B,t,time,l_act,psi_act,theta_act,v_act,psi_dot_act,theta_dot_act,Tl_act,Tr_act)
%This function linearizes the system around the optimized trajectory
    states = interp1(time,[l_act,psi_act,theta_act,v_act,psi_dot_act,theta_dot_act,Tl_act,Tr_act],t); %Obtainig  states and inputs at given point in time 
    l=states(1);psi=states(2);theta=states(3);v=states(4);psi_dot=states(5);theta_dot=states(6);Tl=states(7);Tr=states(8);
    % [l,psi,theta,v,psi_dot,theta_dot] 
     Al = double(subs(A));%,[l,psi,theta,v,psi_dot,theta_dot,Tl,Tr]));%,[l,psi,theta,v,psi_dot,theta_dot,Tl,Tr]));
     Bl = double(subs(B));%,[l,psi,theta,v,psi_dot,theta_dot,Tl,Tr]));%,[l,psi,theta,v,psi_dot,theta_dot,Tl,Tr]));
end 
%% DREQ
function ds = dre(t,s,A,B,Q,R,n)

S = reshape(s,n,n);
% [A,B] = sys(t);

dS = -(A'*S + S*A - S*B*(R\B'*S) + Q);
ds = reshape(dS,n*n,1);

end
%% ODE Solver  
 function dz  = Segway(t,z,A,B,Q,R,Qf,time,l_act,psi_act,theta_act,v_act,psi_dot_act,theta_dot_act,Tl_act,Tr_act)
    disp(t)
    [Al,Bl] =linearize(A,B,t(end),time,l_act,psi_act,theta_act,v_act,psi_dot_act,theta_dot_act,Tl_act,Tr_act);
%     disp(Bl)
    n = 6; % State dimension 
    
    % Solving for P - integrating backwards in time 
    odeFun = @(t,s)dre(t,s,Al,Bl,Q,R,n); 
    s0 = reshape(Qf,n*n,1);
    tSpan = [t(end)+0.001,t(1)];

    options = odeset();
    options.RelTol = 1e-6;
    options.AbsTol = 1e-6;
    P = ode45(odeFun,tSpan,s0);
%     s = deval(S,t);
    p = reshape(deval(P,t),n,n);
%     disp(s)
    K = inv(R)*Bl'*p;
    states = interp1(time,[l_act,psi_act,theta_act,v_act,psi_dot_act,theta_dot_act,Tl_act,Tr_act],t(end));
    l_ref=states(1);psi_ref=states(2);theta_ref=states(3);v_ref=states(4);psi_dot_ref=states(5);theta_dot_ref=states(6);Tl_ref=states(7);Tr_ref=states(8);

    
    dz = zeros(6,1);
    z = num2cell(z);

    [l , psi , theta ,v, psi_dot,theta_dot] = deal(z{:}) ; 
    
    z_new = [l , psi , theta ,v, psi_dot,theta_dot]-[l_ref,psi_ref,theta_ref,v_ref,psi_dot_ref,theta_dot_ref];

    if abs(theta) > 2*pi
        theta = mod(theta, 2*pi);
    end 

    if abs(psi) > 2*pi
        psi = mod(psi, 2*pi);

    end 

    U_ref = [Tl_ref;Tr_ref];

    U = U_ref - K*z_new';
    Tl = U(1);
    Tr = U(2);

    dz(1) = v; 
    dz(2) = psi_dot; 
    dz(3) = theta_dot;
    dz(4) = double((365577050296978040197816623*((2000*Tl)/123 + (2000*Tr)/123 + cos(psi)*((293751*cos(psi)*sin(theta)*psi_dot^2)/500000 + (293751*cos(theta)*sin(psi)*psi_dot*theta_dot)/250000 + (519*v*sin(psi)*psi_dot)/200 + (293751*cos(psi)*sin(theta)*theta_dot^2)/500000) - sin(psi)*((519*psi_dot*v*cos(psi))/200 - (293751*psi_dot^2*sin(psi)*sin(theta))/500000 - (293751*theta_dot^2*sin(psi)*sin(theta))/500000 + (293751*psi_dot*theta_dot*cos(psi)*cos(theta))/250000)))/(281474976710656000000000000000*((189734489104131602862666827337*cos(psi)^2)/56294995342131200000000000000000 + (189734489104131602862666827337*sin(psi)^2)/56294995342131200000000000000000 - (1305476114865129*cos(psi)^4*cos(theta)^2)/1000000000000000000 - (1305476114865129*cos(theta)^2*sin(psi)^4)/1000000000000000000 - (1305476114865129*cos(psi)^2*cos(theta)^2*sin(psi)^2)/500000000000000000 + 24163993013218192887687/140737488355328000000000000)) + (15129*((293751*cos(psi)^2*cos(theta))/500000 + (293751*cos(theta)*sin(psi)^2)/500000)*(Tl + Tr - (288169731*sin(theta))/50000000 - (3010681571494894472289*psi_dot^2*cos(theta)*sin(theta))/8796093022208000000000))/(4000000*((189734489104131602862666827337*cos(psi)^2)/56294995342131200000000000000000 + (189734489104131602862666827337*sin(psi)^2)/56294995342131200000000000000000 - (1305476114865129*cos(psi)^4*cos(theta)^2)/1000000000000000000 - (1305476114865129*cos(theta)^2*sin(psi)^4)/1000000000000000000 - (1305476114865129*cos(psi)^2*cos(theta)^2*sin(psi)^2)/500000000000000000 + 24163993013218192887687/140737488355328000000000000)));
    dz(5) = double(-(15129*((110*Tl)/41 - (110*Tr)/41 + v*((293751*psi_dot*cos(psi)^2*sin(theta))/500000 + (293751*psi_dot*sin(psi)^2*sin(theta))/500000) + (3010681571494894472289*psi_dot*theta_dot*cos(theta)*sin(theta))/4398046511104000000000))/(4000000*((3041890139668318641*cos(theta)^2)/144115188075855872000000 + (77939795071369559382053049*sin(theta)^2)/56294995342131200000000000000 + 3412042821/80000000000000)));
    dz(6) = double(- (((7851951*cos(psi)^2)/800000000 + (7851951*sin(psi)^2)/800000000 + 1/2000)*(Tl + Tr - (288169731*sin(theta))/50000000 - (3010681571494894472289*psi_dot^2*cos(theta)*sin(theta))/8796093022208000000000))/((189734489104131602862666827337*cos(psi)^2)/56294995342131200000000000000000 + (189734489104131602862666827337*sin(psi)^2)/56294995342131200000000000000000 - (1305476114865129*cos(psi)^4*cos(theta)^2)/1000000000000000000 - (1305476114865129*cos(theta)^2*sin(psi)^4)/1000000000000000000 - (1305476114865129*cos(psi)^2*cos(theta)^2*sin(psi)^2)/500000000000000000 + 24163993013218192887687/140737488355328000000000000) - (15129*((293751*cos(psi)^2*cos(theta))/500000 + (293751*cos(theta)*sin(psi)^2)/500000)*((2000*Tl)/123 + (2000*Tr)/123 + cos(psi)*((293751*cos(psi)*sin(theta)*psi_dot^2)/500000 + (293751*cos(theta)*sin(psi)*psi_dot*theta_dot)/250000 + (519*v*sin(psi)*psi_dot)/200 + (293751*cos(psi)*sin(theta)*theta_dot^2)/500000) - sin(psi)*((519*psi_dot*v*cos(psi))/200 - (293751*psi_dot^2*sin(psi)*sin(theta))/500000 - (293751*theta_dot^2*sin(psi)*sin(theta))/500000 + (293751*psi_dot*theta_dot*cos(psi)*cos(theta))/250000)))/(4000000*((189734489104131602862666827337*cos(psi)^2)/56294995342131200000000000000000 + (189734489104131602862666827337*sin(psi)^2)/56294995342131200000000000000000 - (1305476114865129*cos(psi)^4*cos(theta)^2)/1000000000000000000 - (1305476114865129*cos(theta)^2*sin(psi)^4)/1000000000000000000 - (1305476114865129*cos(psi)^2*cos(theta)^2*sin(psi)^2)/500000000000000000 + 24163993013218192887687/140737488355328000000000000)));
 end 