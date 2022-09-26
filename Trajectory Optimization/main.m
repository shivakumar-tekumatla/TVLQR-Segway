clear all;close all;format compact;

[problem,guess]=Segway;          % Fetch the problem definition
options= problem.settings(20);                  % Get options and solver settings 
[solution,MRHistory]=solveMyProblem( problem,guess,options);
[ tv, xv, uv ] = simulateSolution( problem, solution, 'ode113', 0.001 );
csvwrite("OptTrajectory.csv",[ tv, xv, uv ])
%% PLotting 
figure;
subplot(2,1,1);
plot(tv,xv,'linewidth',2);
title('Time vs l, Psi, Theta, V,Psi-dot,Theta-dot')
lgd = legend('l','Psi','Theta','V','Psi-dot','Theta-dot');
lgd.FontSize = 14;

subplot(2,1,2);

plot(tv,uv,'linewidth',2);
title('Generalized Forces - Time vs Left Torque and Right Torque')
lgd = legend('Left Torque','Right Torque');
lgd.FontSize = 14;

t= tv;
y=xv;
figure(2);
subplot(3,2,1)      
plot(t,y(:,1),'linewidth',2);          
title('Time vs l')
xlabel("Time in Seconds")
ylabel("Distance in Meters")

subplot(3,2,2)       
plot(t,y(:,2),'linewidth',2);       
title('Time vs Psi ')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(3,2,3)       
plot(t,y(:,3),'linewidth',2);      
title('Time vs Theta')
xlabel("Time in Seconds")
ylabel("Angle in Radians")

subplot(3,2,4)       
plot(t,y(:,4),'linewidth',2);        
title('Time vs l-Dot')
xlabel("Time in Seconds")
ylabel("Velocity in m/s")

subplot(3,2,5)       
plot(t,y(:,5),'linewidth',2);       
title('Time vs Psi-dot ')
xlabel("Time in Seconds")
ylabel("Anglular Velocity in Radians per second")

subplot(3,2,6)       
plot(t,y(:,6),'linewidth',2);      
title('Time vs Theta-dot')
xlabel("Time in Seconds")
ylabel("Angular Velocity in Radians per second")