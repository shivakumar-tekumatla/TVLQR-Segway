% SFarzan
% Code snippet for integrating differential Riccati equation backward in time
% Fall 2021

% n: state dimension
% sys: function handle for the linearized system
% t: time array
% Q,R,Qf: LQR gains

odeFun = @(t,s)dre(t,s,sys,Q,R,n);
s0 = reshape(Qf,n*n,1);
tSpan = [t(end),t(1)];

options = odeset();
options.RelTol = 1e-6;
options.AbsTol = 1e-6;
S = ode45(odeFun,tSpan,s0);
% s = reshape(deval(S,t),n,n);

function ds = dre(t,s,sys,Q,R,n)

S = reshape(s,n,n);
[A,B] = sys(t);

dS = -(A'*S + S*A - S*B*(R\B'*S) + Q);
ds = reshape(dS,n*n,1);

end
