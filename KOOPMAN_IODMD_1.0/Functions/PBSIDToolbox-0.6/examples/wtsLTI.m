function [sysl,sys_closed] = wtsLTI(h)
%WTSLTI LTI model of coleman transformed wind turbine system

% Parameters of linearised wind turbine model
Rb = 40;            % rotor radius [m]
B = 3;              % number of rotor blades [-]
H = 56;             % tower height [m]
kMx = -3.771e4;     % gain from pitch angle to root lead moment [Nm/deg]
kMz = 1.617e5;      % gain from pitch angle to root flap moment [Nm/deg]
kFx = -6.148e3;     % gain from pitch angle to root lead force [N/deg]
kFz = -1.831e3;     % gain from pitch angle to root flap force [N/deg]
hMx = 8.381e4;      % gain from wind speed to root lead moment [Nm/(m/s)]
hMz = -1.895e5;     % gain from wind speed to root flap moment [Nm/(m/s)]
hFx = 7.202e3;      % gain from wind speed to root lead force [N/(m/s)]
hFz = 4.068e3;      % gain from wind speed to root flap force [N/(m/s)]
Jr = 7.187e6;       % moment of inertia turbine rotor [kgm2]
Jg = 1.067e6;       % moment of inertia generator rotor [kgm2]
ksh = 1.262e8;      % drive-train stiffness [Nm/rad]
dsh = 1.262e5;      % drive-train damping [Nm/rad]
mtw = 156.6e3;      % tower top equivalent mass 1st bend freq. [kg]
ktw = 1.235e6;      % tower top equivalent stiffness 1st bend freq. [kg]
dtw = 2.800e3;      % tower top equivalent damper constant 1st freq. [N/(m/s)]

% System matrices
A0 = [0 0 -3*hMx/Jr 0 0 -ksh/Jr -dsh/Jr;
    0 0 1 0 0 0 0;
    0 -ktw/mtw -dtw/mtw-3*hFx/mtw+hMz/mtw*9/(4*H)*9*Rb/(8*H) 0 0 0 0;
    0 0 0 0 1 0 0;
    0 0 -3/2*9*Rb/(8*H)*hFz/mtw -ktw/mtw -dtw/mtw 0 0;
    0 0 0 0 0 0 1;
    0 0 -3*hMx/Jr 0 0 -(Jr+Jg)/(Jr*Jg)*ksh -(Jr+Jg)/(Jr*Jg)*dsh];

B0 = [3*hMx/Jr 0 0 3*kMx/Jr 0 0 0;
    0 0 0 0 0 0 0;
    -3*hFx/mtw -9*hMz/(4*H*mtw) 0 3*kFx/mtw 9*kMz/(4*H*mtw) 0 0;
    0 0 0 0 0 0 0;
    0 -3*hFz/(2*mtw) 0  0 -3*kFz/(2*mtw) 0 3/(2*H*mtw);
    0 0 0 0 0 0 0;
    3*hMx/Jr 0 0 3*kMx/Jr 0 0 1/Jg];

C0 = [1 0 0 0 0 0 -1;
    0 0 1 0 0 0 0;
    0 0 0 0 1 0 0;
    0 0 -hMz 0 0 0 0;
    0 0 hMz*9/8*Rb/H 0 0 0 0;
    0 0 0 0 0 0 0];

D0 = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    hMz 0 0 kMz 0 0 0;
    0 hMz 0 0 kMz 0 0;
    0 0 hMz 0 0 kMz 0];


%% System
states = {'\Omega_r' 'x_{fa}' '\dot{x}_{fa}' 'x_{sd}' '\dot{x}_{sd}' '\gamma' '\dot{\gamma}'};
inputs = {'v_{cm1}' 'v_{cm2}' 'v_{cm3}' '\theta_{cm1}' '\theta_{cm2}' '\theta_{cm3}' 'T_g'};
outputs = {'\Omega_g' '\dot{x}_{fa}' '\dot{x}_{sd}' 'M_{cm1}' 'M_{cm2}' 'M_{cm3}'};
sysl = c2d(ss(A0,B0,C0,D0,'statename',states,'inputname',inputs,'outputname',outputs),h,'tustin');

% with pitch angles also as output
B0 = [B0 B0(:,[4 7])];
C0 = [C0; zeros(2,size(A0,2))];
D0 = [[D0 D0(:,[4 7])]; 0 0 0 1 0 0 0 1 0; 0 0 0 0 0 0 1 0 1];
inputs = {'v_{cm1}' 'v_{cm2}' 'v_{cm3}' '\theta_{cm1}' '\theta_{cm2}' '\theta_{cm3}' 'T_g' 'C_\theta_{cm1}' 'C_T_g'};
outputs = {'\Omega_g' '\dot{x}_{fa}' '\dot{x}_{sd}' 'M_{cm1}' 'M_{cm2}' 'M_{cm3}' 'F_\theta_{cm1}' 'F_T_g'};
sysc = c2d(ss(A0,B0,C0,D0,'statename',states,'inputname',inputs,'outputname',outputs),h,'tustin');

%% Collective pitch controller (speed regulation)
Delay = tf(1,1,h,'inputdelay',floor(0.14/h));
[b,a] = ellip(4,0.5,14,0.719*(2*h));
Filter = tf(b,a,h);
[b,a] = ellip(2,1,20,[0.3799*(2*h) 0.5140*(2*h)]);
Band = tf(b,a,h);
K = 46.7;
Ti = 10;
PI = c2d(zpk(-1/Ti,0,-K),h,'tustin');
Cpitch = delay2z(Delay*Filter*Band*PI);

%% Torque controller
[b,a] = ellip(4,0.1,30,0.15*(2*h));
Filter = tf(b,a,h);
K = -9.226e5;
Ctorque = delay2z(Filter*K);

%% Closed loop
sys_closed = feedback(sysc,Ctorque,9,1);
sys_closed.inputname{9} = 'Ref_\T_g';
sys_closed = feedback(sys_closed,Cpitch,8,1);
sys_closed.inputname{8} = 'Ref_\Omega_g';




