function [sys]=wtsLPV(Ts)
%WTSLPV LPV model of wind turbine system

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

A1 = zeros(7);

A2 = zeros(7);

A3 = zeros(7);

B0 = [hMx/Jr hMx/Jr hMx/Jr kMx/Jr kMx/Jr kMx/Jr 0;
    0 0 0 0 0 0 0;
    hFx/mtw hFx/mtw hFx/mtw kFx/mtw kFx/mtw kFx/mtw 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 1/mtw*3/(2*H);
    0 0 0 0 0 0 0;
    hMx/Jr hMx/Jr hMx/Jr kMx/Jr kMx/Jr kMx/Jr 1/Jg];

B1 = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    3*hMz/(2*H*mtw) 0 0 3*kMz/(2*H*mtw) 0 0 0;
    0 0 0 0 0 0 0;
    -hFz/mtw 0 0 -kFz/mtw 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];

B2 = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 3*hMz/(2*H*mtw) 0 0 3*kMz/(2*H*mtw) 0 0;
    0 0 0 0 0 0 0;
    0 -hFz/mtw 0 0 -kFz/mtw 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];

B3 = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 3*hMz/(2*H*mtw) 0 0 3*kMz/(2*H*mtw) 0;
    0 0 0 0 0 0 0;
    0 0 -hFz/mtw 0 0 -kFz/mtw 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];

C0 = [1 0 0 0 0 0 -1;
    0 0 1 0 0 0 0;
    0 0 0 0 1 0 0;
    0 0 -hMz 0 0 0 0;
    0 0 -hMz 0 0 0 0;
    0 0 -hMz 0 0 0 0];

C1 = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 9*hMz*Rb/(8*H) 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];

C2 = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 9*hMz*Rb/(8*H) 0 0 0 0;
    0 0 0 0 0 0 0];

C3 = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 9*hMz*Rb/(8*H) 0 0 0 0];

D0 = [0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    hMz 0 0 kMz 0 0 0;
    0 hMz 0 0 kMz 0 0;
    0 0 hMz 0 0 kMz 0];

D1 = zeros(6,7);

D2 = zeros(6,7);

D3 = zeros(6,7);

states = {'\Omega_r' 'x_{fa}' '\dot{x}_{fa}' 'x_{sd}' '\dot{x}_{sd}' '\gamma' '\dot{\gamma}'};
inputs = {'v_{cm1}' 'v_{cm2}' 'v_{cm3}' '\theta_{cm1}' '\theta_{cm2}' '\theta_{cm3}' 'T_g'};
outputs = {'\Omega_g' '\dot{x}_{fa}' '\dot{x}_{sd}' 'M_{cm1}' 'M_{cm2}' 'M_{cm3}'};

sys = idafflpv([A0 zeros(size(A0,1),3*size(A0,1))],[B0 B1 B2 B3],[C0 zeros(size(C0,1),3*size(A0,1))],[D0 D1 D2 D2],'statename',states,'inputname',inputs,'outputname',outputs);
sys = idafflpvA2idss(sys);
sys = c2d(sys,Ts);
sys = idss2idafflpvA(sys,3);

