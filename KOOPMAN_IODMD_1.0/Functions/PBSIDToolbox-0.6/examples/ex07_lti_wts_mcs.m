%% Example 7: Uncertainty bounds on estimation using Monte-Carlo simulations
close all; clear; clc;
addpath D:\Werk_Data\Software\PBSIDToolbox-0.6
%% LTI model of a Coleman tranformed wind turbine system

% LTI system matrices
h = 0.1;             % Sample time
[OL,CL] = wtsLTI(h); % The wind turbine model
n = size(OL.a,1);    % The order of the system
r = size(OL.b,2);    % The number of inputs
l = size(OL.c,1);    % The number of outputs

%% Closed-loop identification experiment
% Simulation of the model in closed loop

% Number of monte carlo simulations
MCS = 100;

% Time sequence
N = 10000;  % number of data points
t = (0:h:h*(N-1))';

% Frequency grid
w = logspace(-2,log10(pi/h),1000);

% Allocate storage vectors
VAF = zeros(MCS,3);
EE = zeros(MCS,7);

for nmcs = 1:MCS
    disp(['Simulation ', num2str(nmcs)])
    
    % Wind disturbance signals
    d = randn(N,3);
    
    % Excitation signal for pitch input
    r_pitch = randn(N,1);
    
    % Excitation signal for Torque input
    r_torque = 1e3.*randn(N,1);

    % Add together for simulation
    r = [r_pitch zeros(N,2) r_torque zeros(N,2)];
    
    % Simulation of the closed-loop system
    y = lsim(CL,[d r],t);

    % Input and output selaction with scaling
    ui = detrend(y(:,7:8),'constant');   % selects input for identification (excitation of pitch + control)
    yi = detrend(y(:,1:3),'constant');   % selects output for identification
    ri = [r_pitch r_torque];
    [us,Du,ys,Dy] = sigscale(ui,yi); % signal scaling

    % Defining a number of constants
    p = 50;     % past window size
    f = 20;     % future window size

    % PBSID-opt
    [S,x] = dordvarx(us,ys,f,p,'tikh','gcv');
    x = dmodx(x,n);
    [Ai,Bi,Ci,Di,Ki] = dx2abcdk(x,us,ys,f,p);
    EE(nmcs,:) = eig(Ai);
    Dat = iddata(ys',us',h);
    Mi = abcdk2idss(Dat,Ai,Bi,Ci,Di,Ki);

    % Variance-accounted-for (by Kalman filter)
    yest = predict(Mi,Dat);
    x0 = findstates(Mi,Dat);
    vfm = vaf(ys,yest.y);

    % store results
    sys = ss(Ai,Bi/Du,Dy*Ci,Dy*Di/Du,h);
    if nmcs == 1
        Hmin = abs(freqresp(sys,w));
        Hmax = abs(freqresp(sys,w));
        Hvaf = abs(freqresp(sys,w));
    else
        Hmin = min(abs(Hmin),abs(freqresp(sys,w)));
        Hmax = max(abs(Hmax),abs(freqresp(sys,w)));
        if mean(vfm) >  max(mean(VAF(:,:),2))
            Hvaf = abs(freqresp(sys,w));
        end
    end
    VAF(nmcs,:) = vfm';
end

%% Identification results
%

% Plot eigenvalues
realeig = eig(minreal(OL(1:5,[4 7])));
figure, subplot(2,2,[1 3]), deigen(EE',realeig);
axis([0 1.1 -1 1]);
subplot(2,2,2), deigen(EE',realeig);
axis([0.948 0.962 0.268 0.282]);
subplot(2,2,4), deigen(EE',realeig);
axis([0.96 1.01 -0.025 0.025]);

% frequency response of identified system
Hr = abs(freqresp(minreal(OL(1:5,[4 7])),w));
figure('Units','normalized','Position',[0 0 1 1]), hold on
dbodemagpatch(Hvaf,Hmin,Hmax,w,h,Hr);
hold off 
