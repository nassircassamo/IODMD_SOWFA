function [Inputs, Outputs, Deterministic,scalingfactors]=preprocessdmdid(beg, rotSpeed,time1,rotorAzimuth,nacelleYaw,pitchmode, pitch,powerGenerator,rho)

 %Def: this function prepares the input/output signals to be used for
 %system identification

    %Input arguments include all relevant singals that can be used for the
    %purposes of input or output
    
    %Output arguments include hte Input and Output signals duly processed
    %to be used for systems identification. Deterministic are the signals
    %which can be used to augment the state vector based on the ideas of
    %the Koopman Operator. Scaling factors can be used if the input/output
    %signals are scaled (for example, by the variance of the signal). These
    %then should be used again for the reconstruction part

    %log:
    %0. first commit October 2020
    %1. function revised and comments added on March 2022


%% Prepepare relevant signals

%the functions detrend and resample were firstly tried to remove means and trends from the
%signals, but it was seen that it had undesired effects near the edges of
%the signals or in the way the signal was being processed, so it was not used later on.
%the function resampleedge effect was created to avoid the undesired
%effects of the resample function and the detrend was achieved by removing
%the mean value from the origina signal, taken manually from the data.

 %X1=detrend(rotSpeed(end-beg*10:1:end,1)');
 %X2=detrend(rotSpeed(end-beg*10:1:end,2)');
 %X3=detrend(rotSpeed(end-beg*10:1:end,1).^2');
 %X4=detrend(rotSpeed(end-beg*10:1:end,2).^2');
 
 %meanX1=7.14; %pitch 0
 meanX1=7.14; %this represents the mean rotor speed of turbine 1 under yaw misalignment of -10 degrees
 X1=(rotSpeed(end-beg*10:1:end,1)-meanX1).^2;
 
 %meanX1=7.14; %pitch 0
 meanX2=4.027; %this represents the mean rotor speed of turbine 2 under yaw misalignment of -10 degrees
 X2=(rotSpeed(end-beg*10:1:end,2)-meanX2).^2;
 
 %rotor speed of turbine 2 as second output
% X1=resample(detrend(rotSpeed(end-750*10:1:end,1)'),1,10);
% X2=resample(detrend(rotSpeed(end-750*10:1:end,2)'),1,10);
% X3=resample(detrend(rotSpeed(end-750*10:1:end,1)'),1,10).^2;
% X4=resample(detrend(rotSpeed(end-750*10:1:end,2)'),1,10).^2;

[X1] = resampleedgeeffect(X1,10); 
[X2] = resampleedgeeffect(X2,10);
%[X3] = resampleedgeeffect(X3,10);
%[X4] = resampleedgeeffect(X4,10);
 
%s4=var(X1);
%s5=var(X2); 
%X1=X1./var(X1);
%X2=X2./var(X2);

% X3=X1.^2;
% X4=X2.^2;
%s6=var(X3); %scaling factor of first deterministc states
%s7=var(X4);
%X3=X3./var(X3);
%X4=X4./var(X4);

% X3=resample(rotSpeed(end-beg*10:1:end,1),1,10).^2;
% X4=resample(rotSpeed(end-beg*10:1:end,2),1,10).^2;
% X3=X3./var(X3);
% X4=X4./var(X4);

%% Peparation of input signals depending on wind farm control strategy

if pitchmode==0
    
    steadyyaw=260; %yaw angle steady state (offset)
    
    %U1=detrend(nacelleYaw(end-beg*10:1:end,1)');
    U1=nacelleYaw(end-beg*10:1:end,1)'; %only taking part of the signal to avoid intial transients
    U1=U1-steadyyaw;
    [U1] = resampleedgeeffect(U1,10);
    %s3=var(U1); %scaling factor to scale input signal. it was seen not to
    %have major impact on numerical robustness and was not used as it made
    %the deisgn of the model predictive controller later on cubmersome
    %U1=U1./var(U1); %scaling
    
    Inputs= U1; %variations with relation to steady state yaw angle are taken as inputs
    
   % transformation to make input more linear
   % Inputs=(1-cosd(Inputs))*6;
  
    %scalingfactors=[s1;s2;s3;s4;s5;s6;s7];
    scalingfactors=[1;1;1;1;1;1;1];
    %scalingfactors=0;
    
elseif pitchmode==1
    
%% MBC: Multi-Blade Coordinate transformation
%A directional thrust force  can be accomplished by implementinf MBC
%transformation, and decoupling/proejcting the blade loads in a non
%-rotating reference frame

% As a result, the measured out-of plane blade root bending moments M(t)
% --> [pitch{ij}(index,1);pitch{ij}(index,2);pitch{ij}(index,3)] are
% projected onto a non rotating reference frame --> PITCH 

    Nturb=2;
    Offset=-8.4*2; 

    for index=1:1:length(time1)
        %for each time instant INDEX get (for each turbine ij below) the 3
        %out-of-plane blade root bending moments, corresponding to the
        %three columns given a certain line
       
        for ij=1:1:Nturb
            
            Azimuth=rotorAzimuth(index,ij);

            PITCH=([1/3 1/3 1/3;
                        2/3*cosd(Azimuth+Offset) 2/3*cosd(Azimuth+120+Offset) 2/3*cosd(Azimuth+240+Offset);  
                        2/3*sind(Azimuth+Offset) 2/3*sind(Azimuth+120+Offset) 2/3*sind(Azimuth+240+Offset);])*...
                    [pitch{ij}(index,1);pitch{ij}(index,2);pitch{ij}(index,3)];         
                
            Pitch1(ij)=PITCH(1);
            Pitch2(ij)=PITCH(2);
            Pitch3(ij)=PITCH(3);
            
            %3 Matrixes containing the different bending moments where each
            %line has the turbine number and each column the time instant
            %INDEX
            PPitch1(ij,index)=PITCH(1);
            PPitch2(ij,index)=PITCH(2);
            PPitch3(ij,index)=PITCH(3);

        end
    end
    
    %U1=resample(detrend(PPitch2(1,end-750*10:1:end)),1,10);
    %s3=var(U1);
    %U1=U1./var(U1);
    %U2=resample(detrend(PPitch3(1,end-750*10:1:end)),1,10);
    %s4=var(U2);
    %U2=U2./var(U2);
    U1=PPitch1(1,end-beg*10:1:end);
    [U1] = resampleedgeeffect(U1,10);
    Inputs= [U1];
    
    scalingfactors=[1;1;1;1;1;1;1];
    
end


%% Preparation of the output signals

%the power production of the two turbines are defined as outputs of the wind turbine system 

%OUTPUT POWER (previously using the in built resample and detrend functions)
% Y1=resample(detrend(powerGenerator(end-750*10:1:end,1)*1e-6'),1,10);
% Y2=resample(detrend(powerGenerator(end-750*10:1:end,2)*1e-6'),1,10);
% Y1=detrend(powerGenerator(end-beg*10:1:end,1)');
% Y2=detrend(powerGenerator(end-beg*10:1:end,2)');
 
meanY1=5.485*10^6;%pitch 0
%meanY1=5.352*10^6; %yaw -10
Y1=powerGenerator(end-beg*10:1:end,1)/rho-meanY1;
   
meanY2=0.7728*10^6; %pitch 0
%meanY2=0.9512*10^6; %yaw -10
Y2=powerGenerator(end-beg*10:1:end,2)/rho-meanY2;
 
[Y1] = resampleedgeeffect(Y1*10^-6,10); %rotor speed of turbine 1 as first output
[Y2] = resampleedgeeffect(Y2*10^-6,10);
  
% TRANSFORMATION
%  Y1=Y1.*Inputs.^(1.88);

 %s1=1;
 %s2=1;
 %s1=var(Y1);
 %s2=var(Y2); 
 
 %Y1=Y1./var(Y1);
 %Y2=Y2./var(Y2);

%OUTPUT ROTOR SPEED
% Y1=X1;
% Y2=X2;

Outputs=[Y1;Y2];

    
%% DETERMINISTIC STATES

Deterministic=[X1; X2]; 
%Deterministic={};

%% Storing of scaling factors

% the scaling factors are stored to be used for validation later, to ensure
% that identification and validation data are in the same order of
% magnitude

%scaling factor for output: S1
%scaling factor for output: S2

%sacling factor for input: S3

%scaling factor for deterministic states: S4
%scaling factor of second deterministc state: S5
    
