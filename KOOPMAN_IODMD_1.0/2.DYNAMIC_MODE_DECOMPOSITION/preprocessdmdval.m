function [Inputs, Outputs, Deterministic]=preprocessdmdval(beg, rotSpeed,time1,rotorAzimuth,nacelleYaw,pitchmode, pitch,scalingfactors,powerGenerator,rho)

%% EVALUATE RELEVANT STATES TO BE UED 
% X1=detrend(rotSpeed(end-beg*10:1:end,1)');
% X2=detrend(rotSpeed(end-beg*10:1:end,2)');
% 
% [X1] = resampleedgeeffect(X1,10); %rotor speed of turbine 1 as first output
% [X2] = resampleedgeeffect(X2,10); %rotor speed of turbine 2 as second output
% X1=resample(detrend(rotSpeed(end-750*10:1:end,1)'),1,10);
% X2=resample(detrend(rotSpeed(end-750*10:1:end,2)'),1,10);
% X3=resample(detrend(rotSpeed(end-750*10:1:end,1)'),1,10).^2;
% X4=resample(detrend(rotSpeed(end-750*10:1:end,2)'),1,10).^2;
 meanX1=7.201;
 X1=rotSpeed(end-beg*10:1:end,1)-meanX1;
 meanX2=3.67;
 X2=rotSpeed(end-beg*10:1:end,2)-meanX2;
 
 [X1] = resampleedgeeffect(X1,10); 
 [X2] = resampleedgeeffect(X2,10);

X1=X1./scalingfactors(4);
X2=X2./scalingfactors(5);

%%
% X3=X1.^2;
% X4=X2.^2;

%X3=X3./scalingfactors(6);
%X4=X4./scalingfactors(7);

% X3=resample(rotSpeed(end-beg*10:1:end,1),1,10).^2;
% X4=resample(rotSpeed(end-beg*10:1:end,2),1,10).^2;
% X3=X3./var(X3);
% X4=X4./var(X4);

%% OUTPUTS
%the rotor speeds of the two turbines are defined as outputs of the wind turbine system 
%POWER OUTPUT
   % meanY1=5.485*10^6;%pitch 0
   meanY1=5.352*10^6; %yaw -10
   %meanY1=mean(powerGenerator(300:500,1));
   Y1=powerGenerator(end-beg*10:1:end,1)/rho-meanY1;
    
  % meanY2=0.7728*10^6; %pitch 0
   meanY2=0.9512*10^6; %yaw -10
   %meanY2=mean(powerGenerator(300:500,2));
   Y2=powerGenerator(end-beg*10:1:end,2)/rho-meanY2;
% 
  [Y1] = resampleedgeeffect(Y1*10^-6,10); %rotor speed of turbine 1 as first output
  [Y2] = resampleedgeeffect(Y2*10^-6,10);

% Y1=resample(detrend(powerGenerator(end-750*10:1:end,1)*1e-6'),1,10);
% Y2=resample(detrend(powerGenerator(end-750*10:1:end,2)*1e-6'),1,10);
% 
 Y1=Y1./scalingfactors(1);
 Y2=Y2./scalingfactors(2);
 
%ROTOR SPEED OUTPUT
% Y1=X1;
% Y2=X2;

Outputs=[Y1;Y2];

%% INPUTS: 
if pitchmode==0

    steadyyaw=260;
    %U1=detrend(nacelleYaw(end-beg*10:1:end,1)');
    U1=nacelleYaw(end-beg*10:1:end,1)';
    U1=U1-steadyyaw;
    [U1] = resampleedgeeffect(U1,10);

    U1=U1./scalingfactors(3);

    Inputs= U1;

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
    
%     U1=resample(detrend(PPitch2(1,end-750*10:1:end)),1,10);
%     U1=U1./scalingfactors(3);
%     U2=resample(detrend(PPitch3(1,end-750*10:1:end)),1,10);
%     U2=U2./scalingfactors(4);
%     Inputs= [U1; U2];

    U1=PPitch1(1,end-beg*10:1:end);
    [U1] = resampleedgeeffect(U1,10);
    Inputs= [U1];
    
end
    

%% DETERMINISTIC STATES
Deterministic=[X1; X2]; 
%Deterministic={};
