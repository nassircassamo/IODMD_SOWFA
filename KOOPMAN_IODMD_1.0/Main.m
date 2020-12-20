%Model Predictive Control For Wake Steering: a Koopman Dynamic Mode
%Decomposition Approach
%Master Thesis Dissertation
%Author: Nassir Rodrigues Cassamo
%Supervisors: Professor Jan-Willem Van Wingerden and Professor João Sousa
%First official commit (1.0): December 2020

%% RELEVAT INFORMATION
%This script has the following goals:
% (1) Assess Simulation data, both qualitatively (animations) and
% quantitatively (graphics)
% (2) Derives a low dimensional model using Dynamica Mode Decomposition
% (variations included to take into account input-output data and other
% known states - deterministic states)
% (3) Validates the models with a set of validaiton data
% (4) Analyses the models, mainly the intrinsic dynamics  modes,
% damping, frequencies, energy) 
% (5) Reconstructs flow based on models an computes deviations from real 
% (6) Designs a Model Predictive Control 

%This script requires the following functions to be added to MATLAB's path
% (1) Functions (LTI toolbox from DCSC, post processing tools from NRL
% (2) cbrewer: color scaling
% (3) altmany-export_fig: to export figures directly to specified directories
% (4) subfolders comprising all parts sub functions

%This script requires data in the following fashion
% (1) 2 folders of name (2.1) steps_yaw and (2.2) steps_yaw_val, that
% contain the data non processed directly from CFD simulation from SOWFA
% (2) Data vectors with post processed information: the post processing is
% performed in the cluster and these vectors contain the resampled
% flowfield and resampled grid points. The grid points were resampled every
% Decimate (variable should be in these data vectors)
    % (2.1) U_data_complete_vec
    % (2.2) U_data_complete_vec_val

%% (0) INITIALISE
    %define type of simulation
    part=0; subpart=1; [f]= MPC_progress(part,subpart,{},{},{});

    clc
    close all
    p=genpath('Functions');
    addpath(p)
    
    maindir='/Volumes/NASSIR/MATLAB/'; %DEFINE MAIN DIRECTORY IN USERS COMPUTER TO STORE ALL RESULTS
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{});
    pitchmode=0; %0 for yaw control and 1 for pith control
    
    if pitchmode==0
        dirName={'steps_yaw_20deg_10offset'}; %directory for identification data, wake redirection control
        dirName_val={'steps_yaw_20deg_10offset_val'}; %directory for validation data, wake redirection control
    elseif pitchmode==1
        dirName={'steps_theta_col_new';}; %directory for identification data, collective pitch control
        dirName_val={'steps_theta_col_new_val'}; %directory for validation data, collective pitch control
    end
    
    detrendingstates=0; %1 to take mean flow and consider turbulent fluctuations
    method=2; %0: DMD ; 1:DMDc; 2:IODMD; 3:EXTIODMD
    koopman=0; %to add deterministic states to flow field data
    videos=0; %generate videos
    snapshots=0; %generate snapshots from simulation data
    
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{});
    % Turbine and flow characteristics to be used 
    rho=1.225; %air density in [kg m^-3]
    D=178; %Rotor Diameter used in simulations= 178 [m]
    % Simulation characterisitc (resampling)
    dt=2; %time sampling
    

%% (1) ASSESS DATA
    part=1; subpart=1; [f]= MPC_progress(part,subpart,f,{},{});

    subpart=2; [f]= MPC_progress(part,subpart,f,{},{});
    if pitchmode==0 %wake redirection control using nacelle yaw angle
         analysis='YAW_MPC_offset'; %name of directory to be created to automatically store results
         filename='U_data_complete_vec_yaw_off.mat'; %directory for matlab file with flow field identification data set
         filenamevalid='U_data_complete_vec_yaw_off_val.mat'; %directory for matlab file with flow field validation data set
         load(filename) ;
         valid=load(filenamevalid);
         
         %for not using all grid points and only part of them (example,
         %only between first and second turbine)
       % [xxx,yyy,zzz,XX,YY,ZZ,QQ_u]=retakepoints(QQ_u,x,y,z,Decimate);
       % [xxx,yyy,zzz,XX,YY,ZZ,valid.QQ_u]=retakepoints(valid.QQ_u,x,y,z,Decimate);
    
         %easy solution to augment u flow field data matricx with other flow
        %field data
        %QQ_u=[double(QQ_u);double(QQ_w)];
        %valid.QQ_u=[double(valid.QQ_u);double(valid.QQ_w)];
     
    elseif pitchmode==1 %collective pitch control was used
        analysis='PULSE_MPC/'; %name of directory to be created to automatically store results
        filename='U_data_complete_vec_pulse.mat'; %directory for matlab file with flow field identification data set
        filenamevalid='U_data_complete_vec_pulse_val.mat'; %directory for matlab file with flow field validation data set
        load(filename) ;
        valid=load(filenamevalid);
        
        %use only some points of the grid 
        %[xxx,yyy,zzz,XX,YY,ZZ,QQ_u]=retakepoints(QQ_u,x,y,z,Decimate);
        %[xxx,yyy,zzz,XX,YY,ZZ,valid.QQ_u]=retakepoints(valid.QQ_u,x,y,z,Decimate);
   
        %easy solution to augment u flow field data matricx with other flow
        %field data instead of using koopmanextension function
        %QQ_u=[double(QQ_u);double(QQ_w)];
        %valid.QQ_u=[double(valid.QQ_u);double(valid.QQ_w)]; %perform same
        %grid resmplaing on validation for comparison purposes later on
         
    end
    
    maindir=strcat(maindir,analysis);    %define main directory
    
    % IDENTIFICATION DATA
    visualisefirstresults(dirName,rho,0,maindir,'SOWFA_id') %0: skip this; %1: see graphs
    % VALIDATION DATA
    visualisefirstresults(dirName_val,rho,0,maindir,'SOWFA_val') %0: skip this; %1: see graphs
 
    % Make movie from results (identification data, steady state (1500-3000)
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    if videos==1
        %second input argument: dirctory to save movies
        [dirpathwake,dirpathveldef]=makeframes(D,strcat(maindir,'MovieWake'),dirName,pitchmode,filename,rho);
        [dirpathfinal]             =makefinalframes(D,strcat(maindir,'Movie_combined'),filename);
        [dirpathpowerinsights]     =makeframespowerinsights(D,rho,strcat(maindir,'Power_insights'),filename);
        [dirpathvelfield]          =plotvecfield(D,strcat(maindir,'power_velfield'),dirName,pitchmode,filename);
        [dirpathcuthubheightvec]   =cuthubheightvec(D,strcat(maindir,'cuthubheight'),filename);
        
        %Movie making. Please specify Frames Per Second (FPS) inside
        makemoviewake('wake_deflection_yawcontrol_sowfa',dirpathwake,rho);
        makemoviewake('wake_deflection_velocitydefifcitslice_yawcontrol_sowfa',dirpathveldef);
        makemoviewake('wake_deflection_combined',dirpathfinal);
        makemoviewake('wake_deflection_powerdelta',dirpathpowerinsights);
        makemoviewake('wake_deflection_velfield',dirpathvelfield);
        makemoviewake('wake_deflection_hubcut',dirpathcuthubheightvec);
        
        movietomp4() %Converts avi video to mp4. Also specify FPS inside
     else
     end
     
    % Make figures to assess wake evolution (snapshots in sequential time
    % instants)
    subpart=4; [f]= MPC_progress(part,subpart,f,{},{}); 
    if snapshots==1
        if pitchmode==0
            wake_vorticity_deflection_yaw(dirName,filename,D,150) %snapshots of wake delfection (three dimensional, vorticity)
            hubheightcut_yaw(dirName,filename,D,150) %snapshots of wake deflection at hub height (velocity)
          
        elseif pitchmode==1
            wake_vorticity_pitch(dirName,filename,D,400) %snapshots of wake contraction (three dimensional, vorticity)
            vefield5D1(dirName,filename,D,400) %snapshots of wake shape taken at several slices downstream
        end
    else 
    end
 
%% (2) DYNAMIC MODE DECOMPOSITION 
part=2; subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
   
    itsf=921; %instant to start from, as certain sample time
    beg=(10001-itsf)/10; %instant to begin defined according to length of data
    %begin=750;
    %beg=750;
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Read and process identification data
    [rotSpeed, nacelleYaw, time1,rotorAzimuth,pitch,powerGenerator]=readdmdinformation(dirName); %read information from simulation
    [Inputs, Outputs, Deterministic,scalingfactors]=preprocessdmdid(beg, rotSpeed,time1,rotorAzimuth,nacelleYaw, pitchmode,pitch,powerGenerator,rho ); %preprocess information (resample, and maintain only relevant data)
    
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Read and process validation data
    [rotSpeed_val, nacelleYaw_val, time1_val,rotorAzimuth_val,pitch_val,powerGenerator_val]=readdmdinformation(dirName_val); %read information from simulation
    [Inputs_val, Outputs_val, Deterministic_val]=preprocessdmdval(beg, rotSpeed_val,time1_val,rotorAzimuth_val,nacelleYaw_val,pitchmode,pitch_val,scalingfactors,powerGenerator_val,rho); %preprocess information (resample and only relevant data)

    % Define states to be used for DMD
    %states=QQ_u(:,(begin-beg)+1:end); % define states: first hypothesis 
    states=QQ_u(:,(itsf-1)*0.1:end); %fluid flow as states, identification data set
    statesvalid=valid.QQ_u(:,(itsf-1)*0.1:end); %fluid flow as states, validaiton data set for comparison
    n=size(states,1);
    
    subpart=4; [f]= MPC_progress(part,subpart,f,{},{}); 
    if detrendingstates
        [states,meansteadystate,scalingfactor]=preprocessstates(states); %remove meanflow or other pre processing techniques to experiment
     %  [statesvalid,meansteadystate,scalingfactor]=preprocessstates(statesvalid);
    else
    end
    
    %include non linear observables - Koopman extensions to better recover
    %non linear dynamics
    if koopman
        [nonlobs]=koopmanstateextension(QQ_u, double(QQ_v), double(QQ_w),rho);
        states=[nonlobs];
    else
    end
    
    subpart=5; [f]= MPC_progress(part,subpart,f,{},{}); 
    r=100; %define truncation level for Singular Value Decomposition 
    
    [sys_red,FITje,U,S,V,method,X,X_p,Xd,dirdmd,xstates]=dynamicmodedecomposition(states,Inputs, Outputs, Deterministic,method,r,maindir,f,dt); 
    save(strcat(dirdmd,'/OPTIONS.mat'),'detrendingstates','method','koopman','rho','D','dt','dirName','dirName_val');

%% (3) DATA VALIDATION 
    part=3; subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    % Validate Models from validation data set
    [FITje_val,dirdmd_val,xstatesvalid]=validatemodels(sys_red,Inputs_val,Outputs_val,r,strcat(dirdmd, '/val'),f,statesvalid,U,Deterministic_val,method);
    save(strcat(dirdmd,'/FIT.mat'),'FITje_val','FITje');
     
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    [modelVAF_val]=idvaloverview(FITje,FITje_val,dirdmd,'VAFidandval');  %overview of models results (identification and validation)
    [a,b]=max(FITje_val(2,1:50)) %best performing model, only analysing first 50
    FITje_val(1,b)
   
    if detrendingstates
        save(strcat(dirdmd,'/RESULTS.mat'),'sys_red',...
            'X','X_p','Xd','Inputs','Outputs','Deterministic','Inputs_val',...
            'Outputs_val','Deterministic_val','U','S','V','dt','r',...
            'dirdmd','method','meansteadystate','scalingfactor','x','y','z',...
            'D','xstates','xstatesvalid','n','Decimate','states','statesvalid');
    else
        save(strcat(dirdmd,'/RESULTS.mat'),'sys_red',...
            'X','X_p','Xd','Inputs','Outputs','Deterministic','Inputs_val',...
            'Outputs_val','Deterministic_val','U','S','V','dt','r',...
            'dirdmd','method','x','y','z','D','xstates',...
            'xstatesvalid','n','Decimate','scalingfactors','states','statesvalid');
    end
    
%% (4) DYNAMICAL ANALYSIS
    part=4; subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    [maxval,modeltouse]=max(FITje_val(2,1:100)); %model to analyse is best performing model out of first 100
    %xo_val=U(1:size(statesvalid,1),1:modeltouse)'*statesvalid(:,1); %initial condition computed based on  low dimension representation 
    %[ysim_val, t, xout]=lsim(sys_red{modeltouse}, [Inputs_val]',[],xo_val);  ysim_val=ysim_val';
   
    % [modalmodels,total,cmodalssd,xmodal]=realdiag(sys_red{modeltouse},Inputs,Outputs,D,9,dt,dirdmd); %NOT SUITABLE FOR USE. Creates state space models for each dynamic mode to evaluate its individual contribution
    %modalanalyse3d(xmodal,[3],x,y,z,Decimate,U,D,modeltouse,modalmodels,dirdmd,Inputs,scalingfactors) %NOT SUITABLE FOR USE. Animation of the different dynamic modes
    % frequencyanalysis(sys_red{modeltouse},modalmodels,D,9,dirdmd)
    [freq,LambdaDiag, P, phi,damping,b]=dynamicalanalysis(sys_red, U, S,V, dt,X_p,X, method,modeltouse,0,D,9,Deterministic,r,dirdmd,n,Xd,modeltouse); %dynamical analysis as proposed by DMD methodology
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{});
    visualisepodmodes(phi,freq, P,x,y,z,Decimate,D,LambdaDiag,damping,method,Xd,dirdmd) %two dimensional POD modes 
    podmodes3dim(x,y,z,Xd, P,phi,fnatural,damping,Decimate,D,dirdmd) %three dimensional POD modes (high dimensional state space eigenvectors)
    %modeanimation([4],x,y,z,P,phi,freq,damping,LambdaDiag, b, Decimate,D,dirdmd,f,1,scalingfactor,meansteadystate,Xd,X) %animation for three dimensional POD mode
%% (5) REBUILD FLOW FIELD AND ASSESS DEVIATIONS
    part=5;subpart=1; [f]= MPC_progress(part,subpart,f,{},{}); 
    
    %project low dimension state to high order representation according to
    %DMD method used
    if method==0
        [statesrebuild]=rebuild(phi,b,LambdaDiag,r,X,Xd); %rebuild states with highest order model when no external forcing is done
    else
        model=length(sys_red);
        if isempty(Xd)
            statesrebuild=U(1:n,1:modeltouse)*xstates{modeltouse}';
            statesrebuildvalid=U(1:n,1:modeltouse)*xstatesvalid{modeltouse}';
        else
            statesrebuild=U(1:n,1:modeltouse)*(xstates{modeltouse}(:,size(Xd,1)+1:end))';
            statesrebuildvalid=U(1:n,1:modeltouse)*(xstatesvalid{modeltouse}(:,size(Xd,1)+1:end))';
        end
    end
    
    %computation of true flow field by taking into account detrending and
    %rescaling that might have been used during identification
    if detrendingstates
        for i=1:size(statesrebuild,2)
            statesrebuild(:,i)=statesrebuild(:,i)*scalingfactor+meansteadystate;
            statesrebuildvalid(:,i)=statesrebuildvalid(:,i)*scalingfactor+meansteadystate;
        end
    else
    end
    
    subpart=2; [f]= MPC_progress(part,subpart,f,{},{}); 
    
    %visual comparison on dmd predictions and true flow field. field should
    %be chaneged according to data ste used
    if pitchmode==0
        comparereconstruction(states(1:size(QQ_u,1),:), statesrebuild(1:size(QQ_u,1),:),D,dirdmd,x,y,z,Decimate,dirName,'reconsversusSOWFA_id',405,Inputs) %comparison at hub height plane
        comparereconstruction(statesvalid(1:size(QQ_u,1),:), statesrebuildvalid(1:size(QQ_u,1),:),D,dirdmd,x,y,z,Decimate,dirName_val,'reconsversusSOWFA_valid',450,Inputs_val)
        comparisonyaw(statesvalid,statesrebuildvalid,Inputs_val,Outputs_val,x,y,z,Decimate,dirdmd,D,ysim_val)
    elseif pitchmode==1  
        comparereconstruction_pitch(statesvalid(1:size(QQ_u,1),:), statesrebuildvalid(1:size(QQ_u,1),:),D,dirdmd,x,y,z,Decimate,dirName_val,'reconsversusSOWFA_valid',450,Inputs_val) %comparison at dosntream rotor plane
        %threeDcomparison(sys_red, modeltouse, Inputs, Outputs, x, y, z, Decimate, D, QQu, flowdmd,dir)
        %comparisonrotorplanes(sys_red, modeltouse, Inputs, Outputs, x, y, z, Decimate, D, QQ_u, flowdmd,dirdmd)
    end
    
    subpart=3; [f]= MPC_progress(part,subpart,f,{},{}); 
    evauatemodelerror(states, statesrebuild,D,dirdmd,filename,dirName,x,y,z,Decimate)
    
    %evaluate time varying error in simualtions according to normalised
    %root mean squared error at each grid lcoation
    [nrmse, nrmsevalid]=evaluatetimevaryingerror(states(1:size(QQ_u,1),:),statesrebuild(1:size(QQ_u,1),:),statesvalid(1:size(QQ_u,1),:), statesrebuildvalid(1:size(QQ_u,1),:),dirdmd,'statestimeerror');
    save(strcat(dirdmd,'/RECONSTRUCTION.mat'),'nrmse','nrmsevalid');
    meannrmse=mean(nrmse)
    meannrmsevalid=mean(nrmsevalid)

    subpart=4; [f]= MPC_progress(part,subpart,f,{},{}); 
        
    close all
    close(f)
%% (6) MODEL PREDICTIVE CONTROL DESIGN 
    freq=1/dt;
    lpf=ss(freq, 1,freq,0,2);
    sys_red_fil=series(lpf, sys_red{modeltouse});
    
   % [u]=power_referencetracking(sys_red_fil,600,600,Inputs_val,Outputs_val,scalingfactors,dirdmd,1000,1);
   % [u]=power_referencetracking_iio(sys_red_fil,500,500,Inputs_val,Outputs_val,scalingfactors,dirdmd);
    %[yf1, yf2]=evaluatepredictionpower(sys_red, modeltouse,Inputs_val, Outputs_val,dt, [1 250],400,200);
    
%% (7) FINAL RESULTS
    Preffnal=generate_reference(1000, 3000,300);
    dirfinalturbine={'/Volumes/NASSIR/MATLAB/FINAL_RESULTS/MPC_Nassir_10_1_higher_setpoint'};
    %dirfinalflow=;
    
    evaluatefinalresults(Preffinal,dirfinalturbine)
   
    