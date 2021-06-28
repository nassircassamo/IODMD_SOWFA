function [Pref]=generate_reference(Tstart, Tend,Hp)
%
%Generate a power reference for the total wind farm. 
%he function should be used for simulation from SOWFA, and it received as
%main arguments:
%Tstart: Time to start the simulation. Will usually b 0
%Tend: Time instantto end the simulation. 
%Hp: Prediction Horizon to be used in MPC, defined in the main script

%For example, given Tstart=0 and Tend=5000, the total time 
%interval will be equally divided into 5 segments of increasing 0.5 MW.
%As the sampling time is 2 seconds, the total number of reference points
%will be equal to Tend/2 plus the prediction Horizon.

%% Generate time axis
ti=Tstart;       		  % starting time
tstep=1;    		      % sample interval
tf = Tend+Hp*tstep;	   	  % final instant
time=ti:tstep:Tend;

l=length(time)-1;
le=l*0.2;

%% Generate reference for tracking
 Pref_wt(:,1) = [0*ones(le,1); 0*ones(le,1);   0*ones(le,1); 0*ones(le,1);  0*ones(le+Hp+1,1)];
 Pref_wt(:,2) = [0*ones(le,1); 0.1*ones(le,1); 0.2*ones(le,1); 0.4*ones(le,1); 0.5*ones(le+Hp+1,1)];
 
 %generate reference for total power in wind farm
 for b=1:length(Pref_wt)
     Pref(b,1)=Pref_wt(b,1)+Pref_wt(b,2);
 end