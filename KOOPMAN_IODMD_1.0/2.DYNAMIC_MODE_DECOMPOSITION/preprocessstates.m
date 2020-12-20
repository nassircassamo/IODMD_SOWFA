function [states,meansteadystate,scalingfactor]=preprocessstates(states)

    %find mean of steady state
    %assuming stady state is from 
     steadystates=states(:,1:5);
 
     for i=1:size(steadystates,1)
         meansteadystate(i)=mean(steadystates(i,:));
     end
     
     meansteadystate=meansteadystate';
 
     %subtract the steady state dynamics from states
 
     for l=1:size(states,2)
         states(:,l)=states(:,l)-meansteadystate;
     end

% states=detrend(states);
%scalingfactor=var(var(states));
scalingfactor=1;
states=states./scalingfactor;

states=states.^2;
end

