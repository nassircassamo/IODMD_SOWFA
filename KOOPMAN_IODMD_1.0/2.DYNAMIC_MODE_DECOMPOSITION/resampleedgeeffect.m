function [resampled_variable] = resampleedgeeffect(variable,scale)


%% resample it according to scaling
count=0;

for i=1:scale:length(variable)
    count=count+1;
    resampled_variable(count)=variable(i);
end




end

