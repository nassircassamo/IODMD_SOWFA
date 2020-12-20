function [nTime,time,timeNames] = getDirectoryTimes(inputDir);


t = dir(inputDir);
t = t(3:end);
nTime = length(t);
for j = 1:nTime
    timeNamesI{j} = t(j).name;
    time(j,1) = str2num(t(j).name);
end
[time,ind] = sortrows(time);
for j = 1:nTime
    timeNames{j} = timeNamesI{ind(j)};
end
