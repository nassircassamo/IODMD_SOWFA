function [time,dt,heights,nTime,nHeights,val] = planarDataRead(dirName,varName)



% Get list of times.
t = dir(dirName);
t = t(3:end);
for j = 1:length(t)
    timesI{j} = t(j).name;
end

for j = 1:length(timesI)
    timeNum(j) = str2num(timesI{j});
end
[timeNum,sortIndex] = sort(timeNum);
for j = 1:length(timesI)
    times{j} = timesI{sortIndex(j)};
end



time = [];
dt = [];
val = [];
printTimeLast = 0.0;
printTimeIncr = 10.0;
for j = 1:length(times)
    disp(['Processing time ',num2str(j),' of ',num2str(length(times)), ':  ',times{j},'s...']);
    
    
    % If on the first time directory, read the list of heights.
    if (j == 1)
        fileName = [dirName,'/',times{j},'/hLevelsCell'];
        fid = fopen(fileName,'r');
        heights = fscanf(fid,'%g',[inf]);
        nHeights = length(heights);
        fclose(fid);
    end
    
    
    % Now read in the real planar data.
    fileName = [dirName,'/',times{j},'/',varName];
    fid = fopen(fileName,'r');
    data = fscanf(fid,'%g',[nHeights+2 inf])';
    fclose(fid);
    if (j == length(times))
        time  = [time; data(:,1)];
        dt = [dt; data(:,2)];
        val = [val; data(:,3:end)];
    else
        t_int = data(:,1);
        [dd tindex] = min(abs(t_int - str2num(times{j+1})));
        time  = [time; data(1:tindex,1)];
        dt = [dt; data(1:tindex,2)];
        val = [val; data(1:tindex,3:end)];
    end
end


nTime = length(time);