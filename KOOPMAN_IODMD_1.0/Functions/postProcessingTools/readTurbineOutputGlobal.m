%% readTurbineOutputGlobal.m
%% FUNCTION TO READ TIME SERIES OF ACTUATOR TURBINE GLOBAL (i.e. POWER, THRUST, etc.) OUTPUT DATA.

% Matthew J. Churchfield
% National Renewable Energy Laboratory (NREL)
% 15013 Denver West Parkway, Golden, CO 80501, USA
% 1.303.384.7080
% matt.churchfield@nrel.gov
% May 20, 2016

function [nTurbine,time,dt,nVal,val] = readTurbineOutputGlobal(dirName,varName)

% Get list of times.
t = dir(dirName);
t = t(3:end);
for j = 1:1%JW length(t)
    times{j} = t(j).name;
end


time = [];
dt = [];
val = [];
printTimeLast = 0.0;
printTimeIncr = 500.0;
for j = 1:1%length(times)
    disp(['Processing time ',num2str(j),' of ',num2str(length(times)), ':  ',times{j},'s...']);
    
    fileName = [dirName,'/',times{j},'/',varName];
    fid = fopen(fileName,'r');
    
    
    if (j == 1)
        % Figure out how many turbines there are.
        fgetl(fid);
        t = fgetl(fid);
        nTurbine = 0;
        while (length(t) > 1)
            t = fgetl(fid);
            nTurbine = nTurbine + 1;
        end
        frewind(fid);
        
        % Figure out how wide the data is.
        fgetl(fid);
        t = fgetl(fid);
        t = str2num(t);
        nVal = length(t) - 3;
        frewind(fid);
    end
    
    % Read the data.
    i = 1;
    fgetl(fid);
    data = fscanf(fid,'%g',[nVal+3 nTurbine])';
    while (isempty(data) == 0)
        timeI(i) = data(1,2);
        dtI(i) = data(1,3);
        valI(i,:,:) = data(:,4:end);
        fgetl(fid);
        data = fscanf(fid,'%g',[nVal+3 nTurbine])';
        if (timeI(min(i,length(timeI))) >= printTimeLast)
            printTimeLast = max((printTimeLast + printTimeIncr),(timeI(i) + printTimeIncr));
            disp(['Time = ',num2str(timeI(i))])
        end
        i = i + 1;
    end
    clear data;
    
    if (j < length(times))
        ind = size(timeI,2);%find(timeI < str2num(times{j+1}),1,'last');
        timeI = timeI(1:ind);
        dtI = dtI(1:ind);
        valI = valI(1:ind,:,:);
    end
    time = [time timeI];
    dt = [dt dtI];
    val = [val; valI];
    
    clear timeI dtI valI data;
    
    % Close the file.
    fclose('all');
end