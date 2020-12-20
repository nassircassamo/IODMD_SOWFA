%% readTurbineOutputSectional.m
%% FUNCTION TO READ TIME SERIES OF ACTUATOR TURBINE SECTIONAL (i.e. lift along the blade span) OUTPUT DATA.

% Matthew J. Churchfield
% National Renewable Energy Laboratory (NREL)
% 15013 Denver West Parkway, Golden, CO 80501, USA
% 1.303.384.7080
% matt.churchfield@nrel.gov

% May 20, 2016

function [nTurbine,nBlade,time,dt,nVal,val] = readTurbineOutputSectional(dirName,varName)

% Get list of times.
t = dir(dirName);
t = t(3:end);
for j = 1:1:1 %JWlength(t)
    times{j} = t(j).name;
end


time = [];
dt = [];
printTimeLast = 0.0;
printTimeIncr = 500.0;
bufferSize = 1E5; %what is this buffer size
for j = 1:length(times)
    disp(['Processing time ',num2str(j),' of ',num2str(length(times)), ':  ',times{j},'s...']);
    
    fileName = [dirName,'/',times{j},'/',varName];
    fid = fopen(fileName,'r');
    
    
    if (j == 1)
        % Figure out how many turbines and blades on each turbine and elements per blade there are.
        fgetl(fid);
        
        t = fgetl(fid);
        nTurbine = 0;
        tNum = -1;
        bNum = -1;
        while (length(t) > 1)
            tNumP = tNum;
            bNumP = bNum;
            data = str2num(t);
            tNum = data(1);
            bNum = data(2);
            if (tNum ~= tNumP)
                nTurbine = nTurbine + 1;
            end
            if (bNum ~= bNumP)
                nBlade(nTurbine) = bNum + 1;
            end
            nVal(nTurbine) = length(data) - 4;
            
            t = fgetl(fid);
        end
        %%edit jW
        
        %why is this necessary ??
        nBlade=nBlade.*ones(nTurbine,1);
        
        %% edit Jw
        frewind(fid);
        for k = 1:nTurbine
            for m = 1:nBlade(k)
                val{k}{m} = [];
                valI{k}{m} = zeros(bufferSize,nVal(k));
            end
        end;
    end
    timeI = zeros(1,bufferSize);
    dtI = zeros(1,bufferSize);
    
    
    % Read the data.
    i = 1;
    flag = 0;
    fgetl(fid);
    while (flag == 0)
        for k = 1:nTurbine
            data = fscanf(fid,'%g',[nVal(k)+4 nBlade(k)])';
            if (isempty(data) == 0)
                timeI(i) = data(1,3);
                dtI(i) = data(1,4);
                for m = 1:nBlade(k)
                    valI{k}{m}(i,1:nVal(k)) = data(m,5:end);
                end
            else
                flag = 1;
                timeI = timeI(1:i-1);
                dtI = dtI(1:i-1);
                valI{k}{m} = valI{k}{m}(1:i-1,:);
            end
            clear data;
        end
        
        if (timeI(min(i,length(timeI))) >= printTimeLast)
            printTimeLast = max((printTimeLast + printTimeIncr),(timeI(i) + printTimeIncr));
            disp(['Time = ',num2str(timeI(i))])
        end
        i = i + 1;
        fgetl(fid);
    end
    
    if (j < length(times))
        ind = find(timeI < str2num(times{j+1}),1,'last');
        timeI = timeI(1:ind);
        dtI = dtI(1:ind);
        for k = 1:nTurbine
            for m = 1:nBlade(k)
                valI{k}{m} = valI{k}{m}(1:ind,:);
            end
        end
    end
    time = [time timeI];
    dt = [dt dtI];
    for k = 1:nTurbine
        for m = 1:nBlade(k)
            val{k}{m} = [val{k}{m}; valI{k}{m}];
        end
    end
    clear timeI dtI valI data;
  
    
    % Close the file.
    fclose('all');
end