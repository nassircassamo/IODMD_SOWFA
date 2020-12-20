%% probeRead.m
%% FUNCTION TO READ OpenFOAM PROBE DATA.

% Matthew J. Churchfield
% National Renewable Energy Laboratory (NREL)
% 15013 Denver West Parkway, Golden, CO 80501, USA
% 1.303.384.7080
% matt.churchfield@nrel.gov

% June 10, 2014

% Modified by Eliot Quon - April 7, 2016


function [t,x,y,z,u,nData,nTime] = probeRead(rootDir,varName,components,headerFormat);


%% Get list of times.
t = dir(rootDir);
if length(t)==0
    fprintf('Directory %s not found!\n',rootDir); 
    return
end
    
t = t(3:end); % skip . and .. 
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




%% Loop through time directories.
for i = 1:length(times)
    disp(['Processing time ',num2str(i),' of ',num2str(length(times)), ':  ',times{i},'s...']);    
    
    %% Open the file.
    file = [rootDir,'/',times{i},'/',varName];
    fid = fopen(file,'r');
    if fid < 0
        fprintf(' - Problem opening %s\n',file);
        return
    end
    
    %% Read the header and get x,y,z locations.
    if (headerFormat == 'old')
        xData = fgetl(fid);
        iStart = findstr('x',xData)+1;
        xData = xData(iStart:end);
        x = str2num(xData);
        clear xData;
        
        yData = fgetl(fid);
        iStart = findstr('y',yData)+1;
        yData = yData(iStart:end);
        y = str2num(yData);
        clear yData;
        
        zData = fgetl(fid);
        iStart = findstr('z',zData)+1;
        zData = zData(iStart:end);
        z = str2num(zData);
        clear zData;
    else
        data = fscanf(fid,'# Probe %d (%g %g %g)\n',[4,inf]);
        assert( all( data(:,end)==0 ) );
        x = data(2,1:end-1)';
        y = data(3,1:end-1)';
        z = data(4,1:end-1)';
    end

    nData = length(x);
    %fprintf('  # probes: %d\n',nData');
    
    %% Read the data.
    fgetl(fid); % header 1
    fgetl(fid); % header 2
    formatString = '%g ';
    for j = 1:nData
        if (components > 1)
            formatString = [formatString,'('];
            for k = 1:components
                formatString = [formatString,'%g '];
            end
            formatString = [formatString(1:end-1),') '];
        else
            for k = 1:components
                formatString = [formatString,'%g '];
            end
        end
    end
    formatString = [formatString,'\n'];
    data = fscanf(fid,formatString,[(nData*components) + 1 Inf])';
    
    %% Split the data into components.
    if (i == 1);
        t = data(:,1);
        for j = 1:components
            u{j} = data(:,j+1:components:end);
        end
        tLast = t(end);
    else
        tI = data(:,1);
        iStart = find(tI > tLast,1,'first');
        t = [t; data(iStart:end,1)];
        for j = 1:components
            u{j} = [u{j}; data(iStart:end,j+1:components:end)];
        end
        tLast = t(end);
        clear tI;
    end
    clear data;
        
    
    %% Close the file.
    fclose('all');
end


nTime = length(t);