%% boundaryDataRead.m
%% FUNCTION TO READ OpenFOAM timeVaryingMappedFixedValue BOUNDARY DATA.

% Matthew J. Churchfield
% National Renewable Energy Laboratory (NREL)
% 15013 Denver West Parkway, Golden, CO 80501, USA
% 1.303.384.7080
% matt.churchfield@nrel.gov

% August 10, 2015




function [xx,yy,zz,t,data] = boundaryDataRead(boundaryDataDir,var);


%% Read in the points data.
fid = fopen([boundaryDataDir,'/points'],'r');
for i = 1:16
    fgetl(fid);
end
nPoints = str2num(fgetl(fid));
fgetl(fid);
points = fscanf(fid,'(%g %g %g)\n',[3 nPoints])';
fclose(fid);
nz = length(unique(points(:,3)));
ny = nPoints/nz;
points = reshape(points,ny,nz,3);

% Strip out the y and z points and get a column of z points, then find the
% z points that lie within a certain range.
xx = points(:,:,1);
yy = points(:,:,2);
zz = points(:,:,3);




%% Read in the flow field data.
d = dir(boundaryDataDir);
d(1) = [];
d(1) = [];
d(end) = [];
nDirs = length(d);


% Read in all times and then come up with the sort index since it does not
% necessarily read in the correct order.
for i = 1:nDirs
    tI(i,1) = str2num(d(i).name);
end
[t,timeInd] = sort(tI);
clear tI;


% Loop over all times.
for i = 1:nDirs
    disp(['Reading data for time ',d(i).name,'...']);
      
    % Open the boundary data
    fid = fopen([boundaryDataDir,'/',d(i).name,'/',var],'r');
    
    % Figure out the data size.
    for j = 1:20
        fgetl(fid);
    end
    nPoints = str2num(fgetl(fid));
    fgetl(fid);
    
    dataTest = fgetl(fid);
    indStart = findstr(dataTest,'(');
    indEnd = findstr(dataTest,')');
    if (isempty(indStart))
        nComponent = 1;
    else
        dataTest = dataTest(indStart+1:indEnd-1);
        dataTest = str2num(dataTest);
        nComponent = length(dataTest);
    end
    clear dataTest;
    
    frewind(fid);
    for j = 1:22
        fgetl(fid);
    end

    if (i == 1)
        data = zeros(ny,nz,nComponent,nDirs);
    end
    
    % Do the actual read of the data and store to a data matrix.
    if (nComponent == 1)
        formatString = '%g\n';
    else
        formatString = '(';
        for j = 1:nComponent-1
            formatString = [formatString,'%g '];
        end
        formatString = [formatString,'%g)\n'];
    end
    dataI = fscanf(fid,formatString,[nComponent nPoints])';
    dataI = reshape(dataI,ny,nz,nComponent);
    dataUnsrt(:,:,:,i) = dataI;
    clear dataI;
        
    % Close the data file.
    fclose(fid);
end



%The data is not in the correct order, so sort it.
data = zeros(size(dataUnsrt));
for i = 1:nDirs
    data(:,:,:,i) = dataUnsrt(:,:,:,timeInd(i));
end
clear dataUnsrt;
