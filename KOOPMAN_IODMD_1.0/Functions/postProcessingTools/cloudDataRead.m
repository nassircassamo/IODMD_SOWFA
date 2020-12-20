%% cloudDataRead.m
%% FUNCTION TO READ DATA SAMPLED WITH THE CLOUD POINT FUNCTION OBJECT.

% Matthew J. Churchfield
% National Renewable Energy Laboratory (NREL)
% 15013 Denver West Parkway, Golden, CO 80501, USA
% 1.303.384.7080
% matt.churchfield@nrel.gov

% May 20, 2016

function [ni,nj,nk,xP,yP,zP,uP] = cloudDataRead(fileIn,nComponents,varNum,rotationAngle,rotationAxis);


% This function reads in regular arrays of cloud data, reshapes them, and rotates them.  It assumes it is given cloud data
% that has points listed in an array order stretched out into one long vector.


% fileIn                    File name of input cloud data file.
% nComponents               Number of components that the data type contains.  1 for scalar, 3 for vector, and so on.
% varNum                    Since a cloud data file can contain more than one variable of a certain type, indicate the 
%                           sequence of the variable you are extracting.  For example, if the variables are 'U' and 'omega',
%                           in that order, and you want 'omega', enter 2.
% rotationAngle             Rotate the data into a local coordinate system.
% rotationAxis              Axis about which to do the rotation.
% ni, nj, nk                Number of points in each direction of the array of data.
% xP, yP, zP                The x'-, y'-, and z'-coordinates of the data (the coordinate system after rotation).
% uP                        The data, u', which is in the rotated coordinate system.



%% READ THE CLOUD DATA
% Open the cloud data file.
fid = fopen(fileIn,'r');

% Read the contents and save.
data = fgetl(fid);
frewind(fid);
data = str2num(data);
nElements = length(data);
clear data;

data = fscanf(fid,'%g',[nElements Inf])';
fclose(fid);
xx = data(:,1);
yy = data(:,2);
zz = data(:,3);
uu = data(:,(nComponents*varNum+1):(nComponents*varNum+nComponents));

% Close the cloud data file.
fclose('all');


%% FIGURE OUT ARRAY SIZE AND RESHAPE
% The cloud data is a vector of data.  Turn it into an array.
dx = xx(2) - xx(1);
dy = yy(2) - yy(1);
dz = zz(2) - zz(1);
vi = [dx dy dz];
i = 1;
flag = 0;
while (flag == 0)
    if (length(xx) >= i + 1)
        dx = xx(i+1) - xx(i);
        dy = yy(i+1) - yy(i);
        dz = zz(i+1) - zz(i);
        v = [dx dy dz];
        vDiff = v - vi;
        vDiffMag = sqrt(vDiff(1)^2 + vDiff(2)^2 + vDiff(3)^2);
        if (vDiffMag > 1.0E-2)
            flag = 1;
        end
        i = i + 1;
    else
        flag = 1;
    end
end
ni = i - 1;
if (length(xx) == ni)
    nj = 1;
    nk = 1;
else (length(xx) > ni)
    dx = xx(ni+1) - xx(1);
    dy = yy(ni+1) - yy(1);
    dz = zz(ni+1) - zz(1);
    vi = [dx dy dz];
    flag = 0;
    j = 1;
    while (flag == 0)
        if (length(xx) >= (j*ni) + 1)
            dx = xx((j*ni)+1) - xx(((j-1)*ni)+1);
            dy = yy((j*ni)+1) - yy(((j-1)*ni)+1);
            dz = zz((j*ni)+1) - zz(((j-1)*ni)+1);
            v = [dx dy dz];
            vDiff = v - vi;
            vDiffMag = sqrt(vDiff(1)^2 + vDiff(2)^2 + vDiff(3)^2);
            if (vDiffMag > 1.0E-2)
                flag = 1;
            end
            j = j + 1;
        else
            flag = 1;
        end
    end
    nj = j;
    if (length(xx) == ni*nj)
        nk = 1;
    else
        dx = xx(ni*nj+1) - xx(1);
        dy = yy(ni*nj+1) - yy(1);
        dz = zz(ni*nj+1) - zz(1);
        vi = [dx dy dz];
        flag = 0;
        k = 1;
        while (flag == 0);
            if (length(xx) >= (j*ni*nj) + 1)
                dx = xx((k*ni*nj)+1) - xx(((k-1)*ni*nj)+1);
                dy = yy((k*ni*nj)+1) - yy(((k-1)*ni*nj)+1);
                dz = zz((k*ni*nj)+1) - zz(((k-1)*ni*nj)+1);
                v = [dx dy dz];
                vDiff = v - vi;
                vDiffMag = sqrt(vDiff(1)^2 + vDiff(2)^2 + vDiff(3)^2);
                if (vDiffMag > 1.0E-2)
                    flag = 1;
                end
                k = k + 1;
            else
                flag = 1;
            end
        end
        nk = k;
    end
end


% Reshape the arrays.
ii = 1;
for k = 1:nk
    for j = 1:nj
        for i = 1:ni
            x(i,j,k) = xx(ii);
            y(i,j,k) = yy(ii);
            z(i,j,k) = zz(ii);
            u(i,j,k,:) = uu(ii,:);
            ii = ii + 1;
        end
    end
end

% Clear data from memory that is no longer needed.
clear data xx yy zz uu;



%% ROTATE THE DATA
% The plane is rotated by some angle with respect to being oriented
% north-south.  Transform the x-y coordinates and velocities into the
% rotated frame.
%R = [cosd(rotationAngle) -sind(rotationAngle);...
%    sind(rotationAngle)  cosd(rotationAngle)];
for i = 1:ni
    for j = 1:nj
        for k = 1:nk
            xP(i,j,k) = x(i,j,k);
            yP(i,j,k) = y(i,j,k);
            zP(i,j,k) = z(i,j,k);
            uP(i,j,k,:) = u(i,j,k,:);
            %xyP = R * [x(j,k); y(j,k)];
            %xP(i,j,k) = xyP(1);
            %yP(i,j,k) = xyP(2);
            %zP(i,j,k) = z(j,k);
            
            %uvP = R * [u(j,k); v(j,k)];
            %uP(i,j,k) = uvP(1);
            %vP(i,j,k) = uvP(2);
            %wP(i,j,k) = w(j,k);
        end
    end
end