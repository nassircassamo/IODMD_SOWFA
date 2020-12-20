function h = plotTurbine3D_yaw(azimuthAngle,yawAngle,xTower,yTower,zHub,D)
if nargin <= 0
    azimuthAngle = 0; % degrees
end
if nargin <= 1
    yawAngle = 10; % degrees
end
if nargin <= 4
    xTower = 300;
    yTower = 500;
    zHub = 90; 
end

yawAngle = yawAngle + 180; % Flip to have hub for yaw=0 pointing towards negative x (assumed incoming flow x-positive)
r_rotor = 90/D; % Rotor dimensions
r_tower = 2/D; % Tower dimensions

% Nacelle dimensions
dxNac = 12/D;
dyNac = 12/D;
dzNac = 4/D;
    
% Tower dimensions
zTower = zHub - dzNac/2;

h = gcf; hold on;
% if length(findobj('Tag','blade'))==0
RotMat = [cosd(yawAngle), -sind(yawAngle); ... % Rotation matrix (z-axis)
          sind(yawAngle) cosd(yawAngle)];

% Plot tower
[X,Y,Z] = cylinder(r_tower);
Z = Z*zTower;
X = X+xTower;
Y = Y+yTower;
surf(X,Y,Z);
hold on;

% Plot nacelle
cNac = [150]; % Color of nacelle
% Bottom & top patches
patch1_Yawed = RotMat*[[-1 -1 1 1]*dxNac/2; [-1 +1 +1 -1]*dyNac/2];
patch(xTower+patch1_Yawed(1,:),yTower+patch1_Yawed(2,:),zTower*ones(1,4),cNac); % Bottom 
patch(xTower+patch1_Yawed(1,:),yTower+patch1_Yawed(2,:),(zTower+dzNac)*ones(1,4),cNac); % Top 

% Front & back patches
patch2_Yawed = RotMat*[[-1 -1 -1 -1]*dxNac/2;[-1 -1 +1 +1]*dyNac/2];
patch3_Yawed = RotMat*[[+1 +1 +1 +1]*dxNac/2;[-1 -1 +1 +1]*dyNac/2];
patch(xTower+patch2_Yawed(1,:),yTower+patch2_Yawed(2,:),zTower+[0 1 1 0]*dzNac,cNac); 
patch(xTower+patch3_Yawed(1,:),yTower+patch3_Yawed(2,:),zTower+[0 1 1 0]*dzNac,cNac); 

% Side patches
patch4_Yawed = RotMat*[[-1 -1 +1 +1]*dxNac/2;[-1 -1 -1 -1]*dyNac/2];
patch5_Yawed = RotMat*[[-1 -1 +1 +1]*dxNac/2;[+1 +1 +1 +1]*dyNac/2];
patch(xTower+patch4_Yawed(1,:),yTower+patch4_Yawed(2,:),zTower+[0 1 1 0]*dzNac,cNac); 
patch(xTower+patch5_Yawed(1,:),yTower+patch5_Yawed(2,:),zTower+[0 1 1 0]*dzNac,cNac);     

%     % Plot rotor circle
theta=-pi:0.01:pi;
x_circ=ones(size(theta))*dxNac/2;
y_circ=r_rotor*cos(theta);
z_circ=r_rotor*sin(theta);
xy_circ_Yaw = RotMat*[x_circ;y_circ];
plot3(xy_circ_Yaw(1,:)+xTower,xy_circ_Yaw(2,:)+yTower,z_circ+zHub,'k','linewidth',3);
% grid on
% axis equal tight
% xlabel('x (m)');
% ylabel('y (m)');
% zlabel('z (m)');
zlim([0 zHub+r_rotor])
% else
% %     No need to replot everything, simply:
%     delete(findobj('Tag','blade')); % Delete blade objects
% end
% Plot rotor blades
nBlades = 3;
nInterp = 40;

bladeAngles = linspace(azimuthAngle,azimuthAngle+360,nBlades+1);
bladeAngles = bladeAngles(1:end-1); % remove last entry

for j = 1:1:length(bladeAngles)
    x_blade = ones(1,nInterp)*dxNac/2;
    y_blade = linspace(0,sind(bladeAngles(j))*r_rotor,nInterp);
    z_blade = linspace(0,cosd(bladeAngles(j))*r_rotor,nInterp);
    xy_blade_Yaw = RotMat*[x_blade;y_blade];
    plot3(xTower+xy_blade_Yaw(1,:),yTower+xy_blade_Yaw(2,:),z_blade+zHub,'k.','markerSize',4,'Tag','blade'); 
    hold on
end
hold off