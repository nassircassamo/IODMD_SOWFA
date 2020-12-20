function []=plotturbines2Dxy(yawAngle,D,yTower,xTower)

 theta=-pi:0.01:pi;
 dxNac = 12/D;
 r_rotor = 90;
 
 x_circ=ones(size(theta))*dxNac/2;
 y_circ=r_rotor*cos(theta);
 z_circ=r_rotor*sin(theta);
 
  RotMat = [cosd(yawAngle), -sind(yawAngle); ... % Rotation matrix (z-axis)
          sind(yawAngle) cosd(yawAngle)];
      
 xy_circ_Yaw = RotMat*[x_circ;y_circ];
    
 p=plot(xy_circ_Yaw(1,:)+xTower, xy_circ_Yaw(2,:)+yTower,'k','linewidth',3);
 p.Color(4)=1;
 
 
 