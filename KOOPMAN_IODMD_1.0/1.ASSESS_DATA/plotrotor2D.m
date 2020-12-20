function [p]=plotrotor2D(yawAngle,D,yTower,zHub)

 theta=-pi:0.01:pi;
 o=length(theta);
 
 
 dxNac = 12/D;
 r_rotor = 90;
 
 x_circ=ones(size(theta))*dxNac/2;
 y_circ=r_rotor*cos(theta);
 z_circ=r_rotor*sin(theta);
 z=20*r_rotor*ones(1,o);
 
 
  RotMat = [cosd(yawAngle), -sind(yawAngle); ... % Rotation matrix (z-axis)
          sind(yawAngle) cosd(yawAngle)];
      
 xy_circ_Yaw = RotMat*[x_circ;y_circ];
    
 p=plot3(xy_circ_Yaw(2,:)+yTower,z_circ+zHub,z,'k','linewidth',3);
 %p.Color(4)=1;
 
end