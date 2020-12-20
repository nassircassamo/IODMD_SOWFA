function [p]=plotturbinefromabove(yawAngle,xcenter, ycenter, D)

coordinateone=[xcenter ycenter+D/2];
coordinatetwo=[xcenter ycenter-D/2];


RotMat = [cosd(yawAngle), -sind(yawAngle); ... % Rotation matrix (z-axis)
          sind(yawAngle) cosd(yawAngle)];
      
coordinateonerot=RotMat*coordinateone';
coordinatetworot=RotMat*coordinatetwo';
      
xcenterum=coordinateonerot(1);
xcenterdois=coordinatetworot(1);

ycenterum=coordinateonerot(2);
ycenterdois=coordinatetworot(2);

p=plot3([xcenterum, xcenterdois],[ycenterum, ycenterdois],[3 3],'k','LineWidth',3);