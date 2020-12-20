function [xxx,yyy,zzz,XX,YY,ZZ,flowreshaped]=retakepoints(flow,x,y,z,Decimate)

[xx,yy,zz]=resamplegrid(x,y,z, Decimate);
X = length(xx);
Y = length(yy);
Z = length(zz);

yyy=yy(1:1:end-1); %yaw
xxx=xx(1:1:end-55); %yaw
zzz=zz(1:1:23); %yaw

% yyy=yy(4:1:end-2); %pitch
% xxx=xx(1:1:end-55); %pitch
% zzz=zz(1:1:23); %pitch

XX=length(xxx);
YY=length(yyy);
ZZ=length(zzz);

%% YAW
 for t=1:size(flow,2)
     UmeanAbs_sh_u{t} = reshape(double(flow(:,t)),Y,X,Z);
     Reconstructedflow=UmeanAbs_sh_u{t}(1:1:end-1,1:1:end-55,1:1:23);
     flowreshaped(:,t)=reshape(Reconstructedflow,[XX*YY*ZZ,1]);
 end

%% PITCH
% for t=1:size(flow,2)
%     UmeanAbs_sh_u{t} = reshape(double(flow(:,t)),Y,X,Z);
%     Reconstructedflow=UmeanAbs_sh_u{t}(4:1:end-2,1:1:end-55,1:1:23);
%     flowreshaped(:,t)=reshape(Reconstructedflow,[XX*YY*ZZ,1]);
% end



