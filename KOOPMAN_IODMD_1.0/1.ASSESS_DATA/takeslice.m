function [dirpath]=takeslice(i, n, Xm_d, Ym_d, Zm_d, UmeanAbs_sh_u , UmeanAbs_sh_v , UmeanAbs_sh_w,rotorPowerrs,Azi,yawanglers,D)
           
%% Open directory for png to be saved

dirpath='/Volumes/NASSIR/MATLAB/Movie_slices_central_velocity_deficit';
if ~exist(dirpath,'dir') 
    mkdir(dirpath);
end

Uups=9;

%% Open figure and set viweing properties
fig2=figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
daspect([1 1 1])
view(3);
view(-25,15)
grid on
%axis tight
%axis([450/D 2000/D 300/D 700/D 0 300/D]) %defining x, y and z limits based on a rotor diameter sacle

%% plot turbine rotors
rotors=plotTurbine3D_yaw(Azi{n}(i,1),yawanglers(i,1)+90,500/D,500/D,115/D,D);
rotors=plotTurbine3D_yaw(Azi{n}(i,2),yawanglers(i,2)+90,500/D+5,500/D,115/D,D);

%% Define axis aspect and length
axis([-0.5+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))+1]);
ax = gca;

ax.XTick = 0:0.985:11;
ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
ax.YTickLabel = {'','',''};
ax.ZTick = 0:0.5:2; 
ax.ZTickLabel = {''};
    
 %% Make slices and plot directly in figure
VEL=sqrt(UmeanAbs_sh_u.^2+UmeanAbs_sh_v.^2+UmeanAbs_sh_w.^2)/Uups;
yslice=500/D;
hold on
%s3=slice(Xm_d,Ym_d,Zm_d,VEL,[],yslice,[]);
warning off
hold on;
VEL=sqrt(UmeanAbs_sh_u.^2+UmeanAbs_sh_v.^2+UmeanAbs_sh_w.^2)/Uups;
xslicetb1=[500/D  ;500/D+1;500/D+2;500/D+3;500/D+4];
xslicetb2=[500/D+5;500/D+6;500/D+7;500/D+8;500/D+9];
yslice=500/D;% [500]/D;
zslice=[115/D];

s1=slice(Xm_d,Ym_d,Zm_d,VEL,xslicetb1,[],[]);
s2=slice(Xm_d,Ym_d,Zm_d,VEL,xslicetb2,[],[]);

shading interp
colormap(jet(4096))
grid off

%% Write wind farm features for time sample i   
STring1=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,2)/1e6));
text(8.0,4,2.0,STring1,'Fontsize',16,'FontName','Tahoma')
STring2=sprintf('%2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6));
text(3.0,4,2.0,STring2,'Fontsize',16,'FontName','Tahoma')
STring3={sprintf('Total Mean Power: %2.1f MW',mean(rotorPowerrs{n}(i,1)/1e6)+mean(rotorPowerrs{n}(i,2)/1e6))};
text(10.0,3.55,2.25,STring3,'Fontsize',16,'FontName','Tahoma')
% STring4={sprintf('Yaw angle: %2.1f º',abs(270-yawanglers(i,1)))};
% text(12.0,3.55,2.3,STring4,'Fontsize',16,'FontName','Tahoma')

%% Time 
k=i+250;
     number=k/30;
     integ=floor(number);
     fract=number-integ;
     
     minutos=integ;
     segundos=60*fract;
     
 
    titlee=suptitle(['First turbine yaw angle of: ', ...
        num2str(round(abs(270-yawanglers(i,1)))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=18;
    
  %% Colorbar
  c = colorbar;
  c.Label.String = '|| (u,v,w) || / U_h ';
  caxis([0 1.25])
  set(gca,'fontsize', 14)
    
%% Export figure
export_fig(fig2,strcat(dirpath,'/image',num2str(20000+i)),'-nocrop','-m2')
warning on

