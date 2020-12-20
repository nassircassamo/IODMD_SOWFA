function [dirpathcuthubheightvec]=cuthubheightvec(D,dirpathcuthubheightvec,cases)

%% Directory 

if ~exist(dirpathcuthubheightvec,'dir') 
    mkdir(dirpathcuthubheightvec);
end

n=1;

[nTurbine,time4,dt,nVal,thrust]        = readTurbineOutputGlobal(cases{n},'rotorAxialForce');
[nTurbine,time4,dt,nVal,thrustv]       = readTurbineOutputGlobal(cases{n},'rotorVerticalForce');
[nTurbine,time4,dt,nVal,thrusth]       = readTurbineOutputGlobal(cases{n},'rotorHorizontalForce');
[nTurbine,time4,dt,nVal,azi]           = readTurbineOutputGlobal(cases{n},'rotorAzimuth');
[nTurbine,time5,dt,nVal,rotorPower{n}] = readTurbineOutputGlobal(cases{n},'rotorPower');
[nTurbine,time6,dt,nVal,yawangle]      = readTurbineOutputGlobal(cases{n},'nacelleYaw');

%% Resample data in order to have according dimensions
beg=750;

thrustrs(:,:)       =resample(thrust(end-beg*10:1:end,1:end),1,10);
thrustvrs(:,:)      =resample(thrustv(end-beg*10:1:end,1:end),1,10);
thrusthrs(:,:)      =resample(thrusth(end-beg*10:1:end,1:end),1,10);
azirs(:,:)          =resample(azi(end-beg*10:1:end,1:end),1,10);
rotorPowerrs{n}(:,:)=resample(rotorPower{n}(end-beg*10:1:end,1:end),1,10);
yawanglers(:,:)     =resample(yawangle(end-beg*10:1:end,1:end),1,10);
timeplot(1,:)       =resample(time6(1,end-beg*10:1:end)',1,10);

%%

load('U_data_complete_vec');

[xx,yy,zz]=resamplegrid(x,y,z, Decimate);
[Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx-500,(yy-500),zz);
X = length(xx);
Y = length(yy);
Z = length(zz);

%% simulation flow field

Uups=9; %[m/s]

[l,c]=size(QQ_u);

for i=10:1:(c-8)
    UmeanAbs_sh_u = reshape(double(QQ_u(:,i)),Y,X,Z);
    UmeanAbs_sh_v = reshape(double(QQ_v(:,i)),Y,X,Z);
    UmeanAbs_sh_w = reshape(double(QQ_w(:,i)),Y,X,Z);
    
    %curl_ circulation density of the fluid and provides mathematical
    %insights into fluid rotation based on fluid velocity field
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    figure1= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    
    %% Get velocity field at hubheight cut
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Vsecv=UmeanAbs_sh_v(:,:,k);
    Wsecv=UmeanAbs_sh_w(:,:,k);
    Usq=squeeze(Usecu);
    Vsq=squeeze(Vsecv);
    Wsq=squeeze(Wsecv);
     
    %% ABSOLUTE VELOCITY 
%     subplot(3,10,1:10)
%     absvel=sqrt(Usq.^2+Vsq.^2+Usq.^2+Wsq.^2)/Uups;
%     a=surf(Xm_shs,Ym_shs,absvel);
%     a.FaceAlpha=0.6;
%     shading interp
%     colormap(jet(4096))
%     view(2)
%     hold on
%     axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
%     ax=gca;
%     ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
%     ylabel('Distance [m]');
%     pbaspect([6 1 1])
%     title('Normalized absolute velocity at hub height');
%     c = colorbar('Box', 'off');
%     c.Label.String = 'Elevation (ft in 1000s)';
%     quiver(Xm_shs,Ym_shs,Usq, Vsq, 'color',[0 0 0]);

    %% STREAMWISE VELOCITY
    subplot(2,10,1:10)
    plotturbinefromabove(yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([4 1 1])
    title('Normalized streamwise velocity at hub height');
    c = colorbar;
    c.Label.String = 'u / U_h';
    set(gca,'fontsize', 14) 
    q=quiver(Xm_shs,Ym_shs,Usq, Vsq, 'color',[0 0 0]);
    q.Color=[0.2 0.2 0.2];
    %q.AutoScaleFactor=0.6;
    caxis([0 1.25])
    hold off 
    
    %% SPANWISE VELOCITY 
    subplot(2,10,11:20)
    plotturbinefromabove(yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Vsq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Vsq)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax=gca;
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([4 1 1])
    title('Normalized spanwise velocity at hub height');
    c = colorbar;
    c.Label.String = 'v / U_h';
    q=quiver(Xm_shs,Ym_shs,Usq, Vsq, 'color',[0 0 0]);
    %q.AutoScaleFactor=0.6;
    set(gca,'fontsize', 14)
    caxis([-0.6 0.5])
    hold off
    
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
    
    
 %%
    warning off
    export_fig(figure1,strcat(dirpathcuthubheightvec,'/image',num2str(10000+i)),'-nocrop','-m2')
    
    warning on
    
    close all
end

