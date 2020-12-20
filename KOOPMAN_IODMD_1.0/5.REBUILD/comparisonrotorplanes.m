function []=comparisonrotorplanes(sys_red, modeltouse, Inputs, Outputs, x, y, z, Decimate, D, QQu, flowdmd,dir)

dircomp='/comparisonrotor';
    dircomp=strcat(dir,dircomp);
    if ~exist(dircomp,'dir') 
        mkdir(dircomp);
    end
    

Uups=9;

%% SET UP 3D BASICS
 Uups=9; %[m/s]
    [xx,yy,zz,X,Y,Z]=retakepoints([],x,y,z,Decimate);
   % [xx,yy,zz]=resamplegrid(x,y,z, Decimate);
    [Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx-500,(yy-500),zz);
   % X = length(xx);
   % Y = length(yy);
   % Z = length(zz);

Xm_d = Xm_sh/D;
Ym_d = Ym_sh/D;
Zm_d = Zm_sh/D;

%% FOR EACH TIME SAMPLE

for t=1:size(QQu,2)

    fig980= figure('Units', 'pixels', 'pos', [75 75 1155 450],'color','white','Visible', 'off');
  %  set(gcf,'color','w','Position', get(0, 'Screensize'));  

     %% TRUE FLOW FIELD
     subplot(1,3,1)
     UmeanAbs_sh_u = reshape(double(QQu(:,t)),Y,X,Z);
     k=72;
     Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
     [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
     Usq=squeeze(Usecu);
%      [M,cc]=contour(Ym_shs',Zm_shs',Usq);
%      hold on
%      cc.LineWidth=1.5;
     a=pcolor(Ym_shs',Zm_shs',Usq);
     set(a,'ZData',-1+zeros(size(Usq)));
     a.FaceAlpha=0.8;
     shading interp
     colormap(jet(4096))
     view(2)
     axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
     xlabel('Distance [m]');
     ylabel('Distance [m]');
     pbaspect([1 1 1])
     hold on
     [p]=plotrotor2D(0,D,0,115);
     daspect([ 1 1 1])
     caxis([0 1.2])
     set(gca,'fontsize', 16)
     c=colorbar('southoutside');
     c.Label.String = 'u_{SOWFA}  / U_\infty ';
     title('SOWFA flow field')
     
     
     %% DMD FLOW FIELD APPROXIMATION
     subplot(1,3,2)
     flow=flowdmd(:,t);
     UmeanAbs_sh_u = reshape(flow,Y,X,Z);
     k=72;
     Usecu=UmeanAbs_sh_u(:,k,:)./Uups;
    [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
    Usq=squeeze(Usecu);
%    [M,cc]=contour(Ym_shs',Zm_shs',Usq);
%     hold on
%     cc.LineWidth=1.5;
     a=pcolor(Ym_shs',Zm_shs',Usq);
     set(a,'ZData',-1+zeros(size(Usq)));
     a.FaceAlpha=0.8;
     shading interp
     colormap(jet(4096))
     view(2)
     axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
     xlabel('Distance [m]');
     ylabel('Distance [m]');
     hold on
     pbaspect([1 1 1])
     hold on
     [p]=plotrotor2D(0,D,0,115);
     daspect([ 1 1 1])
     caxis([0 1.2])
     set(gca,'fontsize', 16)
     c=colorbar('southoutside');
     c.Label.String = 'u_{DMD}  / U_\infty ';
     title('DMD flow field reconstruction')
     
     
     %% PERCENTUAL DEVIATION BETWEEN TRUE FLOW FIELD AND DMD RECONSTRUCTION 
     subplot(1,3,3)
     error=abs( abs((QQu(:,t)-flowdmd(:,t))) ./ abs((QQu(:,t))))*100;
     error3d=reshape(error,Y,X,Z);
     [Ym_shs,Zm_shs] = meshgrid((yy-500),zz);
     err=error3d(:,k,:);
     Usq=squeeze(err);
    % [M,cc]=contour(Ym_shs',Zm_shs',Usq);
    % hold on
    % cc.LineWidth=1.5;
     a=pcolor(Ym_shs',Zm_shs',Usq);
     set(a,'ZData',-1+zeros(size(Usq)));
     a.FaceAlpha=0.9;
     shading interp
     colormap(jet(4096))
     view(2)
     axis([ min(min(min(Ym_sh))) max(max(max(Ym_sh))) min(min(min(Zm_sh))) max(max(max(Zm_sh)))]);
     xlabel('Distance [m]');
     ylabel('Distance [m]');
     hold on
     pbaspect([1 1 1])
     hold on
     [p]=plotrotor2D(0,D,0,115);
     daspect([ 1 1 1])
     caxis([0 100])
     c=colorbar('southoutside');
     c.Label.String = ' [%]';
     set(gca,'fontsize', 16)
      title('DMD relative deviation from SOWFA ')
     %% Major title
     k=t;
     number=k/30;
     integ=floor(number);
     fract=number-integ;
     
     minutos=integ;
     segundos=60*fract;
     
    [ax4,h3]=suplabel(['Flow field reconstruction comparison. \gamma: ', ...
        num2str(round(Inputs(t))-10), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds'] ,'t');
    set(h3,'FontSize',20)
    set(h3,'FontWeight','normal')
    
    %% Exporting
     warning off
     export_fig(fig980,strcat(dircomp,'/image',num2str(10000+t)),'-nocrop','-m2')
     warning on
     close all
   
end
     
     
     