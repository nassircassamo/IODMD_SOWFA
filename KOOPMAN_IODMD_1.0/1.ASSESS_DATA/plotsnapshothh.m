function []=plotsnapshothh(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh,pitchmode)

    UmeanAbs_sh_u = reshape(double(real(states(:,i))),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Usq=squeeze(Usecu);
    
    if pitchmode==0
        plotturbinefromabove(-yawanglers(i,1)+90, 0, 0, D);%yaw
        %plotturbinefromabove(270+90, 0, 0, D);%pitch
        hold on
        plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D); %yaw
    elseif pitchmode==1
        plotturbinefromabove(0, 0, 0, D);%no need to represent yaw misalignment
        %plotturbinefromabove(270+90, 0, 0, D);%pitch
        hold on
        plotturbinefromabove(0, D*5, 0, D); %no need to represent yaw misalignment
    end
        
    %plotturbinefromabove(270+90, D*5, 0, D);%pitch
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
    pbaspect([5 1 1])
    k=i;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    
    
    if pitchmode==0
        titlee=title(['|\delta\gamma| : ', ...
        num2str(round(abs(270-yawanglers(i,1)))), ' degrees. Time: ',num2str(minutos)...
            ,' minutes and ',num2str(segundos),' seconds']);
    elseif pitchmode==1
        titlee=title(['\beta: ', ...
        num2str(round(abs(270-yawanglers(i,1)))), ' degrees. Time: ',num2str(minutos)...
           ,' minutes and ',num2str(segundos),' seconds']); %PITCH
    end
   
    
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u / U_\infty';
    set(gca,'fontsize', 18) 
    caxis([0 1.25])
    set(gca,'Ydir','reverse')
    set(gca,'fontname','times')  % Set it to times
    hold off 
    