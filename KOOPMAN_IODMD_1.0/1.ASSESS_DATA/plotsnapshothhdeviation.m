function []=plotsnapshothhdeviation(states,xx,yy,yawanglers, D, i,X,Y,Z,Uups,Xm_sh,Ym_sh)

    UmeanAbs_sh_u = reshape(double(real(states(:,i))),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu=UmeanAbs_sh_u(:,:,k);
    Usq=squeeze(Usecu);
    plotturbinefromabove(yawanglers(i,1)+90, 0, 0, D);
    hold on
    plotturbinefromabove(yawanglers(i,2)+90, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq/Uups);
    a.FaceAlpha=0.8;
    set(a,'ZData',-1+zeros(size(Usq/Uups)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    k=i+250;
    number=k/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['\gamma : ', ...
    num2str(round(abs(270-yawanglers(i,1)))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds']);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = ' (u_{SOWFA} - u_{DMD} ) / U ';
    set(gca,'fontsize', 12) 
    caxis([-0.5 0.5])
    hold off 