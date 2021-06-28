function []=plotmode(x,y,z,state,mode,Decimate,D,f,phi,P,LambdaDiag,damping)

    u=state;

    [xx,yy,zz]=resamplegrid(x,y,z, Decimate);
    [Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx-500,(yy-500),zz);
    X = length(xx);
    Y = length(yy);
    Z = length(zz);
    U_real = reshape(real(u(:,mode)),Y,X,Z);
    U_imag = reshape(imag(u(:,mode)),Y,X,Z);
    [Xm_shs,Ym_shs] = meshgrid(xx-500,(yy-500));
    k=9; 
    Usecu_real=U_real(:,:,k);
    Usecu_imag=U_imag(:,:,k);
    Usq_real=squeeze(Usecu_real);
    Usq_imag=squeeze(Usecu_imag);
    [p1]=plotturbinefromabove(0, 0, 0, D);
    hold on
    [p2]= plotturbinefromabove(0, D*5, 0, D);
    a=pcolor(Xm_shs,Ym_shs,Usq_real);
    a.FaceAlpha=0.8;
    shading interp
    colormap(jet(4096))
    hold on
    %b=pcolor(Xm_shs,Ym_shs,Usq_imag);
    %b.FaceAlpha=0.8;
    set(a,'ZData',100+zeros(size(Usq_real)))
    shading interp
    colormap(jet(4096))
    view(2)
    axis([ min(min(min(Xm_sh))) max(max(max(Xm_sh))) min(min(min(Ym_sh))) max(max(max(Ym_sh)))]);
    ax=gca;
    set(gca,'XTick',(0:178:max(max(max(Xm_sh)))))
    ax.XTickLabel = {'0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
    ylabel('Distance [m]');
    pbaspect([5 1 1])
    titlee=title(['Mode ', num2str(mode), ' at hub height   |   St: ',num2str(round(f(mode),4))...
        ,' fD/U  |   \xi: '  ,num2str(round(damping(mode),3)) '   |   || \phi ||:  ',num2str(round(P(mode),3))]);
    titlee.FontSize=14;
    titlee.FontWeight='normal';
    c = colorbar;
    c.Label.String = 'u [m/s] ';
    set(gca,'fontsize', 12) 
    %caxis([0 1.25])
    hold off 
    set(gca,'Ydir','reverse')