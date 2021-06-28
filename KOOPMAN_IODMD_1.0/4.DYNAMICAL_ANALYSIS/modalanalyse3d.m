function []=modalanalyse3d(xmodal,modestoanalyse,x,y,z,Decimate,U,D,modeltouse,modalmodels,dir,Inputs,scalingfactors)

    dirvideo='/Modes_animation';
    dirvideo=strcat(dir,dirvideo);
    if ~exist(dirvideo,'dir') 
        mkdir(dirvideo);
    end
   
    
Uups=9;

%% SET UP 3D BASICS
[xx,yy,zz,X,Y,Z]=retakepoints([],x,y,z,Decimate);
%[xx,yy,zz]=resamplegrid(x,y,z, Decimate);
[Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx,yy,zz);
% X = length(xx);
% Y = length(yy);
% Z = length(zz);

Xm_d = Xm_sh/D;
Ym_d = Ym_sh/D;
Zm_d = Zm_sh/D;

%% SEET UP FIGURE
for t=1:size(xmodal{1},2)
    fig980= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
    set(gcf,'color','w','Position', get(0, 'Screensize'));  

    for mm=1:length(modestoanalyse)
        subplot(1,length(modestoanalyse),mm)

        %rotors=plotTurbine3D_yaw(0,Inputs(t)*scalingfactors(3),500/D,500/D,115/D,D);
        rotors=plotTurbine3D_yaw(0,0,500/D,500/D,115/D,D);
        rotors=plotTurbine3D_yaw(0,0,500/D+5,500/D,115/D,D);
        daspect([1 1 1])
        view(3);
        view(285,30)
        grid on
        axis tight
        %axis([450/D 2000/D 300/D 700/D 0 300/D]) %defining x, y and z limits based on a rotor diameter sacle
        axis([-2+min(min(min(Xm_d))) max(max(max(Xm_d))) min(min(min(Ym_d))) max(max(max(Ym_d))) min(min(min(Zm_d))) max(max(max(Zm_d)))]);
        ax = gca;
        grid off
        ax.XTick = 0:0.985:11;
        ax.XTickLabel = {'','','','0D','1D','2D','3D','4D','5D','6D','7D','8D','9D','10D','11D'};
        ax.YTickLabel = {'','',''};
        ax.ZTick = 0:0.5:2; 
        ax.ZTickLabel = {''};

        % PLOT THE MODAL BEHAVIOUR 
        flow=U(:,1:modeltouse)*xmodal{modestoanalyse(mm)}(1:end,t);
        phiunique=reshape(flow, Y,X,Z);
        maxi=max(max(max(phiunique)));
        mini=min(min(min(phiunique)));
        dif=maxi-mini;

        %define want to see
        first=mini+dif/4;
        third=mini+dif/4*2;
        fourth=maxi-dif/3;
        fifth=maxi-dif/4;

        %identify all isosurfaces
        p1 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique,first));
        p3 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique,third));
        p4 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique,fourth));
        p5 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique,fifth));

         %color faces with jet scale by order
        p1.FaceColor = [0 0 0.9];
        p3.FaceColor = [0 0.9 0];
        p4.FaceColor = [0.9 0.4 0];
        p5.FaceColor = [0.9 0 0];

        %take out edges color
        p1.EdgeColor = 'none';
        p3.EdgeColor = 'none';
        p4.EdgeColor = 'none';
        p5.EdgeColor = 'none';

        %reduce mean flow isosurfaces strength and others
        p1.FaceAlpha=0.8;
        p3.FaceAlpha=0.1;
        p4.FaceAlpha=0.5;
        p5.FaceAlpha=0.8;
        %lighting properties
        camlight
        lighting gouraud
        
        %Make title to write properties
        mode=modestoanalyse(mm);
        [wn,zeta,p] = damp(modalmodels{modestoanalyse(mm)});
        title(['Mode ', num2str(mode), '  |   St: ',num2str(round(wn(1)/2/pi*D/Uups,3))...
        ,' fD/U  |   \xi: '  ,num2str(round(zeta(1),3)) ]);
        set(gca, 'FontSize', 14','FontName','Tahoma')
        set(gca, 'TitleFontWeight', 'normal')
        
    end
        % Calculations for time
        k=t+250;
        number=k/30;
        integ=floor(number);
        fract=number-integ;
     
        minutos=integ;
        segundos=60*fract;
        
        [ax4,h3]=suplabel(['Modal analyses. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds'],'t');
        set(h3,'FontSize',16)
        
        warning off
        export_fig(fig980,strcat(dirvideo,'/image',num2str(10000+t)),'-nocrop','-m2')
        warning on
        close all
end

