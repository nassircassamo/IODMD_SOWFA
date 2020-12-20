function []=podmodes3dim(x,y,z,Xd, P,phi,freq, damping, Decimate,D,dirdmd)

dir3dmodes='/3Dmodes';
dirdmd=strcat(dirdmd,dir3dmodes);

if ~exist(dirdmd,'dir') 
        mkdir(dirdmd);
end

%% GRID SPECIFICATIONS
%[xx,yy,zz]=resamplegrid(x,y,z, Decimate);
[xx,yy,zz,X,Y,Z]=retakepoints([],x,y,z,Decimate);
[Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx,yy,zz);
%X = length(xx);
%Y = length(yy);
%Z = length(zz);

Xm_d = Xm_sh/D;
Ym_d = Ym_sh/D;
Zm_d = Zm_sh/D;

%% ORDER MODES BY ENERGY (POD EQUIVALENT)
if isempty(Xd)
     [sortedP,I]=sort(P,'descend');
     sortedphi=phi(:,I);
     sortedf=freq(I);
     sorteddamping=damping(I);
else
%     fromhere=size(Xd,1)+1;
%     phit=phi(fromhere:end,fromhere:end);
%     ft=freq(fromhere:end);
%     Pt=P(fromhere:end);
%     dampingt=damping(fromhere:end);
    [sortedP,I]=sort(P,'descend');
    sortedphi=phi(:,I);
    sortedf=f(I);
    sorteddamping=damping(I);
end

%% RENAMING FOR PLOTTING

 P=sortedP;
 phi=sortedphi;
 freq=sortedf;
 damping=sorteddamping;

%% MODE SELECTION AND RESHPAING
phistates=phi;
close all

    for i=1:size(phistates,2)
       % part=4; subpart=3; [f]= MPC_progress(part,subpart,f,i,size(phistates,2));
        phiunique{i}=reshape(real(phistates(:,i)), Y,X,Z);

        fig580= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
        set(gcf,'color','w','Position', get(0, 'Screensize'));    
        rotors=plotTurbine3D_yaw(0,0,500/D,500/D,115/D,D);
        rotors=plotTurbine3D_yaw(0,0,500/D+5,500/D,115/D,D);
        daspect([1 1 1])
        view(3);
        view(120+180,20)
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
        
        maxi=max(max(max(phiunique{i})));
        mini=min(min(min(phiunique{i})));
        dif=maxi-mini;

        %define want to see
        first=mini+dif/4;
        third=mini+dif/4*2;
        fourth=maxi-dif/3;
        fifth=maxi-dif/4;

        %identify all isosurfaces
        p1 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique{i},first));
        p3 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique{i},third));
        p4 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique{i},fourth));
        p5 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique{i},fifth));

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
        mode=i;
        %insert mode chracterisitcs in title
        title(['POD mode ', num2str(mode), ' at hub height   |   St: ',num2str(round(freq(mode)*178/9,4))...
        ,' fD/U  |   \xi: '  ,num2str(round(damping(mode),3)) '   |   || \phi ||:  ',num2str(round(P(mode),3))]);
        set(gca, 'FontSize', 14','FontName','Tahoma')
        set(gca, 'TitleFontWeight', 'normal')
        legend([p1 p3 p4 p5],['',num2str(first)],...
            ['',num2str(third)],...
             [' ',num2str(fourth)],...
        ['',num2str(fifth)],'Location','bestoutside','Orientation','vertical')
        legend('boxoff')
        
        export_fig(fig580,strcat(dirdmd,'/mode',num2str(i)),'-nocrop','-m2'); 
        close all
    end
    
%     
%     
%       title(['St: ',num2str(round(freq(mode),4))...
%         ,' fD/U  |   \xi: '  ,num2str(round(damping(mode),3)) ]);
%         
%      

     
end



%% FIGURE SPECIFICATIONS
