function []=modeanimation(modes,x,y,z,P,phi,freq,damping,LambdaDiag,b,Decimate,D,dirdmd,f,decision,scalingfactor,meansteadystate,Xd,X)

if decision==0
    
elseif decision==1

%% Grid specifications
    [xx,yy,zz]=resamplegrid(x,y,z, Decimate);
    [Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx,yy,zz);
    XX = length(xx);
    YY = length(yy);
    ZZ = length(zz);

    Xm_d = Xm_sh/D;
    Ym_d = Ym_sh/D;
    Zm_d = Zm_sh/D;

    %% Order by energy to get POD modes
    if isempty(Xd)
        [sortedP,I]=sort(P,'descend');
        sortedphi=phi(:,I);
        sortedf=freq(I);
        sorteddamping=damping(I);
        sortedLd=LambdaDiag(I);
        sortedb=b(I);
    else
        fromhere=size(Xd,1)+1;
        phit=phi(fromhere:end,fromhere:end);
        ft=freq(fromhere:end);
        Pt=P(fromhere:end);
        dampingt=damping(fromhere:end);
        LambdaDiagt=LambdaDiag(fromhere:end);
        bt=b(fromhere:end);
        
        [sortedP,I]=sort(Pt,'descend');
        sortedphi=phit(:,I);
        sortedf=ft(I);
        sorteddamping=dampingt(I);
        sortedLambdaDiag=LambdaDiagt(I);
        b=bt(I);
    end

    %% Rename 
    P=sortedP;
    phi=sortedphi;
    freq=sortedf;
    damping=sorteddamping;
    LambdaDiag=sortedLambdaDiag;

    for mode=1:length(modes)

        %% Rebuild state of interest
        mm1=size(X,2);
        dt=2;
        omega=imag(log(LambdaDiag))/dt/(2*pi);%Hz - continuous time oscill fre
        %omega=omega*2*pi; %rad/s oscillatory frequency 
        t=(0:mm1-1)*dt;  

        for iter=1:length(t)
            time_dynamics(:,iter)=b(modes(mode)).*sin(omega(modes(mode))*2*pi*t(iter));
        end

        %statesrebuild=meansteadystate+phi(:,modes(mode))*time_dynamics*scalingfactor;
        %statesrebuild=meansteadystate;
        statesrebuild=phi(:,modes(mode))*time_dynamics*scalingfactor;
        modalanalysis{mode}=statesrebuild;
    end

    %% Calculate effect on fluid flow
    Uups=9;
    color=jet(5);
    for i=1:length(modes)
        count=0;
        dir=strcat(dirdmd,'/modesanimation/');
        dir=strcat(dir,num2str(modes(i)));
        directo{i}=dir;
        
        if ~exist(dir,'dir') 
            mkdir(dir);
        end
        
        for t=2:150
            count=count+1;
            
            phiunique{count}=reshape(real(modalanalysis{i}(:,t)), YY,XX,ZZ);

            fig580= figure('Units', 'pixels', 'pos', [75 75 1155 650],'color','white','Visible', 'off');
            set(gcf,'color','w','Position', get(0, 'Screensize'));    
            rotors=plotTurbine3D_yaw(0,0,500/D,500/D,115/D,D);
            rotors=plotTurbine3D_yaw(0,0,500/D+5,500/D,115/D,D);
            daspect([1 1 1])
            view(3);
            view(-125,15)
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


            maxi=max(max(max(phiunique{1})));
            mini=min(min(min(phiunique{1})));
            dif=maxi-mini;
       
            %define want to see
            first=mini+dif/4;
            third=mini+dif/4*2;
            fourth=maxi-dif/3;
            fifth=maxi-dif/4;

            %identify all isosurfaces
            p1 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique{count},first));
            p3 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique{count},third));
            p4 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique{count},fourth));
            p5 = patch(isosurface(Xm_d,Ym_d,Zm_d,phiunique{count},fifth));

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
            mode=modes(i);
            %insert mode chracterisitcs in title
            title(['POD mode ', num2str(mode), ' evolution through time   |   St: ',num2str(round(freq(mode),4))...
            ,' fD/U  |   \xi: '  ,num2str(round(damping(mode),3)) '   |   || \phi ||:  ',num2str(round(P(mode),3))]);
            set(gca, 'FontSize', 16','FontName','Tahoma')
            set(gca, 'TitleFontWeight', 'normal')
            legend([p1 p3 p4 p5],['Isosurface of velocity ',num2str(first) ' ms^{-1} '],...
                ['Isosurface of velocity ',num2str(third) ' ms^{-1} '],...
                 ['Isosurface of velocity ',num2str(fourth) ' ms^{-1} '],...
            ['Isosurface of velocity ',num2str(fifth) ' ms^{-1} '],'Location','bestoutside','Orientation','horizontal')
            legend('boxoff')
            export_fig(fig580,strcat(dir,'/image',num2str(t)),'-nocrop','-m2'); 
            close all
        end
    end
    
%     for a=1:length(directo)
%         
%         makemoviewake(strcat('mode',num2str(modes(a))),directo{a});
%         
%     end
    
    
end

    