function []=comparisonyaw(statesvalid,statesrecon,Inputs,Outputs,x,y,z,Decimate,dir,D,Outputs_model)

dircompyaw='/comparisonhubheight';
    dircompyaw=strcat(dir,dircompyaw);
    if ~exist(dircompyaw,'dir') 
        mkdir(dircompyaw);
    end

Inputs=Inputs+260;
[l,c]=size(Inputs);
Inputs(2,1:c)=270;

Inputs=Inputs';
Uups=9;
[xx,yy,zz,X,Y,Z]=retakepoints([],x,y,z,Decimate);
[Xm_sh,Ym_sh,Zm_sh] = meshgrid(xx-500,(yy-500),zz);

for t=1:size(statesvalid,2)
    
    fig980= figure('Units', 'pixels', 'pos', [75 75 1155 450],'color','white','Visible', 'off');
    set(gcf,'color','w','Position', get(0, 'Screensize'));  
    
    % DMD FLOW FIELD
    subplot(2,1,1)
    plotsnapshothh(statesvalid,xx,yy,abs(Inputs), D, t,X,Y,Z,Uups,Xm_sh,Ym_sh)
    set(gca,'fontsize', 20)
    number=t/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
    titlee=title(['SOWFA data: \gamma : -', ...
    num2str(round(abs(Inputs(t)-270))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds. \delta Power: ', num2str(round(Outputs(1,t)+Outputs(2,t),2)),' [MW]']);
    titlee.FontSize=20;
    %titlee.FontWeight='normal';
    
    
    % REAL FLOW FIELD
    subplot(2,1,2)
    plotsnapshothh(statesrecon,xx,yy,abs(Inputs), D, t,X,Y,Z,Uups,Xm_sh,Ym_sh)
    set(gca,'fontsize', 20)
    number=t/30;
    integ=floor(number);
    fract=number-integ;
    minutos=integ;
    segundos=60*fract;
     titlee=title(['DMD reconstruction: \gamma : -', ...
    num2str(round(abs(Inputs(t)-270))), ' degrees. Time: ',num2str(minutos)...
        ,' minutes and ',num2str(segundos),' seconds. \delta Power: ', num2str(round(Outputs_model(1,t)+Outputs_model(2,t),2)),' [MW]']);
    titlee.FontSize=20;
    %titlee.FontWeight='normal';
    
   
    warning off
    export_fig(fig980,strcat(dircompyaw,'/image',num2str(10000+t)),'-nocrop','-m2')
    warning on
    close all
     
end

