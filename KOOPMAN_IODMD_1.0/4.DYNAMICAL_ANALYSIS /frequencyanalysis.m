function []=frequencyanalysis(model,modalmodels,Diam,Uups,maindir)
filename='fre_anal';

  dirfreq='/Modal_analysis';
     dirfreq=strcat(maindir,dirfreq);
    if ~exist(dirfreq,'dir') 
        mkdir(dirfreq);
    end
   
    
%% Change Bode Plot 
bo=bodeoptions;
bo.FreqUnits='Hz';

%% Make file
nomedoficheiro='/freq_analysis';
extensao='.txt';
nomedoficheiro=strcat(nomedoficheiro, extensao);
permissao= 'wt'; 
fid=fopen([maindir nomedoficheiro], permissao); 

   if fid == -1
        disp('Error openning file'); 
   else
   fprintf(fid,['Modal mode number  |  Pole location  |  Frequency [Hz]  |  '...
       'Frequency  [St: fD/U]  |   '...
       'Damping  [ ]   |   ']);
      
    for mm= 1:length(modalmodels)
        [wn,zeta,p] = damp(modalmodels{mm});
        
        fprintf(fid,'\n %10.0f %15.6f %18.6f  %20.6f %20.5f %20.5f %20.4f %20.4f', char(mm),...
            round(p(1),4), wn(1)/2/pi, wn(1)/2/pi/Uups*Diam, round(zeta(1),4));
    end
          
    result = fclose(fid); %fecho do ficheiro
    if result ~=0
        frintf('Erroor closing file');  
    end
        fprintf('\n The File %14s has been created and saved with success \n', nomedoficheiro);
   end

%% Make Bode plot for overall system
fig900=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');
shg

subplot(2,1,1)
[mag,phase,wout] = bode(model,bo);
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=1.4;
p.Color=[0.1 0.1 0.85];
grid on
ylabel('Magnitude [dB]')
xlabel('Frequency [Hz]')
t=title('Frequency response of first turbine reduced order model ');
set(gca,'fontsize', 12)
t.FontSize=12;
xlim([10^-3 0.5])
hold on
p=plot([0.25 0.25],[-20 -5]);
p.Color=[0 0 0];
p.LineWidth=1.5;

subplot(2,1,2)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=1.4;
p.Color=[0.1 0.1 0.85];
grid on
ylabel('Magnitude [dB]')
xlabel('Frequency [Hz]')
set(gca,'fontsize', 12)
t=title('Frequency response of second turbine reduced order model ');
t.FontSize=12;
xlim([10^-3 0.5])
hold on
p=plot([0.25 0.25],[-70 0]);
p.Color=[0 0 0];
p.LineWidth=1.5;

export_fig(fig900,strcat(dirfreq,'/model',strcat(filename,'0')),'-nocrop','-m2'); 
%% Present Bode plot of different modes 1
fig901=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');

%First Modal Model
mm=1;
subplot(5,2,1)
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,2)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,3)
mm=2;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,4)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,5)
mm=3;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,6)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])


subplot(5,2,7)
mm=4;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,8)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,9)
mm=5;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,10)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

export_fig(fig901,strcat(dirfreq,'/modal',strcat(filename,'1')),'-nocrop','-m2'); 

%% Present Bode plot of different modes 2
fig902=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');

%First Modal Model
mm=6;
subplot(5,2,1)
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,2)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,3)
mm=7;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,4)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,5)
mm=8;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,6)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,7)
mm=9;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,8)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,9)
mm=10;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,10)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

export_fig(fig902,strcat(dirfreq,'/modal',strcat(filename,'2')),'-nocrop','-m2'); 
%% Present Bode plot of different modes 3
fig903=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');

%First Modal Model
mm=11;
subplot(5,2,1)
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,2)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,3)
mm=12;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,4)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,5)
mm=13;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,6)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,7)
mm=14;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,8)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,9)
mm=15;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,10)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

export_fig(fig903,strcat(dirfreq,'/modal',strcat(filename,'3')),'-nocrop','-m2'); 

%% Present Bode plot of different modes 4
fig904=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');

%First Modal Model
mm=16;
subplot(5,2,1)
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,2)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])


subplot(5,2,3)
mm=17;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,4)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,5)
mm=18;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
subplot(5,2,6)
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])


subplot(5,2,7)
mm=19;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,8)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,9)
mm=20;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,10)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

export_fig(fig902,strcat(dirfreq,'/modal',strcat(filename,'4')),'-nocrop','-m2'); 
%% Present Bode plot of different modes 5
fig905=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');


%First Modal Model
mm=21;
subplot(5,2,1)
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,2)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])


subplot(5,2,3)
mm=22;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,4)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,5)
mm=23;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,6)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])


subplot(5,2,7)
mm=24;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,8)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,9)
mm=25;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,10)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

export_fig(fig905,strcat(dirfreq,'/modal',strcat(filename,'5')),'-nocrop','-m2'); 

%% Present Bode plot of different modes 6
fig906=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');

%First Modal Model
mm=26;
subplot(5,2,1)
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,2)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,3)
mm=27;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,4)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,5)
mm=28;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,6)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,7)
mm=29;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,8)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,9)
mm=30;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,10)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

export_fig(fig906,strcat(dirfreq,'/modal',strcat(filename,'6')),'-nocrop','-m2'); 

%% Present Bode plot of different modes 6
fig907=figure;
set(gcf,'color','w','Position', get(0, 'Screensize'),'Visible','off');

%First Modal Model
mm=31;
subplot(5,2,1)
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,2)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])


subplot(5,2,3)
mm=32;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,4)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,5)
mm=33;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,6)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,7)
mm=34;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,8)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

subplot(5,2,9)
mm=35;
[mag,phase,wout] = bode(modalmodels{mm});
p=semilogx(wout/2/pi,20*log10(squeeze(mag(1,:,:))));
p.LineWidth=2;
p.Color=[0.8 0.4 0.4];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
[wn,zeta,p] = damp(modalmodels{mm});
t=title(['                                                                                                                                                                                                                           Modal model ', num2str(mm), ' | Frequency: ',num2str(round(wn(1)/2/pi,3)) , ' [Hz] | Frequency: ',...
    num2str(round(wn(1)/2/pi*Diam/Uups,3)) ,' [St] | \xi : ' , num2str(zeta(1)) ,...
    ' | Pole Location : ', num2str(p(1)),'.']);
t.FontSize=12;
t.FontWeight='normal';
xlim([wout(1)/2/pi wout(length(wout))/2/pi])
subplot(5,2,10)
p=semilogx(wout/2/pi,20*log10(squeeze(mag(2,:,:))));
p.LineWidth=2;
p.Color=[0.4 0.4 0.8];
grid on
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
xlim([wout(1)/2/pi wout(length(wout))/2/pi])

export_fig(fig907,strcat(dirfreq,'/modal',strcat(filename,'7')),'-nocrop','-m2');






