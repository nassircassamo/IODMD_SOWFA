% Example 9: PBSIDopt identification of Smart Rotor with periodic effects
close all; clc; clear all; 

%% With Rotor speed of 370 RPM

% Data files
load smartrotor370.mat u y
RPM = 370;
Ts = 0.01;
u = detrend(u);
y = detrend(y);
N = size(y,1);
tt = 0:Ts:Ts*(N-1);
dp = [sin(RPM/60*2*pi*tt); cos(RPM/60*2*pi*tt);...
    sin(2*RPM/60*2*pi*tt); cos(2*RPM/60*2*pi*tt)]';
ud = [u dp];
w = logspace(0.4,log10(1/(2*Ts)),2000)*2*pi;
[us,Du,ys,Dy] = sigscale(u,y); us=us';ys=ys';
[usd,Dud,ysd,Dyd] = sigscale(ud,y); usd=usd';ysd=ysd';

% spectral analysis
[Ga,ws] = spaavf(u,y,Ts,10,[],[],'hamming');
[Gad,wsd] = spaavf(ud,y,Ts,10,[],[],'hamming');
Ga = frd(Ga,ws); 
Gad = frd(Gad,wsd); 

% linear subspace estimate (PBSIDopt)
n = 12;
p = 20;
f = 20;

[S,X] = dordvarx(us,ys,f,p,'tikh','gcv');
figure(2), semilogy(S,'c*')
x = dmodx(X,n);
[Ai,Bi,Ci,Di] = dx2abcdk(x,us,ys,f,p);
Ci=Dy*Ci; Di=Dy*Di/Du; Bi=Bi/Du;
sys_idp=ss(Ai,Bi,Ci,Di,Ts);

[Sd,Xd] = dordvarx(usd,ysd,f,p,'tikh','gcv');
figure(2), hold on, semilogy(Sd,'c+'), hold off;
x = dmodx(Xd,n);
[Aid,Bid,Cid,Did] = dx2abcdk(x,usd,ysd,f,p); 
Cid=Dyd*Cid; Did=Dyd*Did/Dud; Bid=Bid/Dud;
sys_idpd=ss(Aid,Bid,Cid,Did,Ts);

% plotting
figure(100)
[MAG1, PHASE1,ww1]=bode(sys_idp(1:2,1:2),w);
[MAG2, PHASE2,ww2]=bode(sys_idpd(1:2,1:2),w);
[MAGd1, PHASEd1,wwd1]=bode(Ga(1:2,1:2),w);
[MAGd2, PHASEd2,wwd2]=bode(Gad(1:2,1:2),w);

h1=subplot(2,2,1)
loglog(ww1/2/pi,squeeze(MAG1(1,1,:)),'color',[0.8,0.8,0.8]); hold on
loglog(ww2/2/pi,squeeze(MAG2(1,1,:)),'color',[0,0,0]); hold on
loglog(ww1/2/pi,squeeze(MAGd1(1,1,:)),'color',[0.8,0.8,0.8],'linestyle','-.'); hold on
loglog(ww2/2/pi,squeeze(MAGd2(1,1,:)),'color',[0,0,0],'linestyle','-.'); hold on
axis([2.5 50 1e-4 0.03 ])
set(h1,'YTick',[1e-4  1e-3 1e-2 ])
set(h1,'XTick',[ 2 3 4 5 6 7 8 9 10 20 30 40 50])
set(h1,'YTickLabel',[ {'{-4}'}  {'{-3}'} {'{-2}'} {'{2}'}])
set(h1,'XTickLabel',[ {''} {''} 4 {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}])
ylabel('MFC blade 1 [V]')
title('From Flap blade 1 [V]')

h1=subplot(2,2,2)
loglog(ww1/2/pi,squeeze(MAG1(1,2,:)),'color',[0.8,0.8,0.8]); hold on
loglog(ww2/2/pi,squeeze(MAG2(1,2,:)),'color',[0,0,0]); hold on
loglog(ww1/2/pi,squeeze(MAGd1(1,2,:)),'color',[0.8,0.8,0.8],'linestyle','-.'); hold on
loglog(ww2/2/pi,squeeze(MAGd2(1,2,:)),'color',[0,0,0],'linestyle','-.'); hold on
axis([2.5 50 1e-4 0.03 ])
set(h1,'YTick',[1e-4  1e-3 1e-2 ])
set(h1,'XTick',[ 2 3 4 5 6 7 8 9 10 20 30 40 50])
set(h1,'YTickLabel',[ {'{-4}'}  {'{-3}'} {'{-2}'} {'{2}'}])
set(h1,'XTickLabel',[ {''} {''} 4 {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}])
title('From Flap blade 2 [V]')

h1=subplot(2,2,3)
loglog(ww1/2/pi,squeeze(MAG1(2,1,:)),'color',[0.8,0.8,0.8]); hold on
loglog(ww2/2/pi,squeeze(MAG2(2,1,:)),'color',[0,0,0]); hold on
loglog(ww1/2/pi,squeeze(MAGd1(2,1,:)),'color',[0.8,0.8,0.8],'linestyle','-.'); hold on
loglog(ww2/2/pi,squeeze(MAGd2(2,1,:)),'color',[0,0,0],'linestyle','-.'); hold on
axis([2.5 50 1e-4 0.03 ])
set(h1,'YTick',[1e-4  1e-3 1e-2 ])
set(h1,'XTick',[ 2 3 4 5 6 7 8 9 10 20 30 40 50])
set(h1,'YTickLabel',[ {'{-4}'}  {'{-3}'} {'{-2}'} {'{2}'}])
set(h1,'XTickLabel',[ {''} {''} 4 {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}])
ylabel('MFC blade 2 [V]')
xlabel('Frequency [Hz]')

h1=subplot(2,2,4)
loglog(ww1/2/pi,squeeze(MAG1(2,2,:)),'color',[0.8,0.8,0.8]); hold on
loglog(ww2/2/pi,squeeze(MAG2(2,2,:)),'color',[0,0,0]); hold on
loglog(ww1/2/pi,squeeze(MAGd1(2,2,:)),'color',[0.8,0.8,0.8],'linestyle','-.'); hold on
loglog(ww2/2/pi,squeeze(MAGd2(2,2,:)),'color',[0,0,0],'linestyle','-.'); hold on
axis([2.5 50 1e-4 0.03 ])
set(h1,'YTick',[1e-4  1e-3 1e-2 ])
set(h1,'XTick',[ 2 3 4 5 6 7 8 9 10 20 30 40 50])
set(h1,'YTickLabel',[ {'{-4}'}  {'{-3}'} {'{-2}'} {'{2}'}])
set(h1,'XTickLabel',[ {''} {''} 4 {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}])
xlabel('Frequency [Hz]')

%% With Rotor speed of 430 RPM

% Data files
load smartrotor430.mat
RPM = 430;
Ts = 0.01;
u = detrend(u);
y = detrend(y);
N = size(y,1);
tt = 0:Ts:Ts*(N-1);
dp = [sin(RPM/60*2*pi*tt); cos(RPM/60*2*pi*tt);...
    sin(2*RPM/60*2*pi*tt); cos(2*RPM/60*2*pi*tt)]';
ud = [u dp];
w=logspace(0.4,log10(1/(2*Ts)),2000)*2*pi;
[us,Du,ys,Dy] = sigscale(u,y); us=us';ys=ys';
[usd,Dud,ysd,Dyd] = sigscale(ud,y); usd=usd';ysd=ysd';

%% spectral analysis
[Ga,ws] = spaavf(u,y,Ts,10,[],[],'hamming');
[Gad,wsd] = spaavf(ud,y,Ts,10,[],[],'hamming');
Ga = frd(Ga,ws); 
Gad = frd(Gad,wsd); 

%% linear subspace estimate (PBSIDopt)
n = 12;
p = 20;
f = 20;

[S,X] = dordvarx(us,ys,f,p,'tikh','gcv');
figure(2), semilogy(S,'c*')
x = dmodx(X,n);
[Ai,Bi,Ci,Di] = dx2abcdk(x,us,ys,f,p);
Ci=Dy*Ci; Di=Dy*Di/Du; Bi=Bi/Du;
sys_idp=ss(Ai,Bi,Ci,Di,Ts);

[Sd,Xd] = dordvarx(usd,ysd,f,p,'tikh','gcv');
figure(2), hold on, semilogy(Sd,'c+'), hold off;
x = dmodx(Xd,n);
[Aid,Bid,Cid,Did] = dx2abcdk(x,usd,ysd,f,p); 
Cid=Dyd*Cid; Did=Dyd*Did/Dud; Bid=Bid/Dud;
sys_idpd=ss(Aid,Bid,Cid,Did,Ts);

% plotting
figure(100)
[MAG1, PHASE1,ww1]=bode(sys_idp(1:2,1:2),w);
[MAG2, PHASE2,ww2]=bode(sys_idpd(1:2,1:2),w);
[MAGd1, PHASEd1,wwd1]=bode(Ga(1:2,1:2),w);
[MAGd2, PHASEd2,wwd2]=bode(Gad(1:2,1:2),w);

h1=subplot(2,2,1)
loglog(ww1/2/pi,squeeze(MAG1(1,1,:)),'color',[0.8,0.8,0.8]); hold on
loglog(ww2/2/pi,squeeze(MAG2(1,1,:)),'color',[0,0,0]); hold on
loglog(ww1/2/pi,squeeze(MAGd1(1,1,:)),'color',[0.8,0.8,0.8],'linestyle','-.'); hold on
loglog(ww2/2/pi,squeeze(MAGd2(1,1,:)),'color',[0,0,0],'linestyle','-.'); hold on
axis([2.5 50 1e-4 0.03 ])
set(h1,'YTick',[1e-4  1e-3 1e-2 ])
set(h1,'XTick',[ 2 3 4 5 6 7 8 9 10 20 30 40 50])
set(h1,'YTickLabel',[ {'{-4}'}  {'{-3}'} {'{-2}'} {'{2}'}])
set(h1,'XTickLabel',[ {''} {''} 4 {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}])
ylabel('MFC blade 1 [V]')
title('From Flap blade 1 [V]')

h1=subplot(2,2,2)
loglog(ww1/2/pi,squeeze(MAG1(1,2,:)),'color',[0.8,0.8,0.8]); hold on
loglog(ww2/2/pi,squeeze(MAG2(1,2,:)),'color',[0,0,0]); hold on
loglog(ww1/2/pi,squeeze(MAGd1(1,2,:)),'color',[0.8,0.8,0.8],'linestyle','-.'); hold on
loglog(ww2/2/pi,squeeze(MAGd2(1,2,:)),'color',[0,0,0],'linestyle','-.'); hold on
axis([2.5 50 1e-4 0.03 ])
set(h1,'YTick',[1e-4  1e-3 1e-2 ])
set(h1,'XTick',[ 2 3 4 5 6 7 8 9 10 20 30 40 50])
set(h1,'YTickLabel',[ {'{-4}'}  {'{-3}'} {'{-2}'} {'{2}'}])
set(h1,'XTickLabel',[ {''} {''} 4 {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}])
title('From Flap blade 2 [V]')

h1=subplot(2,2,3)
loglog(ww1/2/pi,squeeze(MAG1(2,1,:)),'color',[0.8,0.8,0.8]); hold on
loglog(ww2/2/pi,squeeze(MAG2(2,1,:)),'color',[0,0,0]); hold on
loglog(ww1/2/pi,squeeze(MAGd1(2,1,:)),'color',[0.8,0.8,0.8],'linestyle','-.'); hold on
loglog(ww2/2/pi,squeeze(MAGd2(2,1,:)),'color',[0,0,0],'linestyle','-.'); hold on
axis([2.5 50 1e-4 0.03 ])
set(h1,'YTick',[1e-4  1e-3 1e-2 ])
set(h1,'XTick',[ 2 3 4 5 6 7 8 9 10 20 30 40 50])
set(h1,'YTickLabel',[ {'{-4}'}  {'{-3}'} {'{-2}'} {'{2}'}])
set(h1,'XTickLabel',[ {''} {''} 4 {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}])
ylabel('MFC blade 2 [V]')
xlabel('Frequency [Hz]')

h1=subplot(2,2,4)
loglog(ww1/2/pi,squeeze(MAG1(2,2,:)),'color',[0.8,0.8,0.8]); hold on
loglog(ww2/2/pi,squeeze(MAG2(2,2,:)),'color',[0,0,0]); hold on
loglog(ww1/2/pi,squeeze(MAGd1(2,2,:)),'color',[0.8,0.8,0.8],'linestyle','-.'); hold on
loglog(ww2/2/pi,squeeze(MAGd2(2,2,:)),'color',[0,0,0],'linestyle','-.'); hold on
axis([2.5 50 1e-4 0.03 ])
set(h1,'YTick',[1e-4  1e-3 1e-2 ])
set(h1,'XTick',[ 2 3 4 5 6 7 8 9 10 20 30 40 50])
set(h1,'YTickLabel',[ {'{-4}'}  {'{-3}'} {'{-2}'} {'{2}'}])
set(h1,'XTickLabel',[ {''} {''} 4 {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}])
xlabel('Frequency [Hz]')
